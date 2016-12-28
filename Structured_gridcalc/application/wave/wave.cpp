
/******************************************************************************/
/**
 * \file latticeBoltzmann.cpp
 *
 * \brief Lattice-Boltzmann solution to a problem
 *
 *//*+*************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "LinuxSupport.H"
#include "IntVect.H"
#include "Box.H"
#include "WavePatch.H"
#include "Stopwatch.H"

#ifdef USE_GPU
#include "CudaSupport.H"
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#endif

static const char *const usage =
  "Usage ./wave [-np x] [h [i]]\n"
  "  x : number of threads for OpenMP.  You can also use\n"
  "      'export OMP_NUM_THREADS=x' to use x threads with OpenMP.\n"
  "  h : domain dimensions in y and z (multiple of 32, default=32).\n"
  "  i : number of iterations (i > 0, default=4000*(h/32)).\n"
  "\n  Use 'export OMP_PROC_BIND=TRUE' to lock thread affinity in OpenMP.\n";

//--Prototypes

void writeDotsForBlockIter(const int  a_maxDot,
                           int&       a_numDot,
                           int        a_numDotPrint,
                           const int  a_iteration,
                           const Real a_time);

/*----------------------------------------------------------------------------*/

int main(int argc, const char* argv[])
{
  if (argc > 1 && (std::strcmp(argv[1], "-h") == 0 ||
                   std::strcmp(argv[1], "--help") == 0))
    {
      std::cout << usage;
      return 0;
    }
  Stopwatch<> timerTotal;
  timerTotal.start();

//--Input parameters for the run

  bool badArg = false;
  int iargc = 1;
  while (argc > iargc && argv[iargc][0] == '-')
    {
      if (std::strcmp(argv[iargc], "-np") == 0)
        {
#ifdef _OPENMP
          omp_set_num_threads(std::atoi(argv[iargc+1]));
#else
          std::cout << "OpenMP is not enabled!" << std::endl;
          badArg = true;
#endif
          iargc += 2;
        }
    }

  int h_in = 64;  // 32 is baseline
  CH_assert(h_in % 32 == 0);
  int numIter_in = 1800*(h_in/32);
  {
    int idxArgVal = 0;
    for (int iarg = iargc; iarg < std::min(iargc + 2, argc); ++iarg)
    {
      switch (idxArgVal)
        {
        case 0:
          h_in = std::atoi(argv[iarg]);
          numIter_in = 1800*(h_in/32);
          break;
        case 1:
          numIter_in = std::atoi(argv[iarg]);
          break;
        }
      ++idxArgVal;
    }
  }
  if (badArg)
    {
      std::cout << usage;
      return 1;
    }

//--Other parameters

  const int h = h_in;
  const int numIter = numIter_in;
  const char *const plotDir = "plot";
  const char *const plotFileBase = "plot/plot.";
  // For canonical h = 32
  const int iterBlock = 100*(h/32);
  const int plotFreq = iterBlock;
  const IntVect domainSize(D_DECL(h, h, h));
  const int boxSize = h;
  const Real dx = 1./h;
  const Real c = 1.;
  const Real cfl = 0.125;

//--Check the parameters

  int paramErr = 0;
  if (numIter < 0)
    {
      std::cout << "Number of iterations must be >= 0!" << std::endl;
      ++paramErr;
    }
  if (!System::fileExists(plotDir))
    {
      std::cout << "Directory \'" << plotDir << "\' does not exist!"
                << std::endl;
      ++paramErr;
    }
  if (plotFreq <= 0)
    {
      std::cout << "Plot frequency must be > 0!" << std::endl;
      ++paramErr;
    }
  for (int dir = 0; dir != g_SpaceDim; ++dir)
    {
      if (domainSize[dir] % boxSize != 0)
        {
          std::cout << "Domain size must be a multiple of boxSize "
            "(currently set to " << boxSize << ")!" << std::endl;
          ++paramErr;
        }
    }
  if (numIter < 0)
    {
      std::cout << "Number of iterations must be >= 0!" << std::endl;
      ++paramErr;
    }
  if (paramErr)
    {
      std::cout << usage;
      return 1;
    }

//--Write information about the run

  std::cout << std::left << std::setw(40) << "Box size: " << boxSize
            << std::endl;
  std::cout << std::left << std::setw(40) << "Iterations: " << numIter
            << std::endl;
  std::cout << std::left << std::setw(40) << "Domain: "
            << Box(IntVect::Zero, domainSize-IntVect::Unit)
            << std::endl;
  std::cout << std::left << std::setw(40) << "Plot frequency: "
            << plotFreq << std::endl;
  std::cout << std::left << std::setw(40) << "Precision: "
            << 8*sizeof(Real) << " bits\n";
#ifdef _OPENMP
  std::cout << std::left << std::setw(40) << "Parallel OpenMP threads: "
            << omp_get_max_threads() << std::endl;
#endif
  std::cout << std::endl;

//--Data structures

  const Box domain(IntVect::Zero, domainSize - IntVect::Unit);
  WavePatch patchSolver(domain,
                        boxSize*IntVect::Unit,
                        plotFileBase,
                        c,
                        dx,
                        cfl);

//--Initialize data

  patchSolver.initialData();

//--Write the first plot file if numIter == 0

  if (numIter == 0)
    {
      patchSolver.writePlotFile(patchSolver.currentStepIndex(),
                                patchSolver.iteration());
    }

//--Run the solver

  int numDot = 0;
  const int maxDot = 80;
  bool writeOnFirstSingleIter = true;

//--Run a group of iterations

  int numGroupIter = 0; (void)numGroupIter;
#ifdef USE_GPU
  if (plotFreq % 3 == 0)  // Can copy next solution and write current
    {                     // at same time
      std::cout << "Unable to overlap copy with plot writing." << std::endl;
    }
  Stopwatch<> timerIdleCPU;
  // Timings on GPU
  cudaEvent_t cuEvent_iterGroupStart;
  cudaEvent_t cuEvent_iterGroupEnd;
  CU_SAFE_CALL(cudaEventCreate(&cuEvent_iterGroupStart));
  CU_SAFE_CALL(cudaEventCreate(&cuEvent_iterGroupEnd));
  float gpuTimeAdvance = 0.f;
  cudaEvent_t cuEvent_copyStart;
  cudaEvent_t cuEvent_copyEnd;
  CU_SAFE_CALL(cudaEventCreate(&cuEvent_copyStart));
  CU_SAFE_CALL(cudaEventCreate(&cuEvent_copyEnd));
  float gpuCopy = 0.f;
  // Copy to device
  patchSolver.copyToDeviceAsync(patchSolver.currentStepIndex());
  patchSolver.copyToDeviceAsync(patchSolver.oldStepIndex());
  numGroupIter = numIter/plotFreq;
  for (int iIterGroup = 0; iIterGroup != numGroupIter; ++iIterGroup)
    {
      // Save the current step index and iteration
      const int currentStepIndex = patchSolver.currentStepIndex();
      const int currentIteration = patchSolver.iteration();
      const Real currentTime = patchSolver.time();
      // Launch the kernels
      patchSolver.advanceIterGroup(plotFreq,
                                   cuEvent_iterGroupStart,
                                   cuEvent_iterGroupEnd);
      // Copy results back from GPU
      cudaEventRecord(cuEvent_copyStart, 0);
      // std::cout << "Copying to " << patchSolver.currentStepIndex()
      //           << " at iter " << patchSolver.iteration() << std::endl;
      patchSolver.copyToHostAsync(patchSolver.currentStepIndex());
      cudaEventRecord(cuEvent_copyEnd, 0);
      // Write the plot file (from previous block)
      // std::cout << "Writing from " << currentStepIndex
      //           << " at iter " << currentIteration << std::endl;
      patchSolver.writePlotFile(currentStepIndex, currentIteration);
      // Write progress (complete once plotfile is written)
      if (iIterGroup != 0)
        {
          writeDotsForBlockIter(maxDot,
                                numDot,
                                plotFreq,
                                currentIteration,
                                currentTime);
        }
      // Synchronize
      timerIdleCPU.start();
      CU_SAFE_CALL(cudaDeviceSynchronize());
      // std::cout << "Synchronized" << std::endl;
      timerIdleCPU.stop();
      float elapsed;
      CU_SAFE_CALL(
        cudaEventElapsedTime(&elapsed,
                             cuEvent_iterGroupStart,
                             cuEvent_iterGroupEnd));
      gpuTimeAdvance += elapsed;
      CU_SAFE_CALL(
        cudaEventElapsedTime(&elapsed,
                             cuEvent_copyStart,
                             cuEvent_copyEnd));
      gpuCopy += elapsed;
    }
  if (numGroupIter > 0)
    {
      // Write progress from last iteration group
      writeDotsForBlockIter(maxDot,
                            numDot,
                            plotFreq,
                            patchSolver.iteration(),
                            patchSolver.time());
      // Write the plot file from last iteration group
      // std::cout << "Writing from " << patchSolver.currentStepIndex()
      //           << " at iter " << patchSolver.iteration() << std::endl;
      patchSolver.writePlotFile(patchSolver.currentStepIndex(),
                                patchSolver.iteration());
      writeOnFirstSingleIter = false;
    }
#endif

//--Run single iterations

  const int numSingleIter = numIter - patchSolver.iteration();
  for (int iter = patchSolver.iteration(); iter != numIter; ++iter)
    {
      if (iter % plotFreq == 0 && writeOnFirstSingleIter)
        {
          if (numDot != 0)
            {
              std::cout << std::endl;
              numDot = 0;
            }
          std::cout << "Time step " << std::setw(6)
                    << patchSolver.iteration()
                    << " Old time " << std::scientific
                    << patchSolver.time() << std::endl;
          patchSolver.writePlotFile(patchSolver.currentStepIndex(),
                                    patchSolver.iteration());
        }
      else
        {
          std::cout.put('.');
          std::cout.flush();
          ++numDot;
          if (numDot == maxDot)
            {
              std::cout << std::endl;
              numDot = 0;
            }
        }
      writeOnFirstSingleIter = true;
      patchSolver.advance();
    }

//--Write the final plot file

  if (numDot)
    {
      std::cout << std::endl;
    }
  if (numIter != 0)
    {
      std::cout << "Solution finished " << std::setw(6)
                << patchSolver.iteration() << " iterations at time "
                << std::scientific << patchSolver.time() << std::endl;
      if (numSingleIter != 0)
        {
          patchSolver.writePlotFile(patchSolver.currentStepIndex(),
                                    patchSolver.iteration());
        }
    }

//--Write some times

  timerTotal.stop();
  std::cout << std::endl;
#ifdef USE_GPU
  std::cout << std::left << std::setw(40) << "Time for GPU advance (ms): "
            << gpuTimeAdvance << std::endl;
  CU_SAFE_CALL(cudaEventDestroy(cuEvent_iterGroupStart));
  CU_SAFE_CALL(cudaEventDestroy(cuEvent_iterGroupEnd));
  std::cout << std::left << std::setw(40) << "Time for GPU copy (ms): "
            << gpuCopy << std::endl;
  CU_SAFE_CALL(cudaEventDestroy(cuEvent_copyStart));
  CU_SAFE_CALL(cudaEventDestroy(cuEvent_copyEnd));
  std::cout << std::left << std::setw(40) << "Time CPU idle (ms): "
            << timerIdleCPU.time() << std::endl;
#endif
  std::cout << std::left << std::setw(40) << "Time for CPU advance (ms): "
            << patchSolver.m_timerAdvance.time() << std::endl;
  std::cout << std::left << std::setw(40)
            << "Time for writing plot files (ms): "
            << patchSolver.m_timerWrite.time() << std::endl;
  std::cout << std::left << std::setw(40) << "Total run time (ms): "
            << timerTotal.time() << std::endl;

//--Done

  return 0;
}

/*--------------------------------------------------------------------*/
/// Write dots for a block of iterations performed on GPU
/** \param[in]  maxDot  Maximum dots per row
 *  \param[in]  numDot  Dots written on row at entry
 *  \param[out] numDot  Dots written on row at exit
 *  \param[in]  numDotPrint
 *                      Number of dots to print
 *  \param[in]  a_iteration
 *                      Index of current iteration
 *              a_time  Current time
 *//*-----------------------------------------------------------------*/

void writeDotsForBlockIter(const int  a_maxDot,
                           int&       a_numDot,
                           int        a_numDotPrint,
                           const int  a_iteration,
                           const Real a_time)
{
  while (a_numDotPrint)
    {
      const int numDotPrintLine = std::min(a_maxDot - a_numDot,
                                           a_numDotPrint);
      for (int n = numDotPrintLine; n--;) std::cout.put('.');
      a_numDotPrint -= numDotPrintLine;
      a_numDot += numDotPrintLine;
      if (a_numDot == a_maxDot)
        {
          std::cout << std::endl;
          a_numDot = 0;
        }
    }
  std::cout.flush();
  if (a_numDot != 0)
    {
      std::cout << std::endl;
      a_numDot = 0;
    }
  std::cout << "Time step " << std::setw(6) << a_iteration
            << " Old time " << std::scientific << a_time << std::endl;
}

