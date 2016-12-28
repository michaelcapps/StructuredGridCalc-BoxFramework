
/******************************************************************************/
/**
 * \file WavePatch.cpp
 *
 * \brief Non-inline definitions for classes in WavePatch.H
 *
 *//*+*************************************************************************/

#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>

#include "cgnslib.h"

#include "BaseFabMacros.H"
#include "WavePatch.H"
#ifdef USE_GPU
#include "WavePatch_Cuda.H"
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#endif

#define USE_VEX
#ifdef USE_VEX
#include "VEXTypes.H"           // Vector types
#else
const int VecSz_r = 1;
#endif


/*******************************************************************************
 *
 * Class WavePatch: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/*--------------------------------------------------------------------*/

WavePatch::WavePatch(const Box&        a_domain,
                     const IntVect&    a_maxBoxSize,
                     const char* const a_basePlotName,
                     const Real        a_c,
                     const Real        a_dx,
                     const Real        a_cfl)
  :
  m_boxes(a_domain, a_maxBoxSize),
  m_domain(a_domain),
  m_basePlotName(a_basePlotName),
  m_c(a_c),
  m_dx(a_dx),
  m_dt(a_dx*a_cfl/a_c),
  m_time((Real)0.),
  m_iteration(0),
  m_idxStep(0),
  m_idxStepUpdate(1),
  m_idxStepOld(2)
{
  m_u[0].define(m_boxes, 1, 1);
  m_u[1].define(m_boxes, 1, 1);
  m_u[2].define(m_boxes, 1, 1);
  DataIterator dit(m_boxes);
  m_bidx = *dit;
#ifdef USE_GPU
  // Pack the BaseFabs for each time index into an array
  const BaseFab<Real>* patchData[3];
  patchData[0] = &(m_u[0][m_bidx]);
  patchData[1] = &(m_u[1][m_bidx]);
  patchData[2] = &(m_u[2][m_bidx]);
  WavePatch_Cuda::construct(m_domain,
                            patchData,
                            m_cudaFab_device,
                            m_workBoxesRHS_device,
                            m_workBoxInfoBC_device,
                            m_numBlkRHS,
                            m_numBlkBC);
#endif
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/

WavePatch::~WavePatch()
{
#ifdef USE_GPU
  WavePatch_Cuda::destroy(m_cudaFab_device,
                          m_workBoxesRHS_device,
                          m_workBoxInfoBC_device);
#endif
}

/*--------------------------------------------------------------------*/
//  Set initial data to pulse
/** \param[in]  a_rho   Initial density (default 1)
 *//*-----------------------------------------------------------------*/

void
WavePatch::initialData()
{
  constexpr Real pi = 3.141592653589793;

  MD_ARRAY_RESTRICT(arrun, un());
  MD_BOXLOOP_OMP(m_domain, i)
    {
      D_TERM(const Real x = (i0 + 0.5)*m_dx;,
             const Real y = (i1 + 0.5)*m_dx;,
             const Real z = (i2 + 0.5)*m_dx;);
      const Real r = std::sqrt(D_TERM(  std::pow(x - 0.5, 2),
                                      + std::pow(y - 0.5, 2),
                                      + std::pow(z - 0.5, 2)));
      Real val = 0.;
      if (r <= 0.5)
        {
          val = std::pow(sin(pi*(r + 0.5)), 6);
        }
      arrun[MD_IX(i, 0)] = val;
    }
  unm1().copy(m_domain, un());
}

/*--------------------------------------------------------------------*/
//  Advance one time step
/** Compute unp1() and time n+1 from un() and unm1().
 *//*-----------------------------------------------------------------*/

void
WavePatch::advance()
{
  m_timerAdvance.start();

// Setup for VEX
#ifdef USE_VEX
  const int pencilSize  = m_domain.dimensions()[0];
  const int vecPacked   = pencilSize/VecSz_r;
  const int i0EndPacked = m_domain.loVect(0) + vecPacked*VecSz_r;
  (void)i0EndPacked;
#endif

//--Set BC

#ifdef USE_GPU
  un().copyToDevice();
  WavePatch_Cuda::driverBC(m_numBlkBC, m_idxStep);
  un().copyToHost();
#else
  for (int dir = 0; dir != g_SpaceDim; ++dir)
    {
      for (int side = -1; side < 2; side += 2)
        {
          Box srcBox(m_domain);
          srcBox.adjBox(-1, dir, side);
          Box dstBox(srcBox);
          dstBox.shift(side, dir);
          un().copy(dstBox, 0, un(), srcBox, 0, 1);
        }
    }
#endif

//--Update solution

  const Real factor = std::pow(m_dt*m_c/m_dx, 2)/g_SpaceDim;
#ifdef USE_GPU
  un().copyToDevice();
  unm1().copyToDevice();
  WavePatch_Cuda::driverRHS(m_numBlkRHS,
                            m_idxStep,
                            m_idxStepUpdate,
                            m_idxStepOld,
                            factor);
  unp1().copyToHost();
#else
  MD_ARRAY_RESTRICT(arrunp1, unp1());
  MD_ARRAY_RESTRICT(arrun, un());
  MD_ARRAY_RESTRICT(arrunm1, unm1());

#ifdef USE_VEX
  const __mvr two_vr = _mm_vr(set1)(2.0);
  const __mvr factor_vr = _mm_vr(set1)(factor);
  MD_BOXLOOP_PENCIL_OMP(m_domain, i)
    {
      int i0 = m_domain.loVect(0);
      for (; i0 <= i0EndPacked; i0 += VecSz_r)
        {
          const __mvr unp1_vr =
            two_vr*_mm_vr(loadu)(&arrun[MD_IX(i, 0)]) -
                   _mm_vr(loadu)(&arrunm1[MD_IX(i, 0)]) + factor_vr*
            MD_DIRSUM([=](const int            a_dir,
                          MD_DECLIX(const int, a_o))
              {
                MD_CAPTURE_RESTRICT(arrun);
                return
                         _mm_vr(loadu)(&arrun[MD_OFFSETIX(i,+,a_o, 0)]) -
                  two_vr*_mm_vr(loadu)(&arrun[MD_IX(i, 0)]) +
                         _mm_vr(loadu)(&arrun[MD_OFFSETIX(i,-,a_o, 0)]);
              });
          _mm_vr(storeu)(&arrunp1[MD_IX(i, 0)], unp1_vr);
        }
      // Catch unpacked cells
      for (; i0 <= m_domain.hiVect(0); ++i0)
        {
          arrunp1[MD_IX(i, 0)] =
            2*arrun[MD_IX(i, 0)] - arrunm1[MD_IX(i, 0)] + factor*
            MD_DIRSUM([=](const int a_dir,
                          MD_DECLIX(const int, a_o))
              {
                MD_CAPTURE_RESTRICT(arrun);
                return
                    arrun[MD_OFFSETIX(i,+,a_o, 0)] -
                  2*arrun[MD_IX(i, 0)] +
                    arrun[MD_OFFSETIX(i,-,a_o, 0)];
              });
        }
    }
#else

  // Time terms
  {
    MD_BOXLOOP_OMP(m_domain, i)
      {
        arrunp1[MD_IX(i, 0)] = 2*arrun[MD_IX(i, 0)] - arrunm1[MD_IX(i, 0)];
      }
  }

  // Laplacian
  for (int dir = 0; dir != g_SpaceDim; ++dir)
    {
      const int MD_ID(o, dir);
      MD_BOXLOOP_OMP(m_domain, i)
        {
          arrunp1[MD_IX(i, 0)] += factor*(arrun[MD_OFFSETIX(i,+,o, 0)] -
                                          2*arrun[MD_IX(i, 0)] +
                                          arrun[MD_OFFSETIX(i,-,o, 0)]);
        }         
    }
#endif  /* !VEX */
#endif  /* !GPU */

//--Swap indices (unp1->un, un->unm1)

  advanceStepIndex();
  ++m_iteration;
  m_time += m_dt;
  m_timerAdvance.stop();
}

/*--------------------------------------------------------------------*/
//  Advance a group of time steps using GPU
/** 
 *//*-----------------------------------------------------------------*/

#ifdef USE_GPU
void
WavePatch::advanceIterGroup(const int   a_numIter,
                            cudaEvent_t a_cuEvent_iterGroupStart,
                            cudaEvent_t a_cuEvent_iterGroupEnd)
{
  m_timerAdvance.start();
  const Real factor = std::pow(m_dt*m_c/m_dx, 2)/g_SpaceDim;
  WavePatch_Cuda::driverAdvanceIterGroup(a_numIter,
                                         m_numBlkBC,
                                         m_numBlkRHS,
                                         m_idxStep,
                                         m_idxStepUpdate,
                                         m_idxStepOld,
                                         factor,
                                         a_cuEvent_iterGroupStart,
                                         a_cuEvent_iterGroupEnd);
  m_iteration += a_numIter;
  m_time += a_numIter*m_dt;
  m_timerAdvance.stop();
}
#endif

/*--------------------------------------------------------------------*/
//  Write the plot file
/** Write un().
 *  \param[in]  a_idxStep
 *                      Index of solution in time to write (current is
 *                      given by m_idxStep)
 *  \param[in]  a_iteration
 *                      Index of iteration to write
 *  \return             -1 Error
 *                       0 Success
 *                      >0 CGNS error
 *//*-----------------------------------------------------------------*/

int
WavePatch::writePlotFile(const int a_idxStep, const int a_iteration) const
{
  m_timerWrite.start();
#ifndef NO_CGNS
  int cgerr;
  std::ostringstream fileName;
  fileName << m_basePlotName << std::setw(6) << std::setfill('0')
           << a_iteration << ".cgns";

  // Open the CGNS file
  int indexFile;
  cgerr = cg_open(fileName.str().c_str(), CG_MODE_WRITE, &indexFile);
  if (cgerr)
    {
      cg_error_print();
      m_timerWrite.stop();
      return cgerr;
    }

  // Create the base
  int indexBase;
  int iCellDim = g_SpaceDim;
  int iPhysDim = g_SpaceDim;
  cgerr = cg_base_write(indexFile, "Base", iCellDim, iPhysDim, &indexBase);
  if (cgerr)
    {
      cg_error_print();
      m_timerWrite.stop();
      return cgerr;
    }

//--Write the zones and grids

  // Write the grid coordinates
  int indexZoneOffset;  // The difference between CGNS indexZone and the
                        // globalBoxIndex
  cgerr = m_boxes.writeCGNSZoneGrid(indexFile,
                                    indexBase,
                                    indexZoneOffset,
                                    m_domain.loVect(),
                                    m_dx);
  if (cgerr)
    {
      std::cout << "EE Failed to write zone and grid for box " << cgerr-1 << '!'
                << std::endl;
      m_timerWrite.stop();
      return -1;
    }
  
  // Write the solution data
  static const char *const stateNames[] = { "displacement" };
  cgerr = m_u[a_idxStep].writeCGNSSolData(indexFile,
                                          indexBase,
                                          indexZoneOffset,
                                          stateNames);
  if (cgerr)
    {
      std::cout << "EE Failed to write solution for box " << cgerr-1 << '!'
                << std::endl;
      m_timerWrite.stop();
      return -1;
    }

  // Close the CGNS file
  cgerr = cg_close(indexFile);

#endif  /* CGNS */
  m_timerWrite.stop();
  return 0;
}
