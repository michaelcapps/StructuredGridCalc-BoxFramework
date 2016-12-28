#include <iostream>
#include <iomanip>
#include <sstream>

#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"

#define STUDENTSETUP

int main(int argc, const char* argv[])
{
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
  int status = 0;
  // Use this for output so I/O from processes do not get mixed.
  std::ostringstream ost;

//--Initialize MPI

  DisjointBoxLayout::initMPI(argc, argv);
  int numProc = DisjointBoxLayout::numProc();
  int procID = DisjointBoxLayout::procID();
  const bool masterProc = (procID == 0);

  if (numProc != 2)
    {
       if (masterProc)
         {
           std::cout << "Error: this test must be run with 2 processes!\n";
         }
       MPI_Abort(MPI_COMM_WORLD, 1);
    }
  if (masterProc)
    {
      if (verbose) std::cout << "Using " << numProc << " processors\n";
    }

//--Tests

  Box domain(IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(7, 3, 3)));
  // 2 boxes in x direction
  DisjointBoxLayout dbl(domain, 4*IntVect::Unit);

#ifdef STUDENTSETUP
  // Box 0 should be on process 0, and Box 1 on process 1
  int c = 0;
  for (LayoutIterator lit(dbl); lit.ok(); ++lit, ++c)
    {
      if (c == 0 && dbl.proc(lit) != 0) ++status;
      if (c == 1 && dbl.proc(lit) != 1) ++status;
#else
  // Box 1 should be on process 0, and Box 0 on process 1
  int c = 0;
  for (LayoutIterator lit(dbl); lit.ok(); ++lit, ++c)
    {
      if (c == 0 && dbl.proc(lit) != 1) ++status;
      if (c == 1 && dbl.proc(lit) != 0) ++status;
#endif
      if (masterProc && verbose)
        {
          std::cout << "Box " << c << " is on process " << dbl.proc(lit)
                    << std::endl;
        }
    }

  LevelData<BaseFab<Real> > lvldata(dbl, 1, 1);

  // Initialize
  lvldata.setVal(0.);
  for (DataIterator dit(dbl); dit.ok(); ++dit, ++c)
    {
      lvldata[dit].setVal(procID + 0.5);
    }
  if (verbose)
    {
      for (int iProc = 0; iProc != numProc; ++iProc)
        {
          if (iProc == procID)
            {
              for (DataIterator dit(dbl); dit.ok(); ++dit)
                {
                  BaseFab<Real>& fab = lvldata[dit];
#ifdef STUDENTSETUP
                  int interiorCell = (procID == 0) ? 3 : 4;
                  int ghostCell = (procID == 0) ? 4 : 3;
#else
                  int interiorCell = (procID == 0) ? 4 : 3;
                  int ghostCell = (procID == 0) ? 3 : 4;
#endif
                  ost << "Proc " << procID
                      << "\n  Interior: "
                      << fab(IntVect(D_DECL(interiorCell, 0, 0)), 0)
                      << "\n  Ghost: "
                      << fab(IntVect(D_DECL(ghostCell, 0, 0)), 0)
                      << std::endl;
                  std::cout << ost.str();
                  ost.str("");
                }
            }
          MPI_Barrier(MPI_COMM_WORLD);
        }
    }

#if 1
  // Exchange ghosts
  Copier copier;
  copier.defineExchangeLD(lvldata);
  lvldata.exchangeBegin(copier);
  // Do some computational work during communication
  double dtmp = 2.5*2;
  lvldata.exchangeEnd(copier);

  // Check results
  Box regionSend;
  Box regionRecv;
  c = 0;
  for (DataIterator dit(dbl); dit.ok(); ++dit)
    {
      const Box box = dbl[dit];
#ifdef STUDENTSETUP
      // Process 0 works on box 0 (right or max side is of concern)
      // Process 1 works on box 1 (left or min side is of concern)
      const int side = 1 - 2*procID;
#else
      // Process 0 works on box 1 (left or min side is of concern)
      // Process 1 works on box 0 (right or max side is of concern)
      const int side = 2*procID - 1;
#endif
      regionSend = box;
      regionSend.adjBox(-1, 0, side);
      regionRecv = box;
      regionRecv.adjBox(1, 0, side);
    }
  if (verbose)
    {
      for (int iProc = 0; iProc != numProc; ++iProc)
        {
          if (iProc == procID)
            {
              ost << "Proc " << procID
                  << "\n  RegionSend" << regionSend
                  << "\n  RegionRecv" << regionRecv << std::endl;
              std::cout << ost.str();
              ost.str("");
            }
          MPI_Barrier(MPI_COMM_WORLD);
        }
    }
  for (int iProc = 0; iProc != numProc; ++iProc)
    {
      if (iProc == procID)
        {
          for (DataIterator dit(dbl); dit.ok(); ++dit)
            {
              BaseFab<Real>& fab = lvldata[dit];
              if (verbose)
                {
#ifdef STUDENTSETUP
                  int interiorCell = (procID == 0) ? 3 : 4;
                  int ghostCell = (procID == 0) ? 4 : 3;
#else
                  int interiorCell = (procID == 0) ? 4 : 3;
                  int ghostCell = (procID == 0) ? 3 : 4;
#endif
                  ost << "Proc " << procID
                      << "\n  Interior: "
                      << fab(IntVect(D_DECL(interiorCell, 0, 0)), 0)
                      << "\n  Ghost: "
                      << fab(IntVect(D_DECL(ghostCell, 0, 0)), 0)
                      << std::endl;
                  std::cout<< ost.str();
                  ost.str("");
                }
              // Expect value from other processor
              const Real val = (procID == 0) ? 1.5 : 0.5;
              for (BoxIterator bit(regionRecv); bit.ok(); ++bit)
                {
                  if (fab(*bit, 0) != val) ++status;
                }
            }
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

  // Get sum of all status into master process
  int allStatus;
  MPI_Reduce(&status, &allStatus, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

//--Output status

  if (masterProc)
    {
      if (verbose)
        {
          std::cout << "Status: " << allStatus << std::endl;
        }
      const char* const testName = "testMPIExchange";
      const char* const statLbl[] = {
        "failed",
        "passed"
      };
      std::cout << std::left << std::setw(40) << testName
                << statLbl[(status == 0)] << std::endl;
    }

  // Finalize MPI
  DisjointBoxLayout::finalizeMPI();
  return status;
}
