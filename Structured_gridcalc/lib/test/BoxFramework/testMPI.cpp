#include <iostream>
#include <iomanip>
#include <sstream>
#include <unistd.h>

#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"

int main(int argc, const char* argv[])
{
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
  int status = 0;
  int mpierr;
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
  // if (procID == 1)
  //   {
  //     int i = 0;
  //     char hostname[256];
  //     gethostname(hostname, sizeof(hostname));
  //     printf("PID %d on %s ready for attach\n", getpid(), hostname);
  //     fflush(stdout);
  //     while (0 == i)
  //       sleep(5);
  //   }

//--Tests

  Box domain(IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(7, 3, 3)));
  // 2 boxes in x direction
  DisjointBoxLayout dbl(domain, 4*IntVect::Unit);

  // Box 0 should be on process 0, and Box 1 on process 1
  int c = 0;
  for (LayoutIterator lit(dbl); lit.ok(); ++lit, ++c)
    {
      if (c == 0 && dbl.proc(lit) != 0) ++status;
      if (c == 1 && dbl.proc(lit) != 1) ++status;
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
  MPI_Barrier(MPI_COMM_WORLD);
  if (verbose)
    {
      for (int iProc = 0; iProc != numProc; ++iProc)
        {
          if (iProc == procID)
            {
              for (DataIterator dit(dbl); dit.ok(); ++dit)
                {
                  BaseFab<Real>& fab = lvldata[dit];
                  int interiorCell = (procID == 0) ? 3 : 4;
                  int ghostCell = (procID == 0) ? 4 : 3;
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

  // Locations
  Box regionSend;
  Box regionRecv;
  c = 0;
  for (DataIterator dit(dbl); dit.ok(); ++dit)
    {
      const Box box = dbl[dit];
      // Process 0 works on box 0 (right or max side is of concern)
      // Process 1 works on box 1 (left or min side is of concern)
      const int side = 1 - 2*procID;
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

  // Buffers
  Real *sendBuf = new Real[regionSend.size()];
  Real *recvBuf = new Real[regionRecv.size()];

  // Pack sending
  for (DataIterator dit(dbl); dit.ok(); ++dit)
    {
      CH_assert(dbl[dit].contains(regionSend));
      const BaseFab<Real>& fab = lvldata[dit];
      int bufidx = 0;
      for (BoxIterator bit(regionSend); bit.ok(); ++bit, ++bufidx)
        {
          sendBuf[bufidx] = fab(*bit, 0);
        }
    }

  // Handles
  MPI_Request mpiRequest[2];
  // Status
  MPI_Status mpiStatus[2];

  // Post sends from this process
  int tag = 0;
  int destProc = !procID;
  mpierr = MPI_Isend(sendBuf,
                     regionSend.size(),
                     BXFR_MPI_REAL,
                     destProc,
                     tag,
                     MPI_COMM_WORLD,
                     &mpiRequest[0]);

  // Post receives to this process
  int srcProc = !procID;
  mpierr = MPI_Irecv(recvBuf,
                     regionRecv.size(),
                     BXFR_MPI_REAL,
                     srcProc,
                     tag,
                     MPI_COMM_WORLD,
                     &mpiRequest[1]);

  // Wait for messages
  mpierr = MPI_Waitall(2, mpiRequest, mpiStatus);
  if (mpierr) ++status;

  // Unpack receives
  for (DataIterator dit(dbl); dit.ok(); ++dit)
    {
      BaseFab<Real>& fab = lvldata[dit];
      CH_assert(fab.box().contains(regionRecv));
      int bufidx = 0;
      for (BoxIterator bit(regionRecv); bit.ok(); ++bit, ++bufidx)
        {
          fab(*bit, 0) = recvBuf[bufidx];
        }
    }

  // Check results
  for (int iProc = 0; iProc != numProc; ++iProc)
    {
      if (iProc == procID)
        {
          for (DataIterator dit(dbl); dit.ok(); ++dit)
            {
              BaseFab<Real>& fab = lvldata[dit];
              if (verbose)
                {
                  int interiorCell = (procID == 0) ? 3 : 4;
                  int ghostCell = (procID == 0) ? 4 : 3;
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

  // Get sum of all status into master process
  int allStatus;
  MPI_Reduce(&status, &allStatus, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  // Cleanup
  delete[] sendBuf;
  delete[] recvBuf;

//--Output status

  if (masterProc)
    {
      if (verbose)
        {
          std::cout << "Status: " << allStatus << std::endl;
        }
      const char* const testName = "testMPI";
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
