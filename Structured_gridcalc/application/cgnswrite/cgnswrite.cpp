
/******************************************************************************/
/**
 * \file cgnswrite.cpp
 *
 * \brief Test writing a linear solution to a cgns file
 *
 *//*+*************************************************************************/

#ifdef USE_MPI
#include "pcgnslib.h"
#else
#include "cgnslib.h"
#endif

#include "BoxIterator.H"
#include "LevelData.H"

int main(int argc, const char* argv[])
{
  DisjointBoxLayout::initMPI(argc, argv);
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
//--Problem Domain from (0,0,0) to (63,31,31)

  Box problemDomain(IntVect::Zero, IntVect(D_DECL(63, 31, 31)));
  // Let the origin of the coordinate system be 0
  IntVect origin = IntVect::Zero;
  // Let the mesh spacing be 1
  Real dx = 1.;

//--Create a disjointBoxLayout with boxes of size 16^3

  DisjointBoxLayout boxes(problemDomain, 16*IntVect::Unit);

//--Create a LevelData with 4 components (Density, VelocityX, VelocityY,
//--and VelocityZ) and 1 ghost cell.

  const char *const variableNames[] =
    {
      "Density",
      "VelocityX",
      "VelocityY",
      "VelocityZ"
    };
  LevelData<FArrayBox> U(boxes, 4, 1);

//--Initialize U: Give density a tri-linear profile and set velocity to zero

  U.setVal(0.);
  // IntVect (0,0,0) has density 1 and density gradient = (0.01,0.01.0.01)
  const Real rhoOrigin = 1.;
  const Real rhoGradient[g_SpaceDim] = { D_DECL(0.01, 0.01, 0.01) };
  
  for (DataIterator dit(boxes); dit.ok(); ++dit)
    {
      FArrayBox& fabU = U[dit];
      for (BoxIterator bit(fabU.box()); bit.ok(); ++bit)
        {
          const IntVect iv = *bit;
          // Assumes dx = 1
          const Real density = rhoOrigin + D_TERM(  iv[0]*dx*rhoGradient[0],
                                                  + iv[1]*dx*rhoGradient[1],
                                                  + iv[2]*dx*rhoGradient[2]);
          fabU(iv, 0) = density;
        }
    }

//--Write the data

  int cgerr;

  // Open the CGNS file
  int indexFile;
  const char* const fileName = "solution0000.cgns";
#ifdef USE_MPI
  cgerr = cgp_open(fileName, CG_MODE_WRITE, &indexFile);
#else
  cgerr = cg_open(fileName, CG_MODE_WRITE, &indexFile);
#endif
  if (cgerr)
    {
      if (verbose) cg_error_print();
      return cgerr;
    }

  // Create the base node
  int indexBase;
  int iCellDim = g_SpaceDim;
  int iPhysDim = g_SpaceDim;
  cgerr = cg_base_write(indexFile, "Base", iCellDim, iPhysDim, &indexBase);
  if (cgerr)
    {
      if (verbose) cg_error_print();
      return cgerr;
    }

  // Write the zones and coordinates (indexZoneOffset needs to be determined
  // so we can use 'indexCGNSzone = globalBoxIndex + indexZoneOffset'.  This
  // value is almost certainly 1).
  if (verbose)
    {
      std::cout << "Writing Grid Coordinates\n";
    }
  int indexZoneOffset;  // The difference between CGNS indexZone and the
                        // globalBoxIndex.
  cgerr = boxes.writeCGNSZoneGrid(indexFile,
                                  indexBase,
                                  indexZoneOffset,
                                  origin,
                                  dx);
  cg_error_print();
  if (cgerr)
    {
      if (verbose) cg_error_print();
      return cgerr;
    }

  // Write the solution data
  if (verbose)
    {
      std::cout << "Writing Solution Data\n";
    }
  cgerr = U.writeCGNSSolData(indexFile,
                             indexBase,
                             indexZoneOffset,
                             variableNames);
  if (cgerr)
    {
      if (verbose) cg_error_print();
      return cgerr;
    }

  // Close the CGNS file
#ifdef USE_MPI
  cgerr = cgp_close(indexFile);
#else
  cgerr = cg_close(indexFile);
#endif
  if (cgerr)
    {
      if (verbose) cg_error_print();
      return cgerr;
    }

  if (verbose)
    {
      std::cout << "Successfully wrote " << fileName << std::endl;
    }
  DisjointBoxLayout::finalizeMPI();
}
