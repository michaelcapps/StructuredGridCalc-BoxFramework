
/******************************************************************************/
/**
 * \file LevelData.cpp
 *
 * \brief Specializations for classes in LevelData.H
 *
 *//*+*************************************************************************/

#include <vector>

#ifndef NO_CGNS
#ifdef USE_MPI
#include "pcgnslib.h"
#else
#include "cgnslib.h"
#endif
#endif

#include "LevelData.H"
#include "BaseFab.H"


/*******************************************************************************
 *
 * Class LevelData: member specializations
 *
 ******************************************************************************/

#ifndef NO_CGNS
/*--------------------------------------------------------------------*/
//  Write CGNS solution data to a file (specialized for BaseFab<Real>)
/** The CGNS file must be open
 *  \param[in] a_indexFile
 *                      CGNS index of file
 *  \param[in] a_indexBase
 *                      CGNS index of base node
 *  \param[in] a_indexZoneOffset
 *                      The difference between CGNS indexZone and the
 *                      globalBoxIndex
 *  \param[in] a_varNames
 *                      Array of variable names
 *  \return             0  Success
 *                      >0 1+ the global index of the box that failed
 *  \note
 *  <ul>
 *    <li> Zones are numbered based on their global index (+1) in the
 *         layout
 *  </ul>
 *//*-----------------------------------------------------------------*/

// Local storage of indices for a Box
namespace{
struct CGNSIndices
{
  CGNSIndices()
    {
      indexField.reserve(2 + g_SpaceDim);  // Reserve estimate of space
    }
  int indexSol;
  std::vector<int> indexField;
};
}

template<>
int
LevelData<BaseFab<Real> >::writeCGNSSolData(
  const int                a_indexFile,
  const int                a_indexBase,
  const int                a_indexZoneOffset,
  const char *const *const a_varNames) const
{
  int cgerr;
  std::vector<CGNSIndices> localCGNSIndices(m_disjointBoxLayout.localSize());

//--These variables describe the shape of the array in the CGNS file.  We only
//--want to write the core grid.  The min corner of the core grid has index
//--1.

  cgsize_t rmin[g_SpaceDim];
  cgsize_t rmax[g_SpaceDim];

//--These variables describe the shape of the array in memory.  The min corner
//--of the array in memory has index 1.

  int memnumdim = g_SpaceDim;    // Rank of array in memory
  cgsize_t memdim[g_SpaceDim];   // Size of array in memory
  cgsize_t memrmin[g_SpaceDim];  // Lower limit of range to write
  cgsize_t memrmax[g_SpaceDim];  // Upper limit of range to write

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Bulk writing of the solution meta-data -- all processors write the
 * same information and note use of LayoutIterator
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  std::vector<int> indexField(ncomp());
  for (LayoutIterator lit(m_disjointBoxLayout); lit.ok(); ++lit)
    {
      const int globalBoxIndex = (*lit).globalIndex();
      const int localBoxIndex  = (*lit).localIndex();
      int indexZone = globalBoxIndex + a_indexZoneOffset;

      // Write flow solution nodes
      int indexSol;
      cgerr = cg_sol_write(a_indexFile,
                           a_indexBase,
                           indexZone,
                           "Solution",
                           CellCenter,
                           &indexSol);
      if (cgerr) return indexZone;
      // Write the meta-data to a field solution node for each component
      // (user must use SIDS-standard names here)
#ifdef USE_MPI
      for (int iComp = 0; iComp != ncomp(); ++iComp)
        {
          cgp_field_write(a_indexFile,
                          a_indexBase,
                          indexZone,
                          indexSol,
                          CGNS_REAL,
                          a_varNames[iComp],
                          &indexField[iComp]);
          if (cgerr) return indexZone;
        }
#endif
      if (DisjointBoxLayout::procID() == m_disjointBoxLayout.proc(lit))
        {
          CGNSIndices& thisCGNSIndices = localCGNSIndices[localBoxIndex];
          thisCGNSIndices.indexSol = indexSol;
#ifdef USE_MPI
          for (int iComp = 0; iComp != ncomp(); ++iComp)
            {
              thisCGNSIndices.indexField.push_back(indexField[iComp]);
            }
#endif
        }
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Collective writing of the solution data -- each processor writes
 * its own information (note use of DataIterator)
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      const int globalBoxIndex = (*dit).globalIndex();
      const int localBoxIndex  = (*dit).localIndex();
      int indexZone = globalBoxIndex + a_indexZoneOffset;
      CGNSIndices& thisCGNSIndices = localCGNSIndices[localBoxIndex];
      // Only write the core grid (not ghost cells).  Indexing in CGNS starts
      // at 1, not 0.
      const IntVect boxdim = m_disjointBoxLayout[dit].dimensions();
      const BaseFab<Real>& fab = this->operator[](dit);
      const IntVect fabdim = fab.box().dimensions();
      for (int dir = 0; dir != g_SpaceDim; ++dir)
        {
          // The range of the data in the CGNS file (only contains core grid)
          rmin[dir]    = 1;
          rmax[dir]    = boxdim[dir];
          // The size and range of data in memory
          memdim[dir]  = fabdim[dir];
          memrmin[dir] = 1 + m_nghost;
          memrmax[dir] = memdim[dir] - m_nghost;
        }
      for (int iComp = 0; iComp != ncomp(); ++iComp)
        {
#if USE_MPI
          cgerr = cgp_field_general_write_data(
            a_indexFile,
            a_indexBase,
            indexZone,
            thisCGNSIndices.indexSol,
            thisCGNSIndices.indexField[iComp],
            rmin, rmax,
            memnumdim, memdim, memrmin, memrmax,
            fab.dataPtr(iComp));
#else
          cgerr = cg_field_general_write(
            a_indexFile,
            a_indexBase,
            indexZone,
            thisCGNSIndices.indexSol,
            CGNS_REAL,
            a_varNames[iComp],
            rmin, rmax,
            memnumdim, memdim, memrmin, memrmax,
            fab.dataPtr(iComp),
            &indexField[iComp]);
#endif
          if (cgerr) return indexZone;
        }
    }

  return 0;
}
#endif  /* CGNS */
