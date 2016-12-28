
#ifndef _DISJOINTBOXLAYOUT_H_
#define _DISJOINTBOXLAYOUT_H_


/******************************************************************************/
/**
 * \file DisjointBoxLayout.H
 *
 * \brief Class describing a layout of non-overlapping boxes in a domain
 *
 *//*+*************************************************************************/

#include <memory>
#include <vector>

#include "Parameters.H"
#include "BoxIndex.H"
#include "Box.H"

//--Forward declarations

class LayoutIterator;


/*******************************************************************************
 */
///  Disjoint (non-overlapping) layout of boxes
/**
 *   \note
 *   <ul>
 *     <li> Most copying and assignment only performs a shallow copy of the
 *          array of boxes
 *   </ul>
 *
 *//*+*************************************************************************/

class DisjointBoxLayout
{

//--Types

  struct BoxEntry
  {
    Box box;
    int proc;
  };

//--Friends

  friend class NeighborIterator;
  friend class PeriodicIterator;


/*====================================================================*
 * Public constructors and destructors
 *====================================================================*/

public:

  /// Default constructor
  DisjointBoxLayout();

  /// Constructor
  DisjointBoxLayout(const Box& a_domain, const IntVect& a_maxBoxSize);

  /// Copy constructor
  DisjointBoxLayout(const DisjointBoxLayout&) = default;

  /// Move constructor
  DisjointBoxLayout(DisjointBoxLayout&&) noexcept = default;

  /// Assignment constructor
  DisjointBoxLayout& operator=(const DisjointBoxLayout&) = default;

  /// Move assignment constructor
  DisjointBoxLayout& operator=(DisjointBoxLayout&&) noexcept = default;

  /// Destructor
  ~DisjointBoxLayout() = default;

  /// Define (weak construction)
  void define(const Box& a_domain, const IntVect& a_maxBoxSize);

  /// Define with deep copy
  void defineDeepCopy(const DisjointBoxLayout& a_dbl);


/*====================================================================*
 * Members functions
 *====================================================================*/

public:

  /// Get the problem domain
  const Box& problemDomain() const;

  /// Index with a LayoutIterator
  const Box& operator[](const LayoutIterator& a_layit) const;

  /// Index with a BoxIndex
  const Box& operator[](const BoxIndex& a_bidx) const;

  /// Find local process ID of a box with a LayoutIterator
  int proc(const LayoutIterator& a_layit) const;

  /// Find local process ID of a box with a BoxIndex
  int proc(const BoxIndex& a_bidx) const;

  /// Get box and local process ID with a LayoutIterator
  const Box& box(const LayoutIterator& a_layit, int& a_proc) const;

  /// Number of boxes
  int size() const;

  /// Number of boxes on this process
  int localSize() const;

  /// Return a BoxIndex from a linear index (starting at 0)
  BoxIndex dataIndex(const int a_idx) const;

  /// Number of boxes in each direction
  const IntVect& dimensions() const;

  /// Unique identifying tag for the DBL
  size_t tag() const;

#ifndef NO_CGNS
  /// Write CGNS zone and grid to a file
  int writeCGNSZoneGrid(const int      a_indexFile,
                        const int      a_indexBase,
                        int&           a_indexZoneOffset,
                        const IntVect& a_origin,
                        const Real     a_dx) const;
#endif

  /// Const access to a BoxEntry with a linear index
  //  FOR INTERNAL USE AND TESTING ONLY
  const BoxEntry& getLinear(const int a_idx) const;

  /// Access a BoxEntry with a linear index
  //  FOR INTERNAL USE AND TESTING ONLY
  BoxEntry& getLinear(const int a_idx);

  /// Linear offset to a neighbour based on an IntVect
  int linearNbrOffset(const IntVect& a_nbrOffset) const;

  /// Begin linear index into local boxes
  int localIdxBegin() const;

  /// End linear index into local boxes (gives one past last local box)
  int localIdxEnd() const;

  /// Initialize MPI
  static void initMPI(int argc, const char* argv[]);

  /// Finalize MPI
  static void finalizeMPI();

  /// Number of processes
  static int numProc();

  /// ID of this process
  static int procID();


/*====================================================================*
 * Data members
 *====================================================================*/

protected:

  Box m_domain;                       ///< Box describing the domain
  IntVect m_stride;                   ///< Stride for finding neighbour boxes
  IntVect m_numBox;                   ///< Number of boxes in each direction
  int m_size;                         ///< Total number of boxes
  std::shared_ptr<std::vector<BoxEntry> > m_boxes;
                                      ///< Array of boxes
  int m_localIdxBeg;                  ///< Begin index of boxes local to this
                                      ///< processes in m_boxes
  int m_numLocalBox;                  ///< Number of boxes local to this process

  static int s_numProc;               ///< Total number of processes
  static int s_procID;                ///< ID for this process
};


/*******************************************************************************
 *
 * Class DisjointBoxLayout: inline member definitions
 *
 * Note: a few inline definitions that need a complete specification of
 *   LayoutIterator are in LayoutIterator.H:
 *   - const Box& DisjointBoxLayout::operator[](const LayoutIterator& a_layit)
 *     const;
 *   - int DisjointBoxLayout::proc(const LayoutIterator& a_layit) const;
 *   - const Box& DisjointBoxLayout::box(const LayoutIterator& a_layit,
 *                                       int&                  a_proc) const
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get the problem domain
/*--------------------------------------------------------------------*/

inline const Box&
DisjointBoxLayout::problemDomain() const
{
  return m_domain;
}

/*--------------------------------------------------------------------*/
//  Index with a BoxIndex
/** This is less safe and is discouraged because it is not possible to
 *  check that the Iterator from which the BoxIndex was obtained
 *  matches the Box layout.
 *//*-----------------------------------------------------------------*/

inline const Box&
DisjointBoxLayout::operator[](const BoxIndex& a_bidx) const
{
  return (*m_boxes)[a_bidx.globalIndex()].box;
}

/*--------------------------------------------------------------------*/
//  Find local process ID of a box with a BoxIndex
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::proc(const BoxIndex& a_bidx) const
{
  return (*m_boxes)[a_bidx.globalIndex()].proc;
}

/*--------------------------------------------------------------------*/
//  Number of boxes
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::size() const
{
  return m_size;
}

/*--------------------------------------------------------------------*/
//  Number of boxes on this process
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::localSize() const
{
  return m_numLocalBox;
}

/*--------------------------------------------------------------------*/
//  Return a BoxIndex from a linear index (starting at 0)
/** This function should only be used to access local boxes with
 *  OpenMP
 *  \param[in]  a_idx   Index in range (0 : localSize()-1)
 *  \return             A valid BoxIndex for this range
 *//*-----------------------------------------------------------------*/

inline BoxIndex
DisjointBoxLayout::dataIndex(const int a_idx) const
{
  return BoxIndex(m_localIdxBeg + a_idx, a_idx);
}

/*--------------------------------------------------------------------*/
//  Number of boxes in each direction
/*--------------------------------------------------------------------*/

inline const IntVect&
DisjointBoxLayout::dimensions() const
{
  return m_numBox;
}

/*--------------------------------------------------------------------*/
//  Unique identifying tag for the DBL
/*--------------------------------------------------------------------*/

inline size_t
DisjointBoxLayout::tag() const
{
  CH_assert(m_boxes);
  return reinterpret_cast<size_t>(m_boxes.get());
}

/*--------------------------------------------------------------------*/
//  Const access to a BoxEntry with a linear index
/** FOR INTERNAL USE AND TESTING ONLY
 *  \param[in] a_idx    Linear index into the array of boxes
 *  \return             BoxEntry at that index
 *//*-----------------------------------------------------------------*/

inline const DisjointBoxLayout::BoxEntry&
DisjointBoxLayout::getLinear(const int a_idx) const
{
  CH_assert(a_idx >= 0 && a_idx < m_size);
  CH_assert(m_boxes);
  return m_boxes->operator[](a_idx);
}

/*--------------------------------------------------------------------*/
//  Access a BoxEntry with a linear index
/** FOR INTERNAL USE AND TESTING ONLY
 *  \param[in] a_idx    Linear index into the array of boxes
 *  \return             BoxEntry at that index
 *//*-----------------------------------------------------------------*/

inline DisjointBoxLayout::BoxEntry&
DisjointBoxLayout::getLinear(const int a_idx)
{
  CH_assert(a_idx >= 0 && a_idx < m_size);
  CH_assert(m_boxes);
  return m_boxes->operator[](a_idx);
}

/*--------------------------------------------------------------------*/
//  Linear offset to a neighbour based on an IntVect
/** \param[in] a_nbrOffset
 *                      Offset to the neighbour described by an
 *                      IntVect
 *  \return             Linear offset to the neighbour
 *//*-----------------------------------------------------------------*/

inline int
DisjointBoxLayout::linearNbrOffset(const IntVect& a_nbrOffset) const
{
  return D_TERM(  a_nbrOffset[0]*m_stride[0],
                + a_nbrOffset[1]*m_stride[1],
                + a_nbrOffset[2]*m_stride[2]);
}

/*--------------------------------------------------------------------*/
//  Begin linear index into local boxes
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::localIdxBegin() const
{
  return m_localIdxBeg;
}

/*--------------------------------------------------------------------*/
//  End linear index into local boxes (gives one past last local box)
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::localIdxEnd() const
{
  return m_localIdxBeg + m_numLocalBox;
}

/*--------------------------------------------------------------------*/
//  Number of processes
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::numProc()
{
  return s_numProc;
}

/*--------------------------------------------------------------------*/
//  ID of thie process
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::procID()
{
  return s_procID;
}

#endif  /* ! defined _DISJOINTBOXLAYOUT_H_ */
