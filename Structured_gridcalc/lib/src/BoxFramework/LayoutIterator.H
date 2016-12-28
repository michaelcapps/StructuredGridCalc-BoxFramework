
#ifndef _LAYOUTITERATOR_H_
#define _LAYOUTITERATOR_H_


/******************************************************************************/
/**
 * \file LayoutIterator.H
 *
 * \brief Iterator for indexing all boxes in a domain
 *
 *//*+*************************************************************************/

#include <cstring>
#include <ostream>
#include <iterator>

#include "Parameters.H"
#include "BoxIndex.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "BoxIterator.H"

//--Parameters for configuring the iterators

/// Trimming for neighbor iterators
enum
{
  TrimCenter = (1<<0),
  TrimFace   = (1<<1),
  TrimEdge   = (1<<2),
  TrimCorner = (1<<3)
};

/// Periodic directions for periodic iterators
enum
{
  PeriodicX = (1<<0),
  PeriodicY = (1<<1),
  PeriodicZ = (1<<2)
};


/*******************************************************************************
 */
///  Iterate over boxes in a DisjointBoxLayout
/**
 *//*+*************************************************************************/

class LayoutIterator
{

/*=============================================================================*
 * Types
 *============================================================================*/

public:

  using value_type = BoxIndex;

  using iterator_category = std::bidirectional_iterator_tag;
  using difference_type = int;

  using Self = LayoutIterator;

//--Friends

  friend class NeighborIterator;
  friend class PeriodicIterator;


/*=============================================================================*
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Construct from a DisjointBoxLayout
  LayoutIterator(const DisjointBoxLayout& a_dbl);

//--Use synthesized copy, copy assignment, move, move assignment, and destructor


/*=============================================================================*
 * Members functions
 *============================================================================*/

public:

  /// Reference to the BoxIndex
  value_type operator*() const;

  /// Prefix increment
  Self& operator++();

  /// Postfix increment
  Self operator++(int);

  /// Prefix decrement
  Self& operator--();

  /// Postfix decrement
  Self operator--(int);

  /// Increment by difference
  Self& operator+=(const difference_type& a_delta);

  /// Decrement by difference
  Self& operator-=(const difference_type& a_delta);

  /// Not Equal
  bool operator!=(const Self& a_iter) const;

  /// Still valid
  bool ok() const;

  /// Unique identifying tag for the DBL on which the iterator is constructed
  size_t tag() const;


/*=============================================================================*
 * Data members
 *============================================================================*/

protected:

  DisjointBoxLayout m_disjointBoxLayout;
                                      ///< Box layout in the domain
  int m_current;                      ///< A linear index describing the 
                                      ///< current box
  int m_size;                         ///< Number of boxes to iterate over
};


/*******************************************************************************
 */
///  Iterate only over boxes local to this processor
/**
 *   Requires only a simply modification to the construction of LayoutIterator
 *
 *//*+*************************************************************************/

class DataIterator : public LayoutIterator
{

public:

  /// Construct from a DisjointBoxLayout
  DataIterator(const DisjointBoxLayout& a_dbl)
    :
    LayoutIterator(a_dbl)
    {
      m_size = a_dbl.localIdxEnd();
      m_current = a_dbl.localIdxBegin();
    }

//--Use synthesized copy, copy assignment, move, move assignment, and destructor

  /// Reset to initial state
  void reset()
    {
      m_current = m_disjointBoxLayout.localIdxBegin();
    }
};


/*******************************************************************************
 */
///  Iterate over neighboring boxes to that pointed at by a LayoutIterator
/**
 *   This class is built with a LayoutIterator and iterates over the boxes
 *   neighboring that pointed at by the LayoutIterator.  A box describes
 *   the extent of the neighbors.  Only neighbors within the problem domain
 *   are considered
 *
 *//*+*************************************************************************/

class NeighborIterator : public LayoutIterator
{

/*=============================================================================*
 * Types
 *============================================================================*/

public:

  using value_type = BoxIndex;

  using iterator_category = std::forward_iterator_tag;
  using difference_type = IntVect;

  using Self = NeighborIterator;


/*=============================================================================*
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Construction
  NeighborIterator(
    const LayoutIterator& a_lit,
    const unsigned        a_trim = 0,
    Box                   a_nbr = Box(-IntVect::Unit, IntVect::Unit));

//--Use synthesized copy, copy assignment, move, move assignment, and destructor


/*=============================================================================*
 * Members functions
 *============================================================================*/

  /// Prefix increment
  Self& operator++();

  /// Postfix increment
  Self operator++(int);

  /// Still valid
  bool ok() const;

  /// Neighbor direction
  const IntVect& nbrDir() const;

//--Restrict some member functions from the base

  /// Prefix decrement
  Self& operator--() = delete;

  /// Postfix decrement
  Self operator--(int) = delete;

  /// Increment by difference
  Self& operator+=(const difference_type& a_delta) = delete;

  /// Decrement by difference
  Self& operator-=(const difference_type& a_delta) = delete;


/*=============================================================================*
 * Data members
 *============================================================================*/

protected:

  BoxIterator m_nbrOffset;            ///< An iterator over IntVects marking
                                      ///< neighbor boxes
  int m_base;                         ///< A linear index describing the base
                                      ///< box from the LayoutIterator
  int m_trim;                         ///< Codimensions to trim
};


/*******************************************************************************
 */
///  Iterate over periodic neighbor boxes to that pointed at by a LayoutIterator
/**
 *   This class is built with a LayoutIterator and iterates over the boxes
 *   neighboring that pointed at by the LayoutIterator.  The neighboring
 *   boxes are all on periodic boundaries.  You should test beforehand that
 *   the box pointed at by the LayoutIterator is indeed adjacent to a periodic
 *   boundary before creating this iterator.
 *
 *//*+*************************************************************************/

class PeriodicIterator : public LayoutIterator
{

/*=============================================================================*
 * Types
 *============================================================================*/

public:

  using value_type = BoxIndex;

  using iterator_category = std::forward_iterator_tag;
  using difference_type = IntVect;

  using Self = PeriodicIterator;


/*=============================================================================*
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Construction
  PeriodicIterator(
    const LayoutIterator& a_lit,
    const unsigned        a_trim = 0,
    const unsigned        a_periodic = 7);

//--Use synthesized copy, copy assignment, move, move assignment, and destructor


/*=============================================================================*
 * Members functions
 *============================================================================*/

  /// Prefix increment
  Self& operator++();

  /// Postfix increment
  Self operator++(int);

  /// Still valid
  bool ok() const;

  /// Neighbor direction
  const IntVect& nbrDir() const;

protected:

  /// Set current with considerations for periodic boundaries
  void setCurrent();

//--Restrict some member functions from the base

public:

  /// Prefix decrement
  Self& operator--() = delete;

  /// Postfix decrement
  Self operator--(int) = delete;

  /// Increment by difference
  Self& operator+=(const difference_type& a_delta) = delete;

  /// Decrement by difference
  Self& operator-=(const difference_type& a_delta) = delete;


/*=============================================================================*
 * Data members
 *============================================================================*/

protected:

  BoxIterator m_nbrOffset;            ///< An iterator over IntVects marking
                                      ///< neighbor boxes
  Box m_ivDomain;                     ///< A box describing the domain where
                                      ///< boxes are represented as IntVects
  Box m_ivPeriodicDomain;             ///< ivDomain, grown in periodic
                                      ///< directions
  Box m_ivPeriodicDomainSide[g_SpaceDim][2];
                                      ///< ivDomain, grown in a single periodic
                                      ///< direction on one side (0=low, 1=high)
  IntVect m_ivBase;                   ///< An IntVect describing the base box
                                      ///< in ivDomain
  int m_base;                         ///< A linear index describing the base
                                      ///< box from the LayoutIterator
  int m_trim;                         ///< Codimensions to trim
  int m_periodic;                     ///< Periodic directions
};


/*******************************************************************************
 *
 * Class LayoutIterator: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Construct from a DisjointBoxLayout
/** Starts at lower corner of box
 *//*-----------------------------------------------------------------*/

inline
LayoutIterator::LayoutIterator(const DisjointBoxLayout& a_dbl)
  :
  m_disjointBoxLayout(a_dbl),
  m_current(0),
  //m_size(a_dbl.size())
  m_size(a_dbl.size())
{
}

/*--------------------------------------------------------------------*/
//  Reference to the BoxIndex
/*--------------------------------------------------------------------*/

inline auto
LayoutIterator::operator*() const
  -> value_type
{
  return BoxIndex(m_current, m_current - m_disjointBoxLayout.localIdxBegin());
}

/*--------------------------------------------------------------------*/
//  Prefix increment
/*--------------------------------------------------------------------*/

inline auto
LayoutIterator::operator++()
  -> Self&
{
  ++m_current;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Postfix increment
/*--------------------------------------------------------------------*/

inline auto
LayoutIterator::operator++(int)
  -> Self
{
  Self tmp = *this;
  this->operator++();
  return tmp;
}

/*--------------------------------------------------------------------*/
//  Prefix decrement
/*--------------------------------------------------------------------*/

inline auto
LayoutIterator::operator--()
  -> Self&
{
  --m_current;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Postfix decrement
/*--------------------------------------------------------------------*/

inline auto
LayoutIterator::operator--(int)
  -> Self
{
  Self tmp = *this;
  this->operator--();
  return tmp;
}

/*--------------------------------------------------------------------*/
//  Increment by difference
/*--------------------------------------------------------------------*/

inline auto
LayoutIterator::operator+=(const difference_type& a_delta)
  -> Self&
{
  m_current += a_delta;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Decrement by difference
/*--------------------------------------------------------------------*/

inline auto
LayoutIterator::operator-=(const difference_type& a_delta)
  -> Self&
{
  m_current -= a_delta;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Not equal
/*--------------------------------------------------------------------*/

inline bool
LayoutIterator::operator!=(const Self& a_lit) const
{
  // You should only compare iterators built on the same DBL
  CH_assert(a_lit.tag() == tag());
  return (m_current != a_lit.m_current);
}

/*--------------------------------------------------------------------*/
//  Still valid
/*--------------------------------------------------------------------*/

inline bool
LayoutIterator::ok() const
{
  return (m_current < m_size);
}

/*--------------------------------------------------------------------*/
//  Unique identifying tag for the DBL on which the iterator is
//  constructed
/*--------------------------------------------------------------------*/

inline size_t
LayoutIterator::tag() const
{
  return m_disjointBoxLayout.tag();
}


/*******************************************************************************
 *
 * Class NeighborIterator: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Construction
/** Only neighbours inside the problem domain are considered.
 *  \param[in]  a_lit   Layout iterator to find neighbours too
 *  \param[in]  a_trim  Trimmed sections are not included as
 *                      neighbors.  (0,0,0) is always trimmed since
 *                      that is the base IV.  To trim edges and
 *                      corners, you would pass the value
 *                      TrimEdge | TrimCorner as an argument.  Default
 *                      is no trimming (aside from (0,0,0)).
 *  \param[in]  a_nbr   A box describing the neighbours to consider.
 *                      Checking neighbours more than 1 away (the
 *                      default) would be unusual.
 *//*-----------------------------------------------------------------*/

inline
NeighborIterator::NeighborIterator(const LayoutIterator &a_lit,
                                   const unsigned        a_trim,
                                   Box                   a_nbr)
  :
  LayoutIterator(a_lit),
  m_nbrOffset(),
  m_base(a_lit.m_current),
  m_trim(a_trim | TrimCenter)
{
  // Assume each box is a single IV.  Construct a box representing the domain
  Box ivDomain(IntVect::Zero, m_disjointBoxLayout.m_numBox - IntVect::Unit);

  // From the current index in the given LayoutIterator, find the corresponding
  // IV in our domain.
  const IntVect& stride = m_disjointBoxLayout.m_stride;
  int linearIdx = (*a_lit).globalIndex();
  IntVect ivBase;
  D_INVTERM(ivBase[0] = linearIdx;,
            ivBase[1] = linearIdx/stride[1];
            linearIdx -= stride[1]*ivBase[1];,
            ivBase[2] = linearIdx/stride[2];
            linearIdx -= stride[2]*ivBase[2];)

  // Shift the ivDomain so that 0,0,0 is instead centered on ivBase.  This is
  // required since the box a_nbr is also centered on (0,0,0)
  ivDomain.shift(-ivBase);
  // Intersect to crop the selection of neighbours by the domain
  a_nbr &= ivDomain;
  m_nbrOffset = BoxIterator(a_nbr);
  // Find first valid IntVect, avoiding trimmed values
  while (m_nbrOffset.ok() && ((1 << m_nbrOffset->norm1()) & m_trim))
    {
      ++m_nbrOffset;
    }
  m_current = m_base + m_disjointBoxLayout.linearNbrOffset(*m_nbrOffset);
}

/*--------------------------------------------------------------------*/
//  Prefix increment
/*--------------------------------------------------------------------*/

inline auto
NeighborIterator::operator++()
  -> Self&
{
  ++m_nbrOffset;
  // Find first valid IntVect, avoiding trimmed values
  while (m_nbrOffset.ok() && ((1 << m_nbrOffset->norm1()) & m_trim))
    {
      ++m_nbrOffset;
    }
  m_current = m_base + m_disjointBoxLayout.linearNbrOffset(*m_nbrOffset);
  return *this;
}

/*--------------------------------------------------------------------*/
//  Postfix increment
/*--------------------------------------------------------------------*/

inline auto
NeighborIterator::operator++(int)
  -> Self
{
  Self tmp = *this;
  this->operator++();
  return tmp;
}

/*--------------------------------------------------------------------*/
//  Still valid
/*--------------------------------------------------------------------*/

inline bool
NeighborIterator::ok() const
{
  return m_nbrOffset.ok();
}

/*--------------------------------------------------------------------*/
//  Neighbor direction
/** \return             An Intvect with components 1, -1, or 0, giving
 *                      the direction to the neighbor
 *//*-----------------------------------------------------------------*/

inline const IntVect&
NeighborIterator::nbrDir() const
{
  return *m_nbrOffset;
}


/*******************************************************************************
 *
 * Class PeriodicIterator: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Construction
/** Only neighbours inside the problem domain are considered.
 *  \param[in]  a_lit   Layout iterator to find neighbours too
 *  \param[in]  a_trim  Trimmed sections are not included as
 *                      neighbors.  (0,0,0) is always trimmed since
 *                      that is the base IV.  To trim edges and
 *                      corners, you would pass the value
 *                      TrimEdge | TrimCorner as an argument.  Default
 *                      is no trimming (aside from (0,0,0)).
 *  \param[in]  a_periodic
 *                      Which directions are periodic.  To select X
 *                      and Y, you would pass the value
 *                      PeriodicX | PeriodicY as an argument.  Default
 *                      is X, Y, and Z are periodic.
 *//*-----------------------------------------------------------------*/

inline
PeriodicIterator::PeriodicIterator(const LayoutIterator &a_lit,
                                   const unsigned        a_trim,
                                   const unsigned        a_periodic)
  :
  LayoutIterator(a_lit),
  m_nbrOffset(),
  m_base(a_lit.m_current),
  m_trim(a_trim | TrimCenter),
  m_periodic(a_periodic)
{
  Box nbr(-IntVect::Unit, IntVect::Unit);

  // Assume each box is a single IV.  Construct a box representing the domain
  m_ivDomain.define(IntVect::Zero,
                    m_disjointBoxLayout.m_numBox - IntVect::Unit);

  // From the current index in the given LayoutIterator, find the corresponding
  // IV in our domain.
  const IntVect& stride = m_disjointBoxLayout.m_stride;
  int linearIdx = (*a_lit).globalIndex();
  IntVect m_ivBase;
  D_INVTERM(m_ivBase[0] = linearIdx;,
            m_ivBase[1] = linearIdx/stride[1];
            linearIdx -= stride[1]*m_ivBase[1];,
            m_ivBase[2] = linearIdx/stride[2];
            linearIdx -= stride[2]*m_ivBase[2];)

  // Shift the m_ivDomain so that 0,0,0 is instead centered on m_ivBase.  This
  // is required since the box nbr is also centered on (0,0,0)
  m_ivDomain.shift(-m_ivBase);
  // Periodic domains are grown by 1 in periodic directions
  m_ivPeriodicDomain = m_ivDomain;
  for (int dir = 0; dir != g_SpaceDim; ++dir)
    {
      if (a_periodic & (1<<dir))
        {
          m_ivPeriodicDomain.grow(1, dir);
        }
    }
  for (int dir = 0; dir != g_SpaceDim; ++dir)
    {
      if (a_periodic & (1<<dir))
        {
          m_ivPeriodicDomainSide[dir][0] = m_ivPeriodicDomain;
          m_ivPeriodicDomainSide[dir][0].grow(-1, dir);
          m_ivPeriodicDomainSide[dir][1] = m_ivPeriodicDomainSide[dir][0];
          m_ivPeriodicDomainSide[dir][0].adjBox(1, dir, -1);
          m_ivPeriodicDomainSide[dir][1].adjBox(1, dir, 1);
        }
    }
  // Intersect to crop the selection of neighbours by the domain
  nbr &= m_ivPeriodicDomain;
  // Nothing to do if not near a periodic boundary (set to empty box)
  if (m_ivDomain.contains(nbr))
    {
      nbr = Box();
    }
  m_nbrOffset = BoxIterator(nbr);
  setCurrent();
}

/*--------------------------------------------------------------------*/
//  Prefix increment
/*--------------------------------------------------------------------*/

inline auto
PeriodicIterator::operator++()
  -> Self&
{
  ++m_nbrOffset;
  setCurrent();
  return *this;
}

/*--------------------------------------------------------------------*/
//  Postfix increment
/*--------------------------------------------------------------------*/

inline auto
PeriodicIterator::operator++(int)
  -> Self
{
  Self tmp = *this;
  this->operator++();
  return tmp;
}

/*--------------------------------------------------------------------*/
//  Still valid
/*--------------------------------------------------------------------*/

inline bool
PeriodicIterator::ok() const
{
  return m_nbrOffset.ok();
}

/*--------------------------------------------------------------------*/
//  Neighbor direction
/** \return             An Intvect with components 1, -1, or 0, giving
 *                      the direction to the periodic image of the
 *                      neighbour (i.e., this will point outside the
 *                      domain).
 *//*-----------------------------------------------------------------*/

inline const IntVect&
PeriodicIterator::nbrDir() const
{
  return *m_nbrOffset;
}

/*--------------------------------------------------------------------*/
//  Set current with considerations for periodic boundaries
/*--------------------------------------------------------------------*/

inline void
PeriodicIterator::setCurrent()
{
  // Find first valid IntVect, avoiding trimmed values and those inside domain
  while (m_nbrOffset.ok() && (
           ((1 << m_nbrOffset->norm1()) & m_trim) ||
           m_ivDomain.contains(*m_nbrOffset)))  // Has to be outside domain
    {
      ++m_nbrOffset;
    }
  // Normal offset
  IntVect offset = *m_nbrOffset;
  // Add offsets based on periodicity
  const IntVect& domainDimensions = m_disjointBoxLayout.dimensions();
  const IntVect nbr = m_ivBase + offset;
  for (int dir = 0; dir != g_SpaceDim; ++dir)
    {
      if (m_periodic & (1<<dir))
        {
          if (m_ivPeriodicDomainSide[dir][0].contains(nbr))
            {
              offset[dir] += domainDimensions[dir];
            }
          if (m_ivPeriodicDomainSide[dir][1].contains(nbr))
            {
              offset[dir] -= domainDimensions[dir];
            }
        }
    }
  m_current = m_base + m_disjointBoxLayout.linearNbrOffset(offset);
}


/*******************************************************************************
 *
 * Class DisjointBoxLayout: inline member definitions
 *
 * These are included here because of a dependency on LayoutIterator
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Index with a LayoutIterator
/*--------------------------------------------------------------------*/

inline const Box&
DisjointBoxLayout::operator[](const LayoutIterator& a_layit) const
{
  // Make sure we are indexing with an iterator built on this DBL
  CH_assert(a_layit.tag() == tag());

  return (*m_boxes)[(*a_layit).globalIndex()].box;  
}

/*--------------------------------------------------------------------*/
//  Find local process ID of a box with a LayoutIterator
/*--------------------------------------------------------------------*/

inline int
DisjointBoxLayout::proc(const LayoutIterator& a_layit) const
{
  // Make sure we are indexing with an iterator built on this DBL
  CH_assert(a_layit.tag() == tag());
  return (*m_boxes)[(*a_layit).globalIndex()].proc;
  //return 0;
}

/*--------------------------------------------------------------------*/
//  Get box and local process ID with a LayoutIterator
/*--------------------------------------------------------------------*/

inline const Box&
DisjointBoxLayout::box(const LayoutIterator& a_layit, int& a_proc) const
{
  // Make sure we are indexing with an iterator built on this DBL
  CH_assert(a_layit.tag() == tag());

  a_proc = (*m_boxes)[(*a_layit).globalIndex()].proc;
  return (*m_boxes)[(*a_layit).globalIndex()].box;
  //return Box();
}

#endif
