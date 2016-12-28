
#ifndef _BOX_H_
#define _BOX_H_


/******************************************************************************/
/**
 * \mainpage Box framework for parallel calculations on structured grids
 *
 * \file Box.H
 *
 * \brief Class describing a multidimensional region in space
 *
 *//*+*************************************************************************/

#include <iosfwd>

#include "Parameters.H"
#include "IntVect.H"

/*
  Notes:
    - tag __CUDACC__ means the file is process with nvcc.  The class is
      modified whether on the cpu or gpu.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef HOSTDEVICE
#undef DEVICE
#ifdef __CUDACC__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define HOSTDEVICE
#define DEVICE
#endif

//--Forward declarations

template<typename T, typename Ref, typename Ptr>
class BoxIteratorImpl;

//--BoxIterator Types

using BoxIterator = BoxIteratorImpl<IntVect, const IntVect&, const IntVect*>;


/*******************************************************************************
 */
///  Box in space
/**
 ******************************************************************************/

class Box
{

/*==============================================================================
 * Types
 *============================================================================*/

public:

  using const_iterator =  BoxIterator;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  HOSTDEVICE Box();

  /// Construct with values
  HOSTDEVICE Box(const IntVect& a_lo, const IntVect& a_hi);

  /// Copy constructor
  HOSTDEVICE Box(const Box&) = default;

  /// Move constructor
  HOSTDEVICE Box(Box&&) = default;

  /// Assignment constructor
  HOSTDEVICE Box& operator=(const Box&) = default;

  /// Move assignment constructor
  HOSTDEVICE Box& operator=(Box&&) = default;

  /// Destructor
  HOSTDEVICE ~Box() = default;

  /// Weak construction
  HOSTDEVICE void define(const IntVect& a_lo, const IntVect& a_hi);


/*==============================================================================
 * Members functions
 *============================================================================*/

public:

  /// Return the lower IntVect
  HOSTDEVICE const IntVect& loVect() const;

  /// Return an integer direction in the lower IntVect
  HOSTDEVICE int loVect(const int a_dir) const;

  /// Return a modifiable integer direction in the lower IntVect
  HOSTDEVICE int& loVect(const int a_dir);

  /// Return the upper IntVect
  HOSTDEVICE const IntVect& hiVect() const;

  /// Return an integer direction in the upper IntVect
  HOSTDEVICE int hiVect(const int a_dir) const;

  /// Return a modifiable integer direction in the upper IntVect
  HOSTDEVICE int& hiVect(const int a_dir);

  /// Grow the box in all directions
  HOSTDEVICE Box& grow(const int a_i);

  /// Grow the box in a specific direction
  HOSTDEVICE Box& grow(const int a_i, const int a_dir);

  /// Grow the upper corner in all directions (switch to vertex)
  HOSTDEVICE Box& growHi(const int a_i);

  /// Grow the lower side in a specific direction
  HOSTDEVICE Box& growLo(const int a_i, const int a_dir);

  /// Grow the upper side in a specific direction
  HOSTDEVICE Box& growHi(const int a_i, const int a_dir);

  /// Shift by an IntVect
  HOSTDEVICE Box& shift(const IntVect& a_shift);

  /// Shift by a scalar in a single direction
  HOSTDEVICE Box& shift(const int a_i, const int a_dir);

  /// Adjacent cells on one side of the box
  HOSTDEVICE Box& adjBox(int a_ncell, const int a_dir, const int a_side);

  /// Intersect with another box
  HOSTDEVICE Box& operator&=(const Box& a_box);

  /// Contains an IntVect
  HOSTDEVICE bool contains(const IntVect& a_iv) const;

  /// Contains another box
  HOSTDEVICE bool contains(const Box& a_box) const;

  /// Return the size
  HOSTDEVICE int size() const;

  /// Return the dimensions
  HOSTDEVICE IntVect dimensions() const;

  /// Is the box empty?
  HOSTDEVICE bool isEmpty() const;

  /// Starting iterator
  const_iterator begin() const;

  /// Ending iterator
  const_iterator end() const;

  /// Equal
  HOSTDEVICE bool operator==(const Box& a_box) const;

  /// Not Equal
  HOSTDEVICE bool operator!=(const Box& a_box) const;

#ifdef __CUDACC__
//--These are for easy copying on the gpu
  /// Constant access to data
  DEVICE const int& operator[](const int a_idx) const;

  /// Access to an index
  DEVICE int& operator[](const int a_idx);

  /// Get the strides for an array sized to this box
  DEVICE void getStride(IntVect& a_stride) const;

  /// Get the strides for an array sized to this box
  DEVICE void getStride(int* const a_stride) const;

  /// Get the offset
  HOSTDEVICE int getOffset(const IntVect& a_stride) const;

  /// Get the offset
  HOSTDEVICE int getOffset(const int* const a_stride) const;

  /// Get the offset (strides are computed on the fly)
  HOSTDEVICE int getOffset() const;

  /// Convert a vector index (in the box) to a zero-based linear index
  HOSTDEVICE int vecToLin0(const IntVect& a_vec, const IntVect& a_stride) const;

  /// Convert a vector index (in the box) to a zero-based linear index
  HOSTDEVICE int vecToLin0(const IntVect&   a_vec,
                           const int* const a_stride) const;

  /// Convert a box-based linear index (in the box) to a vector index
  DEVICE void linToVec(int a_lin, const IntVect& a_stride,IntVect& a_vec) const;

  /// Convert a box-based linear index (in the box) to a vector index
  DEVICE void linToVec(int              a_lin,
                       const int* const a_stride,
                       IntVect&         a_vec) const;

  /// Convert a box-based linear index (in the box) to a vector index
  DEVICE void linToVec(int a_lin, IntVect& a_vec) const;
#endif


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  IntVect m_lo;                       ///< Lower corner of box
  IntVect m_hi;                       ///< Upper corner of box
};


/*******************************************************************************
 *
 * Class Box: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/** Creates an empty box
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline
Box::Box()
#ifndef __CUDACC__
  :
  m_lo(IntVect::Unit),
  m_hi(IntVect::Zero)
#endif
{
}

/*--------------------------------------------------------------------*/
//  Constructor with values
/** \param[in] a_i#     Initial values
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline
Box::Box(const IntVect& a_lo, const IntVect& a_hi)
  :
  m_lo(a_lo),
  m_hi(a_hi)
{
  CH_assert(m_lo <= m_hi);
}

/*--------------------------------------------------------------------*/
//  Weak construction
/** \param[in] a_i#     Initial values
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline void
Box::define(const IntVect& a_lo, const IntVect& a_hi)
{
  CH_assert(a_lo <= a_hi);
  m_lo = a_lo;
  m_hi = a_hi;
}

/*--------------------------------------------------------------------*/
//  Return the lower IntVect
/*--------------------------------------------------------------------*/

HOSTDEVICE inline const IntVect&
Box::loVect() const
{
  return m_lo;
}

/*--------------------------------------------------------------------*/
//  Return an integer direction in the lower IntVect
/*--------------------------------------------------------------------*/

HOSTDEVICE inline int
Box::loVect(const int a_dir) const
{
  return m_lo[a_dir];
}

/*--------------------------------------------------------------------*/
//  Return a modifiable integer direction in the lower IntVect
/*--------------------------------------------------------------------*/

HOSTDEVICE inline int&
Box::loVect(const int a_dir)
{
  return m_lo[a_dir];
}

/*--------------------------------------------------------------------*/
//  Return the upper IntVect
/*--------------------------------------------------------------------*/

HOSTDEVICE inline const IntVect&
Box::hiVect() const
{
  return m_hi;
}

/*
inline int
Box::loVect(const int a_dir) const
{
  return m_lo[a_dir];
}

inline int&
Box::loVect(const int a_dir)
{
  return m_lo[a_dir];
}

*/
HOSTDEVICE inline int
Box::hiVect(const int a_dir) const
{
  return m_hi[a_dir];
}

/*--------------------------------------------------------------------*/
//  Return a modifiable integer direction in the upper IntVect
/*--------------------------------------------------------------------*/

HOSTDEVICE inline int&
Box::hiVect(const int a_dir)
{
  return m_hi[a_dir];
}

/*--------------------------------------------------------------------*/
//  Grow the box in all directions
/** \param[in]  a_i     Amount to grow by
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::grow(const int a_i)
{
  m_lo = m_lo - a_i*IntVect(1,1,1);
  m_hi = m_hi + a_i*IntVect(1,1,1);
  return *this;
}

/*--------------------------------------------------------------------*/
//  Grow the box in a specific direction
/** \param[in]  a_i     Amount to grow by
 *  \param[in]  a_dir   Direction to grow in
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::grow(const int a_i, const int a_dir)
{
  CH_assert(a_dir >= 0 && a_dir <= g_SpaceDim);
  m_lo[a_dir] -= a_i;
  m_hi[a_dir] += a_i;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Grow the upper corner in all directions (switch to vertex)
/** \param[in]  a_i     Amount to grow by
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::growHi(const int a_i)
{
  m_hi += a_i;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Grow the lower side in a specific direction
/** \param[in]  a_i     Amount to grow by
 *  \param[in]  a_dir   Direction to grow in
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::growLo(const int a_i, const int a_dir)
{
  m_lo[a_dir] -= a_i;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Grow the upper side in a specific direction
/** \param[in]  a_i     Amount to grow by
 *  \param[in]  a_dir   Direction to grow in
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::growHi(const int a_i, const int a_dir)
{
  m_hi[a_dir] += a_i;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Grow the box in a specific direction
/** \param[in]  a_i     Amount to grow by
 *  \param[in]  a_dir   Direction to grow in
 *//*-----------------------------------------------------------------*/
/*
inline Box&
Box::grow(const int a_i, const int a_dir)
{
  m_hi[a_dir] += a_i;
  m_lo[a_dir] += -a_i;
  return *this;
}
*/

HOSTDEVICE inline Box&
Box::shift(const IntVect& a_shift)
{
  m_lo += a_shift;
  m_hi += a_shift;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Shift by a scalar in a single direction
/** \param[in]  a_i     Number of cells to shift
 *  \param[in]  a_dir   Direction to shift
 *  \return             Shifted box
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::shift(const int a_i, const int a_dir)
{
  m_lo[a_dir] += a_i;
  m_hi[a_dir] += a_i;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Adjacent cells on one side of the box
/** Return a box adjacent to this one on a given side
 *  \param[in]  a_i     Number of cells adjacent to this box.  If
 *                      negative, adjacent cells are within the
 *                      original box
 *  \param[in]  a_dir   Direction of the side
 *  \param[in]  a_side  Sign of the side
 *                      <= 0 - min
 *                      >  0 - max
 *  \return             Adjacent box
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::adjBox(int a_ncell, const int a_dir, const int a_side)
{
  if(a_ncell > 0 && a_side > 0)
  {
    m_hi[a_dir] += a_ncell;
    m_lo[a_dir] = m_hi[a_dir] - a_ncell + 1;
  }
  else if(a_ncell  >0 && a_side <= 0)
  {
    m_lo[a_dir] -= a_ncell;
    m_hi[a_dir] = m_lo[a_dir] + a_ncell -1;
  }
  else if(a_ncell < 0 && a_side > 0)
  {
    m_lo[a_dir] = m_hi[a_dir] + a_ncell + 1;
  }
  else if(a_ncell < 0 && a_side <= 0)
  {
    m_hi[a_dir] = m_lo[a_dir]-a_ncell -1;
  }
  return *this;
}

/*--------------------------------------------------------------------*/
//  Intersect with another box
/** \param[in]  a_box   Box to intersect with
 *  \return             Intersection of the two boxes
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline Box&
Box::operator&=(const Box& a_box)
{
  	m_lo.max(a_box.m_lo);
	m_hi.min(a_box.m_hi);
	return *this;
}

/*--------------------------------------------------------------------*/
//  Contains an IntVect
/** \param[in]  a_iv    IntVect to check
 *  \return             T - contains IntVect
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline bool
Box::contains(const IntVect& a_iv) const
{
  return (m_lo <= a_iv && a_iv <= m_hi);
}

/*--------------------------------------------------------------------*/
//  Contains another box
/** \param[in]  a_box   Box to check
 *  \return             T - contains Box
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline bool
Box::contains(const Box& a_box) const
{
  return (m_lo <= a_box.m_lo && a_box.m_hi <= m_hi);
}

/*--------------------------------------------------------------------*/
//  Return the size
/*--------------------------------------------------------------------*/

HOSTDEVICE inline int
Box::size() const
{
  return D_TERM( (m_hi[0] - m_lo[0] + 1),
                *(m_hi[1] - m_lo[1] + 1),
                *(m_hi[2] - m_lo[2] + 1));
}

/*
inline bool 
Box::isEmpty() const
{
  return D_TERM((m_hi[0] < m_lo[0]),||(m_hi[1] < m_lo[1]),||(m_hi[2]<m_lo[2]));
}
*/
/*--------------------------------------------------------------------*/
//  Return the dimensions
/*--------------------------------------------------------------------*/

HOSTDEVICE inline IntVect
Box::dimensions() const
{
  return IntVect(D_DECL((m_hi[0] - m_lo[0] + 1),
                        (m_hi[1] - m_lo[1] + 1),
                        (m_hi[2] - m_lo[2] + 1)));
}

/*--------------------------------------------------------------------*/
//  Is the box empty?
/*--------------------------------------------------------------------*/

HOSTDEVICE inline bool
Box::isEmpty() const
{
  return D_TERM(   (m_hi[0] < m_lo[0]),
                || (m_hi[1] < m_lo[1]),
                || (m_hi[2] < m_lo[2]));
}

/*--------------------------------------------------------------------*/
//  Are the boxes equivalent?
/*--------------------------------------------------------------------*/

HOSTDEVICE inline bool
Box::operator==(const Box& a_box) const
{
  return (m_lo == a_box.m_lo && m_hi == a_box.m_hi);
}

/*--------------------------------------------------------------------*/
//  Are the boxes not equivalent?
/*--------------------------------------------------------------------*/

HOSTDEVICE inline bool
Box::operator!=(const Box& a_box) const
{
  return (m_lo != a_box.m_lo || m_hi != a_box.m_hi);
}

#ifdef __CUDACC__
/*--------------------------------------------------------------------*/
//  Constant access to an index
/** \param[in]  a_idx   Linear index across lo and hi
 *  \return             Value at requested index
 *//*-----------------------------------------------------------------*/

DEVICE inline const int&
Box::operator[](const int a_idx) const
{
  CH_assert(a_idx >= 0 && a_idx < 2*g_SpaceDim);
  return m_lo.dataPtr()[a_idx];
}

/*--------------------------------------------------------------------*/
//  Access to an index
/** \param[in]  a_idx   Linear index across lo and hi
 *  \return             Value at requested index
 *//*-----------------------------------------------------------------*/

DEVICE inline int&
Box::operator[](const int a_idx)
{
  CH_assert(a_idx >= 0 && a_idx < 2*g_SpaceDim);
  return m_lo.dataPtr()[a_idx];
}

/*--------------------------------------------------------------------*/
//  Get the strides for an array sized to this box
/** \param[out] a_stride
 *                      Strides in an IntVect
 *//*-----------------------------------------------------------------*/

DEVICE inline void
Box::getStride(IntVect& a_stride) const
{
  getStride(a_stride.dataPtr());
}

/*--------------------------------------------------------------------*/
//  Get the strides for an array sized to this box
/** \param[out] a_stride
 *                      Strides in an array
 *//*-----------------------------------------------------------------*/

DEVICE inline void
Box::getStride(int* const a_stride) const
{
  D_TERM(a_stride[0] = 1;,
         a_stride[1] = a_stride[0]*(hiVect()[0] - loVect()[0] + 1);,
         a_stride[2] = a_stride[1]*(hiVect()[1] - loVect()[1] + 1);)
}

/*--------------------------------------------------------------------*/
//  Get the offset
/** Add the offset to the pointer to get how far the pointer needs to
 *  be offset to account for a non-zero smallEnd.  If loVect is zero,
 *  this is a zero-based array and no offset is necessary.
 *  \param[in]  a_stride
 *                      Strides in an array
 *  \return             Offset
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline int
Box::getOffset(const IntVect& a_stride) const
{
  return getOffset(a_stride.dataPtr());
}

/*--------------------------------------------------------------------*/
//  Get the offset
/** Add the offset to the pointer to get how far the pointer needs to
 *  be offset to account for a non-zero smallEnd.  If loVect is zero,
 *  this is a zero-based array and no offset is necessary.
 *  \param[in]  a_stride
 *                      Strides in an IntVect
 *  \return             Offset
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline int
Box::getOffset(const int* const a_stride) const
{
  return vecToLin0(-loVect(), a_stride);
}

/*--------------------------------------------------------------------*/
//  Get the offset (strides are computed on the fly)
/** Add the offset to the pointer to get how far the pointer needs to
 *  be offset to account for a non-zero smallEnd.  If loVect is zero,
 *  this is a zero-based array and no offset is necessary.
 *  \return             Offset
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline int
Box::getOffset() const
{
  IntVect vec = -loVect();
  D_TERM(int stride = 1;
         int offset = vec[0];,
         stride *= (hiVect()[0] - loVect()[0] + 1);
         offset += vec[1]*stride;,
         stride *= (hiVect()[1] - loVect()[1] + 1);
         offset += vec[2]*stride;)
    return offset;
}

/*--------------------------------------------------------------------*/
//  Convert a vector index (in the box) to a zero-based linear index
/** No consideration is made for a non-zero smallEnd.  This must
 *  instead be accounted for by the offset.  An offset data pointer
 *  (fab.dataPtr() + offset) plus lin0 yields a valid index into the
 *  data.  Alternatively a direct data pointer indexed with
 *  [lin + offset] yields a valid index into the data
 *  (fab.dataPtr() + lin0 + offset).
 *  \param[in]  a_vec   Vector index (contained in the box) represented as
 *                      an array
 *  \param[in]  a_stride
 *                      Strides in an IntVect
 *  \return             Linear index (used with offset data pointer)
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline int
Box::vecToLin0(const IntVect& a_vec, const IntVect& a_stride) const
{
  return vecToLin0(a_vec, a_stride.dataPtr());
}

/*--------------------------------------------------------------------*/
//  Convert a vector index (in the box) to a zero-based linear index
/** No consideration is made for a non-zero smallEnd.  This must
 *  instead be accounted for by the offset.  An offset data pointer
 *  (fab.dataPtr() + offset) plus lin0 yields a valid index into the
 *  data.  Alternatively a direct data pointer indexed with
 *  [lin + offset] yields a valid index into the data
 *  (fab.dataPtr() + lin0 + offset).
 *  \param[in]  a_vec   Vector index (contained in the box) represented as
 *                      an array
 *  \param[in]  a_stride
 *                      Strides in an array
 *  \return             Linear index (used with offset data pointer)
 *//*-----------------------------------------------------------------*/

HOSTDEVICE inline int
Box::vecToLin0(const IntVect& a_vec, const int* const a_stride) const
{
  return D_TERM(  a_vec[0],
                + a_vec[1]*a_stride[1],
                + a_vec[2]*a_stride[2]);
}

/*--------------------------------------------------------------------*/
//  Convert a box-based linear index (in the box) to a vector index
/** If you have a zero-based linear index (see vecToLin0), you must
 *  add the offset before calling this routine.  E.g.
 *  linToVec(lin0 + offset, strides, vec).
 *  The linear index must resolve to a location in the box.
 *  \param[in]  a_lin   A box-based linear index (add offset if you have
 *                      a zero-based index)
 *  \param[in]  a_stride
 *                      Strides in an IntVect
 *  \param[out] a_vec   Vector index (contained in the box) represented as
 *                      an IntVect
 *//*-----------------------------------------------------------------*/

DEVICE inline void
Box::linToVec(int a_lin, const IntVect& a_stride, IntVect& a_vec) const
{
  linToVec(a_lin, a_stride.dataPtr(), a_vec);
}

/*--------------------------------------------------------------------*/
//  Convert a box-based linear index (in the box) to a vector index
/** If you have a zero-based linear index (see vecToLin0), you must
 *  add the offset before calling this routine.  E.g.
 *  linToVec(lin0 + offset, strides, vec).
 *  The linear index must resolve to a location in the box.
 *  \param[in]  a_lin   A box-based linear index (add offset if you have
 *                      a zero-based index)
 *  \param[in]  a_stride
 *                      Strides in an array
 *  \param[out] a_vec   Vector index (contained in the box) represented as
 *                      an IntVect
 *//*-----------------------------------------------------------------*/

DEVICE inline void
Box::linToVec(int a_lin, const int* const a_stride, IntVect& a_vec) const
{
  int tmp;
  D_INVTERM(a_vec[0] = a_lin + loVect()[0];,

            tmp = a_lin/a_stride[1];
            a_vec[1] = tmp + loVect()[1];
            a_lin -= tmp*a_stride[1];,

            tmp = a_lin/a_stride[2];
            a_vec[2] = tmp + loVect()[2];
            a_lin -= tmp*a_stride[2];)
}

/*--------------------------------------------------------------------*/
//  Convert a box-based linear index (in the box) to a vector index
/** If you have a zero-based linear index (see vecToLin0), you must
 *  add the offset before calling this routine.  E.g.
 *  linToVec(lin0 + offset, strides, vec).
 *  The linear index must resolve to a location in the box.  Strides
 *  are computed in this version.
 *  \param[in]  a_lin   A box-based linear index (add offset if you have
 *                      a zero-based index)
 *  \param[out] a_vec   Vector index (contained in the box) represented as
 *                      an IntVect
 *//*-----------------------------------------------------------------*/

DEVICE inline void
Box::linToVec(int a_lin, IntVect& a_vec) const
{
  IntVect stride;
  getStride(stride);
  linToVec(a_lin, stride.dataPtr(), a_vec);
}
#endif

/*******************************************************************************
 *
 * Class Box: external related functions
 *
 ******************************************************************************/

/// Output the Box
std::ostream &operator<<(std::ostream &a_os, const Box &a_box);

#endif  /* ! defined _BOX_H_ */
