
#ifndef _BASEFAB_H_
#define _BASEFAB_H_


/******************************************************************************/
/**
 * \file BaseFab.H
 *
 * \brief Class providing data for a box
 *
 *//*+*************************************************************************/

#include <cstring>
#include <ostream>
#include <limits>
#include <memory>

#include "Parameters.H"
#include "Box.H"

#ifdef USE_GPU
#include "CudaSupport.H"
#endif


/*******************************************************************************
 */
///  Data for a box (fortran array box)
/**
 ******************************************************************************/

template <typename T>
class BaseFab
{

/*==============================================================================
 * Types
 *============================================================================*/

public:

  using value_type = T;

  enum class AllocBy
  {
    none,                             ///< Undefined
    array,                            ///< Data allocated by new[]
    alias                             ///< Data aliased
  };


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor (no allocation)
  BaseFab();

  /// Constructor with sizes
  BaseFab(const Box& a_box, const int a_ncomp, T *const a_alias = nullptr);

  // Constructor with sizes and a default value
  BaseFab(const Box& a_box,
          const int  a_ncomp,
          const T&   a_val,
          T *const   a_alias = nullptr);

  /// Copy constructor
  BaseFab(const BaseFab&) = delete;

  /// Move constructor
  BaseFab(BaseFab&&) noexcept;

  /// Assignment constructor
  BaseFab& operator=(const BaseFab&) = delete;

  /// Move assignment constructor
  BaseFab& operator=(BaseFab&&) noexcept;

  /// Weak construction with sizes
  void define(const Box& a_box, const int a_ncomp, T *const a_alias = nullptr);

  /// Weak construction with sizes and a default value
  void define(const Box& a_box,
              const int  a_ncomp,
              const T&   a_val,
              T *const   a_alias = nullptr);

  /// Destructor
  ~BaseFab();


/*==============================================================================
 * Members functions
 *============================================================================*/

public:

  /// Return the box
  const Box& box() const;

  /// Return the number of components
  int ncomp() const;

  /// Return the total number of elements
  int size() const;

  /// Return the total number of bytes used
  size_t sizeBytes() const;

  /// Constant access to an element
  const T& operator()(const IntVect& a_iv, const int a_icomp) const;

  /// Access to an element
  T& operator()(const IntVect& a_iv, const int a_icomp);

  /// Assign a constant to all components
  void setVal(const T& a_val);

  /// Assign a constant to a single component
  void setVal(const int a_icomp, const T& a_val);

  /// Copy a portion of another BaseFab, same region and all components
  void copy(const Box&     a_box,
            const BaseFab& a_src);

  /// Copy a portion of another BaseFab
  void copy(const Box&     a_dstBox,
            const int      a_dstComp, 
            const BaseFab& a_src,
            const Box&     a_srcBox,
            const int      a_srcComp,
            const int      a_numComp,
            const unsigned a_compFlags = std::numeric_limits<unsigned>::max());

  /// Linearize data in a region and place in a buffer
  void linearOut(
    void *const    a_buffer,
    const Box&     a_region,
    const int      a_startComp,
    const int      a_endComp,
    const unsigned a_compFlags = std::numeric_limits<unsigned>::max()) const;

  /// Replace data in a region from a linear buffer
  void linearIn(
    const void* const a_buffer,
    const Box&        a_region,
    const int         a_startComp,
    const int         a_endComp,
    const unsigned    a_compFlags = std::numeric_limits<unsigned>::max());

  /// Obtain a linear index (internal and testing use only)
  int index(IntVect a_iv) const;

  /// Start of data for a component (internal use only)
  const T* dataPtr(const int a_icomp = 0) const;

  /// Start of data for a component (internal use only)
  T* dataPtr(const int a_icomp = 0);

  /// Get spatial strides (internal use only)
  const IntVect& getStride() const;

  /// Get component stride (internal use only)
  int getComponentStride() const;

#ifdef USE_GPU
  /// Copy array to device
  void copyToDevice() const;

  /// Asynchronous copy array to device
  void copyToDeviceAsync(cudaStream_t a_stream = 0) const;

  /// Copy array to host
  void copyToHost();

  /// Asynchronous copy array to host
  void copyToHostAsync(cudaStream_t a_stream = 0);
#endif


/*==============================================================================
 * Private members functions
 *============================================================================*/

private:

  /// Set strides
  void setStride();

  /// Allocate memory
  void allocate();

  /// Deallocate memory
  void deallocate();


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Box m_box;                          ///< Box defining data
  IntVect m_stride;                   ///< Stride for indexing
  int m_ncomp;                        ///< Number of components
  int m_size;                         ///< Size of the box
  T* m_data;                          ///< Data
  AllocBy m_allocBy;                  ///< Method of allocation
#ifdef USE_GPU
public:
  SymbolPair<T> m_dataSymbol;         ///< Pointers to data on host and device
#endif
};


/*******************************************************************************
 *
 * Class BaseFab: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Return the box
/*--------------------------------------------------------------------*/

template <typename T>
inline const Box&
BaseFab<T>::box() const
{
  return m_box;
}

/*--------------------------------------------------------------------*/
//  Return the number of components
/*--------------------------------------------------------------------*/

template <typename T>
inline int
BaseFab<T>::ncomp() const
{
  return m_ncomp;
}

/*--------------------------------------------------------------------*/
//  Return the total number of elements
/*--------------------------------------------------------------------*/

template <typename T>
inline int
BaseFab<T>::size() const
{
  return m_ncomp*m_size;
}

/*--------------------------------------------------------------------*/
//  Return the total number of bytes used
/*--------------------------------------------------------------------*/

template <typename T>
inline size_t
BaseFab<T>::sizeBytes() const
{
  return (*this).size()*sizeof(typeid(T).name());
}

/*--------------------------------------------------------------------*/
//  Constant access to an element
/** \param[in]  a_iv    IntVect location
 *  \param[in]  a_icomp Component index
 *  \return             Constant element
 *//*-----------------------------------------------------------------*/

template <typename T>
inline const T&
BaseFab<T>::operator()(const IntVect& a_iv, const int a_icomp) const
{
  return m_data[a_icomp*m_size + index(a_iv)];
}

/*--------------------------------------------------------------------*/
//  Access to an element
/** \param[in]  a_iv    IntVect location
 *  \param[in]  a_icomp Component index
 *  \return             Modifiable element
 *//*-----------------------------------------------------------------*/

template <typename T>
inline T&
BaseFab<T>::operator()(const IntVect& a_iv, const int a_icomp)
{
  return m_data[a_icomp*m_size + index(a_iv)];
}

/*--------------------------------------------------------------------*/
//  Obtain a linear index
/** \param[in]  a_iv    IntVect to index
 *  \return             Linear index
 *//*-----------------------------------------------------------------*/

template <typename T>
inline int
BaseFab<T>::index(IntVect a_iv) const
{
  CH_assert(m_box.contains(a_iv));
  a_iv -= m_box.loVect();  // Relative to lower corner
  return D_TERM(  a_iv[0]*m_stride[0],
                + a_iv[1]*m_stride[1],
                + a_iv[2]*m_stride[2]);
}

/*--------------------------------------------------------------------*/
//  Start of data for a component (internal use only)
/** \param[in]  a_icomp Component
 *  \return             Pointer to start of data
 *//*-----------------------------------------------------------------*/

template <typename T>
inline const T*
BaseFab<T>::dataPtr(const int a_icomp) const
{
  return &(m_data[a_icomp*m_size]);
}

/*--------------------------------------------------------------------*/
//  Start of data for a component (internal use only)
/** \param[in]  a_icomp Component
 *  \return             Pointer to start of data
 *//*-----------------------------------------------------------------*/

template <typename T>
inline T*
BaseFab<T>::dataPtr(const int a_icomp)
{
  return &(m_data[a_icomp*m_size]);
}

/*--------------------------------------------------------------------*/
//  Get spatial strides (internal use only)
/** \return             IntVect of spatial strides in the FAB
 *//*-----------------------------------------------------------------*/

template <typename T>
inline const IntVect&
BaseFab<T>::getStride() const
{
  return m_stride;
}

/*--------------------------------------------------------------------*/
//  Get component stride (internal use only)
/** \return             Stride from one component to another
 *//*-----------------------------------------------------------------*/

template <typename T>
inline int
BaseFab<T>::getComponentStride() const
{
  return m_size;
}


/*******************************************************************************
 *
 * Type definitions
 *
 ******************************************************************************/

// Real valued BaseFabs are called FArrayBox
using FArrayBox = BaseFab<Real>;


/*******************************************************************************
 *
 * Temporary BaseFab
 *
 ******************************************************************************/

#ifdef USE_STACK

// Defines FAB on the stack
#define FABSTACKTEMP(fabname, box, ncomp)                               \
  Real fabname ## _smem[box.size()*ncomp];                              \
  FArrayBox fabname(box, ncomp, fabname ## _smem)

#else

// Allocate FAB on the heap
#define FABSTACKTEMP(fabname, box, ncomp)                               \
  FArrayBox fabname(box, ncomp)

#endif  /* defined USE_STACK */

#endif  /* ! defined _BASEFAB_H_ */
