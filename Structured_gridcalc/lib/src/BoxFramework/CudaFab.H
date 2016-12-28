
#ifndef _CUDAFAB_H_
#define _CUDAFAB_H_
#ifdef __CUDACC__  /* Can only be compiled with nvcc */


/******************************************************************************/
/**
 * \file CudaFab.H
 *
 * \brief Class providing data for a box on GPU
 *
 *//*+*************************************************************************/

#include "Parameters.H"
#include "Box.H"
#include "CudaSupport.H"

#define HOSTDEVICE __host__ __device__
#define DEVICE __device__


/*******************************************************************************
 */
///  Data layout for a BaseFab on the CPU (used in Cuda)
/**
 ******************************************************************************/

template <typename T>
class BaseFabData
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

//--Friends

template <typename S>
friend class CudaFab;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor (no allocation)
  DEVICE BaseFabData() = delete;

  /// Copy constructor
  DEVICE BaseFabData(const BaseFabData&) = delete;

  /// Move constructor
  DEVICE BaseFabData(BaseFabData&&) noexcept = delete;

  /// Assignment constructor
  DEVICE BaseFabData& operator=(const BaseFabData&) = delete;

  /// Move assignment constructor
  DEVICE BaseFabData& operator=(BaseFabData&&) noexcept = delete;

  /// Destructor
  DEVICE ~BaseFabData() = delete;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Box m_box;                          ///< Box defining data
  IntVect m_stride;                   ///< Stride for indexing
  int m_ncomp;                        ///< Number of components
  int m_size;                         ///< Size of the box
  T* m_data;                          ///< Data on CPU
  AllocBy m_allocBy;                  ///< Method of allocation
#ifdef USE_GPU
public:
  SymbolPair<T> m_dataSymbol;         ///< Pointers to data on host and device
#endif
};


/*******************************************************************************
 */
///  Data for a box (fortran array box) aliased on GPU
/**
 ******************************************************************************/

template <typename T>
class CudaFab
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor (no allocation)
  HOSTDEVICE CudaFab();

  /// Constructor with sizes
  HOSTDEVICE CudaFab(T *const a_alias, const Box& a_box, const int a_ncomp);

  /// Copy construct from BaseFabData
  HOSTDEVICE CudaFab(const BaseFabData<T>& a_baseFab);

  /// Copy constructor (shallow)
  DEVICE CudaFab(const CudaFab&) = default;

  /// Move constructor
  DEVICE CudaFab(CudaFab&&) noexcept = default;

  /// Assignment constructor (shallow)
  DEVICE CudaFab& operator=(const CudaFab&) = default;

  /// Move assignment constructor
  DEVICE CudaFab& operator=(CudaFab&&) noexcept = default;

  /// Weak construction with sizes
  HOSTDEVICE void define(T *const a_alias, const Box& a_box, const int a_ncomp);

  /// Weak construction from BaseFabData
  HOSTDEVICE void define(const BaseFabData<T>& a_baseFab);

  /// Weak construction from CudaFab using multiple threads to copy
  DEVICE void define(const CudaFab<T>& a_cudaFab,
                     const int         a_idxThr0,
                     const int         a_numThr);

  /// Destructor
  DEVICE ~CudaFab() = default;


/*==============================================================================
 * Members functions
 *============================================================================*/

public:

  /// Number of threads required to perform a cooperative copy
  DEVICE static int numThrCopy();
  
  /// Return the box
  DEVICE const Box& box() const;

  /// Return the number of components
  DEVICE int ncomp() const;

  /// Return the total number of elements
  DEVICE int size() const;

  /// Return the total number of bytes used
  DEVICE size_t sizeBytes() const;

  /// Constant access to an element
  DEVICE const T& operator()(const IntVect& a_iv, const int a_icomp) const;

  /// Access to an element
  DEVICE T& operator()(const IntVect& a_iv, const int a_icomp);

  /// Obtain a linear index (internal and testing use only)
  DEVICE int index(IntVect a_iv) const;

  /// Shift the Fab by a_i cells in a specific direction
  DEVICE void shift(const int a_i, const int a_dir);
  
  /// Start of data for a component (internal use only)
  DEVICE const T* dataPtr(const int a_icomp = 0) const;

  /// Start of data for a component (internal use only)
  DEVICE T* dataPtr(const int a_icomp = 0);

  /// Get spatial strides (internal use only)
  DEVICE const IntVect& getStride() const;

  /// Get component stride (internal use only)
  DEVICE int getComponentStride() const;


/*==============================================================================
 * Protected members functions
 *============================================================================*/

protected:

  /// Set strides
  HOSTDEVICE void setStride();


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Box m_box;                          ///< Box defining data
  IntVect m_stride;                   ///< Stride for indexing
  int m_ncomp;                        ///< Number of components
  int m_size;                         ///< Size of the box
  T* m_data;                          ///< Data on GPU with offset
};


/*******************************************************************************
 *
 * Class CudaFab: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor (no allocation)
/*--------------------------------------------------------------------*/

template <typename T>
HOSTDEVICE inline
CudaFab<T>::CudaFab()
{
}

/*--------------------------------------------------------------------*/
//  Constructor with sizes
/** \param[in] a_alias  Memory location to reference
 *  \param[in] a_box    Box defining array dimensions
 *  \param[in] a_ncomp  Number of components
 *//*-----------------------------------------------------------------*/

template <typename T>
HOSTDEVICE inline
CudaFab<T>::CudaFab(T *const a_alias, const Box& a_box, const int a_ncomp)
  :
  m_box(a_box),
  m_ncomp(a_ncomp)
{
  setStride();
  m_data = a_alias + m_box.getOffset(m_stride);
}

/*--------------------------------------------------------------------*/
//  Copy construct from BaseFabData
/** BaseFabData has same layout as full BaseFab.  A single thread does
 *  the construction if on GPU.
 *  \param[in] a_baseFab
 *                      BaseFab data
 *//*-----------------------------------------------------------------*/

template <typename T>
HOSTDEVICE inline
CudaFab<T>::CudaFab(const BaseFabData<T>& a_baseFab)
  :
  m_box(a_baseFab.m_box),
  m_stride(a_baseFab.m_stride),
  m_data(static_cast<T*>(a_baseFab.m_dataSymbol.device) +
         m_box.getOffset(m_stride)),
  m_ncomp(a_baseFab.m_ncomp),
  m_size(a_baseFab.m_size)
{
}

/*--------------------------------------------------------------------*/
//  Weak construction with sizes
/** \param[in] a_alias  Memory location to reference
 *  \param[in] a_box    Box defining array dimensions
 *  \param[in] a_ncomp  Number of components
 *//*-----------------------------------------------------------------*/

template <typename T>
HOSTDEVICE inline
void
CudaFab<T>::define(T *const a_alias, const Box& a_box, const int a_ncomp)
{
  m_box = a_box;
  m_ncomp = a_ncomp;
  setStride();
  m_data = a_alias + m_box.getOffset(m_stride);
}

/*--------------------------------------------------------------------*/
//  Weak construction from BaseFabData
/** BaseFabData has same layout as full BaseFab.  A single thread does
 *  the construction if on GPU.
 *  \param[in] a_baseFab
 *                      BaseFab data
 *//*-----------------------------------------------------------------*/

template <typename T>
HOSTDEVICE inline
void
CudaFab<T>::define(const BaseFabData<T>& a_baseFab)
{
  m_box = a_baseFab.m_box;
  m_stride = a_baseFab.m_stride;
  m_data = static_cast<T*>(a_baseFab.m_dataSymbol.device) +
    m_box.getOffset(m_stride);
  m_ncomp = a_baseFab.m_ncomp;
  m_size = a_baseFab.m_size;
}

/*--------------------------------------------------------------------*/
//  Weak construction from CudaFab using multiple threads to copy
/** Multiple threads participate in the copy, each copying an 32-bit
 *  piece of data.
 *  \param[in] a_cudaFab
 *                      Source CudaFab
 *  \param[in] a_idxThr0
 *                      Index in block of the first thread to use.
 *                      All threads must be in blockDim.x.
 *  \param[in] a_numThr Number of threads to use in the copy.  This
 *                      must be >= numThrCopy()
 *//*-----------------------------------------------------------------*/

template <typename T>
DEVICE inline
void
CudaFab<T>::define(const CudaFab<T>& a_cudaFab,
                   const int         a_idxThr0,
                   const int         a_numThr)
{
  const int idxThr = threadIdx.x - a_idxThr0;
  if (idxThr >= 0 && idxThr < numThrCopy())
    {
      CH_assert(a_numThr <= blockDim.x);
      CH_assert(a_numThr >= numThrCopy());
      CH_assert(static_cast<void*>(this) == static_cast<void*>(&m_box));
      const int* src = reinterpret_cast<const int*>(&a_cudaFab);
      int* dst = reinterpret_cast<int*>(this);
      dst[idxThr] = src[idxThr];
    }
}

/*--------------------------------------------------------------------*/
//  Number of threads required to perform a cooperative copy
/*--------------------------------------------------------------------*/

template <typename T>
DEVICE inline int
CudaFab<T>::numThrCopy()
{
  // Assumes integer copy of everything
  return sizeof(CudaFab<T>)/4;
}

/*--------------------------------------------------------------------*/
//  Return the box
/*--------------------------------------------------------------------*/

template <typename T>
DEVICE inline const Box&
CudaFab<T>::box() const
{
  return m_box;
}

/*--------------------------------------------------------------------*/
//  Return the number of components
/*--------------------------------------------------------------------*/

template <typename T>
DEVICE inline int
CudaFab<T>::ncomp() const
{
  return m_ncomp;
}

/*--------------------------------------------------------------------*/
//  Return the total number of elements
/*--------------------------------------------------------------------*/

template <typename T>
DEVICE inline int
CudaFab<T>::size() const
{
  return m_ncomp*m_size;
}

/*--------------------------------------------------------------------*/
//  Return the total number of bytes used
/*--------------------------------------------------------------------*/

template <typename T>
DEVICE inline size_t
CudaFab<T>::sizeBytes() const
{
  return ((size_t)size())*sizeof(T);
}

/*--------------------------------------------------------------------*/
//  Constant access to an element
/*--------------------------------------------------------------------*/

template <typename T>
DEVICE inline const T&
CudaFab<T>::operator()(const IntVect& a_iv, const int a_icomp) const
{
  return m_data[a_icomp*m_size + index(a_iv)];
}

/*--------------------------------------------------------------------*/
//  Access to an element
/*--------------------------------------------------------------------*/

template <typename T>
DEVICE inline T&
CudaFab<T>::operator()(const IntVect& a_iv, const int a_icomp)
{
  return m_data[a_icomp*m_size + index(a_iv)];
}

/*--------------------------------------------------------------------*/
//  Obtain a linear index
/** \param[in]  a_iv    IntVect to index
 *  \return             Linear index
 *//*-----------------------------------------------------------------*/

template <typename T>
DEVICE inline
int
CudaFab<T>::index(IntVect a_iv) const
{
  CH_assert(m_box.contains(a_iv));
  return D_TERM(  a_iv[0],
                + a_iv[1]*m_stride[1],
                + a_iv[2]*m_stride[2]);
}

/*--------------------------------------------------------------------*/
//  Shift the Fab by a_i cells in a specific direction
/** \param[in]  a_i     Number of cells to shift
 *  \param[in]  a_dir   Direction to shift
 *//*-----------------------------------------------------------------*/

template <typename T>
DEVICE inline
void
CudaFab<T>::shift(const int a_i, const int a_dir)
{
  // Substract the old offset
  m_data -= m_box.getOffset(m_stride);
  // Shift
  m_box.shift(a_i, a_dir);
  // Add the new offset
  m_data += m_box.getOffset(m_stride);
}

/*--------------------------------------------------------------------*/
//  Start of data for a component (internal use only)
/** \param[in]  a_icomp Component
 *  \return             Pointer to start of data
 *//*-----------------------------------------------------------------*/

template <typename T>
DEVICE inline const T*
CudaFab<T>::dataPtr(const int a_icomp) const
{
  // return &(m_data[a_icomp*m_size]);
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Start of data for a component (internal use only)
/** \param[in]  a_icomp Component
 *  \return             Pointer to start of data
 *//*-----------------------------------------------------------------*/

template <typename T>
DEVICE inline T*
CudaFab<T>::dataPtr(const int a_icomp)
{
  // return &(m_data[a_icomp*m_size]);
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Get spatial strides (internal use only)
/** \return             IntVect of spatial strides in the FAB
 *//*-----------------------------------------------------------------*/

template <typename T>
DEVICE inline const IntVect&
CudaFab<T>::getStride() const
{
  return m_stride;
}

/*--------------------------------------------------------------------*/
//  Get component stride (internal use only)
/** \return             Stride from one component to another
 *//*-----------------------------------------------------------------*/

template <typename T>
DEVICE inline int
CudaFab<T>::getComponentStride() const
{
  return m_size;
}

/*--------------------------------------------------------------------*/
//  Set strides
/*--------------------------------------------------------------------*/

template <typename T>
HOSTDEVICE inline
void
CudaFab<T>::setStride()
{
  const IntVect& lo = m_box.loVect();
  const IntVect& hi = m_box.hiVect();
  CH_assert(lo <= hi);
  // Set strides
  D_TERM(m_stride[0] = 1;,
         m_stride[1] = m_stride[0]*(hi[0] - lo[0] + 1);,
         m_stride[2] = m_stride[1]*(hi[1] - lo[1] + 1);)
  // Set size
  m_size = m_stride[g_SpaceDim-1]*(hi[g_SpaceDim-1] - lo[g_SpaceDim-1] + 1);
} 


/*******************************************************************************
 */
///  Data for a moving window in cache implemented with slabs
/**
 ******************************************************************************/

template <typename T, int MxSz>
class SlabFab : public CudaFab<T>
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor (no allocation)
  DEVICE SlabFab()
    :
    CudaFab<T>()
    { }

  /// Weak construction with sizes
  //  Does not synchronize
  DEVICE void define(T *const          a_alias,
                     const Box&        a_box,
                     const int         a_ncomp,
                     const int         a_nrmDir,
                     const int         a_locNDBeg,
                     const int         a_locNDEnd,
                     IntVect&          a_vecLS,
                     const int         a_compBegSrc,
                     const CudaFab<T>& a_src,
                     const int         a_numThrLS);

  /// Copy constructor (shallow)
  DEVICE SlabFab(const SlabFab&) = default;

  /// Move constructor
  DEVICE SlabFab(SlabFab&&) noexcept = default;

  /// Assignment constructor (shallow)
  DEVICE SlabFab& operator=(const SlabFab&) = default;

  /// Move assignment constructor
  DEVICE SlabFab& operator=(SlabFab&&) noexcept = default;

  /// Destructor
  DEVICE ~SlabFab() = default;


/*==============================================================================
 * Members functions
 *============================================================================*/

  /// Number of slabs
  DEVICE int numSlab() const;

  /// Constant access to an element
  DEVICE const T& operator()(const IntVect& a_iv, const int a_icomp) const;

  /// Access to an element
  DEVICE T& operator()(const IntVect& a_iv, const int a_icomp);

  /// Obtain a linear index (internal and testing use only)
  DEVICE int index(IntVect a_iv) const;

  /// Shift the Fab by a_i cells in a specific direction (all threads)
  // Synchronizes
  DEVICE void shift(int a_i);

  /// Shift the Fab by a_i cells in a specific direction and load new data
  // Synchronizes
  DEVICE void shift(int a_i, IntVect& a_vecThrLS);

  /// Select loading method
  DEVICE void selectLoadMethod(const int a_numThr);

  /// Load a slab
  DEVICE void loadSlab(const IntVect& a_vecThr,
                       const int      a_idxThr0,
                       const int      a_numThr);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:
  using CudaFab<T>::m_box;
  using CudaFab<T>::m_stride;
  using CudaFab<T>::m_ncomp;
  using CudaFab<T>::m_size;
  using CudaFab<T>::m_data;
  int m_pTable[MxSz];                 ///< Permuation table
  int m_nrmDir;                       ///< Orthogonal direction to slab
  int m_numPntSlab;                   ///< Number of points in single layer slab
  int m_numSubBox;                    ///< This is a key to the type of loader
                                      ///<  0 : undefined
                                      ///< -3 : m_numPntSlab = blockDim.x - the
                                      ///<      sizes are optimized for loading
                                      ///<      data
                                      ///< -2 : desired numThr <= m_numPntSlab
                                      ///< -1 : Use linear indices to schedule
                                      ///<      threads.  This has worst case
                                      ///<      performance
                                      ///< >0 : Number of subboxes.  The slab
                                      ///<      is split into subboxes for
                                      ///<      loading.  Must be a power of 2.
  int m_compBegSrc;                   ///< First component from source data
  const CudaFab<T>* m_src;            ///< Source data for the cache window
};


/*******************************************************************************
 *
 * Class SlabFab: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Weak construction with sizes
/** Does not synchronize
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
void
SlabFab<T, MxSz>::define(T *const          a_alias,
                         const Box&        a_box,
                         const int         a_ncomp,
                         const int         a_nrmDir,
                         const int         a_locNDBeg,
                         const int         a_locNDEnd,
                         IntVect&          a_vecLS,
                         const int         a_compBegSrc,
                         const CudaFab<T>& a_src,
                         const int         a_numThrLS)
{
  if (threadIdx.x == 0)
    {
      CudaFab<T>::m_box = a_box;
      CudaFab<T>::m_ncomp = a_ncomp;
      CudaFab<T>::setStride();
      IntVect lo = a_box.loVect();
      lo[a_nrmDir] = 0;  // No offset for permuted direction
      CudaFab<T>::m_data =
        a_alias + CudaFab<T>::m_box.vecToLin0(-lo, CudaFab<T>::m_stride);
      #pragma unroll
      for (int i = 0; i != MxSz; ++i)
        {
          m_pTable[i] = i;
        }
      m_nrmDir = a_nrmDir;
      IntVect dims = a_box.dimensions();
      const int numSlab = dims[a_nrmDir];
      CH_assert(numSlab > 1 && numSlab <= MxSz);
      dims[a_nrmDir] = 1;
      m_numPntSlab = dims.product();
      CH_assert(a_numThrLS <= blockDim.x);
      selectLoadMethod(a_numThrLS);  // Sets m_numSubBox
      m_compBegSrc = a_compBegSrc;
      m_src = &a_src;
    }
  __syncthreads();
  // Set vector index from linear index
  if (threadIdx.x < m_numPntSlab)
    {
      m_box.linToVec(threadIdx.x, a_vecLS);
    }
  // Load initial data
  for (int locND = a_locNDBeg; locND <= a_locNDEnd; ++locND)
    {
      a_vecLS[m_nrmDir] = locND;
      loadSlab(a_vecLS, 0, a_numThrLS);
    }
}

/*--------------------------------------------------------------------*/
//  Number of slabs
/** 
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
int
SlabFab<T, MxSz>::numSlab() const
{
  return m_box.dimensions()[m_nrmDir];
}

/*--------------------------------------------------------------------*/
//  Constant access to an element
/** 
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
const T&
SlabFab<T, MxSz>::operator()(const IntVect& a_iv, const int a_icomp) const
{
  return m_data[a_icomp*m_size + index(a_iv)];
}

/*--------------------------------------------------------------------*/
//  Access to an element
/** 
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
T&
SlabFab<T, MxSz>::operator()(const IntVect& a_iv, const int a_icomp)
{
  return m_data[a_icomp*m_size + index(a_iv)];
}

/*--------------------------------------------------------------------*/
//  Obtain a linear index (internal and testing use only)
/** 
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
int
SlabFab<T, MxSz>::index(IntVect a_iv) const
{
  CH_assert(m_box.contains(a_iv));
  // Indirection for index in direction m_nrmDir through permutation table
  a_iv[m_nrmDir] = m_pTable[a_iv[m_nrmDir] - m_box.loVect(m_nrmDir)];
  return D_TERM(  a_iv[0],
                + a_iv[1]*m_stride[1],
                + a_iv[2]*m_stride[2]);
}

/*--------------------------------------------------------------------*/
//  Shift the Fab by a_i cells in a specific direction (all threads)
/** Synchronizes
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
void
SlabFab<T, MxSz>::shift(int a_i)
{
  if (a_i == 0) return;
  const int nSlab = numSlab();
  CH_assert(blockDim.x >= nSlab);
  __syncthreads();  // Safe to move box
  if (threadIdx.x < nSlab)
    {
      if (threadIdx.x == 0)
        {
          m_box.shift(a_i, m_nrmDir);
        }
      // Get a_i >= 0 by adding chunks of size m_numSlab
      a_i += ((abs(a_i) - 1)/nSlab + 1)*nSlab;
      m_pTable[threadIdx.x] = (m_pTable[threadIdx.x] + a_i) % nSlab;
    }
  __syncthreads();  // Safe to use box
}

/*--------------------------------------------------------------------*/
//  Shift the Fab by a_i cells in a specific direction and load new
//  data
/** Synchronizes
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
void
SlabFab<T, MxSz>::shift(int a_i, IntVect& a_vecThrLS)
{
  if (a_i == 0) return;
  shift(a_i);
  // Number of slabs to load
  const int numSlabLoad = abs(a_i);
  int locNDBeg;  // Normal direction index of first slab to load
  int locNDEnd;  // Normal direction index of last slab to load
  if (a_i > 0)
    {
      locNDEnd = m_box.hiVect(m_nrmDir);
      locNDBeg = locNDEnd - numSlabLoad + 1;
    }
  else
    {
      locNDBeg = m_box.loVect(m_nrmDir);
      locNDEnd = locNDBeg + numSlabLoad - 1;
    }

  for (int locND = locNDBeg; locND <= locNDEnd; ++locND)
    {
      a_vecThrLS[m_nrmDir] = locND;
      loadSlab(a_vecThrLS, 0, blockDim.x);
    }
  __syncthreads();
}

/*--------------------------------------------------------------------*/
//  Select loading method
/** Call with a single thread if 'this' is in shared memory.
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
void
SlabFab<T, MxSz>::selectLoadMethod(const int a_numThr)
{
  if (blockDim.x == a_numThr && blockDim.x == m_numPntSlab)
    {
      m_numSubBox = -3;
    }
  else if (m_numPntSlab <= a_numThr)
    {
      m_numSubBox = -2;
    }
  else
    {
      CH_assert(false);
    }
}

/*--------------------------------------------------------------------*/
//  Load a slab
/** The loading method has already been selected by calling
 *  selectLoadMethod
 *//*-----------------------------------------------------------------*/

template <typename T, int MxSz>
DEVICE inline
void
SlabFab<T, MxSz>::loadSlab(const IntVect& a_vecThr,
                           const int      a_idxThr0,
                           const int      a_numThr)
{
  CH_assert(m_numSubBox != 0);
  CH_assert(a_vecThr[m_nrmDir] >= m_box.loVect(m_nrmDir) &&
            a_vecThr[m_nrmDir] <= m_box.hiVect(m_nrmDir));
  CH_assert(a_numThr <= blockDim.x);

  const int compEndSrc = m_compBegSrc + m_ncomp;

  switch (m_numSubBox)
    {
    case -3:

//--Number of threads in block is optimized for loading the slab

      CH_assert(blockDim.x == a_numThr);
      for (int iC = m_compBegSrc; iC != compEndSrc; ++iC)
        {
          this->operator()(a_vecThr, iC) = (*m_src)(a_vecThr, iC);
        }
      break;

    case -2:

//--Ideal loading, there are sufficient threads to load the entire slab at once

    {
      const int idxThr = threadIdx.x - a_idxThr0;
      CH_assert(m_numPntSlab <= a_numThr);
      if (idxThr >= 0 && idxThr < m_numPntSlab)
        {
          for (int iC = m_compBegSrc; iC != compEndSrc; ++iC)
            {
              this->operator()(a_vecThr, iC) = (*m_src)(a_vecThr, iC);
            }
        }
      break;
    }
    }  // Switch
}

#endif  /* __CUDACC__ */
#endif  /* ! defined _CUDAFAB_H_ */
