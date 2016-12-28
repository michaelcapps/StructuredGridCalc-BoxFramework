
/******************************************************************************/
/**
 * \file BaseFab.cpp
 *
 * \brief Non-inline definitions for classes in BaseFab.H
 *
 *//*+*************************************************************************/

// #define DEBUGFAB
#ifdef DEBUGFAB
#include <iostream>
#include <iomanip>
#endif

#include "BaseFab.H"
#include "BaseFabMacros.H"

#ifdef DEBUGFAB
  #define FABDBG(x) x
#else
  #define FABDBG(x) (void)0
#endif


/*******************************************************************************
 *
 * Class BaseFab: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor (no allocation)
/*--------------------------------------------------------------------*/

template <typename T>
BaseFab<T>::BaseFab()
  :
  m_box(),
  m_stride(IntVect::Zero),
  m_ncomp(0),
  m_size(0),
  m_data(nullptr),
  m_allocBy(AllocBy::none)
{
  FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): default construction\n");
}

/*--------------------------------------------------------------------*/
//  Constructor with sizes
/** \param[in]  a_box   Box defining array dimensions
 *  \param[in]  a_ncomp Number of components
 *  \param[in]  a_alias nullptr forces allocation.  Otherwise, memory
 *                      is aliased to this address.  Default parameter
 *                      is nullptr
 *//*-----------------------------------------------------------------*/

template <typename T>
BaseFab<T>::BaseFab(const Box& a_box, const int a_ncomp, T *const a_alias)
{
  m_box = a_box;
  m_ncomp = a_ncomp;
  m_data = a_alias;
  m_allocBy = (a_alias == nullptr) ? AllocBy::array : AllocBy::alias;
  setStride();
  allocate();
  FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): construction with sizes\n");
}

/*--------------------------------------------------------------------*/
//  Constructor with sizes and a default value
/** \param[in]  a_box   Box defining array dimensions
 *  \param[in]  a_ncomp Number of components
 *  \param[in]  a_val   Value to assign
 *  \param[in]  a_alias nullptr forces allocation.  Otherwise, memory
 *                      is aliased to this address.  Default parameter
 *                      is nullptr
 *//*-----------------------------------------------------------------*/

template <typename T>
BaseFab<T>::BaseFab(const Box& a_box, const int a_ncomp, const T& a_val, T *const a_alias)
{
  m_box = a_box;
  m_ncomp = a_ncomp;
  m_data = a_alias;
  m_allocBy = (a_alias == nullptr) ? AllocBy::array : AllocBy::alias;
  setStride();
  allocate();
  setVal(a_val);
  FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): construction with sizes and default value\n");
}

/*--------------------------------------------------------------------*/
//  Move constructor
/** Moving BaseFabs built as an alias will cause an error
 *  \param[in]  a_fab   Rvalue RHS
 *//*-----------------------------------------------------------------*/

template <typename T>
BaseFab<T>::BaseFab(BaseFab&& a_fab) noexcept
  :
  m_box(std::move(a_fab.m_box)),
  m_stride(std::move(a_fab.m_stride)),
  m_ncomp(a_fab.m_ncomp),
  m_size(a_fab.m_size),
  m_data(a_fab.m_data),
  m_allocBy(a_fab.m_allocBy)
#ifdef USE_GPU
  ,m_dataSymbol(a_fab.m_dataSymbol)
#endif
{
  FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): move construction\n");
  CH_assert(a_fab.m_allocBy != AllocBy::alias);
  // Ensure a_fab will not delete
  a_fab.m_data = nullptr;
#ifdef USE_GPU
      a_fab.m_dataSymbol.host = nullptr;
      a_fab.m_dataSymbol.device = nullptr;
#endif
}

/*--------------------------------------------------------------------*/
//  Move assignment constructor
/** Moving BaseFabs built as an alias will cause an error
 * \param[in]  a_fab    Rvalue RHS
 *//*-----------------------------------------------------------------*/

template <typename T>
BaseFab<T>&
BaseFab<T>::operator=(BaseFab&& a_fab) noexcept
{
  //if(this != std::move(&a_fab))
  if(this != &a_fab)
  {
    deallocate();
    m_box = a_fab.m_box;
    m_stride = a_fab.m_stride;
    m_ncomp= a_fab.m_ncomp;
    m_data = a_fab.m_data;
    m_allocBy = a_fab.m_allocBy;
    FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): move construction\n");
    CH_assert(a_fab.m_allocBy != AllocBy::alias);
    // // Ensure a_fab will not delete
    a_fab.m_data = nullptr;
  }
  return *this;
}

/*--------------------------------------------------------------------*/
//  Weak construction with sizes
/** \param[in]  a_box   Box defining array dimensions
 *  \param[in]  a_ncomp Number of components
 *  \param[in]  a_alias nullptr forces allocation.  Otherwise, memory
 *                      is aliased to this address.  Default parameter
 *                      is nullptr
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::define(const Box& a_box, const int a_ncomp, T *const a_alias)
{
  FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): define\n");
  deallocate();
  m_box = a_box;
  m_ncomp = a_ncomp;
  m_data = a_alias;
  m_allocBy = (a_alias == nullptr) ? AllocBy::array : AllocBy::alias;
  allocate();
}

/*--------------------------------------------------------------------*/
//  Weak construction with sizes and a default value
/** \param[in]  a_box   Box defining array dimensions
 *  \param[in]  a_ncomp Number of components
 *  \param[in]  a_val   Value to assign
 *  \param[in]  a_alias nullptr forces allocation.  Otherwise, memory
 *                      is aliased to this address.  Default parameter
 *                      is nullptr
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::define(const Box& a_box, const int a_ncomp,const T& a_val, T *const a_alias)
{
  FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): define\n");
  deallocate();
  m_box = a_box;
  m_ncomp = a_ncomp;
  m_data = a_alias;
  m_allocBy = (a_alias == nullptr) ? AllocBy::array : AllocBy::alias;
  allocate();
  setVal(a_val);
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/

template <typename T>
BaseFab<T>::~BaseFab()
{
  FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
         << "): destructor\n");
  deallocate();
}

/*--------------------------------------------------------------------*/
//  Assign a constant to all components
/** \param[in]  a_val   Value to assign
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::setVal(const T& a_val)
{
  T* p = dataPtr(0);
  for (int n = size(); n--;)
    {
      *p++ = a_val;
    }
}

/*--------------------------------------------------------------------*/
//  Assign a constant to a single component
/** \param[in]  a_icomp Component index
 *  \param[in]  a_val   Value to assign
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::setVal(const int a_icomp, const T& a_val)
{
  CH_assert(a_icomp >= 0 && a_icomp < m_ncomp);
  T* p = dataPtr(a_icomp);
  for (int n = m_box.size(); n--;)
    {
      *p++ = a_val;
    }
}

/*--------------------------------------------------------------------*/
//  Copy a portion of another BaseFab, same region and all components
/** \param[in]  a_box   Region to copy
 *  \param[in]  a_src   Source BaseFab
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::copy(const Box&     a_box,
                 const BaseFab& a_src)
{
  CH_assert(a_src.ncomp() == m_ncomp);
  
  MD_ARRAY_RESTRICT(arrSrc, a_src);
  MD_ARRAY_RESTRICT(arrDst, *this);

  for (int ic = 0; ic != m_ncomp; ++ic)
    {
      MD_BOXLOOP_OMP(a_box, i)
        {
          arrDst[MD_IX(i, ic)] = arrSrc[MD_IX(i,ic)];
        }
    }
}


/*--------------------------------------------------------------------*/
//  Copy a portion of another BaseFab
/** \param[in]  a_dstBox
 *                      Region to copy to in this BaseFab
 *  \param[in]  a_dstComp
 *                      Start index for destination components
 *  \param[in]  a_src   Source BaseFab
 *  \param[in]  a_srcBox
 *                      Region to copy from in source BaseFab
 *  \param[in]  a_srcComp
 *                      Start index for source components
 *  \param[in]  a_numComp
 *                      Number of components to copy
 *  \param[in]  a_compFlags
 *                      Components are selectively copied based on the
 *                      bits in a_compFlags.  This is only used for
 *                      components < (number of bits in unsigned,
 *                      normally 32).  The bits flags start at bit 0 =
 *                      dstComp/srcComp.  Default, all components are
 *                      used.
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::copy(const Box&     a_dstBox,
                 const int      a_dstComp,
                 const BaseFab& a_src,
                 const Box&     a_srcBox,
                 const int      a_srcComp,
                 const int      a_numComp,
                 const unsigned a_compFlags)
{
  const IntVect len = a_dstBox.dimensions();
  CH_assert(this != &a_src);
  CH_assert(len == a_srcBox.dimensions());
  CH_assert(m_box.contains(a_dstBox));
  CH_assert(a_src.box().contains(a_srcBox));
  CH_assert(a_dstComp >= 0 && (a_dstComp + a_numComp) <= m_ncomp);
  CH_assert(a_srcComp >= 0 && (a_srcComp + a_numComp) <= a_src.ncomp());

  MD_ARRAY_RESTRICT(arrSrc, a_src);
  MD_ARRAY_RESTRICT(arrDst, *this);
  IntVect offset = a_srcBox.loVect() - a_dstBox.loVect();
  for (int ic = 0; ic != a_numComp; ++ic)
    {
      const int iDstC = ic + a_dstComp;
      if ((iDstC >= (int)(8*sizeof(unsigned))) || (a_compFlags & (1 << iDstC)))
        {
          const int iSrcC = ic + a_srcComp;
          MD_BOXLOOP_OMP(a_dstBox, i)
            {
              arrDst[MD_IX(i, iDstC)] = arrSrc[MD_OFFSETIV(i,+,offset, iSrcC)];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Linearize data in a region and place in a buffer
/** \param[out] a_buffer
 *                      Linear buffer filled with data
 *  \param[in]  a_region
 *                      Box describing region to place in buffer
 *  \param[in]  a_startComp
 *                      Start of components to add to buffer
 *  \param[in]  a_endComp
 *                      One past last component to replace with buffer
 *  \param[in]  a_compFlags
 *                      Components are selectively copied based on the
 *                      bits in a_compFlags.  This is only used for
 *                      components < (number of bits in unsigned,
 *                      normally 32).  The bits flags start at bit 0 =
 *                      a_startComp.  Default, all components are
 *                      used.
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::linearOut(void *const    a_buffer,
                      const Box&     a_region,
                      const int      a_startComp,
                      const int      a_endComp,
                      const unsigned a_compFlags) const
{
  CH_assert(a_buffer != NULL);
  CH_assert(m_box.contains(a_region));
  CH_assert(a_startComp >= 0);
  CH_assert(a_endComp >= a_startComp && a_endComp <= m_ncomp);

  MD_ARRAY_RESTRICT(arr, *this);
  T* p = static_cast<T*>(a_buffer);
  for (int ic = a_startComp; ic != a_endComp; ++ic)
    {
      if ((ic >= (int)(8*sizeof(unsigned))) || (a_compFlags & (1 << ic)))
        {
          MD_BOXLOOP_OMP(a_region, i)
            {
              *p++ = arr[MD_IX(i, ic)];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Replace data in a region from a linear buffer
/** \param[int] a_buffer
 *                      Linear buffer filled with data
 *  \param[in]  a_region
 *                      Box describing region to replace with buffer
 *  \param[in]  a_startComp
 *                      Start of components to replace with buffer
 *  \param[in]  a_endComp
 *                      One past last component to replace with buffer
 *  \param[in]  a_compFlags
 *                      Components are selectively copied based on the
 *                      bits in a_compFlags.  This is only used for
 *                      components < (number of bits in unsigned,
 *                      a_startComp.  Default, all components are
 *                      used.
 *//*-----------------------------------------------------------------*/
template <typename T>
void
BaseFab<T>::linearIn(const void* const a_buffer,
                     const Box&     a_region,
                     const int      a_startComp,
                     const int      a_endComp,
                     const unsigned a_compFlags)
{
  CH_assert(a_buffer != NULL);
  CH_assert(m_box.contains(a_region));
  CH_assert(a_startComp >= 0);
  CH_assert(a_endComp >= a_startComp && a_endComp <= m_ncomp);

  MD_ARRAY_RESTRICT(arr, *this);
  T* p = (T*)a_buffer;
  for (int ic = a_startComp; ic != a_endComp; ++ic)
    {
      if ((ic >= (int)(8*sizeof(unsigned))) || (a_compFlags & (1 << ic)))
        {
          MD_BOXLOOP_OMP(a_region, i)
            {
              arr[MD_IX(i, ic)] = *p++;
            }
        }
    }
}


#ifdef USE_GPU
/*--------------------------------------------------------------------*/
//  Copy array to device
/*--------------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::copyToDevice() const
{
  CU_SAFE_CALL(cudaMemcpy(m_dataSymbol.device,
                          m_dataSymbol.host,
                          sizeBytes(),
                          cudaMemcpyHostToDevice));
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy array to device
/** \param[in]  a_stream
 *                      Stream index (defaults to default stream)
 *//*-----------------------------------------------------------------*/

template <typename T>
inline void
BaseFab<T>::copyToDeviceAsync(cudaStream_t a_stream) const
{
  CU_SAFE_CALL(cudaMemcpyAsync(m_dataSymbol.device,
                               m_dataSymbol.host,
                               sizeBytes(),
                               cudaMemcpyHostToDevice,
                               a_stream));
}

/*--------------------------------------------------------------------*/
//  Copy array to host
/*--------------------------------------------------------------------*/

template <typename T>
inline void
BaseFab<T>::copyToHost()
{
  CU_SAFE_CALL(cudaMemcpy(m_dataSymbol.host,
                          m_dataSymbol.device,
                          sizeBytes(),
                          cudaMemcpyDeviceToHost));
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy array to host
/** \param[in]  a_stream
 *                      Stream index (defaults to default stream)
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::copyToHostAsync(cudaStream_t a_stream)
{
  CU_SAFE_CALL(cudaMemcpyAsync(m_dataSymbol.host,
                               m_dataSymbol.device,
                               sizeBytes(),
                               cudaMemcpyDeviceToHost,
                               a_stream));
}
#endif  /* CUDA */

/*--------------------------------------------------------------------*/
//  Set strides
/*--------------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::setStride()
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

/*--------------------------------------------------------------------*/
//  Allocate memory
/*--------------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::allocate()
{
  setStride();
  if (m_allocBy == AllocBy::array)
    {
      deallocate();
#ifdef USE_GPU
      const size_t numBytes = sizeBytes();
      CU_SAFE_CALL(cudaMallocHost(&(m_dataSymbol.host), numBytes));
      m_data = m_dataSymbol.host;
      CU_SAFE_CALL(cudaMalloc(&(m_dataSymbol.device), numBytes));
#else
      m_data = new T[size()];
#endif
      FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
             << "): new\n");
    }
}

/*--------------------------------------------------------------------*/
//  Deallocate memory
/*--------------------------------------------------------------------*/

template <typename T>
void
BaseFab<T>::deallocate()
{
  if (m_allocBy == AllocBy::array && m_data != nullptr)
    {
      FABDBG(std::cout << "BaseFab (" << std::setw(14) << m_data
             << "): delete\n");
#ifdef USE_GPU
      CU_SAFE_CALL(cudaFreeHost(m_dataSymbol.host));
      m_dataSymbol.host = nullptr;
      CU_SAFE_CALL(cudaFree(m_dataSymbol.device));
#else
      delete[] m_data;
#endif
      m_data = nullptr;
    }
}


/*******************************************************************************
 *
 * Explicit instantiations of class BaseFab
 *
 ******************************************************************************/

template class BaseFab<bool>;
template class BaseFab<char>;
template class BaseFab<int>;
template class BaseFab<unsigned>;
template class BaseFab<Real>;
