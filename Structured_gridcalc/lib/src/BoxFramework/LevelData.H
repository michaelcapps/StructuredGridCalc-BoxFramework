
#ifndef _LEVELDATA_H_
#define _LEVELDATA_H_


/******************************************************************************/
/**
 * \file LevelData.H
 *
 * \brief Class storing data on a DisjointBoxLayout
 *
 *//*+*************************************************************************/

#include <iostream>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Parameters.H"
#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "Copier.H"

#ifdef USE_GPU
#include "CudaSupport.H"
#endif

#define USE_MPIWAITALL  // Use Waitany if commented out


/*******************************************************************************
 */
///  Data for a layout of boxes
/**
 ******************************************************************************/

template <typename T>
class LevelData
{


/*====================================================================*
 * Public constructors and destructors
 *====================================================================*/

public:

  /// Default constructor
  LevelData();

  /// Constructor with dbl, ncomp, and nghost
  LevelData(DisjointBoxLayout& a_dbl, int a_ncomp, int a_nghost);
  
  /// Const constructor
  LevelData(const DisjointBoxLayout& a_dbl,
                        const int                a_ncomp,
                        const int                a_nghost);

  /// Define (weak construction)
  void define(const DisjointBoxLayout& a_dbl,
              const int                a_ncomp,
              const int                a_nghost);


/*====================================================================*
 * Members functions
 *====================================================================*/

public:

  /// Index with a DataIterator
  T& operator[](const LayoutIterator& a_lit);

  /// Constant index with a DataIterator
  const T& operator[](const LayoutIterator& a_lit) const;

  /// Index with a BoxIndex
  T& operator[](const BoxIndex& a_bidx);

  /// Constant index with a BoxIndex
  const T& operator[](const BoxIndex& a_bidx) const;

  /// Assign a constant to all components
  void setVal(const typename T::value_type& a_val);

  /// Assign a constant to a single component
  void setVal(const int a_icomp, const typename T::value_type& a_val);

  /// Unique identifying tag (from the DBL)
  size_t tag() const;

  /// Constant access to elements with a linear index
  //  FOR INTERNAL USE AND TESTING ONLY
  const T& getLinear(const int a_idx) const;

  /// Number of boxes
  int size() const;

  /// Number of components
  int ncomp() const;

  /// Number of ghosts
  int nghost() const;

  /// The layout of boxes
  const DisjointBoxLayout& disjointBoxLayout() const;

  /// Exchange to fill ghost cells
  void exchange(Copier& a_copier);

  /// Begin exchange to fill ghost cells
  void exchangeBegin(Copier& a_copier);

  /// End exchange to fill ghost cells
  void exchangeEnd(Copier& a_copier);

  /// Write CGNS solution data to a file (specialized for BaseFab<Real>)
#ifndef NO_CGNS
  int writeCGNSSolData(const int                a_indexFile,
                       const int                a_indexBase,
                       const int                a_indexZoneOffset,
                       const char* const *const a_varNames) const;
#endif

#ifdef USE_GPU
  /// Copy all elements to device
  void copyToDevice() const;

  /// Asynchronous copy all elements to device
  void copyToDeviceAsync(cudaStream_t a_stream = 0) const;

  /// Copy all elements to host
  void copyToHost();

  /// Asynchronous copy all elements to host
  void copyToHostAsync(cudaStream_t a_stream = 0);
#endif


/*====================================================================*
 * Data members
 *====================================================================*/

protected:

  DisjointBoxLayout m_disjointBoxLayout;
                                      ///< Disjoint box layout this level data
                                      ///< is built on
  std::vector<T> m_data;              ///< The data (usually BaseFabs)
  int m_ncomp;                        ///< Number of components
  int m_nghost;                       ///< Number of ghosts
};


/*******************************************************************************
 periodic_box*
 * Class LevelData: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/

template <typename T>
inline
LevelData<T>::LevelData()
  :
  m_disjointBoxLayout(),
  m_data(),
  m_ncomp(0),
  m_nghost(0)
{
}

/// Constructor with DisjointBoxLayout, ncomp, and nghost
template <typename T>
inline
LevelData<T>::LevelData(DisjointBoxLayout& a_dbl, int a_ncomp, int a_nghost)
  :
  m_disjointBoxLayout(a_dbl),
  m_data(),
  m_ncomp(a_ncomp),
  m_nghost(a_nghost)
{
  m_data.resize(a_dbl.localSize());
  typename std::vector<T>::iterator it = m_data.begin();
  for(DataIterator dit(a_dbl); dit.ok(); ++dit)
  {
    Box temp_box = a_dbl[dit];
    temp_box.grow(m_nghost);
    (*it).define(temp_box,m_ncomp);
    ++it;
  }
}

/*
/// Weak constructor
template <typename T>
inline void
LevelData<T>::define(const DisjointBoxLayout& a_dbl,const int a_ncomp,const int a_nghost)
{
  m_disjointBoxLayout = a_dbl;
  m_data.resize(a_dbl.localSize());
  m_ncomp = a_ncomp;
  m_nghost = a_nghost;
  typename std::vector<T>::iterator it = m_data.begin();
  for(DataIterator dit(a_dbl); dit.ok(); ++dit)
    {
      Box temp_box = a_dbl[dit];
      temp_box.grow(m_nghost);
      (*it).define(temp_box,m_ncomp);
      ++it;
    }
}
*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_dbl   The disjoint box layout
 *  \param[in]  a_ncomp Number of components
 *  \param[in]  a_nghost
 *                      Number of ghost cells
 *//*-----------------------------------------------------------------*/

template <typename T>
LevelData<T>::LevelData(const DisjointBoxLayout& a_dbl,
                        const int                a_ncomp,
                        const int                a_nghost)
  :
  m_disjointBoxLayout(a_dbl),
  m_data(),
  m_ncomp(a_ncomp),
  m_nghost(a_nghost)
{
  m_data.resize(size());
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      Box box = m_disjointBoxLayout[dit];
      box.grow(a_nghost);
      this->operator[](dit).define(box, a_ncomp);
    }
}

/*--------------------------------------------------------------------*/
//  Define (weak construction)
/** \param[in]  a_dbl   The disjoint box layout
 *  \param[in]  a_ncomp Number of components
 *  \param[in]  a_nghost
 *                      Number of ghost cells
 *//*-----------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::define(const DisjointBoxLayout& a_dbl,
                     const int                a_ncomp,
                     const int                a_nghost)
{
  m_disjointBoxLayout = a_dbl;
  m_ncomp = a_ncomp;
  m_nghost = a_nghost;
  m_data.resize(size());
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      Box box = m_disjointBoxLayout[dit];
      box.grow(a_nghost);
      this->operator[](dit).define(box, a_ncomp);
    }
}

/*--------------------------------------------------------------------*/
//  Index with a LayoutIterator
/** \param[in]  a_lit   Layout iterator
 *  \return             Element at the index of the data iterator
 *//*-----------------------------------------------------------------*/

template <typename T>
inline T&
LevelData<T>::operator[](const LayoutIterator& a_lit)
{
  // Make sure we are indexing with an iterator built on this DBL
  CH_assert(a_lit.tag() == tag());
  const BoxIndex& bidx = (*a_lit);
  CH_assert(bidx.localIndex() < m_disjointBoxLayout.localSize());
  CH_assert(bidx.globalIndex() >= m_disjointBoxLayout.localIdxBegin() &&
            bidx.globalIndex() <  m_disjointBoxLayout.localIdxEnd());
  return m_data[bidx.localIndex()];
}

/*--------------------------------------------------------------------*/
//  Constant index with a LayoutIterator
/** \param[in]  a_lit   Layout iterator
 *  \return             Element at the index of the data iterator
 *//*-----------------------------------------------------------------*/

template <typename T>
inline const T&
LevelData<T>::operator[](const LayoutIterator& a_lit) const
{
  // Make sure we are indexing with an iterator built on this DBL
  CH_assert(a_lit.tag() == tag());
  const BoxIndex& bidx = (*a_lit);
  CH_assert(bidx.localIndex() < m_disjointBoxLayout.localSize());
  CH_assert(bidx.globalIndex() >= m_disjointBoxLayout.localIdxBegin() &&
            bidx.globalIndex() <  m_disjointBoxLayout.localIdxEnd());
  return m_data[bidx.localIndex()];
}

/*--------------------------------------------------------------------*/
//  Index with a BoxIndex
/** This is less safe and is discouraged because it is not possible to
 *  check that the Iterator matches the Box layout.
 *  \param[in]  a_bidx  BoxIndex
 *  \return             Element at the index
 *//*-----------------------------------------------------------------*/

template <typename T>
inline T&
LevelData<T>::operator[](const BoxIndex& a_bidx)
{
  CH_assert(a_bidx.localIndex() < m_disjointBoxLayout.localSize());
  CH_assert(a_bidx.globalIndex() >= m_disjointBoxLayout.localIdxBegin() &&
            a_bidx.globalIndex() <  m_disjointBoxLayout.localIdxEnd());
  return m_data[a_bidx.localIndex()];
}

/*--------------------------------------------------------------------*/
//  Constant index with a BoxIndex
/** This is less safe and is discouraged because it is not possible to
 *  check that the Iterator matches the Box layout.
 *  \param[in]  a_bidx  BoxIndex
 *  \return             Element at the index
  *//*-----------------------------------------------------------------*/

template <typename T>
inline const T&
LevelData<T>::operator[](const BoxIndex& a_bidx) const
{
  CH_assert(a_bidx.localIndex() < m_disjointBoxLayout.localSize());
  CH_assert(a_bidx.globalIndex() >= m_disjointBoxLayout.localIdxBegin() &&
            a_bidx.globalIndex() <  m_disjointBoxLayout.localIdxEnd());
  return m_data[a_bidx.localIndex()];
}

/*--------------------------------------------------------------------*/
//  Assign a constant to all components
/** \param[in] a_val    Value to assign
 *//*-----------------------------------------------------------------*/

template <typename T>
inline void
LevelData<T>::setVal(const typename T::value_type& a_val)
{
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      this->operator[](dit).setVal(a_val);
    }
}

/*--------------------------------------------------------------------*/
//  Assign a constant to a single component
/** \param[in] a_icomp  Component index
 *  \param[in] a_val    Value to assign
 *//*-----------------------------------------------------------------*/

template <typename T>
inline void
LevelData<T>::setVal(const int a_icomp, const typename T::value_type& a_val)
{
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      this->operator[](dit).setVal(a_icomp, a_val);
    }
}

/*--------------------------------------------------------------------*/
//  Unique identifying tag (from the DBL)
/*--------------------------------------------------------------------*/

template <typename T>
inline size_t
LevelData<T>::tag() const
{
  return m_disjointBoxLayout.tag();
}

/*--------------------------------------------------------------------*/
//  Constant access elements with a linear index
/** FOR INTERNAL USE AND TESTING ONLY
 *  \param[in]  a_idx   Linear index to the data
 *  \return             Element at the index
 *//*-----------------------------------------------------------------*/

template <typename T>
inline const T&
LevelData<T>::getLinear(const int a_idx) const
{
  CH_assert(a_idx >= 0 && a_idx < size());
  return m_data[a_idx];
}

/*--------------------------------------------------------------------*/
//  Number of boxes
/*--------------------------------------------------------------------*/

template <typename T>
inline int
LevelData<T>::size() const
{
  return m_disjointBoxLayout.localSize();
}

/*--------------------------------------------------------------------*/
//  Number of components
/*--------------------------------------------------------------------*/

template <typename T>
inline int
LevelData<T>::ncomp() const
{
  return m_ncomp;
}

/*--------------------------------------------------------------------*/
//  Number of ghosts
/*--------------------------------------------------------------------*/

template <typename T>
inline int
LevelData<T>::nghost() const
{
  return m_nghost;
}

/*--------------------------------------------------------------------*/
//  The layout of boxes
/*--------------------------------------------------------------------*/

template <typename T>
inline const DisjointBoxLayout&
LevelData<T>::disjointBoxLayout() const
{
  return m_disjointBoxLayout;
}

/*--------------------------------------------------------------------*/
//  Exchange to fill ghost cells
/** \param[in]  a_copier
 *                      A copier that caches data motion patterns
 *//*-----------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::exchange(Copier& a_copier)
{
  if (m_nghost > 0)
    {
      const int startComp = a_copier.startComp();
      const int numComp   = a_copier.numComp();
#ifdef USE_MPI
      const int endComp   = a_copier.endComp();
      MPI_Request* requests = a_copier.requests();
      int idxReq = 0;
#endif

      int midx;
      const int nmitem = a_copier.numMotionItem();
      for (midx = 0; midx < nmitem; ++midx)
        {
          Motion2Way& motion = a_copier[midx];
#ifdef USE_MPI
          if (motion.isLocal())
#endif
            {
              CH_assert(motion.isLocal());
              
              m_data[motion.bidxRecv().localIndex()].copy(motion.regionRecv(),startComp,
                                                          m_data[motion.bidxSend().localIndex()],
                                                          motion.regionSend(),startComp, numComp,
                                                          motion.compRecvFlags());
             
            }
#ifdef USE_MPI
          else
            {
              this->operator[](motion.m_bidxLocal).linearOut(
                motion.m_sendBuffer.get(),
                motion.m_regionSend,
                startComp,
                endComp);
              motion.postMessages(a_copier.bytesPerCell(),
                                  requests + idxReq,
                                  requests + idxReq + 1);
              idxReq += 2;
            }
#endif
        }
      

#ifdef USE_MPI
      if (DisjointBoxLayout::numProc() > 1)
        {
          const int nReq = a_copier.numRequest();
#ifndef USE_MPIWAITALL
          
          for (int iReq = 0; iReq != nReq; ++iReq)
            {
              int ridx;  
              int mpierr = MPI_Waitany(nReq, requests, &ridx, MPI_STATUS_IGNORE);
              if (mpierr)
                {
                  std::cout << "Error waiting on one message on process "
                        << DisjointBoxLayout::procID() << std::endl;
                  abort();
                }
              CH_assert(mpierr == 0);
                           
              if (ridx & 1)  
                {
                               
                  Motion2Way& motion =a_copier[a_copier.motionItemIndex(ridx)];
                  if (!motion.isLocal())
                    {
                     
                      this->operator[](motion.m_bidxLocal).linearIn(
                        motion.m_recvBuffer.get(),
                        motion.m_regionRecv,
                        startComp,
                        endComp);
                    }
                 
                }
            }
#else
          
         
          int mpierr = MPI_Waitall(nReq, requests, MPI_STATUSES_IGNORE);
          if (mpierr)
            {
              std::cout << "Error waiting for all messages on process "
                        << DisjointBoxLayout::procID() << std::endl;
              abort();
            }
         
          CH_assert(mpierr == 0);
          for (int midx = 0; midx != nmitem; ++midx)
            {
              Motion2Way& motion = a_copier[midx];
              if (!motion.isLocal())
                {
                  this->operator[](motion.m_bidxLocal).linearIn(
                    motion.m_recvBuffer.get(),
                    motion.m_regionRecv,
                    startComp,
                    endComp);
                 
                }
            }
#endif
        }
#endif
    }
}



/*--------------------------------------------------------------------*/
//  Begin exchange to fill ghost cells
/** Use with exchangeEnd to overlap computation with communication
 *  \param[in]  a_copier
 *                      A copier that caches data motion patterns
 *//*-----------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::exchangeBegin(Copier& a_copier)
{
  //**FIXME
	if(m_nghost > 0)
	{
		const int startComp = a_copier.startComp();
       		const int numComp   = a_copier.numComp();
	#ifdef USE_MPI
        	const int endComp   = a_copier.endComp();
    		MPI_Request* requests = a_copier.requests();
		int idxReq = 0;
 	#endif
       		int midx;
       		const int nmitem = a_copier.numMotionItem();
       		for (midx = 0; midx < nmitem; ++midx)
        	{
              	Motion2Way& motion = a_copier[midx];
 	#ifdef USE_MPI
         	if (motion.isLocal())
 	#endif
            	{
                	(m_data[(motion.bidxRecv()).localIndex()]).copy(motion.regionSend(),m_data[(motion.bidxSend()).localIndex()]);
             	}
 	#ifdef USE_MPI
          	else //goes with if(motion.isLocal())
             	{
               		this->operator[](motion.m_bidxLocal).linearOut(
                 	motion.m_sendBuffer.get(),
                 	motion.m_regionSend,
                 	startComp,
                 	endComp);
              	 	motion.postMessages(a_copier.bytesPerCell(),
                        	requests + idxReq,
                                requests + idxReq + 1);
               		idxReq += 2;
             	}
 	#endif
         	}//end for				
	}//end if(m_nghost) 
}

/*--------------------------------------------------------------------*/
//  End exchange to fill ghost cells/** Use with exchangeBegin to overlap computation with communication
/** Use with exchangeBegin to overlap computation with communication
 *  \param[in]  a_copier
 *                      A copier that caches data motion patterns
 *//*-----------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::exchangeEnd(Copier& a_copier)
{
  //**FIXME
	#ifdef USE_MPI
	if (DisjointBoxLayout::numProc() > 1)
        {
        	const int nReq = a_copier.numRequest();
 	#ifndef USE_MPIWAITALL
        // Wait for first message, unpack as soon as received
        	for (int iReq = 0; iReq != nReq; ++iReq)
             	{
               		int ridx;  // Request index
               //**FIXME
               		int mpierr = MPI_Waitany(nReq,a_copier.requests(),&ridx,MPI_STATUS_IGNORE);
               		if (ridx & 1)  // This is a receive (has odd request index)
                 	{
                   //**FIXME
                   		Motion2Way &motion = a_copier[ridx/2];
                  			(m_data[(motion.bidxRecv()).localIndex()]). //this is the BaseFab
                     				linearIn((motion.m_recvBuffer).get(),motion.m_regionRecv,a_copier.startComp(),a_copier.endComp(),motion.m_compRecvFlags);

			}
		}//end for
 	#else
        // Full barrier wait
        int mpierr = MPI_Waitall(nReq, requests, MPI_STATUSES_IGNORE);
        if (mpierr)
   	{
        	<< DisjointBoxLayout::procID() << std::endl;
               	abort();
        }
        CH_assert(mpierr == 0);
        for (int midx = 0; midx != nmitem; ++midx)
        {
        	Motion2Way& motion = a_copier[midx];
        	if (!motion.isLocal())
        	{
        	this->operator[](motion.m_bidxLocal).linearIn(
                     motion.m_recvBuffer.get(),
                     motion.m_regionRecv,
                     a_copier.startComp,
                     a_copier.endComp);
        	}
        }
 #endif
   } 
 #endif

}//end exchangeEnd

#ifndef NO_CGNS
/*--------------------------------------------------------------------*/
//  Write CGNS solution data to a file (specialized for BaseFab<Real>)
//  Not available in general.  See LevelData.cpp for specialization.
/*--------------------------------------------------------------------*/

template <typename T>
int
LevelData<T>::writeCGNSSolData(const int                a_indexFile,
                               const int                a_indexBase,
                               const int                a_indexZoneOffset,
                               const char *const *const a_varNames) const
{
  assert(false);
  return 0;
}

// Specialized for BaseFab<Real>
template<>
int
LevelData<BaseFab<Real> >::writeCGNSSolData(
  const int                a_indexFile,
  const int                a_indexBase,
  const int                a_indexZoneOffset,
  const char *const *const a_varNames) const;
#endif  /* CGNS */

#ifdef USE_GPU
/*--------------------------------------------------------------------*/
//  Copy all elements to device
/*--------------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::copyToDevice() const
{
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      this->operator[](dit).copyToDevice();
    }
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy all elements to device
/** \param[in]  a_stream
 *                      Stream index (defaults to default stream)
 *//*-----------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::copyToDeviceAsync(cudaStream_t a_stream) const
{
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      this->operator[](dit).copyToDeviceAsync(a_stream);
    }
}

/*--------------------------------------------------------------------*/
//  Copy all elements to host
/*--------------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::copyToHost()
{
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      this->operator[](dit).copyToHost();
    }
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy all elements to host
/** \param[in]  a_stream
 *                      Stream index (defaults to default stream)
 *//*-----------------------------------------------------------------*/

template <typename T>
void
LevelData<T>::copyToHostAsync(cudaStream_t a_stream)
{
  for (DataIterator dit(m_disjointBoxLayout); dit.ok(); ++dit)
    {
      this->operator[](dit).copyToHostAsync(a_stream);
    }
}
#endif  /* CUDA */

#endif  /* ! defined _LEVELDATA_H_ */

