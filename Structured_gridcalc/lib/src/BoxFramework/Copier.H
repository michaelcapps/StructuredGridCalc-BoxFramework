
#ifndef _COPIER_H_
#define _COPIER_H_


/******************************************************************************/
/**
 * \file Copier.H
 *
 * \brief Data motion caches and patterns for copying information between boxes
 *
 *//*+*************************************************************************/

#include <cstdlib>
#include <memory>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "BoxIndex.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"

//--Forward declarations

template <typename T>
class LevelData;


/*******************************************************************************
 */
///  2-way motion involves the sending and reception of a message
/**
 ******************************************************************************/

class Motion2Way
{

//--Deleter type for buffers

  struct DelBuffer
  {
    void operator()(void* addr)
      {
        free(addr);
      }
  };

//--Friends

  template <typename T>
  friend class LevelData;


/*====================================================================*
 * Public constructors and destructors
 *====================================================================*/

public:

  /// Default constructor
  Motion2Way();

  /// Constructor
  Motion2Way(const int               a_bytesPerCell,
             const DisjointBoxLayout a_disjointBoxLayout,
             const BoxIndex&         a_bidxLocal,
             const BoxIndex&         a_bidxRemote,
             const Box&              a_regionRecv,
             const Box&              a_regionSend,
             const Box&              a_regionSendRemote,
             const IntVect&          a_sendDir);

  // Use synthesized copy, move, copy assignment, move assignment, and
  // destructor.


/*====================================================================*
 * Members functions
 *====================================================================*/

public:

  /// Are operations local (both boxes on same process)?
  bool isLocal() const;

  /// Generate a unique tag (based on sending process) for this motion item
  int uniqueTag(const BoxIndex& a_bidxSend, const IntVect& a_sendDir) const;

#ifdef USE_MPI
  /// Post messages from this motion item
  void postMessages(const int          a_bytesPerCell,
                    MPI_Request *const a_sendRequest,
                    MPI_Request *const a_recvRequest) const;
#endif

//--Access for local operations

  /// Index to the local recieving box
  const BoxIndex& bidxRecv() const { return m_bidxLocal; }

  /// Index to the local sending box
  const BoxIndex& bidxSend() const { return m_bidxRemote; }

  /// Receiving box region
  const Box& regionRecv() const { return m_regionRecv; }

  /// Sending box region
  const Box& regionSend() const { return m_regionSendRemote; }

  /// Receive direction
  IntVect recvDir() const { return (-m_sendDir); }

  /// Send direction
  const IntVect& sendDir() const { return m_sendDir; }

  /// Component receive flags
  unsigned compRecvFlags() const { return m_compRecvFlags; }

  /// Component send flags
  unsigned compSendFlags() const { return m_compSendFlags; }

  /// Modify component receive flags
  void setCompRecvFlags(const unsigned a_flags);

  /// Modify component send flags
  void setCompSendFlags(const unsigned a_flags);


/*====================================================================*
 * Data members
 *====================================================================*/

protected:

  BoxIndex m_bidxLocal;               ///< Index to the local box
  BoxIndex m_bidxRemote;              ///< Index to the remote or neighbor box
  Box m_regionRecv;                   ///< Region to receive a message into
                                      ///< this box
  Box m_regionSend;                   ///< Region to send a message from this
                                      ///< box.  This is only required for MPI.
  Box m_regionSendRemote;             ///< The region to send a message from the
                                      ///< remote box.  This is the same as
                                      ///< regionRecv except for periodic
                                      ///< boundaries
  int m_localProcID;                  ///< ID of the local process
  int m_remoteProcID;                 ///< ID of the remote process
  int m_tagSend;                      ///< Unique tag for sending messages
  int m_tagRecv;                      ///< Unique tag for receiving messages
  IntVect m_sendDir;                  ///< Direction to send information
  unsigned m_compRecvFlags;           ///< Bit flags describing components
                                      ///< to transfer in receive direction
  unsigned m_compSendFlags;           ///< Bit flags describing components
                                      ///< to transfer in send direction
  std::unique_ptr<void, DelBuffer> m_recvBuffer;
                                      ///< Buffer for receiving messages
  std::unique_ptr<void, DelBuffer> m_sendBuffer;
                                      ///< Buffer for sending messages
};


/*******************************************************************************
 */
///  Class to cache all data motion in a layout of boxes
/**
 ******************************************************************************/

class Copier
{

/*====================================================================*
 * Public constructors and destructors
 *====================================================================*/

public:

  /// Default constructor
  Copier();

  /// Copy constructor not allowed
  Copier(const Copier&) = delete;

  /// Move constructor
  Copier(Copier&&) noexcept = default;

  /// Assignment constructor not allowed
  Copier& operator=(const Copier&) = delete;

  /// Move assignment constructor
  Copier& operator=(Copier&&) noexcept = default;

  /// Destructor
  ~Copier() = default;

  /// Weak construction of an exchange copier for all components of a LevelData
  template <typename S>
  void defineExchangeLD(const LevelData<S>&      a_lvlData,
                        const unsigned           a_periodic = 0u,
                        const unsigned           a_trim = 0u);

  /// Weak construction of an exchange copier from a DBL
  template <typename T>
  void defineExchangeDBL(const DisjointBoxLayout& a_disjointBoxLayout,
                         const int                a_numGhost,
                         const int                a_startComp,
                         const int                a_numComp,
                         const unsigned           a_periodic = 0u,
                         const unsigned           a_trim = 0u);


/*====================================================================*
 * Members functions
 *====================================================================*/

public:

  /// Unique tag identifying the DisjointBoxLayout this Copier is valid for
  size_t tag() const;

  /// Number of bytes per cell to copy (includes all components)
  int bytesPerCell() const;

  /// Number of motion items
  int numMotionItem() const;

  /// Const access to a motion item
  const Motion2Way& operator[](const int a_idx) const;

  /// Access to a motion item
  Motion2Way& operator[](const int a_idx);

  /// Start of component range
  int startComp() const;

  /// One past end of component range
  int endComp() const;

  /// Number of components to copy
  int numComp() const;

  /// Calculate a binomial coefficient
  static int binomial(const int n, int k);

#ifdef USE_MPI
  /// Number of MPI requests
  int numRequest() const;

  /// The MPI requests
  MPI_Request* requests();

  /// Get the motion item index for a specific request index
  int motionItemIndex(const int a_idxReq) const;
#endif


/*====================================================================*
 * Data members
 *====================================================================*/

protected:

  size_t m_tag;                       ///< A unique tag identifying the
                                      ///< DisjointBoxLayout for which this
                                      ///< Copier was built
  int m_bytesPerCell;                 ///< Number of bytes of data per cell
                                      ///< in a BaseFab (for all components)
  int m_startComp;                    ///< Start for a range of components
  int m_endComp;                      ///< One past last component in range
  std::vector<Motion2Way> m_motionItem;
                                      ///< An array of items describing 2-way
                                      ///< exchanges of data between boxes
#ifdef USE_MPI
  std::vector<MPI_Request> m_mpiRequest;
                                      ///< MPI handles for non-blocking calls
  std::vector<int> m_midxForReq;      ///< Motion item index for a given
                                      ///< request (only required for Waitany)
#endif
  int m_numReq;                       ///< Number of messages
};


/*******************************************************************************
 *
 * Class Motion2Way: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/

inline
Motion2Way::Motion2Way()
  :
  m_bidxLocal(),
  m_bidxRemote(),
  m_regionRecv(),
  m_regionSend(),
  m_regionSendRemote(),
  m_localProcID(-1),
  m_remoteProcID(-1),
  m_tagSend(-1),
  m_tagRecv(-1),
  m_compRecvFlags(std::numeric_limits<unsigned>::max()),
  m_compSendFlags(std::numeric_limits<unsigned>::max()),
  m_recvBuffer(nullptr, DelBuffer()),
  m_sendBuffer(nullptr, DelBuffer())
{ }

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_bytesPerCell
 *                      Bytes of date per cell for exchange
 *  \param[in]  a_disjointBoxLayout
 *                      Layout of boxes
 *  \param[in]  a_bidxLocal
 *                      BoxIndex of local box
 *  \param[in]  a_bidxRemote
 *                      BoxIndex of neighbor (remote) box
 *  \param[in]  a_regionRecv
 *                      Region to receive into in local box
 *  \param[in]  a_regionSend
 *                      Region to send from local box (only required
 *                      for MPI)
 *  \param[in]  a_regionSendRemote
 *                      Region to send from remote box (only required
 *                      for local copies)
 *  \param[in]  a_sendDir
 *                      Direction to send information
 *//*-----------------------------------------------------------------*/

inline
Motion2Way::Motion2Way(const int               a_bytesPerCell,
                       const DisjointBoxLayout a_disjointBoxLayout,
                       const BoxIndex&         a_bidxLocal,
                       const BoxIndex&         a_bidxRemote,
                       const Box&              a_regionRecv,
                       const Box&              a_regionSend,
                       const Box&              a_regionSendRemote,
                       const IntVect&          a_sendDir)
  :
  m_bidxLocal(a_bidxLocal),
  m_bidxRemote(a_bidxRemote),
  m_regionRecv(a_regionRecv),
  m_regionSend(a_regionSend),
  m_regionSendRemote(a_regionSendRemote),
  m_localProcID(a_disjointBoxLayout.proc(a_bidxLocal)),
  m_remoteProcID(a_disjointBoxLayout.proc(a_bidxRemote)),
  m_tagSend(uniqueTag(a_bidxLocal, a_sendDir)),
  m_tagRecv(uniqueTag(a_bidxRemote, -a_sendDir)),
  m_sendDir(a_sendDir),
  m_compRecvFlags(std::numeric_limits<unsigned>::max()),
  m_compSendFlags(std::numeric_limits<unsigned>::max()),
  m_recvBuffer(nullptr, DelBuffer()),
  m_sendBuffer(nullptr, DelBuffer())
{
  if (!isLocal())
    {
      m_recvBuffer.reset(std::malloc(a_bytesPerCell*m_regionRecv.size()));
      m_sendBuffer.reset(std::malloc(a_bytesPerCell*m_regionSend.size()));
    }
}

/*--------------------------------------------------------------------*/
//  Are operations local (both boxes on same process)?
/*--------------------------------------------------------------------*/

inline bool
Motion2Way::isLocal() const
{
  return (m_localProcID == m_remoteProcID);
}

/*--------------------------------------------------------------------*/
//  Generate a unique tag (based on sending process) for this motion
//  item
/** \note
 *  <ul>
 *    <li> This could be done better.  It is currently a completely
 *         unique index for the entire program rather than just a
 *         unique index between two processes.  Doing better generally
 *         requires the local index of the sending process, which we
 *         don't have for receives.
 *  </ul>
 *//*-----------------------------------------------------------------*/

inline int
Motion2Way::uniqueTag(const BoxIndex& a_bidxSend,
                      const IntVect&  a_sendDir) const
{
  return 27*a_bidxSend.globalIndex() + D_TERM(    (a_sendDir[0]+1),
                                              + 3*(a_sendDir[1]+1),
                                              + 9*(a_sendDir[2]+1));
}

/*--------------------------------------------------------------------*/
//  Post messages from this motion item
/*--------------------------------------------------------------------*/

#ifdef USE_MPI
inline void
Motion2Way::postMessages(const int          a_bytesPerCell,
                         MPI_Request *const a_sendRequest,
                         MPI_Request *const a_recvRequest) const
{ //FIXMES WERE HERE
  CH_assert(m_sendBuffer != NULL);
  MPI_Isend(m_sendBuffer.get(),a_bytesPerCell*m_regionSend.size(), MPI_BYTE, m_remoteProcID, m_tagSend, MPI_COMM_WORLD, a_sendRequest);
  CH_assert(m_recvBuffer != NULL);
  MPI_Irecv(m_recvBuffer.get(),a_bytesPerCell*m_regionRecv.size(), MPI_BYTE, m_remoteProcID, m_tagRecv, MPI_COMM_WORLD, a_recvRequest);
}
#endif

/*--------------------------------------------------------------------*/
//  Modify component receive flags
/*--------------------------------------------------------------------*/

inline void
Motion2Way::setCompRecvFlags(const unsigned a_flags)
{
  m_compRecvFlags = a_flags;
}

/*--------------------------------------------------------------------*/
//  Modify component send flags
/*--------------------------------------------------------------------*/

inline void
Motion2Way::setCompSendFlags(const unsigned a_flags)
{
  m_compSendFlags = a_flags;
}


/*******************************************************************************
 *
 * Class Copier: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/

inline
Copier::Copier()
  :
  m_tag(0),
  m_bytesPerCell(-1),
  m_motionItem(),
#ifdef USE_MPI
  m_mpiRequest(),
  m_midxForReq(),
#endif
  m_numReq(0)
{
}

/*--------------------------------------------------------------------*/
//  Weak construction of an exchange copier for all components of a
//  LevelData
/** This does not mean the copier can only be used with this
 *  LevelData, any similar one built on the same DisjointBoxLayout
 *  would work.
 *  \tparam T           Type of data in a LevelData (BaseFab?)
 *  \param[in]  a_lvlData
 *                      LevelData to build the copier for
 *  \param[in]  a_periodic
 *                      Which directions are periodic.  To select X
 *                      and Y, you would pass the value
 *                      PeriodicX | PeriodicY as an argument.  Default
 *                      is no periodic.
 *  \param[in]  a_trim  Trimmed sections are not included as
 *                      neighbors.  (0,0,0) is always trimmed since
 *                      that is the base IV.  To trim edges and
 *                      corners, you would pass the value
 *                      TrimEdge | TrimCorner as an argument.  Default
 *                      is no trimming (aside from (0,0,0)).
 *//*-----------------------------------------------------------------*/

template <typename S>
inline void
Copier::defineExchangeLD(const LevelData<S>&      a_lvlData,
                         const unsigned           a_periodic,
                         const unsigned           a_trim)
{
  typedef typename S::value_type T;
  defineExchangeDBL<T>(a_lvlData.disjointBoxLayout(),
                       a_lvlData.nghost(),
                       0,
                       a_lvlData.ncomp(),
                       a_periodic,
                       a_trim);
}

/*--------------------------------------------------------------------*/
//  Weak construction of an exchange copier from a DBL
/** \tparam T           Type of data in a cell
 *  \param[in]  a_disjointBoxLayout
 *                      The disjoint box layout to build the copier
 *                      for
 *  \param[in]  a_numGhost
 *                      Number of ghosts to copy
 *  \param[in]  a_startComp
 *                      Start of range of components to copy
 *  \param[in]  a_numComp
 *                      Total number of components to copy
 *  \param[in]  a_periodic
 *                      Which directions are periodic.  To select X
 *                      and Y, you would pass the value
 *                      PeriodicX | PeriodicY as an argument.  Default
 *                      is no periodic.
 *  \param[in]  a_trim  Trimmed sections are not included as
 *                      neighbors.  (0,0,0) is always trimmed since
 *                      that is the base IV.  To trim edges and
 *                      corners, you would pass the value
 *                      TrimEdge | TrimCorner as an argument.  Default
 *                      is no trimming (aside from (0,0,0)).
 *//*-----------------------------------------------------------------*/

template <typename T>
inline void
Copier::defineExchangeDBL(const DisjointBoxLayout& a_disjointBoxLayout,
                          const int                a_numGhost,
                          const int                a_startComp,
                          const int                a_numComp,
                          const unsigned           a_periodic,
                          const unsigned           a_trim)
{
  CH_assert(a_startComp >= 0);
  CH_assert(a_numComp > 0);
  m_tag = a_disjointBoxLayout.tag();
  m_bytesPerCell = sizeof(T)*a_numComp;
  m_startComp = a_startComp;
  m_endComp = a_startComp + a_numComp;
  m_motionItem.clear();
#ifdef USE_MPI
  m_mpiRequest.clear();
  m_midxForReq.clear();
#endif
  m_numReq = 0;
  if (a_numGhost > 0)
    {
      Box periodicTestDomain = a_disjointBoxLayout.problemDomain();
      for (int dir = 0; dir != g_SpaceDim; ++dir)
        {
          if (a_periodic & (1<<dir))
            {
              periodicTestDomain.grow(-1, dir);
            }
        }
      // Try to predict the number of motion items
      int predNumMotionItem = 0;
      for (int m = 0; m < g_SpaceDim; ++m)
        {
          // Get the number of m-dimension hypercubes on the surfaces of the box
          if (!(a_trim & (1<<(m+1))))
            {
              predNumMotionItem += (1<<(g_SpaceDim-m))*binomial(g_SpaceDim, m);
            }
        }
      predNumMotionItem *= a_disjointBoxLayout.localSize();
      m_motionItem.reserve(predNumMotionItem);

//--Iterate over boxes on the process

      for (DataIterator dit(a_disjointBoxLayout); dit.ok(); ++dit)
        {
          // The fab and box at the current position of the DataIterator
          int localProcID;
          const Box localBox = a_disjointBoxLayout.box(dit, localProcID);
//--Interior neighbors

          for(NeighborIterator nbrit(dit, a_trim); nbrit.ok(); ++nbrit)
            {
              Box remoteBox = a_disjointBoxLayout[nbrit];
	      Box regionRecv(localBox);
              regionRecv.grow(a_numGhost);
	      regionRecv &= remoteBox;
#ifdef USE_MPI
              remoteBox.grow(a_numGhost);
              Box regionSend(localBox);
              regionSend &= remoteBox;
#else
              Box regionSend;
#endif
              m_motionItem.emplace_back(m_bytesPerCell,
                                        a_disjointBoxLayout,
                                        *dit,
                                        *nbrit,
                                        regionRecv,
                                        regionSend,
                                        regionRecv,
                                        nbrit.nbrDir());
		if (!m_motionItem.back().isLocal())
                {
			  m_numReq += 2;
                }
            }

//--Periodic neighbors

          if (!periodicTestDomain.contains(localBox))
            {
              for (PeriodicIterator perit(dit, a_trim, a_periodic); perit.ok();
                   ++perit)
                {
                  Box remoteBox = a_disjointBoxLayout[perit];
                  // We need to shift the remoteBox (which is inside the domain)
                  // to its periodic location outside the domain.
                  const IntVect& shiftDir = perit.nbrDir();
                  IntVect shiftBy =
                    localBox.loVect() - remoteBox.loVect()  // Shift to local
                    + shiftDir*localBox.dimensions();       // Shift past local
                  remoteBox.shift(shiftBy);
                  Box regionRecv(localBox);
                  regionRecv.grow(a_numGhost);
                  regionRecv &= remoteBox;
#ifdef USE_MPI
                  remoteBox.grow(a_numGhost);
                  Box regionSend(localBox);
                  regionSend &= remoteBox;
#else
                  Box regionSend;
#endif
                  Box regionSendRemote(regionRecv);
                  regionSendRemote.shift(-shiftBy);
                  m_motionItem.emplace_back(m_bytesPerCell,
                                            a_disjointBoxLayout,
                                            *dit,
                                            *perit,
                                            regionRecv,
                                            regionSend,
                                            regionSendRemote,
                                            perit.nbrDir());
                  if (!m_motionItem.back().isLocal())
                    {
                      m_numReq += 2;
                    }
                }
            }
        }

      // Allocate MPI constructs if required
#ifdef USE_MPI
      m_mpiRequest.resize(m_numReq);
      m_midxForReq.resize(m_numReq/2);
      int cRecvReq = 0;
      const int nMotionItem = numMotionItem();
      for (int i = 0; i != nMotionItem; ++i)
        {
          if (!m_motionItem[i].isLocal())
            {
              m_midxForReq[cRecvReq++] = i;
            }
        }
#endif
    }
}

/*--------------------------------------------------------------------*/
//  Unique tag identifying the DisjointBoxLayout this Copier is valid
//  for
/*--------------------------------------------------------------------*/

inline size_t
Copier::tag() const
{
  return m_tag;
}

/*--------------------------------------------------------------------*/
//  Number of bytes per cell to copy (includes all components)
/*--------------------------------------------------------------------*/

inline int
Copier::bytesPerCell() const
{
  return m_bytesPerCell;
}

/*--------------------------------------------------------------------*/
//  Number of motion items
/*--------------------------------------------------------------------*/

inline int
Copier::numMotionItem() const
{
  return m_motionItem.size();
}

/*--------------------------------------------------------------------*/
//  Const access to a motion item
/*--------------------------------------------------------------------*/

inline const Motion2Way&
Copier::operator[](const int a_idx) const
{
  CH_assert(a_idx >= 0 && a_idx < numMotionItem());
  return m_motionItem[a_idx];
}

/*--------------------------------------------------------------------*/
//  Access to a motion item
/*--------------------------------------------------------------------*/

inline Motion2Way&
Copier::operator[](const int a_idx)
{
  CH_assert(a_idx >= 0 && a_idx < numMotionItem());
  return m_motionItem[a_idx];
}

/*--------------------------------------------------------------------*/
//  Start of component range
/*--------------------------------------------------------------------*/

inline int
Copier::startComp() const
{
  return m_startComp;
}

/*--------------------------------------------------------------------*/
//  One past end of component range
/*--------------------------------------------------------------------*/

inline int
Copier::endComp() const
{
  return m_endComp;
}

/*--------------------------------------------------------------------*/
//  Number of components to copy
/*--------------------------------------------------------------------*/

inline int
Copier::numComp() const
{
  return m_endComp - m_startComp;
}

/*--------------------------------------------------------------------*/
/// Calculate a binomial coefficient
/** \param[in]          n
 *  \param[in]          k
 *  \return             Binomial coefficient
 *//*-----------------------------------------------------------------*/

inline int
Copier::binomial(const int n, int k)
{
  CH_assert(k >= 0);
  CH_assert(n >= k);
  if (k == 0) return 1;
  int cnum = n - k + 1;
  int cden = k;
  while (--k)
    {
      cnum *= n - k + 1;
      cden *= k;
    }
  return cnum/cden;
}

#ifdef USE_MPI
/*--------------------------------------------------------------------*/
//  Number of MPI requests
/*--------------------------------------------------------------------*/

inline int
Copier::numRequest() const
{
  return m_numReq;
}

/*--------------------------------------------------------------------*/
//  The MPI requests
/*--------------------------------------------------------------------*/

inline MPI_Request*
Copier::requests()
{
  return m_mpiRequest.data();
}

/*--------------------------------------------------------------------*/
//  Get the motion item index for a specific request index
/*--------------------------------------------------------------------*/

inline int
Copier::motionItemIndex(const int a_idxReq) const
{
  return m_midxForReq[a_idxReq/2];  // Since there are 2 requests per item
}
#endif

#endif  /* ! defined _COPIER_H_ */
