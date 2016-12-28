
#ifndef _WAVEPATCH_H_
#define _WAVEPATCH_H_


/******************************************************************************/
/**
 * \file wavePatch.H
 *
 * \brief Implements wave equation on a patch
 *
 *//*+*************************************************************************/

#include "Parameters.H"
#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "Stopwatch.H"

#ifdef USE_GPU
#include "CudaSupport.H"
#endif


/*******************************************************************************
 */
///  Solution on a patch
/**
 ******************************************************************************/

class WavePatch
{

/*====================================================================*
 * Types
 *====================================================================*/

public:

  using PatchSolData = BaseFab<Real>;
  using LevelSolData = LevelData<BaseFab<Real> >;


/*====================================================================*
 * Public constructors and destructors
 *====================================================================*/

public:

  /// Constructor
  WavePatch(const Box&        a_domain,
            const IntVect&    a_maxBoxSize,
            const char* const a_basePlotName,
            const Real        a_c,
            const Real        a_dx,
            const Real        a_cfl);

  /// Copy constructor not permitted
  WavePatch(const WavePatch&) = delete;

  /// Assignment constructor not permitted
  WavePatch& operator=(const WavePatch&) = delete;
  
  /// Destructor
  ~WavePatch();


/*====================================================================*
 * Members functions
 *====================================================================*/

public:

  /// Set initial data to rest velocity
  void initialData();

  /// Advance one time step
  void advance();

#ifdef USE_GPU
  /// Advance a group of time steps using GPU
  void advanceIterGroup(const int a_numIter,
                        cudaEvent_t a_cuEvent_iterGroupStart,
                        cudaEvent_t a_cuEvent_iterGroupEnd);
#endif

  /// Write the plot file
  int writePlotFile(const int a_idxStep, const int a_iteration) const;

  /// Access u (given time index)
  PatchSolData& u(const int a_idxStep);

  /// Const access u (given time index)
  const PatchSolData& u(const int a_idxStep) const;

  /// Access u (time level n)
  PatchSolData& un();

  /// Const access u (time level n)
  const PatchSolData& un() const;

  /// Access u update (time level n+1)
  PatchSolData& unp1();

  /// Const access u update (time level n+1)
  const PatchSolData& unp1() const;

  /// Access older u (time level n-1)
  PatchSolData& unm1();

  /// Const access older u (time level n-1)
  const PatchSolData& unm1() const;

  /// Advance indices (unp1->un, un->unm1)
  void advanceStepIndex();

  /// Current step index
  int currentStepIndex() const;

  /// Old step index
  int oldStepIndex() const;

  /// Current time
  Real time() const;

  /// Current iteration
  int iteration() const;

#ifdef USE_GPU
  /// Copy data to host
  void copyToHostAsync(const int a_idxStep);

  /// Copy data to device
  void copyToDeviceAsync(const int a_idxStep);
#endif


/*====================================================================*
 * Data members
 *====================================================================*/

protected:
  DisjointBoxLayout m_boxes;          ///< Layout of boxes in the domain
  LevelSolData m_u[3];                ///< Scalar displacement
  Box m_domain;                       ///< Problem domain
  const char* m_basePlotName;         ///< Name for plot files
  Real m_c;                           ///< Wave speed
  Real m_dx;                          ///< Mesh spacing
  Real m_dt;                          ///< Time step
  Real m_time;                        ///< Current time
  int m_iteration;                    ///< Current iteration
  int m_idxStep;                      ///< Index of \f$u^n\f$
  int m_idxStepUpdate;                ///< Index of \f$u^{n+1}\f$
  int m_idxStepOld;                   ///< Index of \f$u^{n-1}\f$
  BoxIndex m_bidx;                    ///< Since we only have a single box,
                                      ///< store the index to it.
public:
  Stopwatch<> m_timerAdvance;         ///< Timer for advance function
  mutable Stopwatch<> m_timerWrite;   ///< Timer for plot writing

#ifdef USE_GPU
protected:
  AccelPointer m_cudaFab_device;      ///< Pointer to CudaFabs on device
  AccelPointer m_workBoxesRHS_device; ///< Pointer to work boxes for computing
                                      ///< RHS on device
  AccelPointer m_workBoxInfoBC_device;
                                      ///< Pointer to work box information for
                                      ///< updated BC on device
  int m_numBlkRHS;                    ///< Number of blocks when computing RHS
  int m_numBlkBC;                     ///< Number of blocks when updating BC
#endif
};


/*******************************************************************************
 *
 * Class WavePatch: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Access u (given time index)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::u(const int a_idxStep)
  -> PatchSolData&
{
  return m_u[a_idxStep][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Const access u (given time index)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::u(const int a_idxStep) const
  -> const PatchSolData&
{
  return m_u[a_idxStep][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Access un (time level n)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::un()
  -> PatchSolData&
{
  return m_u[m_idxStep][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Constant access un (time level n)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::un() const
  -> const PatchSolData&
{
  return m_u[m_idxStep][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Access unp1 (time level n+1)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::unp1()
  -> PatchSolData&
{
  return m_u[m_idxStepUpdate][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Constant access unp1 (time level n+1)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::unp1() const
  -> const PatchSolData&
{
  return m_u[m_idxStepUpdate][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Access unm1 (time level n-1)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::unm1()
  -> PatchSolData&
{
  return m_u[m_idxStepOld][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Constant access unm1 (time level n-1)
/*--------------------------------------------------------------------*/

inline auto
WavePatch::unm1() const
  -> const PatchSolData&
{
  return m_u[m_idxStepOld][m_bidx];
}

/*--------------------------------------------------------------------*/
//  Advance indices (unp1->un, un->unm1)
/*--------------------------------------------------------------------*/

inline void
WavePatch::advanceStepIndex()
{
  const int tmp   = m_idxStepOld;
  m_idxStepOld    = m_idxStep;
  m_idxStep       = m_idxStepUpdate;
  m_idxStepUpdate = tmp;
}

/*--------------------------------------------------------------------*/
//  Current step index
/*--------------------------------------------------------------------*/

inline int
WavePatch::currentStepIndex() const
{
  return m_idxStep;
}

/*--------------------------------------------------------------------*/
//  Old step index
/*--------------------------------------------------------------------*/

inline int
WavePatch::oldStepIndex() const
{
  return m_idxStepOld;
}

/*--------------------------------------------------------------------*/
//  Current time
/*--------------------------------------------------------------------*/

inline Real
WavePatch::time() const
{
  return m_time;
}

/*--------------------------------------------------------------------*/
//  Current iteration
/*--------------------------------------------------------------------*/

inline int
WavePatch::iteration() const
{
  return m_iteration;
}

#ifdef USE_GPU
/*--------------------------------------------------------------------*/
//  Copy data to host
/*--------------------------------------------------------------------*/

inline void
WavePatch::copyToHostAsync(const int a_idxStep)
{
  u(a_idxStep).copyToHostAsync(0);
}

/*--------------------------------------------------------------------*/
//  Copy data to device
/*--------------------------------------------------------------------*/

inline void
WavePatch::copyToDeviceAsync(const int a_idxStep)
{
  u(a_idxStep).copyToDeviceAsync(0);
}
#endif

#endif  /* ! defined _WAVEPATCH_H_ */
