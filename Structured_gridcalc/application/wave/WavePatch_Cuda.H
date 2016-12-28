
#ifndef _WAVEPATCH_CUDA_H_
#define _WAVEPATCH_CUDA_H_


/******************************************************************************/
/**
 * \file WavePatchCuda.H
 *
 * \brief Isolated Cuda drivers and kernels for class WavePatch
 *
 *//*+*************************************************************************/

//--Forward declarations

template <typename T>
class BaseFab;
typedef void* AccelPointer;

namespace WavePatch_Cuda
{

  // Standard co-dimension 1 number of threads
  constexpr int g_blkSize       = 16; ///< Size in one dimension of a canonical
                                      ///< compute block
  constexpr int g_blkSizeGhost  = g_blkSize + 2;
                                      ///< Size to load all ghosts (for RHS
                                      ///< kernel)
  constexpr int g_numThrRHSAr   = g_blkSize*g_blkSize;
                                      ///< Number of threads for compute in RHS
                                      ///< kernel
  constexpr int g_numThrRHSLS   = g_blkSizeGhost*g_blkSizeGhost;
                                      ///< Number of threads for loading or
                                      ///< storing data in RHS kernel
  constexpr int g_numThrBC      = g_blkSize*g_blkSize;
                                      ///< Number of threads for BC kernel
  constexpr int g_RHSNrmDirSlab = g_SpaceDim-1;
                                      ///< Direction normal to the slab for RHS
                                      ///< kernel

  /// Information required for each work box for the BC kernel
  struct WorkBoxInfoBC
  {
    Box m_workBox;
    int m_dir;
    int m_side;  // Needs to be either -1 or +1
  };

  /// Construct necessary data for executing the GPU kernels
  void construct(const Box&                        a_domain,
                 const BaseFab<Real> *const *const a_patchData,
                 AccelPointer&                     a_cudaFab_device,
                 AccelPointer&                     a_workBoxesRHS_device,
                 AccelPointer&                     a_workBoxInfoBC_device,
                 int&                              a_numBlkRHS,
                 int&                              a_numThrBC);

  /// Destroy data on the GPU created during construction
  void destroy(AccelPointer& a_cudaFab_device,
               AccelPointer& a_workBoxesRHS_device,
               AccelPointer& a_workBoxInfoBC_device);

  /// Driver to run BC kernel on GPU
  void driverBC(int a_numBlkBC, const int a_idxStep);

  /// Driver to run RHS kernel on GPU
  void driverRHS(const int  a_numBlkRHS,
                 const int  a_idxStep,
                 const int  a_idxStepUpdate,
                 const int  a_idxStepOld,
                 const Real a_factor);

  /// Driver to advance a group of iterations
  void driverAdvanceIterGroup(const int   a_numIter,
                              const int   a_numBlkBC,
                              const int   a_numBlkRHS,
                              int&        a_idxStep,
                              int&        a_idxStepUpdate,
                              int&        a_idxStepOld,
                              const Real  a_factor,
                              cudaEvent_t a_cuEvent_iterGroupStart,
                              cudaEvent_t a_cuEvent_iterGroupEnd);

}  /* namespace WavePatch_Cuda */

#endif  /* ! defined _WAVEPATCH_CUDA_H_ */
