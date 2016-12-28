
/******************************************************************************/
/**
 * \file WavePatchCuda.cu
 *
 * \brief Isolated Cuda drivers and kernels for class WavePatch
 *
 *//*+*************************************************************************/

#include "CudaFab.H"
#include "WavePatch_Cuda.H"

template <typename T>
class BaseFab;

//--Constant memory

/* All of these structures live in global memory.  In constant memory, there
 * is only a pointer to the location in global memory.
 */

__constant__ CudaFab<Real>* c_fabs;
__constant__ Box* c_workBoxesRHS;
__constant__ WavePatch_Cuda::WorkBoxInfoBC* c_workBoxInfoBC;

__global__ void kernelBC(const int a_idxStep);
__global__ void kernelRHS(const int a_timeN,
                          const int a_timeNp1,
                          const int a_timeNm1,
                          const Real a_factor);


/*******************************************************************************
 *
 * Class wavePatch: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Construct necessary data for executing the GPU kernels
/** All data is stored in global memory and only pointers to this
 *  data is stored in constant memory.
 *  <ul>
 *    <li> The fabs for times n-1, n, n+1 are stored
 *    <li> The work boxes, directions, and sides are stored for the
 *         BC kernel
 *    <li> The work boxes are stored for the RHS kernel.
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
WavePatch_Cuda::construct(
  const Box&                        a_domain,
  const BaseFab<Real> *const *const a_patchData,
  AccelPointer&                     a_cudaFab_device,
  AccelPointer&                     a_workBoxesRHS_device,
  AccelPointer&                     a_workBoxInfoBC_device,
  int&                              a_numBlkRHS,
  int&                              a_numBlkBC)
{
  a_numBlkRHS = a_domain.dimensions().product()/
    (g_blkSize*g_blkSize*g_blkSize);
  a_numBlkBC =
    2*(  a_domain.dimensions()[0]*a_domain.dimensions()[1]
       + a_domain.dimensions()[1]*a_domain.dimensions()[2]
       + a_domain.dimensions()[2]*a_domain.dimensions()[0])/
    (g_blkSize*g_blkSize);

//--Setup CudaFabs (note that these are aliases and the array data is managed
//--separately via the BaseFab class).

  const BaseFabData<Real> *const *const patchData =
    reinterpret_cast<const BaseFabData<Real> *const *>(a_patchData);
  // Fabs on host
  const int numBytesCudaFab = 3*sizeof(CudaFab<Real>);
  CudaFab<Real>* cudaFab_host;
  CU_SAFE_CALL(cudaHostAlloc(&cudaFab_host,
                             numBytesCudaFab,
                             cudaHostAllocWriteCombined));
  cudaFab_host[0].define(*patchData[0]);
  cudaFab_host[1].define(*patchData[1]);
  cudaFab_host[2].define(*patchData[2]);

  // Fabs on device
  CU_SAFE_CALL(cudaMalloc(&a_cudaFab_device, numBytesCudaFab));
  // Copy to device
  CU_SAFE_CALL(cudaMemcpy(a_cudaFab_device,
                          cudaFab_host,
                          numBytesCudaFab,
                          cudaMemcpyHostToDevice));

  // Set pointer in constant memory
  CU_SAFE_CALL(cudaMemcpyToSymbol(c_fabs,
                                  &a_cudaFab_device,
                                  sizeof(AccelPointer)));

  // Memory on host can be released
  CU_SAFE_CALL(cudaFreeHost(cudaFab_host));
  
//--Setup work boxes for RHS

  // RHS on host
  const int numBytesRHS = a_numBlkRHS*sizeof(Box);
  Box* workBoxesRHS_host;
  CU_SAFE_CALL(cudaHostAlloc(&workBoxesRHS_host,
                             numBytesRHS,
                             cudaHostAllocWriteCombined));
  int c = 0;
  D_INVTERM(for (int i0 = a_domain.loVect()[0]; i0 <= a_domain.hiVect()[0];
                 i0 += g_blkSize),
            for (int i1 = a_domain.loVect()[1]; i1 <= a_domain.hiVect()[1];
                 i1 += g_blkSize),
            for (int i2 = a_domain.loVect()[2]; i2 <= a_domain.hiVect()[2];
                 i2 += g_blkSize))
    {
      Box computeBox(IntVect(D_DECL(i0, i1, i2)),
                     IntVect(D_DECL(i0 + g_blkSize - 1,
                                    i1 + g_blkSize - 1,
                                    i2 + g_blkSize - 1)));
      // The work boxes are the core cells in the box, but there is one layer
      // of cells in the last spatial dimension.  The work box can then be
      // shifted through the cells in this direction during kernel execution.
      computeBox.hiVect(g_SpaceDim-1) = computeBox.loVect(g_SpaceDim-1);
      CH_assert(c < a_numBlkRHS);
      workBoxesRHS_host[c++] = computeBox;

      // Testing
      CH_assert(computeBox.size() == g_numThrRHSAr);
      Box loadstoreBox(computeBox);
      loadstoreBox.grow(1);
      CH_assert(loadstoreBox.size() == 3*g_numThrRHSLS);
      loadstoreBox.hiVect(g_SpaceDim-1) = loadstoreBox.loVect(g_SpaceDim-1);
      CH_assert(loadstoreBox.size() == g_numThrRHSLS);
    }
  CH_assert(c == a_numBlkRHS);

  // RHS on device
  CU_SAFE_CALL(cudaMalloc(&a_workBoxesRHS_device, numBytesRHS));
  // Copy to device
  CU_SAFE_CALL(cudaMemcpy(a_workBoxesRHS_device,
                          workBoxesRHS_host,
                          numBytesRHS,
                          cudaMemcpyHostToDevice));

  // Set pointer in constant memory
  CU_SAFE_CALL(cudaMemcpyToSymbol(c_workBoxesRHS,
                                  &a_workBoxesRHS_device,
                                  sizeof(AccelPointer)));

  // Memory on host can be released
  CU_SAFE_CALL(cudaFreeHost(workBoxesRHS_host));

//--Setup work boxes for BC

  // BC on host
  const int numBytesBC = a_numBlkBC*sizeof(WorkBoxInfoBC);
  WorkBoxInfoBC* workBoxInfoBC_host;
  CU_SAFE_CALL(cudaHostAlloc(&workBoxInfoBC_host,
                             numBytesBC,
                             cudaHostAllocWriteCombined));
  Box domainTest(a_domain);
  domainTest.grow(-1);
  c = 0;
  D_INVTERM(for (int i0 = a_domain.loVect()[0]; i0 <= a_domain.hiVect()[0];
                 i0 += g_blkSize),
            for (int i1 = a_domain.loVect()[1]; i1 <= a_domain.hiVect()[1];
                 i1 += g_blkSize),
            for (int i2 = a_domain.loVect()[2]; i2 <= a_domain.hiVect()[2];
                 i2 += g_blkSize))
    {
      Box computeBox(IntVect(D_DECL(i0, i1, i2)),
                     IntVect(D_DECL(i0 + g_blkSize - 1,
                                    i1 + g_blkSize - 1,
                                    i2 + g_blkSize - 1)));
      if (!domainTest.contains(computeBox))
        {
          for (int dir = 0; dir != g_SpaceDim; ++dir)
            {
              if (computeBox.loVect()[dir] == a_domain.loVect()[dir])
                {
                  // The work boxes are the interior cells adjacent to the
                  // domain in this direction.
                  Box workBox(computeBox);
                  workBox.adjBox(-1, dir, -1);
                  CH_assert(c < a_numBlkBC);
                  WorkBoxInfoBC& workBoxInfoBC = workBoxInfoBC_host[c++];
                  workBoxInfoBC.m_workBox = workBox;
                  workBoxInfoBC.m_dir = dir;
                  workBoxInfoBC.m_side = -1;
                }
              if (computeBox.hiVect()[dir] == a_domain.hiVect()[dir])
                {
                  Box workBox(computeBox);
                  workBox.adjBox(-1, dir,  1);
                  CH_assert(c < a_numBlkBC);
                  WorkBoxInfoBC& workBoxInfoBC = workBoxInfoBC_host[c++];
                  workBoxInfoBC.m_workBox = workBox;
                  workBoxInfoBC.m_dir = dir;
                  workBoxInfoBC.m_side =  1;
                }
            }
        }
    }
  CH_assert(c == a_numBlkBC);

  // BC on device
  CU_SAFE_CALL(cudaMalloc(&a_workBoxInfoBC_device, numBytesBC));
  // Copy to device
  CU_SAFE_CALL(cudaMemcpy(a_workBoxInfoBC_device,
                          workBoxInfoBC_host,
                          numBytesBC,
                          cudaMemcpyHostToDevice));

  // Set pointer in constant memory
  CU_SAFE_CALL(cudaMemcpyToSymbol(c_workBoxInfoBC,
                                  &a_workBoxInfoBC_device,
                                  sizeof(AccelPointer)));

  // Memory on host can be released
  CU_SAFE_CALL(cudaFreeHost(workBoxInfoBC_host));

//--Configure cache for kernels

  CU_SAFE_CALL(cudaFuncSetCacheConfig(kernelBC, cudaFuncCachePreferL1));
  CU_SAFE_CALL(cudaFuncSetCacheConfig(kernelRHS, cudaFuncCachePreferShared));
}

/*--------------------------------------------------------------------*/
//  Destroy data on the GPU created during construction
/**
 *//*-----------------------------------------------------------------*/

void
WavePatch_Cuda::destroy(AccelPointer& a_cudaFab_device,
                        AccelPointer& a_workBoxesRHS_device,
                        AccelPointer& a_workBoxInfoBC_device)
{
  CU_SAFE_CALL(cudaFree(a_cudaFab_device));
  a_cudaFab_device = nullptr;
  CU_SAFE_CALL(cudaFree(a_workBoxesRHS_device));
  a_workBoxesRHS_device = nullptr;
  CU_SAFE_CALL(cudaFree(a_workBoxInfoBC_device));
  a_workBoxInfoBC_device = nullptr;
}

/*--------------------------------------------------------------------*/
/// Sets reflection BC for wave solution
/** \param[in]  a_idxStep
 *                      Index of BaseFab at time 'n'
 *//*-----------------------------------------------------------------*/

__global__ void kernelBC(const int a_idxStep)
{
  WavePatch_Cuda::WorkBoxInfoBC& r_workBoxInfo = c_workBoxInfoBC[blockIdx.x];

  // Get the source cell
  IntVect r_idxVecSrc;
  r_workBoxInfo.m_workBox.linToVec(threadIdx.x, r_idxVecSrc);

  // Get the destination cell
  IntVect r_idxVecDst(r_idxVecSrc);
  r_idxVecDst[r_workBoxInfo.m_dir] += r_workBoxInfo.m_side;
  c_fabs[a_idxStep](r_idxVecDst, 0) = c_fabs[a_idxStep](r_idxVecSrc, 0);
}

/*--------------------------------------------------------------------*/
//  Driver to run BC kernel on GPU
/** \param[in]  a_numBlkBC
 *                      Number of blocks for the kernel.  Each block
 *                      is a face of a block on a boundary.
 *  \param[in]  a_idxStep
 *                      Index of BaseFab at time 'n'
 *//*-----------------------------------------------------------------*/

void
WavePatch_Cuda::driverBC(const int a_numBlkBC, const int a_idxStep)
{
  kernelBC<<<a_numBlkBC, g_numThrBC, 0, 0>>>(a_idxStep);
}

/*--------------------------------------------------------------------*/
/// RHS kernel on GPU
/** \param[in]  a_timeN Index for fab at time \f$u^n\f$,
 *  \param[in]  a_timeNp1
 *                      Index for fab at time \f$u^{n+1}\f$
 *  \param[in]  a_timeNp1
 *                      Index for fab at time \f$u^{n-1}\f$
 *  \param[in]  a_factor
 *                      \f$(\frac{\Delta t c}{\Delta x})^2\frac{1}{D}\f$
 *//*-----------------------------------------------------------------*/

__global__ void
__launch_bounds__(WavePatch_Cuda::g_numThrRHSAr, 2)
kernelRHS(const int a_timeN,
          const int a_timeNp1,
          const int a_timeNm1,
          const Real a_factor)
{
  //**FIXME Implement advance by updating u^{n+1}

}

/*--------------------------------------------------------------------*/
//  Driver to run RHS kernel on GPU
/** \param[in]  a_numBlkRHS
 *                      Number of CUDA blocks for RHS kernel
 *  \param[in]  a_idxStep
 *                      Index for fab at time \f$u^n\f$,
 *  \param[in]  a_idxStepUpdate
 *                      Index for fab at time \f$u^{n+1}\f$
 *  \param[in]  a_idxStepOld
 *                      Index for fab at time \f$u^{n-1}\f$
 *  \param[in]  a_factor
 *                      \f$(\frac{\Delta t c}{\Delta x})^2\frac{1}{D}\f$
 *//*-----------------------------------------------------------------*/

void
WavePatch_Cuda::driverRHS(const int  a_numBlkRHS,
                          const int  a_idxStep,
                          const int  a_idxStepUpdate,
                          const int  a_idxStepOld,
                          const Real a_factor)
{
  kernelRHS<<<a_numBlkRHS, g_numThrRHSLS, 0, 0>>>(a_idxStep,
                                                  a_idxStepUpdate,
                                                  a_idxStepOld,
                                                  a_factor);
}

/*--------------------------------------------------------------------*/
/** Advances a group of iterations on the GPU
 *  \param[in]  a_numIter
 *                      Number of iterations in the group
 *  \param[in]  a_numBlkBC
 *                      Number of CUDA blocks for BC kernel
 *  \param[in]  a_numBlkRHS
 *                      Number of CUDA blocks for RHS kernel
 *  \param[in]  a_idxStep
 *                      Index for fab at time \f$u^n\f$,
 *  \param[in]  a_idxStepUpdate
 *                      Index for fab at time \f$u^{n+1}\f$
 *  \param[in]  a_idxStepOld
 *                      Index for fab at time \f$u^{n-1}\f$
 *  \param[in]  a_factor
 *                      \f$(\frac{\Delta t c}{\Delta x})^2\frac{1}{D}\f$
 *  \param[out] a_cuEvent_iterGroupStart
 *                      Timer at start of iteration group
 *  \param[out] a_cuEvent_iterGroupEnd
 *                      Timer at end of iteration group
 *//*-----------------------------------------------------------------*/

void
WavePatch_Cuda::driverAdvanceIterGroup(const int   a_numIter,
                                       const int   a_numBlkBC,
                                       const int   a_numBlkRHS,
                                       int&        a_idxStep,
                                       int&        a_idxStepUpdate,
                                       int&        a_idxStepOld,
                                       const Real  a_factor,
                                       cudaEvent_t a_cuEvent_iterGroupStart,
                                       cudaEvent_t a_cuEvent_iterGroupEnd)
{
  cudaEventRecord(a_cuEvent_iterGroupStart, 0);
  for (int iter = 0; iter != a_numIter; ++iter)
    {
      kernelBC<<<a_numBlkBC, g_numThrBC, 0, 0>>>(a_idxStep);
      kernelRHS<<<a_numBlkRHS, g_numThrRHSAr, 0, 0>>>(a_idxStep,
                                                      a_idxStepUpdate,
                                                      a_idxStepOld,
                                                      a_factor);
      // Advance step indices
      const int tmp   = a_idxStepOld;
      a_idxStepOld    = a_idxStep;
      a_idxStep       = a_idxStepUpdate;
      a_idxStepUpdate = tmp;
    }
  cudaEventRecord(a_cuEvent_iterGroupEnd, 0);
}
