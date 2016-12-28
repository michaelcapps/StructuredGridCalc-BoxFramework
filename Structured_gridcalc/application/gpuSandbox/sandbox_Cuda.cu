#include <cstdio>

#include "sandbox_Cuda.H"
#include "Config.H"
#include "IntVect.H"
#include "Box.H"
#include "CudaFab.H"

__global__ void
testKernel1(Real* a_fab, IntVect a_ivA)
{
  __shared__ IntVect s_iv;
  if (threadIdx.x < g_SpaceDim)
    {
      s_iv[threadIdx.x] = a_ivA[threadIdx.x];
    }
  IntVect r_iv(s_iv);
  r_iv += s_iv;
  if (threadIdx.x < g_SpaceDim)
    {
      a_fab[threadIdx.x] = r_iv[threadIdx.x];
      printf("r_iv[%d]: %d ", threadIdx.x, r_iv[threadIdx.x]);
    }
}

void
testCuda1(SymbolPair<Real> a_fab)
{
  IntVect ivA(D_DECL(2, 3, 4));
  IntVect ivB(D_DECL(1, 2, 3));
  std::cout << ivA << std::endl;
  ivB.max(ivA);
  testKernel1<<<1, 16>>>(static_cast<Real*>(a_fab.device), ivA);
}

/*----------------------------------------------------------------------------*/

__global__ void
testKernel2(Real* a_fab, Box a_box)
{
  __shared__ Box s_box;
  if (threadIdx.x < 2*g_SpaceDim)
    {
      s_box[threadIdx.x] = a_box[threadIdx.x];
    }
  Box r_box(s_box);
  r_box.grow(1);
  // __syncthreads();
  if (threadIdx.x < 2*g_SpaceDim)
    {
      a_fab[threadIdx.x] = r_box[threadIdx.x];
      printf("r_box[%d]: %d ", threadIdx.x, r_box[threadIdx.x]);
    }
  
}

void
testCuda2(SymbolPair<Real> a_fab)
{
  Box box(IntVect(D_DECL(-1, 0, 0)), IntVect(D_DECL(0, 1, 1)));
  box.shift(1, 0);
  std::cout << box << std::endl;
  testKernel2<<<1, 16>>>(static_cast<Real*>(a_fab.device), box);
}

/*----------------------------------------------------------------------------*/

__global__ void
testKernel3(Real* a_fab, Box a_fabBox)
{
  __shared__ int stride[g_SpaceDim];
  __shared__ int offset;
  __shared__ int cstride;
  if (threadIdx.x == 0)
    {
      a_fabBox.getStride(stride);
      offset = a_fabBox.getOffset(stride);
      cstride = a_fabBox.size();
    }
  __syncthreads();
  IntVect idxVec;
  a_fabBox.linToVec(threadIdx.x, stride, idxVec);
  int idxLin0 = a_fabBox.vecToLin0(idxVec, stride);
  if (idxLin0 + offset != threadIdx.x)
    {
      printf("Conversion failed for thread %2d: vec: (%2d,%2d,%2d) lin: %2d\n",
             threadIdx.x, idxVec[0], idxVec[1], idxVec[2], idxLin0);
    }
  a_fab[idxLin0 + offset          ] = (Real)-1.0;
  a_fab[idxLin0 + offset + cstride] = (Real)-2.0;
}

void
testCuda3(SymbolPair<Real> a_fab, const Box& a_fabBox)
{
  CH_assert(a_fabBox.size() <= CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_BLOCK);
  testKernel3<<<1, a_fabBox.size()>>>(static_cast<Real*>(a_fab.device),
                                      a_fabBox);
}

/*----------------------------------------------------------------------------*/

__global__ void
testKernel4(Real* a_fab, Box a_fabBox, Box a_workBox)
{
  __shared__ int fabStride[g_SpaceDim];
  __shared__ int fabCstride;
  if (threadIdx.x == 0)
    {
      a_fabBox.getStride(fabStride);
      fabCstride = a_fabBox.size();
    }
  __syncthreads();

  // Get the cell
  IntVect idxVec;
  {
    a_workBox.linToVec(threadIdx.x, idxVec);
  }

  // Get index into fab
  int idxLin0 = a_fabBox.vecToLin0(idxVec, fabStride);
  a_fab[idxLin0             ] = (Real)1.0;
  a_fab[idxLin0 + fabCstride] = (Real)2.0;
}

void
testCuda4(SymbolPair<Real> a_fab, const Box& a_fabBox, const Box& a_workBox)
{
  CH_assert(a_workBox.size() <= CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_BLOCK);
  testKernel4<<<1, a_workBox.size()>>>(
    // Add the offset into the pointer address
    static_cast<Real*>(a_fab.device) + a_fabBox.getOffset(),
    a_fabBox,
    a_workBox);
}

/*----------------------------------------------------------------------------*/

__global__ void
testKernel5(Real* a_fabData, Box a_fabBox, int a_fabNcomp, Box a_workBox)
{
  __shared__ CudaFab<Real> fab;
  if (threadIdx.x == 0)
    {
      // Add the offset into the pointer address
      fab.define(a_fabData, a_fabBox, a_fabNcomp);
    }
  __syncthreads();

  // Get the cell
  IntVect idxVec;
  {
    a_workBox.linToVec(threadIdx.x, idxVec);
  }

  // Get index into fab
  fab(idxVec, 0) = (Real)3.0;
  fab(idxVec, 1) = (Real)4.0;
}

void
testCuda5(SymbolPair<Real> a_fabData, const Box& a_fabBox, const int a_fabNcomp,
          const Box& a_workBox)
{
  CH_assert(a_workBox.size() <= CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_BLOCK);
  testKernel5<<<1, a_workBox.size()>>>(
    static_cast<Real*>(a_fabData.device),
    a_fabBox,
    a_fabNcomp,
    a_workBox);
}

/*----------------------------------------------------------------------------*/

__global__ void
testKernel6(CudaFab<Real> a_fabA, CudaFab<Real> a_fabB, const Box a_workBox)
{
  __shared__ Real s_slabData[2*3*g_blkSizeGhost*g_blkSizeGhost];
  __shared__ SlabFab<Real, 3> s_slabFab;
  __shared__ CudaFab<Real> s_fabA;
  __shared__ CudaFab<Real> s_fabB;

  // Load fab meta-data for Fabs A and B.
  {
    const int numThrCopy = CudaFab<Real>::numThrCopy();
    CH_assert(numThrCopy < blockDim.x);
    s_fabA.define(a_fabA, 0, numThrCopy);
    s_fabB.define(a_fabB, 0, numThrCopy);
  }

  // Compute (_Ar_ithmetic) index, saved to avoid repeat linear->vector
  // conversion
  IntVect ivecAr;
  a_workBox.linToVec(threadIdx.x, ivecAr);

  // _L_oad/_S_tore index, saved to avoid repeat linear->vector conversion
  IntVect ivecLS;
  // Set up the cache window
  {
    Box LSbox(a_workBox);  // This is the initial cache window
    LSbox.hiVect(2) = LSbox.loVect(2);
    LSbox.grow(1);
    // Initialize the slab Fab cache.  Note that it is shifted one cell towards
    // the low end so it can be shifted back at the beginning of the first
    // iteration.
    LSbox.shift(-1, 2);
    int locNDEnd = LSbox.hiVect(2);  // vec[2] to end loading initial data
    int locNDBeg = locNDEnd - 1;     // vec[2] to begin loading initial data
    s_slabFab.define(s_slabData,          // Data
                     LSbox,               // Initial window
                     2,                   // Number of components
                     2,                   // Normal direction
                     locNDBeg, locNDEnd,  // Initial loading
                     ivecLS,              // Vector index (output)
                     0,                   // Start component in source Fab
                     s_fabA,              // Source data
                     g_blkSizeGhost*g_blkSizeGhost);  // # threads for loading
  }

  // Loop over slabs
  for (int iSlab = 0; iSlab != g_blkSize; ++iSlab)
    {
      // Shift the slab
      s_slabFab.shift(1, ivecLS);

      if (threadIdx.x < g_blkSize*g_blkSize)
        {
          IntVect basis(D_DECL(0, 0, 0));
          basis[2] = 1;
          s_fabB(ivecAr, 0) = (s_slabFab(ivecAr - basis, 0) +
                               s_slabFab(ivecAr, 0) +
                               s_slabFab(ivecAr + basis, 0));
          s_fabB(ivecAr, 1) = (s_slabFab(ivecAr - basis, 1) +
                               s_slabFab(ivecAr, 1) +
                               s_slabFab(ivecAr + basis, 1));

          // Shift the arithmetic IntVect for the next slab iteration
          ivecAr[2] += 1;
        }
    }
}

void
testCuda6(BaseFab<Real> *const a_fabA,
          BaseFab<Real> *const a_fabB,
          const Box&           a_workBox)
{
  BaseFabData<Real> *const fabA =
    reinterpret_cast<BaseFabData<Real>*>(a_fabA);
  CudaFab<Real> cudaFabA;
  cudaFabA.define(*fabA);

  BaseFabData<Real> *const fabB =
    reinterpret_cast<BaseFabData<Real>*>(a_fabB);
  CudaFab<Real> cudaFabB;
  cudaFabB.define(*fabB);

  testKernel6<<<1, g_blkSizeGhost*g_blkSizeGhost>>>(
    cudaFabA,
    cudaFabB,
    a_workBox);
}
