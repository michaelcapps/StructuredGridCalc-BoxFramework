
#ifndef _CUDASUPPORT_H_
#define _CUDASUPPORT_H_

#include <iostream>

#ifdef USE_GPU
#include <cuda_runtime.h>
#endif


/******************************************************************************/
/**
 * \file CudaSupport.H
 *
 * \brief Support for using CUDA
 *
 *//*+*************************************************************************/

/// SDK safe call wrapper.
# define CU_SAFE_CALL( call )                                           \
  {                                                                     \
    cudaError_t cuer = call;                                            \
    if (cuer != cudaSuccess)                                            \
      {                                                                 \
        std::cout << cudaGetErrorString(cuer) << std::endl;             \
      }                                                                 \
    CH_assert(cuer == cudaSuccess);                                     \
  }

//--Type of pointer on device

typedef void* AccelSymbol;
typedef void* AccelPointer;

//--Matching pointers on host and device

template <typename T>
struct SymbolPair
{
  SymbolPair()
    : host(0)
    { }
  SymbolPair(T *const a_host, AccelSymbol a_device)
    : host(a_host), device(a_device)
    { }
  T *host;
  AccelSymbol device;
};

#endif
