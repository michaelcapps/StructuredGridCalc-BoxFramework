
#ifndef _SANDBOX_CUDA_H_
#define _SANDBOX_CUDA_H_

// This file should only be included if USE_GPU is defined

#include "CudaSupport.H"
#include "Parameters.H"

//--Forward declarations

class Box;
template <typename T>
class BaseFab;

void testCuda1(SymbolPair<Real> a_fab);
void testCuda2(SymbolPair<Real> a_fab);
void testCuda3(SymbolPair<Real> a_fab, const Box& a_box);
void testCuda4(SymbolPair<Real> a_fab, const Box& a_fabBox,
               const Box& a_workBox);
void testCuda5(SymbolPair<Real> a_fabData, const Box& a_fabBox,
               const int a_fabNcomp, const Box& a_workBox);

constexpr int g_blkSize      = 4;
constexpr int g_blkSizeGhost = 6;

void testCuda6(BaseFab<Real> *const a_fabA,
               BaseFab<Real> *const a_fabB,
               const Box& a_workBox);

#endif  /* ! defined _SANDBOX_CUDA_H_ */
