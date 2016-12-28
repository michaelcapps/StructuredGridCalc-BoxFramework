#include <iostream>
#include <cmath>

#include "BaseFab.H"
#include "BaseFabMacros.H"
#include "sandbox_Cuda.H"

inline Real
test6SetVal(const IntVect& a_iv, const int a_comp)
{
  switch (a_comp)
    {
    case 0:
      return std::fabs((Real)a_iv[2]) + std::fabs(0.1*a_iv[1])
        + std::fabs(0.01*a_iv[0]);
      break;
    case 1:
      return -99.9;
    default:
      CH_assert(false);
    }
  return 1.2345678E99;
}

int main()
{
  Box coreBox(-IntVect::Unit, 2*IntVect::Unit);
  Box ghostBox(coreBox);
  ghostBox.grow(1);
  FArrayBox fabA(ghostBox, 2);
  MD_ARRAY_RESTRICT(arrA, fabA);

//--Test working with IntVects on the GPU

  std::cout << "--TEST 1--\n\n";
  fabA.copyToDevice();
  testCuda1(fabA.m_dataSymbol);
  fabA.copyToHost();
  cudaDeviceSynchronize();
  std::cout << "  " << fabA.dataPtr()[0]
            << "  " << fabA.dataPtr()[1]
            << "  " << fabA.dataPtr()[2]
            << std::endl;

//--Test working with Boxes on the GPU

  std::cout << "\n--TEST 2--\n\n";
  testCuda2(fabA.m_dataSymbol);
  fabA.copyToHost();
  cudaDeviceSynchronize();
  std::cout << "  " << fabA.dataPtr()[0]
            << "  " << fabA.dataPtr()[1]
            << "  " << fabA.dataPtr()[2]
            << "  " << fabA.dataPtr()[3]
            << "  " << fabA.dataPtr()[4]
            << "  " << fabA.dataPtr()[5]
            << std::endl;

//--Set all values in the fab

  std::cout << "\n--TEST 3--\n\n";
  testCuda3(fabA.m_dataSymbol, fabA.box());
  fabA.copyToHost();
  cudaDeviceSynchronize();

  // Validate
  MD_BOXLOOP(fabA.box(), i)
    {
      CH_assert(arrA[MD_IX(i, 0)] == (Real)-1.0);
      CH_assert(arrA[MD_IX(i, 1)] == (Real)-2.0);
    }

//--Set all core-cell values in the fab

  std::cout << "\n--TEST 4--\n\n";
  testCuda4(fabA.m_dataSymbol, fabA.box(), coreBox);
  fabA.copyToHost();
  cudaDeviceSynchronize();

  // Validate
  MD_BOXLOOP(coreBox, i)
    {
      CH_assert(arrA[MD_IX(i, 0)] == (Real)1.0);
      CH_assert(arrA[MD_IX(i, 1)] == (Real)2.0);
    }
  MD_BOXLOOP(ghostBox, i)
    {
      if (!coreBox.contains(IntVect(D_DECL(i0, i1, i2))))
        {
          CH_assert(arrA[MD_IX(i, 0)] == (Real)-1.0);
          CH_assert(arrA[MD_IX(i, 1)] == (Real)-2.0);
        }
    }

//--Set all core-cell values in the fab

  std::cout << "\n--TEST 5--\n\n";
  testCuda5(fabA.m_dataSymbol, fabA.box(), fabA.ncomp(), coreBox);
  fabA.copyToHost();
  cudaDeviceSynchronize();

  // Validate
  MD_BOXLOOP(coreBox, i)
    {
      CH_assert(arrA[MD_IX(i, 0)] == (Real)3.0);
      CH_assert(arrA[MD_IX(i, 1)] == (Real)4.0);
    }
  MD_BOXLOOP(ghostBox, i)
    {
      if (!coreBox.contains(IntVect(D_DECL(i0, i1, i2))))
        {
          CH_assert(arrA[MD_IX(i, 0)] == (Real)-1.0);
          CH_assert(arrA[MD_IX(i, 1)] == (Real)-2.0);
        }
    }

//--In this test the CudaFabs are passed as arguments and the shared version is
//--constructed from them.  A moving cache window, SlabFab, is used to perform
//--a stencil operation on all cells.

  std::cout << "\n--TEST 6--\n\n";
  CH_assert(g_SpaceDim == 3);
  std::cout << "Core: " << coreBox << std::endl;
  MD_BOXLOOP(ghostBox, i)
    {
      IntVect iv(i0, i1, i2);
      arrA[MD_IX(i, 0)] = i0;
      arrA[MD_IX(i, 1)] = i0 + 1;
    }
  fabA.copyToDevice();
  FArrayBox fabB(ghostBox, 2);
  MD_ARRAY_RESTRICT(arrB, fabB);
  fabB.setVal(1.234567E89);
  fabB.copyToDevice();
  testCuda6(&fabA, &fabB, coreBox);
  fabB.copyToHost();
  cudaDeviceSynchronize();

  // Validate
  {
    const int MD_ID(o, 2);
    MD_BOXLOOP(coreBox, i)
      {
        CH_assert(arrB[MD_IX(i, 0)] == (arrA[MD_OFFSETIX(i,-,o, 0)] +
                                        arrA[MD_IX(i, 0)] +
                                        arrA[MD_OFFSETIX(i,+,o, 0)]));
        CH_assert(arrB[MD_IX(i, 1)] == (arrA[MD_OFFSETIX(i,-,o, 1)] +
                                        arrA[MD_IX(i, 1)] +
                                        arrA[MD_OFFSETIX(i,+,o, 1)]));
      }
  }

  return 0;
}
