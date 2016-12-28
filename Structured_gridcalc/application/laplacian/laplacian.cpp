#include <iostream>
#include <iomanip>

#include "BaseFab.H"
#include "BoxIterator.H"
#include "BaseFabMacros.H"
#include "Stopwatch.H"

int main()
{
  const int n = 64;
  const int niter = 16;
  Box allocBox(IntVect::Unit, n*IntVect::Unit);
  FArrayBox fabA(allocBox, 1);

//--Initialize

  for (BoxIterator bit(allocBox); bit.ok(); ++bit)
    {
      fabA(*bit, 0) = bit->sum();
    }

  Stopwatch<> timer;
  std::cout << std::left << std::setw(40) << "Clock resolution (ns): "
            << timer.resolution() << std::endl;

/*--------------------------------------------------------------------*
 * Lapacian
 *--------------------------------------------------------------------*/

  Box computeBox(allocBox);
  computeBox.grow(-1);
  FArrayBox fabB(computeBox, 1);

//--Warm up cache

  {
    fabB.setVal(0.);
    MD_ARRAY_RESTRICT(arrA, fabA);
    MD_ARRAY_RESTRICT(arrB, fabB);
    for (int iter = 0; iter != niter; ++iter)
      {
        for (int dir = 0; dir != g_SpaceDim; ++dir)
          {
            const int MD_ID(o, dir);
            MD_BOXLOOP(computeBox, i)
              {
                arrB[MD_IX(i, 0)] -= 0.5*(arrA[MD_OFFSETIX(i,+,o, 0)] -
                                          2*arrA[MD_IX(i, 0)] +
                                          arrA[MD_OFFSETIX(i,-,o, 0)]);
              }
          }
      }
  }

//--BoxIterator

  timer.reset();
  {
    timer.start();
    fabB.setVal(0.);
    for (int iter = 0; iter != niter; ++iter)
      {
        for (int dir = 0; dir != g_SpaceDim; ++dir)
          {
            const IntVect ivo((dir == 0), (dir == 1), (dir == 2));
            for (BoxIterator bit(computeBox); bit.ok(); ++bit)
              {
                const IntVect iv = *bit;
                fabB(iv, 0) -= 0.5*(fabA(iv+ivo, 0) -
                                    2*fabA(iv, 0) +
                                    fabA(iv-ivo, 0));
              }
          }
      }
    timer.stop();
  }
  std::cout << std::left << std::setw(40) << "Time for BoxIterator (ms): "
            << timer.time() << std::endl;
  int stat = 0;
  for (BoxIterator bit(computeBox); bit.ok(); ++bit)
    {
      const IntVect iv = *bit;
      if (fabB(iv, 0) != 0.) ++stat;
    }
  CH_assert(stat == 0);

//--VLA macros

  timer.reset();
  {
    timer.start();
    fabB.setVal(0.);
    MD_ARRAY_RESTRICT(arrA, fabA);
    MD_ARRAY_RESTRICT(arrB, fabB);
    for (int iter = 0; iter != niter; ++iter)
      {
        for (int dir = 0; dir != g_SpaceDim; ++dir)
          {
            const int MD_ID(o, dir);
            MD_BOXLOOP(computeBox, i)
              {
                arrB[MD_IX(i, 0)] -= 0.5*(arrA[MD_OFFSETIX(i,+,o, 0)] -
                                          2*arrA[MD_IX(i, 0)] +
                                          arrA[MD_OFFSETIX(i,-,o, 0)]);
              }
          }
      }
    timer.stop();
  }
  std::cout << std::left << std::setw(40) << "Time for VLA macro (ms): "
            << timer.time()
            << " : " << fabB(IntVect(rand() % n, rand() % n, rand() % n), 0)
            << std::endl;
  for (BoxIterator bit(computeBox); bit.ok(); ++bit)
    {
      const IntVect iv = *bit;
      if (fabB(iv, 0) != 0.) ++stat;
    }
  CH_assert(stat == 0);

  return stat;
}
