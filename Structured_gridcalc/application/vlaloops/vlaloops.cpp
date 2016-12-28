#include <iostream>

#include "BaseFab.H"
#include "BoxIterator.H"
#include "BaseFabMacros.H"

int main()
{
  const int n = 4;
  // Note that box is (1,1,1) to (n,n,n)
  Box box(IntVect::Unit, n*IntVect::Unit);
  FArrayBox fab(box, 1);

//--Initialize

  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      fab(*bit, 0) = bit->sum();
    }

//--VLA

  const IntVect lo = fab.box().loVect();
  const IntVect hi = fab.box().hiVect();
  const IntVect dims = fab.box().dimensions();
  const int n0 = dims[0];
  const int n1 = dims[1];
  const int n2 = dims[2];

  // Basic VLA
  {
    Real (*__restrict__ vla)[n2][n1][n0] = // type and symbol
      (Real (*__restrict__)[n2][n1][n0])   // cast
      (fab.dataPtr());                     // memory location
    for (int k = lo[2]; k <= hi[2]; ++k)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int i = lo[0]; i <= hi[0]; ++i)
        {
          CH_assert(&(vla[0][k-1][j-1][i-1]) == &(fab(IntVect(i,j,k), 0)));
          /* We don't want to do this ---^
             Note that this expands to

             0*cstride + (k-1)*stride2 + (j-1)*stride1 + (i-1)*stride0

             and we can rewrite as

             0*cstride + k*stride2 + j*stride1 + i*stride0 -
               (1*stride2 + 1*stride1 + 1*stride0)

             The second row is constant for a fab and gives the offset for
             0-based indexing.  In the next example, we subtract that from
             the data pointer
          */
        }
  }

  // Subtract offset due to 0-based indexing from pointer
  {
    Real (*__restrict__ vla)[n2][n1][n0] = // type and symbol
      (Real (*__restrict__)[n2][n1][n0])   // cast
      // memory location with offset subtracted
      (fab.dataPtr() - (lo[0] + n0*(lo[1] + n1*lo[2])));
    for (int k = lo[2]; k <= hi[2]; ++k)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int i = lo[0]; i <= hi[0]; ++i)
        {
          CH_assert(&(vla[0][k][j][i]) == &(fab(IntVect(i,j,k), 0)));
          /* Yay!  No more offsets ^
          */
        }
  }

  // Add in some C++ syntatic sugar to simplify the definition of the VLA
  // pointer
  {
    auto vla =                                                // type and symbol
      (decltype(fab)::value_type (*__restrict__)[n2][n1][n0]) // cast  
      // memory location: subtract offset
      (fab.dataPtr() - (lo[0] + n0*(lo[1] + n1*lo[2])));
    for (int k = lo[2]; k <= hi[2]; ++k)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int i = lo[0]; i <= hi[0]; ++i)
        {
          CH_assert(&(vla[0][k][j][i]) == &(fab(IntVect(i,j,k), 0)));
        }
  }

  // Use macros to encode the VLAs and loops
  {
    MD_ARRAY_RESTRICT(vla, fab);
    MD_BOXLOOP(box, i)
      {
        CH_assert(&(vla[MD_IX(i, 0)]) == &(fab(IntVect(i0,i1,i2), 0)));
      }
  }

  return 0;
}
