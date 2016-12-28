#include <iostream>
#include <iomanip>

#include "BoxIterator.H"
#include "BaseFab.H"

int main(const int argc, const char* argv[])
{
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
  int status = 0;

//--Tests

  Box boxA(IntVect::Zero, IntVect(D_DECL(2, 2, 2)));

  // Test increments
  {
    if (verbose)
      {
        std::cout << "--- Increments\n";
      }
    BoxIterator boxit(boxA);
    if (*boxit != IntVect::Zero) ++status;
    if (verbose) std::cout << *boxit << std::endl;
#if 1
    if (*++boxit != IntVect(D_DECL(1, 0, 0))) ++status;
    if (verbose)
      {
        std::cout << *boxit << std::endl;
        // Print postfix test
        BoxIterator boxittmp(boxit);
        std::cout << *boxittmp++ << std::endl;
      }
    // Test postfix
    if (*boxit++ != IntVect(D_DECL(1, 0, 0))) ++status;
    if (*boxit != IntVect(D_DECL(2, 0, 0))) ++status;
    if (verbose) std::cout << *boxit << std::endl;
    // Did we loop?
    if (g_SpaceDim > 1)
      {
        if (*++boxit != IntVect(D_DECL(0, 1, 0))) ++status;
        if (verbose) std::cout << *boxit << std::endl;
      }
    boxit += IntVect::Unit;
    if (*boxit != IntVect(D_DECL(1, 2, 1))) ++status;
    if (verbose) std::cout << *boxit << std::endl;
#endif
  }

#if 1
  if (verbose)
    {
      std::cout << "--- The full box\n";
      for (BoxIterator boxit(boxA); boxit.ok(); ++boxit)
        {
          std::cout << *boxit << std::endl;
        }
    }

  // Iterating with BoxIterator should access contiguous elements in a BaseFab
  // with only 1 component.
  {
    if (verbose)
      {
        std::cout << "--- Check the stride between elements accessed with the "
          "iterator\n";
      }
    BaseFab<int> fabA(boxA, 1);
    BoxIterator boxit1(boxA);
    // Also tests assignment
    BoxIterator boxit2 = boxit1;
    ++boxit2;
    if ((&fabA(*boxit2, 0) - &fabA(*boxit1, 0)) != 1) ++status;
    if (verbose)
      {
        std::cout << "Stride0 contig?: "
                  << (&fabA(*boxit2, 0) - &fabA(*boxit1, 0)) << std::endl;
      }
    if (g_SpaceDim >= 2)
      {
        BoxIterator boxit3 = boxit1;
        boxit3 += IntVect(D_DECL(0, 1, 0));
        if ((&fabA(*boxit3, 0) - &fabA(*boxit1, 0)) != 3) ++status;
        if (verbose)
          {
            std::cout << "Stride1?       : "
                      << (&fabA(*boxit3, 0) - &fabA(*boxit1, 0)) << std::endl;
          }
      }
    if (g_SpaceDim >= 3)
      {
        BoxIterator boxit4 = boxit1;
        boxit4 += IntVect(D_DECL(0, 0, 1));
        if ((&fabA(*boxit4, 0) - &fabA(*boxit1, 0)) != 9) ++ status;
        if (verbose)
          {
            std::cout << "Stride2?       : "
                      << (&fabA(*boxit4, 0) - &fabA(*boxit1, 0)) << std::endl;
          }
      }
  }

  // Test assignment to the fab using iterators
  {
    if (verbose)
      {
        std::cout << "--- Test assignment to the fab using iterators\n";
      }
    BaseFab<Real> fabA(boxA, 1, 0.5);
    for (BoxIterator boxit(boxA); boxit.ok(); ++boxit)
      {
        fabA(*boxit, 0) = 1.5;
      }
    Box boxB(IntVect::Zero, IntVect(D_DECL(1, 1, 1)));
    // STL-like iterator usage.  The syntax above is preferred.
    const Box::const_iterator boxitEnd = boxB.end();
    for (Box::const_iterator boxit(boxB.begin()); boxit != boxitEnd; ++boxit)
      {
        fabA(*boxit, 0) = -1.5;
      }
    // Using loops instead of an iterator.  Fast but not as concise.  More error
    // prone.  Here is were we test the above assignments.
    IntVect iv;
    D_INVTERM(
      for (int i = 0; i <= 2; ++i)
        {
          iv[0] = i;,
          for (int j = 0; j <= 2; ++j)
            {
              iv[1] = j;,
              for (int k = 0; k <= 2; ++k)
                {
                  iv[2] = k;)
                  const Real sample = fabA(iv, 0);
                  if (sample != (boxB.contains(iv) ? -1.5 : 1.5)) ++status;
                  if (verbose)
                    {
                      std::cout << iv << ' ' << std::setw(4) << sample
                                << std::endl;
                    }
    D_TERM(},},})
    if (verbose)
      {
        std::cout << std::endl;
      }
  }
#endif

//--Output status

  if (verbose)
    {
      std::cout << "Status: " << status << std::endl;
    }
  const char* const testName = "testBoxIterator";
  const char* const statLbl[] = {
    "failed",
    "passed"
  };
  std::cout << std::left << std::setw(40) << testName
            << statLbl[(status == 0)] << std::endl;
  return status;
}
