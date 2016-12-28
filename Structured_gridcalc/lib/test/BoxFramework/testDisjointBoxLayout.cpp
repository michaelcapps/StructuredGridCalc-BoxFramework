#include <iostream>
#include <iomanip>
#include <cstring>

#include "DisjointBoxLayout.H"

int main(const int argc, const char* argv[])
{
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
  int status = 0;

//--Tests

//#if 0
  const Box domain(IntVect::Zero, 9*IntVect::Unit);

  // Should fit evenly
  DisjointBoxLayout dbl1(domain, 5*IntVect::Unit);
  // ... with this many boxes per dimension
  const IntVect boxDim = 2*IntVect::Unit;
  if (dbl1.size() != ((2*IntVect::Unit).product())) ++status;
  // Last box should have maxsize
  if (dbl1.getLinear(dbl1.size()-1).box.size() != ((5*IntVect::Unit).product()))
   ++status;
  //#if 0
  // Test the size of the first box (general box size testing is performed
  // in testLayoutIterator.cpp)
  {
    const Box& testBox = dbl1.getLinear(0).box;
    if (testBox.loVect() != IntVect::Zero) ++status;
    if (testBox.hiVect() != 4*IntVect::Unit) ++status;
  }
  //#if 0
  // Test the size of the last box
  {
    const Box& testBox = dbl1.getLinear(dbl1.size()-1).box;
    if (testBox.loVect() != 5*IntVect::Unit) ++status;
    if (testBox.hiVect() != 9*IntVect::Unit) ++status;
  }
  //#if 0
  // Misc
  {
    if (dbl1.problemDomain() != domain) ++status;
    if (dbl1.dimensions() != boxDim) ++status;
    if (dbl1.size() != boxDim.product()) ++status;
  }
  //#if 0
  // Test indexing with boxIndex
  {
    int linIdxBox = 0;
    D_INVTERM(
      for (int i = 0; i != 2; ++i)
        {,
          for (int j = 0; j != 2; ++j)
            {,
              for (int k = 0; k != 2; ++k)
                {)
                  IntVect lo(5*IntVect(D_DECL(i, j, k)));
                  IntVect hi(lo + 4*IntVect::Unit);
                  // Here we set the local index to 0 since the global index
                  // should be used for finding boxes
                  const Box& testBox = dbl1[BoxIndex(linIdxBox, 0)];
                  // Access through the linear index should give the same box
                  if (testBox != dbl1.getLinear(linIdxBox).box) ++status;
                  if (testBox.loVect() != lo) ++status;
                  if (testBox.hiVect() != hi) ++status;
                  ++linIdxBox;
    D_TERM(},},})
  }
  #if 1
//--Serial testing

  if (DisjointBoxLayout::numProc() == 1)
    {
      if (dbl1.localSize() != dbl1.size()) ++status;
      if (dbl1.localIdxBegin() != 0) ++status;
      if (dbl1.localIdxEnd() != dbl1.size()) ++status;
      for (int linIdxBox = 0; linIdxBox != boxDim.product(); ++linIdxBox)
        {
          if (dbl1.getLinear(linIdxBox).proc != 0) ++status;
        }
    }
#endif


//--Output status
  if (verbose)
    {
      std::cout << "Status: " << status << std::endl;
    }
  const char* const testName = "testDisjointBoxLayout";
  const char* const statLbl[] = {
    "failed",
    "passed"};
  std::cout << std::left << std::setw(40) << testName
            << statLbl[(status == 0)] << std::endl;
  return status;
}
