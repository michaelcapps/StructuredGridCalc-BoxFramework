#include <cstring>
#include <iostream>
#include <iomanip>

#include "IntVect.H"

int main(const int argc, const char* argv[])
{
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
  int status = 0;

//--Tests

  // Test statics
  IntVect ivStatics = IntVect::Zero;
  if (D_TERM(   ivStatics[0] != 0,
             || ivStatics[1] != 0,
             || ivStatics[2] != 0)) ++status;
  if (verbose)
    {
      std::cout << "IntVect::Zero: " << ivStatics << std::endl;
    }
  ivStatics = IntVect::Unit;
  if (D_TERM(   ivStatics[0] != 1,
             || ivStatics[1] != 1,
             || ivStatics[2] != 1)) ++status;
  if (verbose)
    {
      std::cout << "IntVect::Unit: " << ivStatics << std::endl;
    }
  
  // Test construction
  IntVect ivA(D_DECL(0, 1, 2));
  // Test copy
  const IntVect ivB(ivA);
  // Test default construction
  IntVect ivC;
  if (D_TERM(   ivC[0] != 0,
             || ivC[1] != 0,
             || ivC[2] != 0)) ++status;
  ivC = ivB;
  if (D_TERM(   ivC[0] != 0,
             || ivC[1] != 1,
             || ivC[2] != 2)) ++status;
  if (verbose)
    {
      std::cout << "ivC: " << ivC << std::endl;
    }

  // Test const indexing
  if (ivB[0] != 0) ++status;
  // Test indexing
  ivC[1] = 3;
  if (ivC[1] != 3) ++status;

  // Test some operators
  if ((ivA == ivB) != true) ++status;
  if ((ivA != ivC) != true) ++status;
#if 1
  if ((ivA <= ivB) != true) ++status;
  // Test relational
  ivC = IntVect::Unit;
  IntVect ivD(IntVect::Zero);
  if ((ivD < ivC) != true) ++status;
  if ((ivC <= ivC) != true) ++status;
  if ((ivD <= ivC) != true) ++status;
  ivD[1] = 2;
  if ((ivD < ivC) != false) ++status;
  if ((ivD <= ivC) != false) ++status;
  // Test + with IntVect
  ivC = ivA + ivB;
  if (ivC != IntVect(D_DECL(0, 2, 4))) ++status;
  ivC[0] = 1;
  // Test less than
  if ((ivA <  ivC) != true) ++status;
  // Test add scalar
  ivC += 1;
  if (ivC != IntVect(D_DECL(2, 3, 5))) ++status;
#if SPACEDIM>2
  ivC[2] = -1;
#endif
  // Test max for each component
  ivC.max(ivA);
  if (ivC != IntVect(D_DECL(2, 3, 2))) ++status;
  // Test min for each component
  ivC.min(ivA);
  if (ivC != IntVect(D_DECL(0, 1, 2))) ++status;
#if SPACEDIM>2
  ivC[2] = -2;
#endif
  // Test norm
  if (ivC.norm1() != D_SELECT(0, 1, 3)) ++status;
  ivC[0] = -1;
  // Test sum
  if (ivC.sum() != D_SELECT(-1, 0, -2)) ++status; //PROBLEM?!
  // Test product
  if (ivC.product() != D_SELECT(-1, -1, 2)) ++status;
  // Here, ivC = (-1, 1,-2), ivA = ivB = (0, 1, 2)
  // Test add another IntVect
  ivC += ivA;
  if (ivC != IntVect(D_DECL(-1, 2, 0))) ++status;
  // Test subtract another IntVect
  ivC -= ivA;
  if (ivC != IntVect(D_DECL(-1, 1, -2))) ++status;
  // Test subtract scalar
  ivC -= 2;
  if (ivC != IntVect(D_DECL(-3, -1, -4))) ++status;
  // Test - with IntVect
  ivC = ivA - ivB;
  if (ivC != IntVect(D_DECL(0, 0, 0))) ++status;
  // Test * with IntVect
  ivC = ivA * ivB;
  if (ivC != IntVect(D_DECL(0, 1, 4))) ++status;
  // Test / with IntVect
  ivD = IntVect(D_DECL(2, 2, 4));
  IntVect ivE(D_DECL(2, 1, 2));
  ivC = ivD / ivE;
  if (ivC != IntVect(D_DECL(1, 2, 2))) ++status;
  // Test Unary operator
  ivC = -(ivD);
  if (ivC != IntVect(D_DECL(-2, -2, -4))) ++status;
  ivC = ivA + (-(ivD));
  if (ivC != IntVect(D_DECL(-2, -1, -2))) ++status;
#endif

  // Scaling
  ivC = IntVect(D_DECL(-1, 1, -2));
  if (ivC*2 != IntVect(D_DECL(-2, 2, -4))) ++status;
#if 1
  if (2*ivC != IntVect(D_DECL(-2, 2, -4))) ++status;
#endif

//--Output status

  if (verbose)
    {
      std::cout << "Status: " << status << std::endl;
    }
  const char* const testName = "testIntVect";
  const char* const statLbl[] = {
    "failed",
    "passed"
  };
  std::cout << std::left << std::setw(40) << testName
            << statLbl[(status == 0)] << std::endl;
  return status;
}
