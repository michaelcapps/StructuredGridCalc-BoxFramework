#include <iostream>
#include <iomanip>
#include <vector>

#include "BaseFab.H"


/*==============================================================================
 * How is a fortran function named
 *============================================================================*/

#define CH_FORT_UNDERSCORE

#ifdef CH_FORT_UPPERCASE
  #ifdef CH_FORT_UNDERSCORE
    #define FORTRAN_NAME( NAME ,name ) NAME ## _
  #else
    #define FORTRAN_NAME( NAME ,name ) NAME
  #endif
#else
  #ifdef CH_FORT_UNDERSCORE
    #define FORTRAN_NAME( NAME ,name ) name ## _
  #else
    #define FORTRAN_NAME( NAME ,name ) name
  #endif
#endif


/*==============================================================================
 * BLAS & LAPACK routines
 *============================================================================*/

#define CH_USE_DOUBLE

#ifdef CH_USE_DOUBLE
  #define BLAS2_GEMV FORTRAN_NAME(DGEMV, dgemv)
  #define BLAS3_GEMM FORTRAN_NAME(DGEMM, dgemm)
  #define LAPACK_GETRF FORTRAN_NAME(DGETRF, dgetrf)
  #define LAPACK_GETRI FORTRAN_NAME(DGETRI, dgetri)
#else
  #define BLAS2_GEMV FORTRAN_NAME(SGEMV, sgemv)
  #define BLAS3_GEMM FORTRAN_NAME(SGEMM, sgemm)
  #define LAPACK_GETRF FORTRAN_NAME(SGETRF, sgetrf)
  #define LAPACK_GETRI FORTRAN_NAME(SGETRI, sgetri)
#endif

extern "C" void BLAS2_GEMV(char*, int*, int*, Real*, Real*, int*, Real*, int*,
                           Real*, Real*, int*);
extern "C" void BLAS3_GEMM(char*, char*, int*, int*, int*, Real*, Real*, int*,
                           Real*, int*, Real*, Real*, int*);
extern "C" void LAPACK_GETRF(int*, int*, Real*, int*, int*, int*);
extern "C" void LAPACK_GETRI(int*, Real*, int*, int*, Real*, int*, int*);


/*==============================================================================
 * Matrix macros
 *============================================================================*/

// Set matrix size in an FArrayBox matrix
#define MSZ(x, y) Box(IntVect::Zero, IntVect(D_DECL((x)-1, (y)-1, 0))), 1

// Index an element in an FArrayBox matrix
#define MIX(x, y) IntVect(D_DECL((x), (y), 0)), 0

// Matrix
using Matrix = FArrayBox;

// Vector
using Vector = std::vector<Real>;


/*==============================================================================
 * Functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
/// Pretty print an FArrayBox matrix
/** \param[in]  a_os    Output stream
 *  \param[in]  a_matrix
 *                      Matrix to print
 *  \return             Output stream with matrix printed
 *//*-----------------------------------------------------------------*/

inline std::ostream &operator<<(std::ostream& a_os, const Matrix& a_matrix)
{
  const int prec = 2;
  a_os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  a_os.precision(prec);
  const int iMax = a_matrix.box().dimensions()[0];
  const int jMax = a_matrix.box().dimensions()[1];
  for (int i = 0; i != iMax; ++i)
    {
      if (jMax > 0)
        {
          a_os << std::setw(prec+7) << a_matrix(MIX(i, 0));
          for (int j = 1; j != jMax; ++j)
            {
              a_os << ' ' << std::setw(prec+7) << a_matrix(MIX(i, j));
            }
          a_os << std::endl;
        }
    }
  a_os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
  a_os.precision(6);
  return a_os;
}

/*--------------------------------------------------------------------*/
/// Pretty print a vector
/** \param[in]  a_os    Output stream
 *  \param[in]  a_vec   Vector to print
 *  \return             Output stream with matrix printed
 *//*-----------------------------------------------------------------*/

inline std::ostream &operator<<(std::ostream& a_os, const Vector& a_vec)
{
  const int prec = 2;
  a_os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  a_os.precision(prec);
  const int iMax = a_vec.size();
  for (int i = 0; i != iMax; ++i)
    {
      a_os << std::setw(prec+7) << a_vec[i] << std::endl;
    }
  a_os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
  a_os.precision(6);
  return a_os;
}

/*--------------------------------------------------------------------*/
/// Invert a matrix using Lapack
/** Work and pivot arrays are created using VLA
 *  \param[in]  a_A     Matrix to invert
 *  \param[out] a_A     Inverted matrix
 *  \param[in]  a_lwork Size of work array.  Use appropriate value or
 *                      set to -1 to let Lapack determine value.
 *//*-----------------------------------------------------------------*/

int inverse(Matrix& a_A, int a_lwork)
{
  int ier;
  int lfi_N = a_A.box().dimensions()[0];
  CH_assert(a_A.box().dimensions()[1] >= lfi_N);
  int lfi_LDA = lfi_N;

  // Setup pivots
  int ipiv[lfi_N];

  // Setup work
  if (a_lwork < lfi_N)
    {
      a_lwork = -1;
      Real query;
      // Get the optimal work array size
      LAPACK_GETRI(&lfi_N, a_A.dataPtr(), &lfi_LDA, ipiv, &query, &a_lwork,
                   &ier);
      if (ier != 0)
        {
          std::cout << "GETRI failure status: " << ier << std::endl;
          return 0;
        }
      a_lwork = (int)query;
    }
  Real work[a_lwork];

  LAPACK_GETRF(&lfi_N, &lfi_N, a_A.dataPtr(), &lfi_LDA, ipiv, &ier);
  if (ier != 0)
    {
      std::cout << "GETRF failure status: " << ier << std::endl;
      return 0;
    }

  LAPACK_GETRI(&lfi_N, a_A.dataPtr(), &lfi_LDA, ipiv, work, &a_lwork, &ier);
  if (ier != 0)
    {
      std::cout << "GETRI failure status: " << ier << std::endl;
      return 0;
    }
  return 1;
}

/*--------------------------------------------------------------------*/
/// multiply two matrices
/** 
 *  \param[in]  a_A     left matrix to multiply
 *  \param[in]  a_B     right matrix to multiply
 *  \param[in]  a_C     result of mutiplication (C = A*B)
 *//*-----------------------------------------------------------------*/

int gemm(Matrix& matA, Matrix& matB, Matrix& matC)
{
  char n = 'n';
  int M = matA.box().hiVect(0) + 1; //num rows in A
  int N = matB.box().hiVect(1) + 1; //num cols in B
  int K = matA.box().hiVect(1) + 1; //num cols in A (same as num rows in B)
  double alpha = 1;
  double beta = 0;
  int LDA = std::max(1,M);
  int LDB = std::max(1,K);
  int LDC = std::max(1,M);

  BLAS3_GEMM(&n,&n,&M,&N,&K,&alpha,matA.dataPtr(),&LDA,matB.dataPtr(),&LDB,&beta,matC.dataPtr(),&LDC);  
  return 1;
}

/*--------------------------------------------------------------------*/
/// multiply a matrix and a vetor
/** 
 *  \param[in]  a_A     left matrix to multiply
 *  \param[in]  a_B     right matrix to multiply
 *  \param[in]  a_C     result of mutiplication (C = A*B)
 *//*-----------------------------------------------------------------*/

int gemv(Matrix& matA,Vector& vecA,Vector& vecB)
{
  char n = 'n';
  int M = matA.box().hiVect(0) + 1; //num rows in A
  int N = matA.box().hiVect(1) + 1; //num cols in A
  double alpha = 1;
  int LDA = std::max(1,M);
  int INCX = 1;
  double beta = 0;
  int INCY = 1;

  BLAS2_GEMV(&n,&M,&N,&alpha,matA.dataPtr(),&LDA,&(*vecA.begin()),&INCX,&beta,&(*vecB.begin()),&INCY);
  return 1;
}

/*--------------------------------------------------------------------*/
/// Main routine
/** Construct a FArrayBox matrix and test inversion, matrix multiply,
 *  and matrix vector multiple.  Note that
 *  Matrix is type FArrayBox
 *  Vector is type std::vector<Real>
 *//*-----------------------------------------------------------------*/

int main()
{
  Matrix matA(MSZ(3, 3));

  matA.setVal(1.0);
  matA(MIX(2, 0)) = -1.0;
  matA(MIX(1, 1)) =  2.0;
  matA(MIX(0, 2)) =  3.0;
  std::cout << "matrix:\n" << matA << std::endl;

  int lwork = -1;
  inverse(matA, lwork);  
  std::cout << "inverted matrix:\n" << matA << std::endl;

  Matrix matB(MSZ(3, 4));
  matB(MIX(0, 0)) = -14.0;
  matB(MIX(1, 0)) = -9.50;
  matB(MIX(2, 0)) = -5.0;
  matB(MIX(0, 1)) = -10.0;
  matB(MIX(1, 1)) =  0.0;
  matB(MIX(2, 1)) = -6.0;
  matB(MIX(0, 2)) =  11.0;
  matB(MIX(1, 2)) =  13.5;
  matB(MIX(2, 2)) =  2.0;
  matB(MIX(0, 3)) =  26.2;
  matB(MIX(1, 3)) =  20.8;
  matB(MIX(2, 3)) =  5.4;

  Matrix matC(MSZ(3, 4));
  gemm(matA, matB, matC);
  std::cout << "matmul:\n" << matC << std::endl;

  Vector vecA = { 1.5, -1.5, 2.0, 1.0 };
  std::cout << "vector:\n" << vecA << std::endl;
  Vector vecB(3);
  gemv(matC, vecA, vecB);
  std::cout << "matvecmul:\n" << vecB << std::endl;
}
