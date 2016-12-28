
#ifndef _VEXTYPES_H_
#define _VEXTYPES_H_


/******************************************************************************/
/**
 * \file VEXTypes.H
 *
 * \brief Types for vector extension operations
 *
 *//*+*************************************************************************/


/*******************************************************************************
 *
 * NOTATION
 *
 * Something more similar to GNU notation is used for labeling the typs and
 * definitions defined herein.
 *
 * Intel notation
 * ==============
 *
 *   Floating point functions:
 *     _mm_add_xy()
 *       x : s - single value
 *         : p - packed (or vector)
 *       y : f - 32-bit floating point
 *         : d - 64-bit floating point
 *
 *   Variables:
 *     _m128i - integer vector (signed or unsigned)
 *     _m128  - 32-bit floating point vector
 *     _m128d - 64-bit floating point vector
 *     _m256...
 *
 * GNU notation
 * ============
 *
 *   Floating point functions are similar to Intel
 *
 *   Variables:
 *     vxyz
 *       x : 0-9 - number of scalars in a vector
 *       y :   q -  8-bit
 *         :   h - 16-bit
 *         :   s - 32-bit
 *         :   d - 64-bit
 *       z :   i - integer (signed or unsigned)
 *         :   f - floating point
 *
 * Chombo generalized notation
 * ===========================
 *
 *   v  - (vector) the size of the register on the CPU (i.e., either 128 or 256
 *        bits)
 *   r  - (real) either sf or df (see y for GNU notation) depending on if the
 *        code is compiled for 32-bit or 64-bit floating point numbers
 *
 *   In generalized functions:
 *     _mm_yz(add) - y and z as defined for GNU notation.  All function
 *       generalizations assume packed vectors.  For scalar equivalents, you
 *       must explicitly call the appropriate vector intrinsic.
 *     _mm_r(add)  - Using r to define sf or df
 *
 *   In generalized variables:
 *     __mvr   - a floating-point type, as defined by the build system, that
 *               fits in the vector register
 *     CHvr_t  - union of __mvr to allow access to components
 *
 *     Integers will not be generalized until AVX2 is commonplace
 *   
 ******************************************************************************/

#include <stdint.h>

#include "Config.H"

#ifdef CHDEF_SYSTEM_X86VECEXT_INTEL_H
#include CHDEF_SYSTEM_X86VECEXT_INTEL_H
#elif (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
#include CHDEF_SYSTEM_X86VECEXT_AVX_H
#elif (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_SSE2)
#include CHDEF_SYSTEM_X86VECEXT_SSE2_H
#endif

/*--------------------------------------------------------------------*
 * SSE2
 *--------------------------------------------------------------------*/

#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_SSE2)

//--128 bit

union CH128i32_t
{
  __m128i  m;
  uint32_t u[4];
  int32_t  i[4];
};

union CH128i64_t
{
  __m128i  m;
  uint64_t u[2];
  int64_t  i[2];
};

union CH128sf_t
{
  __m128 m;
  float  f[4];
};

union CH128df_t
{
  __m128d m;
  double  f[2];
};

#ifdef USE_SINGLE_PRECISION
typedef CH128i32_t CH128i_t;  // Same size int as float
typedef CH128sf_t CH128r_t;
typedef __m128 __m128r;
#else
typedef CH128i64_t CH128i_t;  // Same size int as float
typedef CH128df_t CH128r_t;
typedef __m128d __m128r;
#endif
#endif

/*--------------------------------------------------------------------*
 * AVX
 *--------------------------------------------------------------------*/

#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)

//--256 bit

union CH256i32_t
{
  __m256i  m;
  uint32_t u[8];
  int32_t  i[8];
};

union CH256i64_t
{
  __m256i  m;
  uint64_t u[4];
  int64_t  i[4];
};

union CH256sf_t
{
  __m256 m;
  float  f[8];
};

union CH256df_t
{
  __m256d m;
  double  f[4];
};

#ifdef USE_SINGLE_PRECISION
typedef CH256i32_t CH256i_t;  // Same size int as float
typedef CH256sf_t CH256r_t;
typedef __m256 __m256r;
#else
typedef CH256i64_t CH256i_t;  // Same size int as float
typedef CH256df_t CH256r_t;
typedef __m256d __m256r;
#endif
#endif


/*==============================================================================
 *
 * Now define types independent of register size
 *
 *                              *** WARNING *** 
 *
 * Integer types are *not* generalized because arithmetic is not supported
 * unless AVX2 is available.  I suggest using __m128i everywhere or using
 * defines in the code until AVX2 is commonplace.
 *
 *============================================================================*/

//--256 bit registers

#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)

// Alignement requirements for aligned loads/stores
#define CH_VECLS_ALIGN 32

#define VecSz_s 8
#define VecSz_d 4

#ifdef USE_SINGLE_PRECISION
#define VecSz_r 8
#else
#define VecSz_r 4
#endif

typedef CH256i_t CHvi_t;  // Same size int as float
typedef CH256r_t CHvr_t;
typedef __m256r __mvr;
typedef __m256i __mvi;

//--128 bit registers

#elif (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_SSE2)

// Alignement requirements for aligned loads/stores
#define CH_VECLS_ALIGN 16

#define VecSz_s 4
#define VecSz_d 2

#ifdef USE_SINGLE_PRECISION
#define VecSz_r 4
#else
#define VecSz_r 2
#endif

typedef CH128i_t CHvi_t;  // Same size int as float
typedef CH128r_t CHvr_t;
typedef __m128r __mvr;
typedef __m128i __mvi;

#endif


/*==============================================================================
 *
 * Name macros
 *
 *                              *** WARNING *** 
 *
 * Integer types are *not* generalized because arithmetic is not supported
 * unless AVX2 is available.  I suggest using __m128i everywhere or using
 * defines in the code until AVX2 is commonplace.
 *
 *============================================================================*/

#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)

#define _mm_i8(x) _mm256_ ## x ## _epi8
#define _mm_i16(x) _mm256_ ## x ## _epi16
#define _mm_i32(x) _mm256_ ## x ## _epi32
#define _mm_i64(x) _mm256_ ## x ## _epi64
#define _mm_si(x) _mm256_ ## x ## _si256

#define _mm_sf(x) _mm256_ ## x ## _ps
#define _mm_df(x) _mm256_ ## x ## _pd

#ifdef USE_SINGLE_PRECISION
#define _mm_vr(x) _mm256_ ## x ## _ps
#else  // DOUBLE
#define _mm_vr(x) _mm256_ ## x ## _pd
#endif

#else

#define _mm_i8(x) _mm_ ## x ## _epi8
#define _mm_i16(x) _mm_ ## x ## _epi16
#define _mm_i32(x) _mm_ ## x ## _epi32
#define _mm_i64(x) _mm_ ## x ## _epi64
#define _mm_si(x) _mm_ ## x ## _si128

#define _mm_sf(x) _mm_ ## x ## _ps
#define _mm_df(x) _mm_ ## x ## _pd

#ifdef USE_SINGLE_PRECISION
#define _mm_vr(x) _mm_ ## x ## _ps
#else  // DOUBLE
#define _mm_vr(x) _mm_ ## x ## _pd
#endif
#endif

#endif  /* ! defined _VEXTYPES_H_ */
