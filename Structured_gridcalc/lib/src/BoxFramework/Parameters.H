
#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_


/******************************************************************************/
/**
 * \file Parameters.H
 *
 * \brief Parameters for the box framework
 *
 *//*+*************************************************************************/

#include <cassert>

// Define the number of space dimensions
#ifdef SPACEDIM2
#define SPACEDIM 2
#endif

#ifndef SPACEDIM
#define SPACEDIM 3
#endif

// Number of space dimensions
const int g_SpaceDim = SPACEDIM;

/*------------------------------------------------------------------------------
 * SpaceDim Macros
 *----------------------------------------------------------------------------*/

#if ! defined(SPACEDIM)
#error SPACEDIM must be defined
#endif

#ifdef D_EXPR
#undef D_EXPR
#undef D_DECL
#undef D_TERM
#undef D_INVTERM
#undef D_INVTERMPENCIL
#endif

#if SPACEDIM==1
#  define D_EXPR(a,b,c) ((void)((a),0))
#  define D_DECL(a,b,c) a
#  define D_TERM(a,b,c) a
#  define D_INVTERM(a,b,c) a
#  define D_INVTERMPENCIL(a,b,c) (void)0;
#  define D_SELECT(a,b,c) a
#elif SPACEDIM==2
#  define D_EXPR(a,b,c) ((void)((a),(b),0))
#  define D_DECL(a,b,c) a,b
#  define D_TERM(a,b,c) a b
#  define D_INVTERM(a,b,c) b a
#  define D_INVTERMPENCIL(a,b,c) b
#  define D_SELECT(a,b,c) b
#elif SPACEDIM==3
/**
@ingroup macros

  The D_EXPR(a,b,c) macro expands to a comma expression that will
  evaluate SpaceDim of its arguments, counting from the left; i.e.
  CH_SPACEDIM==1 implies only `a' is evaluated, CH_SPACEDIM==2
  implies both `a' and `b', and CH_SPACEDIM==3 implies all three
  arguments are evaluated.  This macro can be used to enhance
  the portability of code that works with structures whose size
  depends on CH_SPACEDIM.

  For example:

     D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);

  More explicitly, this macro evaluates to one of

     ((void)((a),0))

     ((void)((a),(b),0))

     ((void)((a),(b),(c),0))

  corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#  define D_EXPR(a,b,c) ((void)((a),(b),(c),0))

/**
@ingroup macros

  The D_DECL(a,b,c) macro expands to a comma-separated list of
  1, 2, or all 3 of the arguments of the call, depending on the
  value of CH_SPACEDIM.  This can be used to write portable
  function calls that depend on 1, 2, or 3 arguments, corresponding
  to the value of CH_SPACEDIM.

  For example:

     return IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));

  More explicitly, this macro evaluates to one of

     a

     a,b

     a,b,c

  corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#  define D_DECL(a,b,c) a,b,c

/**
@ingroup macros
   define D_TERM(a,b,c) a b c

  More explicitly, this macro evaluates to one of

     a

     a b

     a b c

    corresponding to CH_SPACEDIM values of 1, 2, or 3.

    handy for situations like:

    Real norm = sqrt(D_TERM(u[0]*u[0],+u[1]*u[1],+u[2]*u[2]));
*/
#  define D_TERM(a,b,c) a b c

/**
@ingroup macros
   define D_INVTERM(a,b,c) c b a

  More explicitly, this macro evaluates to one of

     a

     b a

     c b a

    corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#  define D_INVTERM(a,b,c) c b a

/**
@ingroup macros
   define D_INVTERMPENCIL(a,b,c) c b

  More explicitly, this macro evaluates to one of

     (void)0;

     b

     c b

    corresponding to CH_SPACEDIM values of 1, 2, or 3.

    It is the same as D_INVTERM except it excludes the inner-
    most item.
*/
#  define D_INVTERMPENCIL(a,b,c) c b

/**
@ingroup macros
   define D_SELECT(a,b,c) a b c

  More explicitly, this macro evaluates to one of

     a

     b

     c

    corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#  define D_SELECT(a,b,c) c
#endif

/*------------------------------------------------------------------------------
 * Assertions
 *----------------------------------------------------------------------------*/

#ifdef RELEASE
  #define CH_assert(cond) (void)0
#else
  #define CH_assert(cond) assert(cond)
#endif

/*------------------------------------------------------------------------------
 * Precision
 *----------------------------------------------------------------------------*/

#ifdef USE_SINGLE_PRECISION
  // Single
  typedef float Real;
  #define CGNS_REAL RealSingle
  #define BXFR_MPI_REAL MPI_FLOAT 
#else
  // Double
  typedef double Real;
  #define CGNS_REAL RealDouble
  #define BXFR_MPI_REAL MPI_DOUBLE 
#endif

/*------------------------------------------------------------------------------
 * Define USE_STACK to build temporary FArrayBoxs on the stack
 *----------------------------------------------------------------------------*/

// #define USE_STACK

#endif  /* ! defined _PARAMETERS_H_ */
