
#ifndef _BASEFABMACROS_H_
#define _BASEFABMACROS_H_

#include <cassert>

#include "Parameters.H"

// Expands macros and converts to a string (used internally)
#define STRINGIFY(x) #x

/*--------------------------------------------------------------------*
 *  Macro to generate a pointer to VLA from a BaseFab.  'x' is the
 *  name of the multi-dimensional array
 *  Example:
 *    MD_ARRAY(arrA, fabA);
 *  Notes:
 *    - the assert is only to work around what appears to be a
 *      compiler bug in gcc
 *--------------------------------------------------------------------*/

#define MD_ARRAY(x, _fab)                                               \
  D_TERM(                                                               \
    const int _ ## x ## n0 = (_fab).box().dimensions()[0];,             \
    const int _ ## x ## n1 = (_fab).box().dimensions()[1];,             \
    const int _ ## x ## n2 = (_fab).box().dimensions()[2];)             \
  using x ## _value_t = std::conditional_t<                             \
    std::is_const<std::remove_reference_t<decltype(_fab)> >::value,     \
    std::add_const_t<typename std::decay_t<decltype(_fab)>::value_type>, \
    typename std::decay_t<decltype(_fab)>::value_type>;                 \
  x ## _value_t *const _ ## x ## dataPtr =                              \
    ((_fab).dataPtr() - (((_fab).box().loVect()*(_fab).getStride()).sum())); \
  auto x =                                                              \
    (x ## _value_t (*)                                                  \
     D_INVTERM([_ ## x ## n0],[_ ## x ## n1],[_ ## x ## n2]))           \
    (_ ## x ## dataPtr);                                                \
  assert(&((_fab).operator()((_fab).box().hiVect(), (_fab).ncomp()-1)) == \
         &(         x[(_fab).ncomp()-1]                                 \
           D_INVTERM([(_fab).box().hiVect()[0]],                        \
                     [(_fab).box().hiVect()[1]],                        \
                     [(_fab).box().hiVect()[2]])));                     \
  (void)x

/*--------------------------------------------------------------------*
 *  Macro to generate a pointer to VLA from a BaseFab.  The pointer is
 *  annotated with the restrict qualifier.  'x' is the name of the
 *  multi-dimensional array
 *  Example:
 *    MD_ARRAY(arrA, fabA);
 *--------------------------------------------------------------------*/

#define MD_ARRAY_RESTRICT(x, _fab)                                      \
  D_TERM(                                                               \
    const int _ ## x ## n0 = (_fab).box().dimensions()[0];,             \
    const int _ ## x ## n1 = (_fab).box().dimensions()[1];,             \
    const int _ ## x ## n2 = (_fab).box().dimensions()[2];)             \
  using x ## _value_t = std::conditional_t<                             \
    std::is_const<std::remove_reference_t<decltype(_fab)> >::value,     \
    std::add_const_t<typename std::decay_t<decltype(_fab)>::value_type>, \
    typename std::decay_t<decltype(_fab)>::value_type>;                 \
  x ## _value_t *const _ ## x ## dataPtr =                              \
    ((_fab).dataPtr() - (((_fab).box().loVect()*(_fab).getStride()).sum())); \
  auto x =                                                              \
    (x ## _value_t (*__restrict__)                                      \
     D_INVTERM([_ ## x ## n0],[_ ## x ## n1],[_ ## x ## n2]))           \
    (_ ## x ## dataPtr);                                                \
  assert(&((_fab).operator()((_fab).box().hiVect(), (_fab).ncomp()-1)) == \
         &(         x[(_fab).ncomp()-1]                                 \
           D_INVTERM([(_fab).box().hiVect()[0]],                        \
                     [(_fab).box().hiVect()[1]],                        \
                     [(_fab).box().hiVect()[2]])));                     \
  (void)x

/*--------------------------------------------------------------------*
 *  Macro to rebuild a pointer to VLA from constituents which may be
 *  captured via lambda (including the type).
 *  Example:
 *    MD_ARRAY(arrA, fabA);
 *    auto rhs = [=](MD_DECLIX(int, o))
 *      {
 *        MD_CAPTURE(arrA, fabA);
 *      }
 *--------------------------------------------------------------------*/

#define MD_CAPTURE(x)                                                   \
  auto x =                                                              \
    (x ## _value_t (*)                                                  \
     D_INVTERM([_ ## x ## n0],[_ ## x ## n1],[_ ## x ## n2]))           \
    (_ ## x ## dataPtr)

/*--------------------------------------------------------------------*
 *  Macro to rebuild a pointer to VLA from constituents which may be
 *  captured via lambda (including the type).  The pointer is
 *  annotated with the restrict qualifier.
 *  Example:
 *    MD_ARRAY_RESTRICT(arrA, fabA);
 *    auto rhs = [=](MD_DECLIX(int, o))
 *      {
 *        MD_CAPTURE_RESTRICT(arrA);
 *      }
 *--------------------------------------------------------------------*/

#define MD_CAPTURE_RESTRICT(x)                                          \
  auto x =                                                              \
    (x ## _value_t (*__restrict__)                                      \
     D_INVTERM([_ ## x ## n0],[_ ## x ## n1],[_ ## x ## n2]))           \
    (_ ## x ## dataPtr)

/*--------------------------------------------------------------------*
 *  Macro to declare a multidimensional index (useful in callee)
 *--------------------------------------------------------------------*/

#define MD_DECLIX(T, x)                                                 \
  D_DECL(T x##0, T x##1, T x##2)

/*--------------------------------------------------------------------*
 *  Macro to expand a multidimensional index (useful in caller)
 *--------------------------------------------------------------------*/

#define MD_EXPANDIX(x)                                                  \
  D_DECL(x##0, x##1, x##2)

/*--------------------------------------------------------------------*
 *  Macro to generate a nested loop from a box.  The multidimensional
 *  index has values x0, x1, x2, etc.
 *  Example:
 *    MD_BOXLOOP(box, i)
 *      {
 *        std::cout << IntVect(D_DECL(i0, i1, i2)) << std::endl;
 *      }
 *--------------------------------------------------------------------*/

#define MD_BOXLOOP(_box, x)                                             \
  D_INVTERM(                                                            \
    for (int x ## 0 = (_box).loVect()[0]; x ## 0 <= (_box).hiVect()[0]; ++ x ## 0), \
    for (int x ## 1 = (_box).loVect()[1]; x ## 1 <= (_box).hiVect()[1]; ++ x ## 1), \
    for (int x ## 2 = (_box).loVect()[2]; x ## 2 <= (_box).hiVect()[2]; ++ x ## 2))

/*--------------------------------------------------------------------*
 *  Macro to generate a nested loop from a box for all but the first
 *  dimension.  The multidimensional index has values x1, x2, etc.
 *  Example:
 *    MD_BOXLOOP_PENCIL(box, i);
 *      {
 *        for (int i0 = box.smallEnd(0); i0 <= box.bigEnd(0); ++i0)
 *          {
 *            std::cout << IntVect(D_DECL(i0, i1, i2)) << std::endl;
 *          }
 *      }
 *--------------------------------------------------------------------*/

#define MD_BOXLOOP_PENCIL(_box, x)                                      \
  D_INVTERMPENCIL(                                                      \
    (void)0;,                                                           \
    for (int x ## 1 = (_box).loVect()[1]; x ## 1 <= (_box).hiVect()[1]; ++ x ## 1), \
    for (int x ## 2 = (_box).loVect()[2]; x ## 2 <= (_box).hiVect()[2]; ++ x ## 2))

/*--------------------------------------------------------------------*
 *  Macro to generate a nested loop from a box.  The multidimensional
 *  index has values x0, x1, x2, etc.  The outermost loop is
 *  parallelized using OpenMP.
 *  IMPORTANT: The outermost index (e.g., i2) is defined in the
 *  encompassing scope.  If there are conflicts, you will need to
 *  wrap the entire loop in braces as in the example below.
 *  Example:
 *    {
 *      MD_BOXLOOP_OMP(box, i)
 *        {
 *          std::cout << IntVect(D_DECL(i0, i1, i2)) << std::endl;
 *        }
 *    }
 *--------------------------------------------------------------------*/

#define MD_BOXLOOP_OMP(_box, x)                                         \
    int D_SELECT(x ## 0, x ## 1, x ## 2);                               \
    _Pragma( STRINGIFY(omp parallel for default(shared) private(D_SELECT(x ## 0, x ## 1, x ## 2))) ) \
    D_INVTERM(                                                          \
      for (D_SELECT(, int, int) x ## 0 = (_box).loVect()[0]; x ## 0 <= (_box).hiVect()[0]; ++ x ## 0), \
      for (D_SELECT(int, , int) x ## 1 = (_box).loVect()[1]; x ## 1 <= (_box).hiVect()[1]; ++ x ## 1), \
      for (D_SELECT(int, int, ) x ## 2 = (_box).loVect()[2]; x ## 2 <= (_box).hiVect()[2]; ++ x ## 2))

/*--------------------------------------------------------------------*
 *  Macro to generate a nested loop from a box for all but the first
 *  dimension.  The multidimensional index has values x0, x1, x2, etc.
 *  The outermost loop is parallelized using OpenMP.
 *  IMPORTANT: The outermost index (e.g., i2) is defined in the
 *  encompassing scope.  If there are conflicts, you will need to
 *  wrap the entire loop in braces as in the example below.
 *  Example:
 *    {
 *      MD_BOXLOOP_PENCIL_OMP(box, i)
 *        {
 *          for (int i0 = box.smallEnd(0); i0 <= box.bigEnd(0); ++i0)
 *            {
 *              std::cout << IntVect(D_DECL(i0, i1, i2)) << std::endl;
 *            }
 *        }
 *    }
 *--------------------------------------------------------------------*/

#define MD_BOXLOOP_PENCIL_OMP(_box, x)                                  \
    int D_SELECT(x ## 0, x ## 1, x ## 2);                               \
    _Pragma( STRINGIFY(omp parallel for default(shared) private(D_SELECT(x ## 0, x ## 1, x ## 2))) ) \
    D_INVTERMPENCIL(                                                    \
      (void)0;,                                                         \
      for (D_SELECT(int, , int) x ## 1 = (_box).loVect()[1]; x ## 1 <= (_box).hiVect()[1]; ++ x ## 1), \
      for (D_SELECT(int, int, ) x ## 2 = (_box).loVect()[2]; x ## 2 <= (_box).hiVect()[2]; ++ x ## 2))

/*--------------------------------------------------------------------*
 *  Macro to index an array based on a multidimensional index 'x', and
 *  a component index
 *  Example:
 *    MD_BOXLOOP(box, i);
 *      {
 *        arrA[MD_IX(i, 0)] = arrB[MD_IX(i, 1)]
 *      }
 *--------------------------------------------------------------------*/

#define MD_IX(x,c)                                                      \
  (c)] D_INVTERM([ (x ## 0) , [ (x ## 1) ], [ (x ## 2) ])

/*--------------------------------------------------------------------*
 *  Macro to generate a unit vector of components pointing to a
 *  direction.  Most often used with MD_OFFSETIX.
 *  Example:
 *    for (int dir = 0; dir != g_SpaceDim; ++dir)
 *      {
 *        const int MD_ID(o, dir);
 *        // If dir = 0, o0 = 1, o1 = 0, o2 = 0.
 *        // If dir = 1, o0 = 0, o1 = 1, o2 = 0.
 *        // etc
 *      }
 *--------------------------------------------------------------------*/

#define MD_ID(x,_dir)                                                   \
  D_DECL(x ## 0 = ((_dir) == 0),                                        \
         x ## 1 = ((_dir) == 1),                                        \
         x ## 2 = ((_dir) == 2))

/*--------------------------------------------------------------------*
 *  Macro to generate an offset index to an array based on a
 *  multidimensional index 'x', operator 'op', offset index 'o', and
 *  a component index
 *  Example:
 *    for (int dir = 0; dir != g_SpaceDim; ++dir)
 *      {
 *        const int MD_ID(o, dir);
 *        MD_BOXLOOP(box, i);
 *          {
 *            arrA[MD_IX(i, 0)] = arrB[MD_OFFSETIX(i,+,o, 1)]
 *          }
 *      }
 *--------------------------------------------------------------------*/

#define MD_OFFSETIX(x,op,o,c)                                           \
  (c)] D_INVTERM([ (x ## 0 op o ## 0) ,                                 \
                 [ (x ## 1 op o ## 1) ] ,                               \
                 [ (x ## 2 op o ## 2) ])

/*--------------------------------------------------------------------*
 *  Same as above but applies to an IntVect offset 'ov'.
 *  Example:
 *    for (int dir = 0; dir != g_SpaceDim; ++dir)
 *      {
 *        IntVect iv(IntVect::Zero);
 *        iv[dir] = 1;
 *        MD_BOXLOOP(box, i);
 *          {
 *            arrA[MD_IX(i, 0)] = arrB[MD_OFFSETIX(i,+,iv, 1)]
 *          }
 *      }
 *--------------------------------------------------------------------*/

#define MD_OFFSETIV(x,op,ov,c)                                         \
  (c)] D_INVTERM([ (x ## 0 op ov[0]) ,                                 \
                 [ (x ## 1 op ov[1]) ] ,                               \
                 [ (x ## 2 op ov[2]) ])

/*--------------------------------------------------------------------*
 * Template metaprogramming of time-loop.  'DIR' and ii are passed
 * to function F.  Helper templates follow that allow type deduction
 * for F
 *--------------------------------------------------------------------*/

//--General class

template <typename F, int Dir>
class MD_DirLoopFunc
{
public:
  static auto sum(F f)
    {
      constexpr int MD_ID(ii, Dir);
      return f(Dir, MD_EXPANDIX(ii)) + MD_DirLoopFunc<F, Dir-1>::sum(f);
    }
};

//--Specialization for Dir == 0

template <typename F>
class MD_DirLoopFunc<F, 0>
{
public:
  static auto sum(F f)
    {
      constexpr int MD_ID(ii, 0);
      return f(0, MD_EXPANDIX(ii));
    }
};

//--Helper function templates

template <typename F>
inline auto MD_DIRSUM(F f)
{
  return MD_DirLoopFunc<F, g_SpaceDim>::sum(f);
}

/*--------------------------------------------------------------------*
 *  Macro to efficiently iterate over a box, 'b'.
 *  Inside the body of this macro, jR, and kR, are j and k positions
 *  inside the box.  The "pencil" in the i direction starts at iR
 *  and ends at iEndR-1.
 *--------------------------------------------------------------------*/

#define ForAllBPencil(b)                                                \
{                                                                       \
  const IntVect _b_len = (b).dimensions();                              \
  const int iEndR = (b).hiVect()[0] + 1;                                \
  D_TERM(                                                               \
    const int iR = (b).loVect()[0];,                                    \
    const int _j_lo = (b).loVect()[1];,                                 \
    const int _k_lo = (b).loVect()[2];)                                 \
  D_INVTERM(                                                            \
    (void)0;,                                                           \
    for(int _j = 0; _j != _b_len[1]; ++_j) {                            \
      const int jR = _j_lo + _j;,                                       \
      for(int _k = 0; _k != _b_len[2]; ++_k) {                          \
        const int kR = _k_lo + _k;)

/*--------------------------------------------------------------------*
 *  Macro to efficiently iterate over a box, 'bx', and box 'by'.
 *  Both boxes must have the same dimensions.  Inside the body of this
 *  macro, jbxR, and kbxR, are j and k positions inside box 'bx'.
 *  jbyR, and kbyR, are j and k positions inside box 'by'.  The
 *  "pencil" in the i direction starts at ibxR and ends at ibxEndR-1
 *  for box 'bx' and likewise for box 'by'.
 *--------------------------------------------------------------------*/

#define ForAllBBPencil(bx,by)                                           \
{                                                                       \
  const IntVect _b_len = (bx).dimensions();                             \
  const int i##bx##EndR = (bx).hiVect()[0] + 1;                         \
  const int i##by##EndR = (by).hiVect()[0] + 1;                         \
  const IntVect _bx_lo = (bx).loVect();                                 \
  const IntVect _by_lo = (by).loVect();                                 \
  const int i##bx##R = _bx_lo[0];                                       \
  const int i##by##R = _by_lo[0];                                       \
  D_INVTERM(                                                            \
    (void)0;,                                                           \
    for(int _j = 0; _j != _b_len[1]; ++_j) {                            \
      const int j##bx##R = _bx_lo[1] + _j;                              \
      const int j##by##R = _by_lo[1] + _j;,                             \
      for(int _k = 0; _k != _b_len[2]; ++_k) {                          \
      const int k##bx##R = _bx_lo[2] + _k;                              \
      const int k##by##R = _by_lo[2] + _k;)

/*--------------------------------------------------------------------*
 *  Adds closing braces for pencil macros
 *--------------------------------------------------------------------*/

#define EndForPencil D_TERM(},},})

#endif  /* ! defined _BASEFABMACROS_H_ */
