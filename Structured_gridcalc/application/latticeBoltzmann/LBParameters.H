
#ifndef _LBPARAMETERS_H_
#define _LBPARAMETERS_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Lattice-Boltzmann parameters (constant) for the run
 *
 *//*+*************************************************************************/

#include <cstdio>
#include <cmath>

#include "Parameters.H"
#include "IntVect.H"

namespace LBParameters
{

//--Absolute constants for any problem

constexpr int g_numVelDir = 19;       ///< Number of velocity directions
const int g_indexFaceVelBegin = 1;    ///< Start of velocity indexes for faces
const int g_indexFaceVelEnd   = 6;    ///< End of velocity indexes for faces
const int g_indexEdgeVelBegin = 7;    ///< Start of velocity indexes for edges
const int g_indexEdgeVelEnd   = 18;   ///< End of velocity indexes for edges
constexpr int g_numGhost = 1;         ///< Number of ghost cells
const IntVect g_ghostVect = g_numGhost*IntVect::Unit;
                                      ///< Number of ghost cells (this is
                                      ///< constant for all directions but some
                                      ///< routines require a vector)
constexpr int g_numState = 1 + g_SpaceDim;
                                      ///< Number of macroscopic conservative
                                      ///< state variables

constexpr Real g_weight[g_numVelDir] =
                                      ///< Weights defining the distribution
{ 
  1./3.,                                          // Cell
  1./18., 1./18., 1./18., 1./18., 1./18., 1./18., // Faces
  1./36., 1./36., 1./36., 1./36., 1./36., 1./36., // Edges
  1./36., 1./36., 1./36., 1./36., 1./36., 1./36.
};
constexpr int g_verbosity = 1;        ///< Amount of output

constexpr Real g_pi = 3.141592653589793;
constexpr Real g_tau = 0.516;
constexpr Real g_cs2 = 1.0/3;
constexpr Real g_bodyForce = 0.000001042;



/*==============================================================================
 *
 * Functions returning global constants
 *
 *============================================================================*/

/*--------------------------------------------------------------------*/
/// Get the velocity vector from a sequential index
/**  
 *  \param[in]  a_ei   Sequential velocity direction
 *  \return            Pointer to velocity vector
 *//*-----------------------------------------------------------------*/

inline const int* latticeVelocityP(const int a_ei)
{
  CH_assert(a_ei < g_numVelDir);
  static constexpr int idxNeighbour[g_numVelDir][3] =
    {
      {  0,  0,  0 },  // this    0
      { -1,  0,  0 },  // -x      1
      {  1,  0,  0 },  // +x      2
      {  0, -1,  0 },  // -y      3
      {  0,  1,  0 },  // +y      4
      {  0,  0, -1 },  // -z      5
      {  0,  0,  1 },  // +z      6
      { -1, -1,  0 },  // -x,-y   7
      {  1, -1,  0 },  // +x,-y   8
      { -1,  1,  0 },  // -x,+y   9
      {  1,  1,  0 },  // +x,+y  10
      { -1,  0, -1 },  // -x,-z  11
      {  1,  0, -1 },  // +x,-z  12
      { -1,  0,  1 },  // -x,+z  13
      {  1,  0,  1 },  // +x,+z  14
      {  0, -1, -1 },  // -y,-z  15
      {  0,  1, -1 },  // +y,-z  16
      {  0, -1,  1 },  // -y,+z  17
      {  0,  1,  1 }   // +y,+z  18
    };
  return idxNeighbour[a_ei];
}

/*--------------------------------------------------------------------*/
/// Get the velocity vector from a sequential index
/**  
 *  \param[in]  a_ei   Sequential velocity direction
 *  \return            Velocity vector as an IntVect
 *//*-----------------------------------------------------------------*/

inline IntVect latticeVelocity(const int a_ei)
{
  const int *const lVel = latticeVelocityP(a_ei);
  return IntVect(lVel[0], lVel[1], lVel[2]);
}

/*--------------------------------------------------------------------*/
/// Neighbour lattice site
/** \param[in]  a_iv   Current index
 *  \param[in]  a_ei   Sequential velocity direction
 *  \return            Neighbour index
 *//*-----------------------------------------------------------------*/

inline IntVect neighbourLatticeSite(const IntVect& a_iv, const int a_ei)
{
  const int *const lVel = latticeVelocityP(a_ei);
  return IntVect(a_iv[0] + lVel[0],
                 a_iv[1] + lVel[1],
                 a_iv[2] + lVel[2]);
}

/*--------------------------------------------------------------------*/
/// Get opposite velocity direction 
/** \param[in]  a_ei   Sequential velocity direction
 *  \return            Index of opposite vel direction
 *//*-----------------------------------------------------------------*/

inline int oppositeVelDir(const int a_ei)
{
  CH_assert(a_ei >= 0 && a_ei < g_numVelDir);
  static constexpr int idxOpposite[g_numVelDir] =
    { 0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15 };
  return idxOpposite[a_ei];
}

/*--------------------------------------------------------------------*/
/// Get sequential velocity index from a vector
/** \param[in]  a_iv   Velocity vector
 *  \return            Sequential velocity direction
 *//*-----------------------------------------------------------------*/

inline int velIndex(const IntVect& a_iv)
{
  static constexpr int listVelIndex[27] =
    {
      -1,  7, -1, 11,  1, 13, -1,  9, -1,
      15,  3, 17,  5,  0,  6, 16,  4, 18,
      -1,  8, -1, 12,  2, 14, -1, 10, -1
    };
  return listVelIndex[9*a_iv[0] + 3*a_iv[1] + a_iv[2] + 13];
}

/*--------------------------------------------------------------------*/
/// Flags of directions
/** \param[in]  a_ei   Sequential velocity direction
 *  \return            Flag
 *//*-----------------------------------------------------------------*/

inline unsigned streamFillFlags(const int a_ei)
{
  static constexpr unsigned flags[19] =
    {
      0x00000,
      0x02a82,
      0x05504,
      0x28188,
      0x50610,
      0x19820,
      0x66040,
      0x00080,
      0x00100,
      0x00200,
      0x00400,
      0x00800,
      0x01000,
      0x02000,
      0x04000,
      0x08000,
      0x10000,
      0x20000,
      0x40000
    };
  return flags[a_ei];
}

/*--------------------------------------------------------------------*/
/// Get variable state names
/** \param[in]  a_iVar Variable index as in 'U'
 *  \return            Name
 *//*-----------------------------------------------------------------*/

inline const char *const * stateNames()
{
  static const char *const names[] =
  {
    "density",
    "x-velocity",
    "y-velocity",
    "z-velocity"
  };
  return names;
}

/*--------------------------------------------------------------------*/
/// Return the distribution function name for a velocity direction
/** \param[in]  a_iVel Velocity direction
 *  \param[out] a_str  Location to create name
 *  \param[int] a_len  Length of 'a_str'.  Must be >= 6.
 *  \return            a_str
 *//*-----------------------------------------------------------------*/

inline const char* distrFuncNameVel(const int   a_iVel,
                                    char *const a_str,
                                    const int   a_len)
{
  CH_assert(a_len >= 6);
  std::sprintf(a_str, "fi_%02d", a_iVel);
  return a_str;
}

/*--------------------------------------------------------------------*/
/// Comparison with limit tol^2 as x and y -> 0
/** Not really an LB function but often used for verification
 *  \return             T - not equal
 *//*-----------------------------------------------------------------*/

template <typename T>
inline bool compare(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
  return std::fabs(x - y) >
    (std::min(std::fabs(x), std::fabs(y)) + tol)*tol;
}

}  // namespace LBParameters

#endif
