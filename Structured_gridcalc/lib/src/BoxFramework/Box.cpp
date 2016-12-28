
/******************************************************************************/
/**
 * \file Box.cpp
 *
 * \brief Non-inline definitions for classes in Box.H
 *
 *//*+*************************************************************************/

#include <ostream>

#include "Box.H"
#include "BoxIterator.H"


/*******************************************************************************
 *
 * Class Box: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Starting iterator
/*--------------------------------------------------------------------*/

Box::const_iterator
Box::begin() const
{
	return const_iterator(*this,m_lo);
}

/*--------------------------------------------------------------------*/
//  Ending iterator
/*--------------------------------------------------------------------*/

Box::const_iterator
Box::end() const
{
  return const_iterator(
    *this,
    m_lo + IntVect(D_DECL((g_SpaceDim == 1)*(m_hi[0] - m_lo[0] + 1),
                          (g_SpaceDim == 2)*(m_hi[1] - m_lo[1] + 1),
                          (g_SpaceDim == 3)*(m_hi[2] - m_lo[2] + 1))));
}


/*******************************************************************************
 *
 * Class Box: external related functions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Output the Box
/** \param[in] a_os     The output stream
 *  \param[in] a_box    Box to output
 *  \return             The output stream
 *//*-----------------------------------------------------------------*/

std::ostream &operator<<(std::ostream &a_os, const Box& a_box)
{
  a_os << '(' << a_box.loVect() << ", " << a_box.hiVect() << ')';
  return a_os;
}
