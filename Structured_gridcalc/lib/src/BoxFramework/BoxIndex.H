
#ifndef _BOXINDEX_H_
#define _BOXINDEX_H_


/******************************************************************************/
/**
 * \file BoxIndex.H
 *
 * \brief Class providing data for indexing containers based on layouts of boxes
 *
 *//*+*************************************************************************/


/*******************************************************************************
 */
///  Information about a box in a layout
/**
 *//*+*************************************************************************/

class BoxIndex
{
public:

  /// Default constructor
  BoxIndex()
    :
    m_globalIndex(-1),
    m_localIndex(-1)
    {
    }

  /// Constructor
  BoxIndex(const int a_globalIndex, const int a_localIndex)
    :
    m_globalIndex(a_globalIndex),
    m_localIndex(a_localIndex)
    {
    }

//--Use synthesized copy, assignment, destructor

  /// Access global index
  int globalIndex() const
    {
      return m_globalIndex;
    }

  /// Access local index
  int localIndex() const
    {
      return m_localIndex;
    }
  
protected:
  int m_globalIndex;                  ///< Global linear index of the box used
                                      ///< by DisjointBoxLayout
  int m_localIndex;                   ///< Local linear index to this proc used
                                      ///< by LevelData
};

#endif  /* ! defined _BOXINDEX_H_ */
