
#ifndef _BOXITERATOR_H_
#define _BOXITERATOR_H_


/******************************************************************************/
/**
 * \file BoxIterator.H
 *
 * \brief Iterator over all IntVects in a Box
 *
 *//*+*************************************************************************/

#include <cstring>
#include <ostream>
#include <iterator>

#include "Parameters.H"
#include "Box.H"


/*******************************************************************************
 */
///  Iterate over IntVects in a Box
/**
 ******************************************************************************/

template <typename T, typename Ref, typename Ptr>
class BoxIteratorImpl
{


/*==============================================================================
 * Types
 *============================================================================*/

public:

  using value_type = T;
  using reference = Ref;
  using pointer = Ptr;

  using iterator_category = std::random_access_iterator_tag;
  using difference_type = T;

  using Self = BoxIteratorImpl<T, Ref, Ptr>;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  BoxIteratorImpl();

  /// Construct with a box
  BoxIteratorImpl(const Box& a_box);

  /// Construct with a box and set current IV
  BoxIteratorImpl(const Box& a_box, const T& a_iv);

  /// Copy constructor
  BoxIteratorImpl(const BoxIteratorImpl&) = default;

  /// Move constructor
  BoxIteratorImpl(BoxIteratorImpl&&) = default;

  /// Assignment constructor
  BoxIteratorImpl& operator=(const BoxIteratorImpl&) = default;

  /// Move assignment constructor
  BoxIteratorImpl& operator=(BoxIteratorImpl&&) = default;

  /// Destructor
  ~BoxIteratorImpl() = default;


/*==============================================================================
 * Members functions
 *============================================================================*/

public:

  /// Reference to the IntVect
  reference operator*() const;

  /// Pointer to the IntVect
  pointer operator->() const;

  /// Prefix increment
  Self& operator++();

  /// Postfix increment
  Self operator++(int);

  /// Assignment to type T
  void operator=(const T& a_iv);

  /// Increment by difference
  Self& operator+=(const difference_type& a_deltaiv);
  
  /// Equal
  bool operator==(const T& a_iv) const;

  /// Not Equal
  bool operator!=(const T& a_iv) const;

  /// Not Equal
  bool operator!=(const Self& a_iter) const;

  /// Still valid
  bool ok() const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  T m_current;                        ///< Current IntVect
  Box m_box;                          ///< Dimensions
};


/*******************************************************************************
 *
 * Class BoxIteratorImpl: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline
BoxIteratorImpl<T, Ref, Ptr>::BoxIteratorImpl()
  :
  m_current(IntVect::Zero)
{
}

/*--------------------------------------------------------------------*/
//  Construct with a box
/** Starts at lower corner of box
 *//*-----------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline
BoxIteratorImpl<T, Ref, Ptr>::BoxIteratorImpl(const Box& a_box)
  :
  m_current(a_box.loVect()),
  m_box(a_box)
{
}

/*--------------------------------------------------------------------*/
//  Construct with a box and set current IV
/** Starts at specified IntVect.  Mostly used for creating 'end'
 *  iterators
 *//*-----------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline
BoxIteratorImpl<T, Ref, Ptr>::BoxIteratorImpl(const Box& a_box,
                                              const T&   a_iv)
  :
  m_current(a_iv),
  m_box(a_box)
{
}

/*--------------------------------------------------------------------*/
//  Reference to the IntVect
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline typename BoxIteratorImpl<T, Ref, Ptr>::reference
BoxIteratorImpl<T, Ref, Ptr>::operator*() const
{
  return m_current;
}

/*--------------------------------------------------------------------*/
//  Pointer to the IntVect
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline typename BoxIteratorImpl<T, Ref, Ptr>::pointer
BoxIteratorImpl<T, Ref, Ptr>::operator->() const
{
  return &m_current;
}

/*--------------------------------------------------------------------*/
//  Prefix increment
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline typename BoxIteratorImpl<T, Ref, Ptr>::Self&
BoxIteratorImpl<T, Ref, Ptr>::operator++()
{

 //  if(m_current[0] < m_box.hiVect(0))
 // {
 //  m_current[0]++;
 //  return *this;
 // }
 // // if the lowest dim is maxed and we can increment the second dim, and reset the lowest dim.  
 // else if(m_current[0] == m_box.hiVect(0) && m_current[1] < m_box.hiVect(1))
 // {
 //  m_current[0] = m_box.loVect(0);
 //  m_current[1]++;
 //  return *this;
 // }
 // // Checking to see if we have another dimension of freedom.
 // else if (g_SpaceDim == 3)
 // {
 //  //if the first two dimensions are maxed and we have a 3rd, increment the 3rd and reset the first two.
 //  if(m_current[0] == m_box.hiVect(0) && m_current[1] == m_box.hiVect(1))
 //  {
 //   m_current[0] = m_box.loVect(0);
 //   m_current[1] = m_box.loVect(1);
 //   m_current[2]++;
 //   return *this;
 //  }
 //  else // this is the end case where we should return the next value it would be. 
 //  {
 //   m_current[0] = m_box.loVect(0);
 //   m_current[1] = m_box.loVect(1);
 //   m_current[2] = m_box.hiVect(2) + 1; 
 //   return *this;
 //  }
 // }
 // else // this is the case where we have SPACEDIM = 2 and need to return the end case.
 // {
 //  m_current[0] = m_box.loVect(0);
 //  m_current[1] = m_box.hiVect(1) + 1;
 //  return *this;
 // }

  switch (g_SpaceDim)
  {
    case 2:
      if(m_current[0] < m_box.hiVect(0))
      {
        ++m_current[0];
      }
      else
      {
        m_current[0] = m_box.loVect(0);
        ++m_current[1];
      }
      return *this;

    case 3:
      if(m_current[0] < m_box.hiVect(0))
      {
        ++m_current[0];
      }
      else if(m_current[1] < m_box.hiVect(1))
      {
        m_current[0] = m_box.loVect(0);
        ++m_current[1];
      }
      else
      {
        m_current[0] = m_box.loVect(0);
        m_current[1] = m_box.loVect(1);
        ++m_current[2];
      }
      return *this;

    default:
      return *this;
  }//end switch
}

/*--------------------------------------------------------------------*/
//  Postfix increment
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline typename BoxIteratorImpl<T, Ref, Ptr>::Self
BoxIteratorImpl<T, Ref, Ptr>::operator++(int)
{
  Self tmp = *this;
  this->operator++();
  return tmp;
}

/*--------------------------------------------------------------------*/
//  Assign to an IntVect
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline void
BoxIteratorImpl<T, Ref, Ptr>::operator=(const T& a_iv)
{
  m_current = a_iv;
}

/*--------------------------------------------------------------------*/
//  Increment by difference
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline typename BoxIteratorImpl<T, Ref, Ptr>::Self&
BoxIteratorImpl<T, Ref, Ptr>::operator+=(const T& a_deltaiv)
{
  m_current += a_deltaiv;
  return *this;
}

/*--------------------------------------------------------------------*/
//  Equal
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline bool
BoxIteratorImpl<T, Ref, Ptr>::operator==(const T& a_iv) const
{
  return (m_current == a_iv);
}

/*--------------------------------------------------------------------*/
//  Not equal
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline bool
BoxIteratorImpl<T, Ref, Ptr>::operator!=(const T& a_iv) const
{
  return (m_current != a_iv);
}

/*--------------------------------------------------------------------*/
//  Not equal
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline bool
BoxIteratorImpl<T, Ref, Ptr>::operator!=(const Self& a_iter) const
{
  // It doesn't make sense to compare iterators built with different boxes
  CH_assert(m_box == a_iter.m_box);
  return (m_current != a_iter.m_current);
}

/*--------------------------------------------------------------------*/
//  Still valid
/*--------------------------------------------------------------------*/

template <typename T, typename Ref, typename Ptr>
inline bool
BoxIteratorImpl<T, Ref, Ptr>::ok() const
{
  if(D_TERM( (m_current[0] >= m_box.loVect(0) && m_current[0]<=m_box.hiVect(0)),
          && (m_current[0] >= m_box.loVect(0) && m_current[1]<=m_box.hiVect(1)),
          && (m_current[0] >= m_box.loVect(0) && m_current[2]<=m_box.hiVect(2))))
  {
    return true;
  }
  else{return false;}
  
}

#endif
