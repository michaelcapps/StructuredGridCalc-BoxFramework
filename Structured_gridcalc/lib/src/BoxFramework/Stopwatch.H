
#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

/******************************************************************************/
/**
 * \file Stopwatch.H
 *
 * \brief Implements a stopwatch using C++ libraries
 *
 *//*+*************************************************************************/

#include <chrono>


/*******************************************************************************
 */
///  Stopwatch
/*
 ******************************************************************************/

template <typename Clock = std::chrono::steady_clock>
class Stopwatch
{
public:

  // Types
  using duration = typename Clock::duration;
  using time_point = typename Clock::time_point;
  using ticks = typename duration::rep;

  /// Default constructor
  Stopwatch();

  /// Start
  void start();

  /// Stop
  void stop();

  /// Reset time (will continue running if already doing so)
  void reset();

  /// Return duration of latest session or between start and current time
  duration interval();

  /// Total time for all start/stop sessions (not including current)
  template <typename Units = std::milli>
  double time() const;

  /// Average time for a start/stop session (not including current)
  template <typename Units = std::milli>
  double averageTime() const;

  /// Total ticks for all start/stop sessions (not including current)
  ticks tickCount() const;

  /// Average ticks for a start/stop session (not including current)
  ticks averageTickCount() const;

  /// Get the resolution of the clock
  template <typename Units = std::nano>
  double static resolution();

//--Data

  time_point m_startTime;             ///< Time when 'start' called
  duration m_runTime;                 ///< Run time (between start and stop)
  duration m_totalTime;               ///< Total time (all sessions)
  long m_sessions;                    ///< Number of start/stop sessions
  bool m_running;                     ///< T - Timer running
};


/*******************************************************************************
 *
 * Class Stopwatch: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/

template <typename Clock>
inline
Stopwatch<Clock>::Stopwatch()
  : m_running(false)
{
  static_assert(Clock::is_steady, "Clock chosen for stopwatch is not steady");
  reset();
}

/*--------------------------------------------------------------------*/
//  Start
/*--------------------------------------------------------------------*/

template <typename Clock>
inline void
Stopwatch<Clock>::start()
{
  m_startTime = Clock::now();
  m_running = true;
}

/*--------------------------------------------------------------------*/
//  Stop
/*--------------------------------------------------------------------*/

template <typename Clock>
inline void
Stopwatch<Clock>::stop()
{
  m_runTime = Clock::now() - m_startTime;
  if (!m_running) return;
  m_totalTime += m_runTime;
  ++m_sessions;
  m_running = false;
}

/*--------------------------------------------------------------------*/
//  Reset time (will continue running if already doing so)
/*--------------------------------------------------------------------*/

template <typename Clock>
inline void
Stopwatch<Clock>::reset()
{
  m_totalTime = duration::zero();
  m_sessions = 0;
  if (m_running)
    {
      m_startTime = Clock::now();  // Reset start time
    }
}

/*--------------------------------------------------------------------*/
//  Return time (ms) of latest session or between start and current
//  time
/** \return             Last session  in milliseconds by default but
 *                      units can be changed)
 *//*-----------------------------------------------------------------*/

template <typename Clock>
inline typename Stopwatch<Clock>::duration
Stopwatch<Clock>::interval()
{
  duration tmp = Clock::now() - m_startTime;
  if (m_running) m_runTime = tmp;
  return m_runTime;
}

/*--------------------------------------------------------------------*/
//  Total time for all start/stop sessions (not including current)
/** The clock must be stopped to include the latest interval
 *  \return             Total time from all sessions (in milliseconds
 *                      by default but units can be changed)
 *//*-----------------------------------------------------------------*/

template <typename Clock>
template <typename Units>
inline double
Stopwatch<Clock>::time() const
{
  return std::chrono::duration_cast<std::chrono::duration<double, Units> >(
    m_totalTime).count();
}

/*--------------------------------------------------------------------*/
//  Average time for a start/stop session (not including current)
/** The clock must be stopped to include the latest interval
 *  \return             Average time of a session (in milliseconds
 *                      by default but units can be changed)
 *//*-----------------------------------------------------------------*/

template <typename Clock>
template <typename Units>
inline double
Stopwatch<Clock>::averageTime() const
{
  return time<Units>()/m_sessions;
}

/*--------------------------------------------------------------------*/
//  Total ticks for all start/stop sessions (not including current)
/** The clock must be stopped to include the latest interval
 *  \return             Total ticks in units of clock
 *//*-----------------------------------------------------------------*/

template <typename Clock>
inline typename Stopwatch<Clock>::ticks
Stopwatch<Clock>::tickCount() const
{
  return m_totalTime.count();
}

/*--------------------------------------------------------------------*/
//  Average time for a start/stop session (not including current)
/** The clock must be stopped to include the latest interval
 *  \return             Average ticks in units of clock
 *//*-----------------------------------------------------------------*/

template <typename Clock>
inline typename Stopwatch<Clock>::ticks
Stopwatch<Clock>::averageTickCount() const
{
  return tickCount()/m_sessions;
}

/*--------------------------------------------------------------------*/
//  Get the resolution of the clock
/** \return             Resolution of the clock in chosen units.
 *//*-----------------------------------------------------------------*/

template <typename Clock>
template <typename Units>
inline double
Stopwatch<Clock>::resolution()
{
  return std::chrono::duration_cast<std::chrono::duration<double, Units> >(
    duration(1)).count();
}

#endif
