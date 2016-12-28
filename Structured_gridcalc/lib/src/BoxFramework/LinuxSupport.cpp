
/******************************************************************************/
/**
 * \file LinuxSupport.cpp
 *
 * \brief System routines for Linux OS
 *
 *//*+*************************************************************************/

// Feature Test Macro for posix_memalign
#ifdef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN
#define _XOPEN_SOURCE 600
#endif
// Feature Test Macro for nanosleep
#define _POSIX_C_SOURCE 199309L
#include <unistd.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

#include "LinuxSupport.H"


/*============================================================================*/
//  Get the path and name of the currently running executable
/**
 *  \param[out] procname
 *                      Path and name of the executable
 *  \return             > 0     - Successful -- length of string a_procPath
 *                      =-1     - Error
 *                      =-a_len - Likely ran out of space in a_procPath
 *//*=========================================================================*/

int System::getProcessPath(char *const a_procPath, int a_len)
{
  int len = readlink("/proc/self/exe", a_procPath, a_len-1);
  if (len > 0)
    {
      a_procPath[len]='\0';
      if (len == a_len-1)
        {
          len = -a_len;
        }
    }
  return len;
}


/*============================================================================*/
//  Allocate aligned memory
/**
 *  \param[out] a_memptr
 *                      Pointer to allocated memory
 *  \param[in]  a_alignment
 *                      Alignment in bytes.  Must be a multiple of sizeof(void*)
 *                      and a power of 2.
 *  \param[in]  a_size  Number of bytes to allocate
 *  \return             0       - Success
 *                      <posix_memalign>
 *                      EINVAL  - Invalid alignment parameter
 *                      ENOMEM  - Out of memory
 *                      <malloc>
 *                      1       - Out of memory
 *  \note
 *  <ul>
 *    <li> This function returns raw memory.  Use placement new for construction
 *         of objects if required.
 *    <li> Memory allocated with memalign should be deallocated with free()
 *  </ul>
 *//*=========================================================================*/

int System::memalign(void **a_memptr, size_t a_alignment, size_t a_size)
{
  return posix_memalign(a_memptr, a_alignment, a_size);
}


/*============================================================================*/
//  Checks if a file exists
/**
 *  \param[in]  a_fileName
 *                      Name of the file
 *  \return             1 - File exists
 *                      0 - File does not exist
 *//*=========================================================================*/

int System::fileExists(const char *const a_filename)
{
  struct stat buf;
  if (stat(a_filename, &buf) == 0)
    {
      return 1;
    }
  return 0;
}


/*============================================================================*/
//  Sleep for a while
/**
 *  \param[in]  s       Time to sleep in seconds
 *  \return              0 - successful sleep
 *                      -1 - sleep interrupted
 *//*=========================================================================*/

int System::sleep(const double s)
{
  const unsigned sec = static_cast<unsigned int>(std::fabs(s));
  timespec req, rem;
  req.tv_sec = sec;
  req.tv_nsec = static_cast<long>((std::fabs(s) - sec)*1.E9);
  return nanosleep(&req, &rem);
}
