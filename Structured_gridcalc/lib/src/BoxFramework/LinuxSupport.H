
#ifndef _LINUXSUPPORT_H_
#define _LINUXSUPPORT_H_

/******************************************************************************/
/**
 * \file LinuxSupport.H
 *
 * \brief System routines for Linux OS
 *
 *//*+*************************************************************************/

#include <cstring>

namespace System
{

/// Get the path and name of the currently running executable
int getProcessPath(char *const a_procPath, int a_len);

/// Allocate aligned memory
int memalign(void **a_memptr, size_t a_alignment, size_t a_size);

/// Checks if a file exists
int fileExists(const char *const a_filename);

/// Sleep for a while
int sleep(const double s);

}  // Namespace System

#endif
