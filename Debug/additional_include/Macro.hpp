#ifndef __GEDIM_MACRO_H
#define __GEDIM_MACRO_H

/// Activate deprecated class use
/// - 0 false
/// - 1 true
#define USE_DEPRECATED 0

/// Use MPI
/// - 0 false
/// - 1 true
#define USE_MPI 0

/// Verbose Levels
/// - 0 None
/// - 1 Error
/// - 2 Warning
/// - 3 Info
/// - 4 Debug
#define VERBOSE 3

/// Logging Levels
/// - 0 None
/// - 1 Only Console
/// - 2 Only Files
/// - 3 Console and Files
#define LOGGING 3

// the configured options and settings for Tutorial
#define gedim_VERSION_MAJOR 
#define gedim_VERSION_MINOR 

/// @name Code Simplifications
///@{
#ifndef MIN
#define MIN(a,b) (a < b) ? a : b
#endif

#ifndef MAX
#define MAX(a,b) (a > b) ? a : b
#endif

///@}

#endif
