/**************************************************************************************************
 *                                     DISCLAIMER AND LICENSE                                     *
 **************************************************************************************************/

/* MacOS implementations originally taken from
 * http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c (public domain)
 * by David N. Williams
 */

/*
 * Copyright 2017 Stefano Weidmann
 * Author: Weidmann, Stefano
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/**************************************************************************************************
 *                                              TLDR                                              *
 **************************************************************************************************/

 /*
 #include "stefFenv.h"
 // Link with -lm if using in plain C!
 int main() {
    StefFenv_CrashOnFPEs(FE_ALL_EXCEPT & ~FE_INEXACT);
 }
 */


/**************************************************************************************************
 *                                          DESCRIPTION                                           *
 **************************************************************************************************/

/* Ports the GNU floating point exception API to MacOS and Windows:
 * int feenableexcept(int excepts);
 * int fedisableexcept(int excepts);
 * int fegetexcept(void);
 *
 * And add some other niceties
 */

/**************************************************************************************************
 *                                             USAGE                                              *
 **************************************************************************************************/


/* int excepts is a ORed composition of:
 * FE_DIVBYZERO (Pole error: division by zero, or some other asymptotically infinite result (from finite arguments)),
 * FE_INEXACT (Inexact: the result is not exact),
 * FE_INVALID (Domain error: At least one of the arguments is a value for which the function is not defined)
 * FE_OVERFLOW (Overflow range error: The result is too large in magnitude to be represented as a value of the return type)
 * FE_UNDERFLOW (Underflow range error: The result is too small in magnitude to be represented as a value of the return type)
 *
 * FE_ALL_EXCEPT is a combination of all flags above
 *
 * int feenableexcept(int excepts) makes a unixoid OS send the signal SIGFPE to the process if any of the exceptions specified happen.
 * The default action of a process is to abort, at least on MacOS and Linux.
 * On Windows it throws an exception. You must enable SEH (structured exception handling) for it to work.
 *
 * int fedisableexcept(int excepts) disables the specified exceptions (counterpart to feenableexcept)
 *
 * int fegetexcept(void) returns the currently enables exceptions (format like int excepts)
 *
 *
 * Example
 * -----------
 * feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT) // crash on everything except inexact math

 * If you want prettier crashing, do
 */

/**************************************************************************************************
 *                                      INCLUDES AND MACROS                                       *
 **************************************************************************************************/


#pragma once
#ifdef __cplusplus
extern "C" {
#endif

#ifdef _GNU_SOURCE
    #define STEF_FENV_WAS_GNU_SOURCE 1
#else
    #define _GNU_SOURCE
    #define STEF_FENV_WAS_GNU_SOURCE 0
#endif


#if defined(__unix__) || defined(__unix) || \
        (defined(__APPLE__) && defined(__MACH__))
    #if defined(_POSIX_C_SOURCE)
        #if _POSIX_C_SOURCE < 200112L
            #error "This needs _POSIX_C_SOURCE >= 200112L for siginfo_t!"
        #endif
    #else
        #define _POSIX_C_SOURCE 200112L
    #endif
#endif

#include <stdlib.h>
#include <cfenv>
#include <fenv.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef BSD
#pragma STDC FENV_ACCESS ON
#endif


/**************************************************************************************************
 *                                           MACOS PORT                                           *
 **************************************************************************************************/

#if defined(__APPLE__) && defined(__MACH__)

#if !defined(__i386__) && !defined(__x86_64__)
#error "Doesn't work on non-intel macs!"
#endif


static inline
int
fegetexcept(void) {
    static fenv_t fenv;
    return fegetenv(&fenv) ? -1 : (fenv.__control & FE_ALL_EXCEPT);
}

static inline
int
feenableexcept(int excepts) {
    static fenv_t fenv;
    int new_excepts = excepts & FE_ALL_EXCEPT, old_excepts;  // previous masks

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);

    return (fesetenv(&fenv) ? -1 : old_excepts);
}

static inline
int
fedisableexcept(int excepts) {
    static fenv_t fenv;
    int new_excepts = excepts & FE_ALL_EXCEPT,
                      old_excepts;  // all previous masks

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // mask
    fenv.__control |= new_excepts;
    fenv.__mxcsr   |= new_excepts << 7;

    return (fesetenv(&fenv) ? -1 : old_excepts);
}
#endif


/**************************************************************************************************
 *                                          WINDOWS PORT                                          *
 **************************************************************************************************/
// TODO Test it

#if defined(_WIN32) || defined(_WIN64)
#include <float.h>
#include <stdint.h>

// _EM_DENORMAL is left out because it doesn't have a fenv counterpart
const static int StefFenv_WindowsExceptionConstants_C[] = {_EM_INVALID, _EM_ZERODIVIDE, _EM_OVERFLOW, _EM_UNDERFLOW, _EM_INEXACT};
const static int StefFenv_FenvExceptionConstants_C[] =    {FE_INVALID,  FE_DIVBYZERO,   FE_OVERFLOW,  FE_UNDERFLOW,  FE_INEXACT};
const static int StefFenv_NumberOfExceptionConstants_C = sizeof(StefFenv_WindowsExceptionConstants_C) / sizeof(int);

/* windows exception masks work the other way around than fenv masks
 * the exception is enabled if the constant is NOT set
 */

static inline
int
StefFenv_TranslateWindowsExceptionMaskToFenv(const int windowsExceptionMask) {
    int fenvExceptionMask = 0;

    for (int i = 0; i < StefFenv_NumberOfExceptionConstants_C; ++i) {
        if (windowsExceptionMask & ~StefFenv_WindowsExceptionConstants_C[i]) {
            fenvExceptionMask |= StefFenv_FenvExceptionConstants_C[i];
        }
    }

    return fenvExceptionMask;
}


static inline
int
StefFenv_TranslateFenvExceptionMaskToWindows(const int fenvExceptionMask) {
    int windowsExceptionMask = -1;

    for (int i = 0; i < StefFenv_NumberOfExceptionConstants_C; ++i) {
        if (fenvExceptionMask & StefFenv_FenvExceptionConstants_C[i]) {
            windowsExceptionMask &= ~StefFenv_WindowsExceptionConstants_C[i];
        }
    }

    return windowsExceptionMask;
}



static inline
int
fegetexcept(void) {
    _clearfp(); // Clearing enables exceptions afterwards
    const int currentState = _controlfp(0, 0);
    return StefFenv_TranslateWindowsExceptionMaskToFenv(currentState);
}

static inline
int
feenableexcept(int excepts) {
    _clearfp(); // Clearing enables exceptions afterwards
    const int stateBefore = fegetexcept();
    const int windowsExceptionMask = StefFenv_TranslateFenvExceptionMaskToWindows(excepts);
    _controlfp(windowsExceptionMask, _MCW_EM);
    return stateBefore;
}

static inline
int
fedisableexcept(int excepts) {
    _clearfp(); // Clearing enables exceptions afterwards
    const int stateBefore = fegetexcept();
    const int negatedWindowsExceptionMask = ~StefFenv_TranslateFenvExceptionMaskToWindows(excepts);

    _controlfp(negatedWindowsExceptionMask, _MCW_EM);
    return stateBefore;
}

#endif


/**************************************************************************************************
 *                                  PRETTIER CRASHING WITH POSIX                                  *
 **************************************************************************************************/

// Signal handler adapted from http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c

#if defined(__unix__) || defined(__unix) || \
        (defined(__APPLE__) && defined(__MACH__))

#include <signal.h>
#include <unistd.h>

typedef void (*StefFenv_FloatingPointExceptionHandler_t)(
    int signal,
    siginfo_t* signalInfos,
    void* somethingIdontCareAbout
);

static inline
const char*
StefFenv_GetFloatingPointExceptionDescription(
    const int signalExceptionCode
) {
    const static int signalExceptionCodes[] = {
	FPE_INTDIV,
	FPE_INTOVF,
	FPE_FLTDIV,
	FPE_FLTOVF,
	FPE_FLTUND,
	FPE_FLTRES,
	FPE_FLTINV,
	FPE_FLTSUB
    };

    const static char* exceptionDescriptions[] = {
	"integer division by zero",
	"integer overflow",
	"floating point divide by zero",
	"floating point overflow",
	"floating point underflow",
	"floating point inexact result",
	"floating point invalid operation",
	"subscript out of range"
    };
    const int numberOfExceptions = sizeof(signalExceptionCodes) / sizeof(signalExceptionCodes[0]);

    for (int i = 0; i < numberOfExceptions; ++i) {
        if (signalExceptionCode == signalExceptionCodes[i]) {
            return exceptionDescriptions[i];
        }
    }

    return "no floating point exception";
}



static inline
void
StefFenv_DefaultFloatingPointExceptionHandler(
    int signal,
    siginfo_t* signalInfos,
    void* somethingIdontCareAbout
) {
    if (signal == SIGFPE) {
        const char* exceptionDescription = StefFenv_GetFloatingPointExceptionDescription(signalInfos->si_code);
        fputs("**************************************************************************************************\n", stderr);
    	fprintf(stderr, "SIGNAL SIGFPE CAUGHT\n%s\n", exceptionDescription);
        fputs("**************************************************************************************************\n", stderr);
    } else {
        fputs("Should handle a floating point exception signal, but the signal sent wasn't SIGFPE!", stderr);
    }
    fflush(stderr);
    abort();
}


static inline
void
StefFenv_RegisterFloatingPointExceptionHandler(
    const StefFenv_FloatingPointExceptionHandler_t handler
) {

    struct sigaction signalAction;
    sigemptyset(&signalAction.sa_mask);
	/* The SA_SIGINFO flag tells sigaction() to use the sa_sigaction field, not sa_handler. */
    signalAction.sa_flags = SA_SIGINFO;
    signalAction.sa_sigaction = handler;

    if (sigaction(SIGFPE, &signalAction, NULL) < 0) {
		perror("error registering with sigaction");
        abort();
    }
}

static inline
void
StefFenv_ClearFloatingPointExceptionHandler(
) {
    struct sigaction signalAction;
    sigemptyset(&signalAction.sa_mask);
    signalAction.sa_flags = 0;
    signalAction.sa_handler = SIG_DFL;

    if (sigaction(SIGFPE, &signalAction, NULL) < 0) {
		perror("error registering with sigaction");
        abort();
    }
}

static inline
void
StefFenv_CrashOnFPEs(
    const int fpExceptions
) {
    fedisableexcept(FE_ALL_EXCEPT);
    feenableexcept(fpExceptions);
    StefFenv_RegisterFloatingPointExceptionHandler(StefFenv_DefaultFloatingPointExceptionHandler);
}

#else
// not unix
static inline
void
StefFenv_RegisterFloatingPointExceptionHandler(
) {
    fputs("Signal handling not supported on this platform!\n", stderr);
    abort();
}

static inline
void
StefFenv_ClearFloatingPointExceptionHandler(
) {
    fputs("Signal handling not supported on this platform!\n", stderr);
    abort();
}

static inline
void
StefFenv_CrashOnFPEs(
    const int fpExceptions
) {
    fputs("Signal handling not supported on this platform!\n", stderr);
    abort();
}
#endif // not windows


/**************************************************************************************************
 *                                         MACRO CLEANUP                                          *
 **************************************************************************************************/
#ifdef __cplusplus
}
#endif

#if !STEF_FENV_WAS_GNU_SOURCE
#undef _GNU_SOURCE
#endif
#undef STEF_FENV_WAS_GNU_SOURCE
