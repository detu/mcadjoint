#pragma once
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

static int MpiDebug_AttachedSuccessfully_G = 0;

static inline
void
MpiDebug_WaitUntilAttachedImplementation(const char* fileName, const int lineNumber) {
    if (getenv("MPIDEBUG") == NULL) {
        return;
    }
    printf("\nmpiDebugPid=%d mpiDebugLine=%d mpiDebugFile=\"%s\"\n", getpid(), lineNumber, fileName);

    fflush(stdout);

    // the debugger must set MpiDebug_AttachedSuccessfully_G to nonzero after attaching
    while (!MpiDebug_AttachedSuccessfully_G) {
        sleep(5);
    }
}



#define MpiDebug_WaitUntilAttached(ignoredArguments) MpiDebug_WaitUntilAttachedImplementation(__FILE__, __LINE__)

#ifdef __cplusplus
}
#endif
