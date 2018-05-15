#ifdef NDEBUG
#define ASSERT(cond)
#define PRECONDITION(cond)
#define POSTCONDITION(cond)
#define INVARIANT(cond)
#else
#define BASIC_ASSERT_MESSAGE(type, cond, file, line) file ":" #line " " type " " #cond " failed"
#ifdef __cplusplus
#include <iostream>
#define BASIC_ASSERT_FAILED(type, cond, file, line) \
    do { \
        std::cerr << (BASIC_ASSERT_MESSAGE(type, cond, file, line)) << "\n"; \
        std::cerr.flush(); \
        std::abort(); \
    } while (0)
#else
#include <stdio.h>
#include <stdlib.h>
#define BASIC_ASSERT_FAILED(type, cond, file, line) \
do { \
    fputs(BASIC_ASSERT_MESSAGE(type, cond, file, line), stderr); \
    fflush(stderr); \
    abort(); \
} while (0)
#endif

#define BASIC_ASSERT(type, cond) \
    do { \
        if (!(cond)) { \
            BASIC_ASSERT_FAILED(type, cond, __FILE__, __LINE__); \
        } \
    } while (0)

#define ASSERT(cond) BASIC_ASSERT("assert", cond)
#define PRECONDITION(cond) BASIC_ASSERT("precondition", cond)
#define POSTCONDITION(cond) BASIC_ASSERT("postcondition", cond)
#define INVARIANT(cond) BASIC_ASSERT("invariant", cond)

#endif
