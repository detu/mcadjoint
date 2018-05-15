#pragma once
#if !defined(THREAD_LOCAL)
    #if __STDC_VERSION__ < 201112L
        #if defined(__GNUC__)
            #define THREAD_LOCAL __thread
        #elif defined(_MSC_VER)
            #define THREAD_LOCAL __declspec(thread)
        #else
            #error "Don't know how to define a thread local variable!"
        #endif
    #else
        #define THREAD_LOCAL _Thread_local
    #endif
#endif

