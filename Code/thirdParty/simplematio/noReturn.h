#pragma once
#if defined(__cplusplus) && __cplusplus >= 201103L
    #define NO_RETURN [[noreturn]]
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
    #define NO_RETURN _Noreturn
#elif defined(__GNUC__)
    #define NO_RETURN __attribute__((noreturn))
#elif defined(_MSC_VER)
    #define NO_RETURN __declspec(noreturn)
#else
    #define NO_RETURN
    #warning "Didn't find suitable implementation of macro NO_RETURN; defined it to be empty"
#endif
