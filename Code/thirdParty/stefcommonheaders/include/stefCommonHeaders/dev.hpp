#pragma once
#ifndef DEV_HPP
#define DEV_HPP

#include <iostream>
#include <cstdlib>

#if __cplusplus < 201103L
#error "This header needs at least C++ 11!"
#endif
namespace dev {

    [[noreturn]] static inline void stub(const char* function, const char* fileName, const int line) {
        std::cerr << fileName << ":" << line << " " << function << " is just a stub!\n";
        std::abort();
    }

    #define DEV_STUB(ignoredArgument) ::dev::stub(__func__, __FILE__, __LINE__)
}

#endif
