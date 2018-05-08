//
// Created by Stefano Weidmann on 08.05.18.
//

#include "logging.hpp"

#ifdef MULTITHREADING
omp_mutex LOG_MUTEX;
#endif