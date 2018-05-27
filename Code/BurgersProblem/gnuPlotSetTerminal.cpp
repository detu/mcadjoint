//
// Created by Stefano Weidmann on 04.03.18.
//

#define _POSIX_C_SOURCE 200809L
#include <unistd.h>
#include "gnuPlotSetTerminal.hpp"


void setTerminal(Gnuplot& g) {
    g.cmd("set term png");
    return;
    static char hostname[66];
    gethostname(hostname, sizeof(hostname) / sizeof(hostname[0]));
    if (std::strcmp(hostname, "swmb") == 0) {
        g.cmd("set term aqua"); // use aqua on my laptop
    } else {
        g.cmd("set term x11");                                    // setting output terminal to "x11"
    }
}
