//
// Created by Stefano Weidmann on 30.04.18.
//

#pragma once
#include <EigenSimplematio.hpp>
#include "typedefs.hpp"
#include <set>
#include <string>


void dumpInThisMatFile(const std::string& matFileName);
void dumpThis(const char* varName, ConstMatrixRef matrix);
void dumpThis(const char* varName, const Real scalar);
void dumpThis(const char* varName, const SparseMatrix& matrix);
void writeToMatFile();


