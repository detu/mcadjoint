//
// Created by stefano on 22.12.17.
//

#ifndef STEFCOMMONHEADERS_EIGENTOCSV_HPP
#define STEFCOMMONHEADERS_EIGENTOCSV_HPP
#pragma once
#include <Eigen/Core>
#include <string>

namespace stefCommonHeaders {

    static inline void writeToCsv(const std::string& fileName, const Eigen::Ref<const Eigen::MatrixXd>& matrixToWrite) {
        const static Eigen::IOFormat csvFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream outFile(fileName.c_str());
        outFile << matrixToWrite.format(csvFormat);

    }
}
#endif //STEFCOMMONHEADERS_EIGENTOCSV_HPP
