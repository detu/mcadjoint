//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_CELLINDEX_HPP
#define STEFCOMMONHEADERS_CELLINDEX_HPP

#include "typedefs.hpp"
#include <array>

struct CellIndex {
    int i;
    int j;


    enum class Direction {
        NORTH = 0, EAST = 1, WEST = 2, SOUTH = 3
    };

    template <typename Derived>
    Real& operator()(Eigen::MatrixBase<Derived>& expression) const {
        return expression(i, j);
    }

    template <typename Derived>
    Real operator()(const Eigen::MatrixBase<Derived>& expression) const {
        return expression(i, j);
    }

    inline Real& operator()(SparseMatrix& sparseMatrix) const {
        return sparseMatrix.coeffRef(i, j);
    }

    inline Real operator()(const SparseMatrix& sparseMatrix) const {
        return sparseMatrix.coeff(i, j);
    }

    inline int linearIndex(const int numberOfRows) const {
        return i + numberOfRows * j;
    }


    inline bool operator==(const CellIndex& other) {
        return i == other.i && j == other.j;
    }

    inline bool operator!=(const CellIndex& other) {
        return !(*this == other);
    }



    inline CellIndex neighbor(const Direction direction) const {
        switch (direction) {
            case (Direction::NORTH): {
                return {i-1, j};
            }

            case (Direction::SOUTH): {
                return {i+1, j};
            }

            case (Direction::EAST): {
                return {i, j+1};
            }

            case (Direction::WEST): {
                return {i, j-1};
            }
        }
    }

    inline bool hasNeighbor(const Direction direction, const long numberOfRows, const long numberOfCols) const {
        switch (direction) {
            case (Direction::NORTH): {
                return i > 0;
            }

            case (Direction::SOUTH): {
                return i < numberOfRows - 1;
            }

            case (Direction::EAST): {
                return j < numberOfCols - 1;
            }

            case (Direction::WEST): {
                return j > 0;
            }
        }
    }

    inline CellIndex transpose() const {
        return {j, i};
    }



};

#endif //STEFCOMMONHEADERS_CELLINDEX_HPP
