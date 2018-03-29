//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_CELLINDEX_HPP
#define STEFCOMMONHEADERS_CELLINDEX_HPP

#include "typedefs.hpp"
#include <array>

#if !defined(__GNUC__) && !defined(__attribute__)
    #define __attribute__(ignored)
#endif

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


    inline bool operator==(const CellIndex& other) const {
        return i == other.i && j == other.j;
    }

    inline bool operator!=(const CellIndex& other) const {
        return !(*this == other);
    }


    void shiftPoint(Point& point, const int numberOfRows, const int numberOfCols, const Real meshWidth) const {
        point.x += Real(j) / Real(numberOfCols) * meshWidth;
        point.y -= Real(i) / Real(numberOfRows) * meshWidth;
    }

    inline Point toCenterPoint(const int numberOfRows, const int numberOfCols, const Real meshWidth) const {
        Point point = {meshWidth/2.0, 1.0 - meshWidth/2.0};
        shiftPoint(point, numberOfRows, numberOfCols, meshWidth);
        return point;
    }

    inline Point toEastWestBorderPoint(const int numberOfRows, const int numberOfCols, const Real meshWidth) const {
        Point point = {0.0, 1.0 - meshWidth/2.0};

        shiftPoint(point, numberOfRows, numberOfCols, meshWidth);
        return point;
    }

    inline Point toNorthSouthBorderPoint(const int numberOfRows, const int numberOfCols, const Real meshWidth) const {
        Point point = {meshWidth/2.0, 1.0};

        shiftPoint(point, numberOfRows, numberOfCols, meshWidth);
        return point;
    }

    inline Point toBorderPoint(const int numberOfRows, const int numberOfCols, const Real meshWidth, const CellIndex::Direction whichBorder) const {
        switch (whichBorder) {
            case CellIndex::Direction::NORTH:
            case CellIndex::Direction::SOUTH: {
                return toNorthSouthBorderPoint(numberOfRows, numberOfCols, meshWidth);
            }

            case CellIndex::Direction::EAST:
            case CellIndex::Direction::WEST: {
                return toEastWestBorderPoint(numberOfRows, numberOfCols, meshWidth);
            }

        }
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
