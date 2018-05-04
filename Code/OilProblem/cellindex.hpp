//
// Created by Stefano Weidmann on 24.03.18.
//

#pragma once

#include "typedefs.hpp"
#include <array>
#include <string>
#include <stefCommonHeaders/assert.h>
#include <cstdlib>
#include <stdlib.h>

#if !defined(__GNUC__) && !defined(__attribute__)
    #define __attribute__(ignored)
#endif

struct CellIndex {
    int i;
    int j;

    constexpr static inline CellIndex invalidCell() {
        return {-1, -1};
    }

    enum class Direction {
        NORTH = 0, EAST = 1, WEST = 2, SOUTH = 3
    };

    static inline const char* directionToString(const CellIndex::Direction direction) {
        constexpr static std::array<const char*, 4> directionNames = {
              "north", "east", "west", "south"
        };

        return directionNames[static_cast<int>(direction)];
    }


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

    
    inline int linearIndexWithSkip(const int numberOfRows, const int linearIndexToSkip) const {
        int candidate = linearIndex(numberOfRows);
        ASSERT(candidate != linearIndexToSkip);

        if (candidate > linearIndexToSkip) {
            --candidate;
        }

        return candidate;
    }

    
    inline bool operator==(const CellIndex& other) const {
        return i == other.i && j == other.j;
    }

    
    inline bool operator!=(const CellIndex& other) const {
        return !(*this == other);
    }


    inline void shiftPoint(Point& point, const Real meshWidth) const {
        point.x += Real(j) * meshWidth;
        point.y -= Real(i) * meshWidth;
    }

    
    inline Point toCenterPoint(const Real meshWidth) const {
        Point point = {meshWidth/2.0, 1.0 - meshWidth/2.0};
        shiftPoint(point, meshWidth);
        return point;
    }

    
    inline Point toBorderPoint(const Real meshWidth, const CellIndex::Direction whichBorder) const {
        Point point;
        switch (whichBorder) {
            case CellIndex::Direction::NORTH: {
                point = {meshWidth/2.0, 1.0};
                break;
            }
            case CellIndex::Direction::SOUTH: {
                point = {meshWidth/2.0, 1.0 - meshWidth};
                break;
            }

            case CellIndex::Direction::WEST: {
                point = {0.0, 1.0 - meshWidth/2.0};
                break;
            }

            case CellIndex::Direction::EAST: {
                point = {meshWidth, 1.0 - meshWidth/2.0};
                break;
            }
        }

        shiftPoint(point, meshWidth);
        return point;
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


    inline std::vector<CellIndex> neighbors(const int numberOfRows, const int numberOfCols) {
        std::vector<CellIndex> foundNeighbors(4);

        constexpr static std::array<Direction, 4> directionsToCheck = {
              Direction::EAST, Direction::WEST, Direction::NORTH, Direction::SOUTH
        };

        for (const Direction direction: directionsToCheck) {
            if (hasNeighbor(direction, numberOfRows, numberOfCols)) {
                foundNeighbors.push_back(neighbor(direction));
            }
        }

        return foundNeighbors;
    }

    inline CellIndex transpose() const {
        return {j, i};
    }

    static inline CellIndex fromLinearIndex(const int linearIndex, const int numberOfRows) {
        // linearIndex == i + numberOfRows * j
        // i < numberOfRows

        const std::div_t quotientAndRemainder = std::div(linearIndex, numberOfRows);
        return {quotientAndRemainder.rem, quotientAndRemainder.quot};
    }

};

inline std::ostream& operator<<(std::ostream& ostream, const CellIndex& cellindex) {
    ostream << "(" << cellindex.i << ", " << cellindex.j << ")";
    return ostream;
}