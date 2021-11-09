//
// Created by JUNHONG-NIMS on 2021/01/09.
//

#ifndef AGM_COORDINATE_H
#define AGM_COORDINATE_H

#include "matrixrow.h"

namespace AGM {
    class Coordinate {
    private:
        double x{}, y{};
    public:
        Coordinate();

        virtual ~Coordinate();

        Coordinate(double x, double y);

        [[nodiscard]] double getX() const;

        [[nodiscard]] double getY() const;

        void setX(double d);

        void setY(double d);

        [[nodiscard]] double norm() const;

        double &operator[](int i);

        const double &operator[](int i) const;

        Coordinate operator+(const Coordinate &src) const;

        Coordinate operator-(const Coordinate &src) const;

        Coordinate operator*(double d) const;

        bool operator==(const Coordinate &rhs) const;

        bool operator!=(const Coordinate &rhs) const;

        Coordinate &operator=(const Coordinate &rhs);
    };
}


#endif //AGM_COORDINATE_H
