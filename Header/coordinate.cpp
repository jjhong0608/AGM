//
// Created by JUNHONG-NIMS on 2021/01/09.
//

#include "coordinate.h"

AGM::Coordinate::Coordinate() = default;

AGM::Coordinate::Coordinate(double x, double y) : x(x), y(y) {}

double AGM::Coordinate::getX() const {
    return x;
}

void AGM::Coordinate::setX(double d) {
    Coordinate::x = d;
}

double AGM::Coordinate::getY() const {
    return y;
}

void AGM::Coordinate::setY(double d) {
    Coordinate::y = d;
}

double AGM::Coordinate::norm() const {
    return sqrt(x * x + y * y);
}

double &AGM::Coordinate::operator[](int i) {
    if (i == 0) {
        return x;
    } else if (i == 1) {
        return y;
    } else {
        printError("Coordinate::operator[]", "input = %d", i);
    }
    return x;
}

const double &AGM::Coordinate::operator[](int i) const {
    if (i == 0) {
        return x;
    } else if (i == 1) {
        return y;
    } else {
        printError("Coordinate::operator[]", "input = %d", i);
    }
    return x;
}

AGM::Coordinate AGM::Coordinate::operator+(const AGM::Coordinate &src) const {
    auto c = Coordinate();
    c[0] = x + src.x;
    c[1] = y + src.y;
    return c;
}

AGM::Coordinate AGM::Coordinate::operator-(const AGM::Coordinate &src) const {
    auto c = Coordinate();
    c[0] = x - src.x;
    c[1] = y - src.y;
    return c;
}

AGM::Coordinate AGM::Coordinate::operator*(double d) const {
    auto c = Coordinate();
    c[0] = x * d;
    c[1] = y * d;
    return c;
}

AGM::Coordinate::~Coordinate() {
    x = 0.0E0;
    y = 0.0E0;
}

bool AGM::Coordinate::operator==(const AGM::Coordinate &rhs) const {
    return x == rhs.x &&
           y == rhs.y;
}

bool AGM::Coordinate::operator!=(const AGM::Coordinate &rhs) const {
    return !(rhs == *this);
}

AGM::Coordinate &AGM::Coordinate::operator=(const AGM::Coordinate &rhs) {
    if (this != &rhs) {
        x = rhs.x;
        y = rhs.y;
    }
    return *this;
}
