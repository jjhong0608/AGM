//
// Created by NIMS-JUNHONG on 2020/12/24.
//

#include "function.h"

AGM::function::function() = default;

AGM::function::function(const AGM::Point &pt) {
    x = pt[0];
    y = pt[1];
    eps = pt.getMp();
}

AGM::function::function(const AGM::PointHeat &pt) {
    x = pt[0];
    y = pt[1];
    eps = pt.getMp();
    time = AGM::PointHeat::getTime();
}

AGM::function::function(double x, double y) : x(x), y(y) {}

AGM::function::function(double eps, double x, double y) : eps(eps), x(x), y(y) {}

AGM::function::~function() = default;

void AGM::function::setPoint(const AGM::Point &pt) {
    x = pt[0];
    y = pt[1];
    eps = pt.getMp();
}

void AGM::function::setPoint(const AGM::PointHeat &pt) {
    x = pt[0];
    y = pt[1];
    eps = pt.getMp();
//    time = AGM::PointHeat::getTime();
}

double AGM::function::getEps() const {
    return eps;
}

void AGM::function::setEps(double d) {
    function::eps = d;
}

double AGM::function::getX() const {
    return x;
}

void AGM::function::setX(double d) {
    function::x = d;
}

double AGM::function::getY() const {
    return y;
}

void AGM::function::setY(double d) {
    function::y = d;
}

double AGM::function::getTime() const {
    return time;
}

void AGM::function::setTime(double d) {
    function::time = d;
}

double AGM::function::u() const {
//    for IEEE
    double r{std::sqrt(x * x + y * y)};
    return std::log(r);

    if (isclose(x, 15.0)) {
        return (1.5 - y) * y * y;
    } else if (isclose(x, 0.0) && y > 0.5) {
        return -2 * y * (4 * y * y - 9 * y + 6) + 2.5;
    } else if (isclose(y, UNITVALUE)) {
        return HALFVALUE;
    } else {
        return ZEROVALUE;
    }

    return -exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * y) * cos(2 * M_PI * x);
    double Tm = 1e1, T1 = 1e2;
    double W = UNITVALUE, H = UNITVALUE;
    return Tm * std::sinh(M_PI * y / W) / std::sinh(M_PI * H / W) * std::sin(M_PI * x / W) + T1;
}

double AGM::function::phi() const {
    return ZEROVALUE;
}

double AGM::function::f() const {
    return ZEROVALUE;
}

double AGM::function::dx() const {
//    for IEEE
    return x / (x * x + y * y);
    return 2 * M_PI * exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

double AGM::function::dy() const {
//    for IEEE
    return y / (x * x + y * y);
    return -2 * M_PI * exp(-8 * pow(M_PI, 2) * eps * time) * cos(2 * M_PI * x) * cos(2 * M_PI * y);
}

double AGM::function::dxx() const {
    return ZEROVALUE;
}

double AGM::function::dyy() const {
    return ZEROVALUE;
}

double AGM::function::grad() const {
    return std::sqrt(dx() * dx() + dy() * dy());
}
