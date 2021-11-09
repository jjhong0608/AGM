//
// Created by 조준홍 on 2021/03/05.
//

#include "function2D.h"

double AGM::function2D::u() const {
    // Kalman vortex
    double a{HALFVALUE};
    if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
        return ZEROVALUE;
    }
    return UNITVALUE - pow(a, 2) / (pow(x, 2) + pow(y, 2)) + 2 * pow(a * y, 2) / pow(pow(x, 2) + pow(y, 2), 2);
    // BFS
//    if (iszero(x) && y > HALFVALUE) {
//        return 24.0 * (y - HALFVALUE) * (UNITVALUE - y);
//    } else {
//        return ZEROVALUE;
//    }
    // cavity
    if (isclose(y, UNITVALUE)) {
        return UNITVALUE;
    } else {
        return ZEROVALUE;
    }
    // Taylor-Green vortex
    return -exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * y) * cos(2 * M_PI * x);
}

double AGM::function2D::v() const {
    // Kalman vortex
//    double a{HALFVALUE};
//    if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
//        return ZEROVALUE;
//    }
//    return -2 * pow(a, 2) * x * y / pow(pow(x, 2) + pow(y, 2), 2);
    return ZEROVALUE;
    return exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * x) * cos(2 * M_PI * y);
}

double AGM::function2D::p() const {
    return ZEROVALUE;
    return -1.0 / 4.0 * (cos(4 * M_PI * x) + cos(4 * M_PI * y)) * exp(-16.0 * pow(M_PI, 2) * eps * time);
}

double AGM::function2D::ut() const {
    return ZEROVALUE;
    return 8 * pow(M_PI, 2) * eps * exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * y) * cos(2 * M_PI * x);
}

double AGM::function2D::vt() const {
    return ZEROVALUE;
    return -8 * pow(M_PI, 2) * eps * exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * x) * cos(2 * M_PI * y);
}

double AGM::function2D::ux() const {
    return ZEROVALUE;
    return 2 * M_PI * exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

double AGM::function2D::uy() const {
    return ZEROVALUE;
    return -2 * M_PI * exp(-8 * pow(M_PI, 2) * eps * time) * cos(2 * M_PI * x) * cos(2 * M_PI * y);
}

double AGM::function2D::vx() const {
    return ZEROVALUE;
    return 2 * M_PI * exp(-8 * pow(M_PI, 2) * eps * time) * cos(2 * M_PI * x) * cos(2 * M_PI * y);
}

double AGM::function2D::vy() const {
    return ZEROVALUE;
    return -2 * M_PI * exp(-8 * pow(M_PI, 2) * eps * time) * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

double AGM::function2D::px() const {
    return ZEROVALUE;
    return M_PI * exp(-16 * pow(M_PI, 2) * eps * time) * sin(4 * M_PI * x);
}

double AGM::function2D::py() const {
    return ZEROVALUE;
    return M_PI * exp(-16 * pow(M_PI, 2) * eps * time) * sin(4 * M_PI * y);
}

double AGM::function2D::f1() const {
//    return -M_PI * exp(-16 * pow(M_PI, 2) * eps * time) * sin(4 * M_PI * x);
    return ZEROVALUE;
}

double AGM::function2D::f2() const {
//    return -M_PI * exp(-16 * pow(M_PI, 2) * eps * time) * sin(4 * M_PI * y);
    return ZEROVALUE;
}

double AGM::function2D::uxx() const {
    return ZEROVALUE;
}

double AGM::function2D::vyy() const {
    return ZEROVALUE;
}

double AGM::function2D::phi() const {
    return ZEROVALUE;
    return (1.0 / 2.0) * M_PI *
           (sin(4 * M_PI * x) + (1.0 / 2.0) * sin(M_PI * (4 * x - 4 * y)) + (1.0 / 2.0) * sin(M_PI * (4 * x + 4 * y))) *
           exp(-16 * pow(M_PI, 2) * eps * time);
}

double AGM::function2D::psi() const {
    return ZEROVALUE;
    return (1.0 / 2.0) * M_PI *
           (-sin(4 * M_PI * y) + (1.0 / 2.0) * sin(M_PI * (4 * x - 4 * y)) - 1.0 / 2.0 * sin(M_PI * (4 * x + 4 * y))) *
           exp(-16 * pow(M_PI, 2) * eps * time);
}


