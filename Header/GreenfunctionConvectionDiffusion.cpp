//
// Created by NIMS-JUNHONG on 2021/02/08.
//

#include "GreenfunctionConvectionDiffusion.h"

AGM::GreenfunctionConvectionDiffusion::GreenfunctionConvectionDiffusion(double tm, double tau, double tp, double mpl,
                                                                        double mpr) : Greenfunction(tm, tau, tp, mpl,
                                                                                                    mpr) {}

AGM::GreenfunctionConvectionDiffusion::GreenfunctionConvectionDiffusion(double tm, double tau, double tp, double mpl,
                                                                        double mpr, double a)
        : Greenfunction(tm, tau, tp, mpl, mpr), a(a) {}

double AGM::GreenfunctionConvectionDiffusion::g(double u, double v) const {
    double b = std::fabs(a);
    if (v > tau - NEARZERO) return (UNITVALUE - exp(-b * (u - v) / mpr));
    else if (u < tau + NEARZERO) return (UNITVALUE - exp(-b * (u - v) / mpl));
    else return (UNITVALUE - exp(-b * (u - tau) / mpr) * exp(-b * (tau - v) / mpl));
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_square_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return pow(a * t / mpl, 2) + 2.0E0 * a * t / mpl + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            return g(tp, tau) * (b * (pow(tau, 3) - pow(tm, 3)) / 3.0E0 -
                                 pow(mpl / a, 3) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm)))) / g(tp, tm) / a;
        case 'r':
            f = [&](double t) -> double {
                return pow(a * t / mpr, 2) + 2.0E0 * a * t / mpr + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            return g(tau, tm) * (b * (pow(tau, 3) - pow(tp, 3)) / 3.0E0 -
                                 pow(mpr / a, 3) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau)))) / g(tp, tm) / a;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_linear_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return a * t / mpl + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            return g(tp, tau) * (b * (pow(tau, 2) - pow(tm, 2)) / 2.0E0 -
                                 pow(mpl / a, 2) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm)))) / g(tp, tm) / a;
        case 'r':
            f = [&](double t) -> double {
                return a * t / mpr + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            return g(tau, tm) * (b * (pow(tau, 2) - pow(tp, 2)) / 2.0E0 -
                                 pow(mpr / a, 2) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau)))) / g(tp, tm) / a;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_const_zeroconvection(i);
    }
    double b{};
    switch (i) {
        case 'l':
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            return g(tp, tau) * (b * (tau - tm) - mpl / a * (sgn(a) * g(tau, tm))) / g(tp, tm) / a;
        case 'r':
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            return g(tau, tm) * (b * (tau - tp) + mpr / a * (sgn(a) * g(tp, tau))) / g(tp, tm) / a;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square_zeroconvection(char i) const {
    double f{};
    switch (i) {
        case 'l':
            if (a > ZEROVALUE) {
                f = (1.0 / 60.0) *
                    (-a * (6 * pow(tau, 5) - 15 * pow(tau, 4) * tm + 10 * pow(tau, 3) * pow(tm, 2) - pow(tm, 5)) +
                     5 * mpl * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tm + pow(tm, 4))) / pow(mpl, 2);
            } else {
                f = (1.0 / 60.0) * (a * (9 * pow(tau, 5) - 20 * pow(tau, 4) * tm + 10 * pow(tau, 3) * pow(tm, 2) +
                                         5 * tau * pow(tm, 4) - 4 * pow(tm, 5)) +
                                    5 * mpl * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tm + pow(tm, 4))) / pow(mpl, 2);
            }
            return g(tp, tau) * f / g(tp, tm);
        case 'r':
            if (a > ZEROVALUE) {
                f = (1.0 / 60.0) * (a * (9 * pow(tau, 5) - 20 * pow(tau, 4) * tp + 10 * pow(tau, 3) * pow(tp, 2) +
                                         5 * tau * pow(tp, 4) - 4 * pow(tp, 5)) +
                                    5 * mpr * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tp + pow(tp, 4))) / pow(mpr, 2);
            } else {
                f = (1.0 / 60.0) *
                    (-a * (6 * pow(tau, 5) - 15 * pow(tau, 4) * tp + 10 * pow(tau, 3) * pow(tp, 2) - pow(tp, 5)) +
                     5 * mpr * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tp + pow(tp, 4))) / pow(mpr, 2);
            }
            return g(tau, tm) * f / g(tp, tm);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear_zeroconvection(char i) const {
    double f{};
    switch (i) {
        case 'l':
            if (a > ZEROVALUE) {
//                f = (1.0 / 24.0) *
//                    (-a * (3 * pow(tau, 4) - 8 * pow(tau, 3) * tm + 6 * pow(tau, 2) * pow(tm, 2) - pow(tm, 4)) +
//                     4 * mpl * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3))) / pow(mpl, 2);
                f = (1.0 / 120.0) * (pow(a, 2) *
                                     (4 * pow(tau, 5) - 15 * pow(tau, 4) * tm + 20 * pow(tau, 3) * pow(tm, 2) -
                                      10 * pow(tau, 2) * pow(tm, 3) + pow(tm, 5)) - 5 * a * mpl * (3 * pow(tau, 4) -
                                                                                                   8 * pow(tau, 3) *
                                                                                                   tm +
                                                                                                   6 * pow(tau, 2) *
                                                                                                   pow(tm, 2) -
                                                                                                   pow(tm, 4)) +
                                     20 * pow(mpl, 2) * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3))) /
                    pow(mpl, 3);
            } else {
//                f = (1.0 / 24.0) * (a * (5 * pow(tau, 4) - 12 * pow(tau, 3) * tm + 6 * pow(tau, 2) * pow(tm, 2) +
//                                         4 * tau * pow(tm, 3) - 3 * pow(tm, 4)) +
//                                    4 * mpl * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3))) / pow(mpl, 2);
                f = (1.0 / 120.0) * (3 * pow(a, 2) *
                                     (3 * pow(tau, 5) - 10 * pow(tau, 4) * tm + 10 * pow(tau, 3) * pow(tm, 2) -
                                      5 * tau * pow(tm, 4) + 2 * pow(tm, 5)) + 5 * a * mpl * (5 * pow(tau, 4) -
                                                                                              12 * pow(tau, 3) * tm +
                                                                                              6 * pow(tau, 2) *
                                                                                              pow(tm, 2) +
                                                                                              4 * tau * pow(tm, 3) -
                                                                                              3 * pow(tm, 4)) +
                                     20 * pow(mpl, 2) * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3))) /
                    pow(mpl, 3);
            }
            return g(tp, tau) * f / g(tp, tm);
        case 'r':
            if (a > ZEROVALUE) {
//                f = (1.0 / 24.0) * (a * (5 * pow(tau, 4) - 12 * pow(tau, 3) * tp + 6 * pow(tau, 2) * pow(tp, 2) +
//                                         4 * tau * pow(tp, 3) - 3 * pow(tp, 4)) +
//                                    4 * mpr * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3))) / pow(mpr, 2);
                f = (1.0 / 120.0) * (3 * pow(a, 2) *
                                     (3 * pow(tau, 5) - 10 * pow(tau, 4) * tp + 10 * pow(tau, 3) * pow(tp, 2) -
                                      5 * tau * pow(tp, 4) + 2 * pow(tp, 5)) + 5 * a * mpr * (5 * pow(tau, 4) -
                                                                                              12 * pow(tau, 3) * tp +
                                                                                              6 * pow(tau, 2) *
                                                                                              pow(tp, 2) +
                                                                                              4 * tau * pow(tp, 3) -
                                                                                              3 * pow(tp, 4)) +
                                     20 * pow(mpr, 2) * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3))) /
                    pow(mpr, 3);
            } else {
//                f = (1.0 / 24.0) *
//                    (-a * (3 * pow(tau, 4) - 8 * pow(tau, 3) * tp + 6 * pow(tau, 2) * pow(tp, 2) - pow(tp, 4)) +
//                     4 * mpr * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3))) / pow(mpr, 2);
                f = (1.0 / 120.0) * (pow(a, 2) *
                                     (4 * pow(tau, 5) - 15 * pow(tau, 4) * tp + 20 * pow(tau, 3) * pow(tp, 2) -
                                      10 * pow(tau, 2) * pow(tp, 3) + pow(tp, 5)) - 5 * a * mpr * (3 * pow(tau, 4) -
                                                                                                   8 * pow(tau, 3) *
                                                                                                   tp +
                                                                                                   6 * pow(tau, 2) *
                                                                                                   pow(tp, 2) -
                                                                                                   pow(tp, 4)) +
                                     20 * pow(mpr, 2) * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3))) /
                    pow(mpr, 3);
            }
            return g(tau, tm) * f / g(tp, tm);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const_zeroconvection(char i) const {
    double f{};
    switch (i) {
        case 'l':
            if (a > ZEROVALUE) {
//                f = (1.0 / 6.0) * (-a * (pow(tau, 3) - 3 * pow(tau, 2) * tm + 3 * tau * pow(tm, 2) - pow(tm, 3)) +
//                                   3 * mpl * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2))) / pow(mpl, 2);
                f = (1.0 / 24.0) * (pow(a, 2) * (pow(tau, 4) - 4 * pow(tau, 3) * tm + 6 * pow(tau, 2) * pow(tm, 2) -
                                                 4 * tau * pow(tm, 3) + pow(tm, 4)) - 4 * a * mpl * (pow(tau, 3) -
                                                                                                     3 * pow(tau, 2) *
                                                                                                     tm + 3 * tau *
                                                                                                          pow(tm, 2) -
                                                                                                     pow(tm, 3)) +
                                    12 * pow(mpl, 2) * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2))) / pow(mpl, 3);
            } else {
//                f = ((1.0 / 3.0) * a * (pow(tau, 3) - 3 * pow(tau, 2) * tm + 3 * tau * pow(tm, 2) - pow(tm, 3)) +
//                     (1.0 / 2.0) * mpl * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2))) / pow(mpl, 2);
                f = (1.0 / 24.0) * (3 * pow(a, 2) * (pow(tau, 4) - 4 * pow(tau, 3) * tm + 6 * pow(tau, 2) * pow(tm, 2) -
                                                     4 * tau * pow(tm, 3) + pow(tm, 4)) + 8 * a * mpl * (pow(tau, 3) -
                                                                                                         3 *
                                                                                                         pow(tau, 2) *
                                                                                                         tm + 3 * tau *
                                                                                                              pow(tm,
                                                                                                                  2) -
                                                                                                         pow(tm, 3)) +
                                    12 * pow(mpl, 2) * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2))) / pow(mpl, 3);
            }
            return g(tp, tau) * f / g(tp, tm);
        case 'r':
            if (a > ZEROVALUE) {
//                f = ((1.0 / 3.0) * a * (pow(tau, 3) - 3 * pow(tau, 2) * tp + 3 * tau * pow(tp, 2) - pow(tp, 3)) +
//                     (1.0 / 2.0) * mpr * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2))) / pow(mpr, 2);
                f = (1.0 / 24.0) * (3 * pow(a, 2) * (pow(tau, 4) - 4 * pow(tau, 3) * tp + 6 * pow(tau, 2) * pow(tp, 2) -
                                                     4 * tau * pow(tp, 3) + pow(tp, 4)) + 8 * a * mpr * (pow(tau, 3) -
                                                                                                         3 *
                                                                                                         pow(tau, 2) *
                                                                                                         tp + 3 * tau *
                                                                                                              pow(tp,
                                                                                                                  2) -
                                                                                                         pow(tp, 3)) +
                                    12 * pow(mpr, 2) * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2))) / pow(mpr, 3);
            } else {
//                f = (1.0 / 6.0) * (-a * (pow(tau, 3) - 3 * pow(tau, 2) * tp + 3 * tau * pow(tp, 2) - pow(tp, 3)) +
//                                   3 * mpr * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2))) / pow(mpr, 2);
                f = (1.0 / 24.0) * (pow(a, 2) * (pow(tau, 4) - 4 * pow(tau, 3) * tp + 6 * pow(tau, 2) * pow(tp, 2) -
                                                 4 * tau * pow(tp, 3) + pow(tp, 4)) - 4 * a * mpr * (pow(tau, 3) -
                                                                                                     3 * pow(tau, 2) *
                                                                                                     tp + 3 * tau *
                                                                                                          pow(tp, 2) -
                                                                                                     pow(tp, 3)) +
                                    12 * pow(mpr, 2) * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2))) / pow(mpr, 3);
            }
            return g(tau, tm) * f / g(tp, tm);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square_t(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_square_t_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return pow(a * t / mpl, 2) + 2.0E0 * a * t / mpl + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            return g(tp, tau) * pow(mpl / a, 3) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm))) / g(tp, tm) / mpl;
        case 'r':
            f = [&](double t) -> double {
                return pow(a * t / mpr, 2) + 2.0E0 * a * t / mpr + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            return g(tau, tm) * pow(mpr / a, 3) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau))) / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear_t(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_linear_t_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return a * t / mpl + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            return g(tp, tau) * pow(mpl / a, 2) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm))) / g(tp, tm) / mpl;
        case 'r':
            f = [&](double t) -> double {
                return a * t / mpr + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            return g(tau, tm) * pow(mpr / a, 2) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau))) / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const_t(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_const_t_zeroconvection(i);
    }
    std::function<double(double)> f;
    switch (i) {
        case 'l':
            return g(tp, tau) * mpl / a * sgn(a) * g(tau, tm) / g(tp, tm) / mpl;
        case 'r':
            return -g(tau, tm) * mpr / a * sgn(a) * g(tp, tau) / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square_t_zeroconvection(char i) const {
    double f{};
    switch (i) {
        case 'l':
            if (a > ZEROVALUE) {
                f = (1.0 / 12.0) * (-a * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tm + pow(tm, 4)) +
                                    4 * mpl * (pow(tau, 3) - pow(tm, 3))) / mpl;
            } else {
                f = (1.0 / 12.0) *
                    (a * (pow(tau, 4) - 4 * tau * pow(tm, 3) + 3 * pow(tm, 4)) + 4 * mpl * (pow(tau, 3) - pow(tm, 3))) /
                    mpl;
            }
            return g(tp, tau) * f / g(tp, tm) / mpl;
        case 'r':
            if (a > ZEROVALUE) {
                f = (1.0 / 12.0) *
                    (a * (pow(tau, 4) - 4 * tau * pow(tp, 3) + 3 * pow(tp, 4)) + 4 * mpr * (pow(tau, 3) - pow(tp, 3))) /
                    mpr;
            } else {
                f = (1.0 / 12.0) * (-a * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tp + pow(tp, 4)) +
                                    4 * mpr * (pow(tau, 3) - pow(tp, 3))) / mpr;
            }
            return g(tau, tm) * f / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear_t_zeroconvection(char i) const {
    double f{};
    switch (i) {
        case 'l':
            if (a > ZEROVALUE) {
                f = (1.0 / 6.0) * (-a * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3)) +
                                   3 * mpl * (pow(tau, 2) - pow(tm, 2))) / mpl;
            } else {
                f = (1.0 / 6.0) *
                    (a * (pow(tau, 3) - 3 * tau * pow(tm, 2) + 2 * pow(tm, 3)) + 3 * mpl * (pow(tau, 2) - pow(tm, 2))) /
                    mpl;
            }
            return g(tp, tau) * f / g(tp, tm) / mpl;
        case 'r':
            if (a > ZEROVALUE) {
                f = (1.0 / 6.0) *
                    (a * (pow(tau, 3) - 3 * tau * pow(tp, 2) + 2 * pow(tp, 3)) + 3 * mpr * (pow(tau, 2) - pow(tp, 2))) /
                    mpr;
            } else {
                f = (1.0 / 6.0) * (-a * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3)) +
                                   3 * mpr * (pow(tau, 2) - pow(tp, 2))) / mpr;
            }
            return g(tau, tm) * f / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const_t_zeroconvection(char i) const {
    double f{};
    switch (i) {
        case 'l':
            if (a > ZEROVALUE) {
                f = (-1.0 / 2.0 * a * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2)) + mpl * (tau - tm)) / mpl;
            } else {
                f = ((1.0 / 2.0) * a * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2)) + mpl * (tau - tm)) / mpl;
            }
            return g(tp, tau) * f / g(tp, tm) / mpl;
        case 'r':
            if (a > ZEROVALUE) {
                f = ((1.0 / 2.0) * a * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2)) + mpr * (tau - tp)) / mpr;
            } else {
                f = (-1.0 / 2.0 * a * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2)) + mpr * (tau - tp)) / mpr;
            }
            return g(tau, tm) * f / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square_tau(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_square_tau_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{}, c{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return pow(a * t / mpl, 2) + 2.0E0 * a * t / mpl + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            return c * (b * (pow(tau, 3) - pow(tm, 3)) / 3.0E0 -
                        pow(mpl / a, 3) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm)))) / g(tp, tm) / mpl;
        case 'r':
            f = [&](double t) -> double {
                return pow(a * t / mpr, 2) + 2.0E0 * a * t / mpr + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            return c * (b * (pow(tau, 3) - pow(tp, 3)) / 3.0E0 -
                        pow(mpr / a, 3) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau)))) / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear_tau(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_linear_tau_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{}, c{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return a * t / mpl + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            return c * (b * (pow(tau, 2) - pow(tm, 2)) / 2.0E0 -
                        pow(mpl / a, 2) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm)))) / g(tp, tm) / mpl;
        case 'r':
            f = [&](double t) -> double {
                return a * t / mpr + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            return c * (b * (pow(tau, 2) - pow(tp, 2)) / 2.0E0 -
                        pow(mpr / a, 2) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau)))) / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const_tau(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_const_tau_zeroconvection(i);
    }
    double b{}, c{};
    switch (i) {
        case 'l':
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            return c * (b * (tau - tm) - mpl / a * (sgn(a) * g(tau, tm))) / g(tp, tm) / mpl;
        case 'r':
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            return c * (b * (tau - tp) + mpr / a * (sgn(a) * g(tp, tau))) / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square_tau_zeroconvection(char i) const {
    double f{}, c{};
    switch (i) {
        case 'l':
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            f = (1.0 / 12.0) * a * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tm + pow(tm, 4)) / mpl;
            return c * f / g(tp, tm) / mpl;
        case 'r':
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            f = (1.0 / 12.0) * a * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tp + pow(tp, 4)) / mpr;
            return c * f / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear_tau_zeroconvection(char i) const {
    double f{}, c{};
    switch (i) {
        case 'l':
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
//            f = (1.0 / 6.0) * a * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3)) / mpl;
            if (ispositive(a)) {
                f = (1.0 / 24.0) * a *
                    (-a * (3 * pow(tau, 4) - 8 * pow(tau, 3) * tm + 6 * pow(tau, 2) * pow(tm, 2) - pow(tm, 4)) +
                     4 * mpl * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3))) / pow(mpl, 2);
            } else {
                f = (1.0 / 24.0) * a * (a * (5 * pow(tau, 4) - 12 * pow(tau, 3) * tm + 6 * pow(tau, 2) * pow(tm, 2) +
                                             4 * tau * pow(tm, 3) - 3 * pow(tm, 4)) +
                                        4 * mpl * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3))) / pow(mpl, 2);
            }
            return c * f / g(tp, tm) / mpl;
        case 'r':
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
//            f = (1.0 / 6.0) * a * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3)) / mpr;
            if (ispositive(a)) {
                f = (1.0 / 24.0) * a * (a * (5 * pow(tau, 4) - 12 * pow(tau, 3) * tp + 6 * pow(tau, 2) * pow(tp, 2) +
                                             4 * tau * pow(tp, 3) - 3 * pow(tp, 4)) +
                                        4 * mpr * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3))) / pow(mpr, 2);
            } else {
                f = (1.0 / 24.0) * a *
                    (-a * (3 * pow(tau, 4) - 8 * pow(tau, 3) * tp + 6 * pow(tau, 2) * pow(tp, 2) - pow(tp, 4)) +
                     4 * mpr * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3))) / pow(mpr, 2);
            }
            return c * f / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const_tau_zeroconvection(char i) const {
    double f{}, c{};
    switch (i) {
        case 'l':
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
//            f = (1.0 / 2.0) * a * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2)) / mpl;
            if (ispositive(a)) {
                f = (1.0 / 6.0) * a * (-a * (pow(tau, 3) - 3 * pow(tau, 2) * tm + 3 * tau * pow(tm, 2) - pow(tm, 3)) +
                                       3 * mpl * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2))) / pow(mpl, 2);
            } else {
                f = (1.0 / 6.0) * a *
                    (2 * a * (pow(tau, 3) - 3 * pow(tau, 2) * tm + 3 * tau * pow(tm, 2) - pow(tm, 3)) +
                     3 * mpl * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2))) / pow(mpl, 2);
            }
            return c * f / g(tp, tm) / mpl;
        case 'r':
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
//            f = (1.0 / 2.0) * a * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2)) / mpr;
            if (ispositive(a)) {
                f = (1.0 / 6.0) * a *
                    (2 * a * (pow(tau, 3) - 3 * pow(tau, 2) * tp + 3 * tau * pow(tp, 2) - pow(tp, 3)) +
                     3 * mpr * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2))) / pow(mpr, 2);
            } else {
                f = (1.0 / 6.0) * a * (-a * (pow(tau, 3) - 3 * pow(tau, 2) * tp + 3 * tau * pow(tp, 2) - pow(tp, 3)) +
                                       3 * mpr * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2))) / pow(mpr, 2);
            }
            return c * f / g(tp, tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square_ttau(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_square_ttau_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{}, c{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return pow(a * t / mpl, 2) + 2.0E0 * a * t / mpl + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            return a * c * pow(mpl / a, 3) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm))) / g(tp, tm) / pow(mpl, 2);
        case 'r':
            f = [&](double t) -> double {
                return pow(a * t / mpr, 2) + 2.0E0 * a * t / mpr + 2.0E0;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            return a * c * pow(mpr / a, 3) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau))) / g(tp, tm) / pow(mpr, 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear_ttau(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_linear_ttau_zeroconvection(i);
    }
    std::function<double(double)> f;
    double b{}, c{};
    switch (i) {
        case 'l':
            f = [&](double t) -> double {
                return a * t / mpl + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE - g(tau, tm) : UNITVALUE;
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            return a * c * pow(mpl / a, 2) * (b * f(tm) - f(tau) * (b - sgn(a) * g(tau, tm))) / g(tp, tm) / pow(mpl, 2);
        case 'r':
            f = [&](double t) -> double {
                return a * t / mpr + UNITVALUE;
            };
            b = std::signbit(a) ? UNITVALUE : UNITVALUE - g(tp, tau);
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            return a * c * pow(mpr / a, 2) * (b * f(tp) - f(tau) * (b + sgn(a) * g(tp, tau))) / g(tp, tm) / pow(mpr, 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const_ttau(char i) const {
    if (fabs(a / mpl) < 1e-2) {
        return integrate_const_ttau_zeroconvection(i);
    }
    std::function<double(double)> f;
    double c{};
    switch (i) {
        case 'l':
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            return c * mpl * (sgn(a) * g(tau, tm)) / g(tp, tm) / pow(mpl, 2);
        case 'r':
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            return -c * mpr * (sgn(a) * g(tp, tau)) / g(tp, tm) / pow(mpr, 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_square_ttau_zeroconvection(char i) const {
    double f{}, c{};
    switch (i) {
        case 'l':
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            if (a > ZEROVALUE) {
                f = (1.0 / 12.0) * (-a * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tm + pow(tm, 4)) +
                                    4 * mpl * (pow(tau, 3) - pow(tm, 3))) / mpl;
            } else {
                f = (1.0 / 12.0) *
                    (a * (pow(tau, 4) - 4 * tau * pow(tm, 3) + 3 * pow(tm, 4)) + 4 * mpl * (pow(tau, 3) - pow(tm, 3))) /
                    mpl;
            }
            return a * c * f / g(tp, tm) / pow(mpl, 2);
        case 'r':
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            if (a > ZEROVALUE) {
                f = (1.0 / 12.0) *
                    (a * (pow(tau, 4) - 4 * tau * pow(tp, 3) + 3 * pow(tp, 4)) + 4 * mpr * (pow(tau, 3) - pow(tp, 3))) /
                    mpr;
            } else {
                f = (1.0 / 12.0) * (-a * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tp + pow(tp, 4)) +
                                    4 * mpr * (pow(tau, 3) - pow(tp, 3))) / mpr;
            }
            return a * c * f / g(tp, tm) / pow(mpr, 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_linear_ttau_zeroconvection(char i) const {
    double f{}, c{};
    switch (i) {
        case 'l':
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            if (a > ZEROVALUE) {
                f = (1.0 / 6.0) * (-a * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3)) +
                                   3 * mpl * (pow(tau, 2) - pow(tm, 2))) / mpl;
            } else {
                f = (1.0 / 6.0) *
                    (a * (pow(tau, 3) - 3 * tau * pow(tm, 2) + 2 * pow(tm, 3)) + 3 * mpl * (pow(tau, 2) - pow(tm, 2))) /
                    mpl;
            }
            return a * c * f / g(tp, tm) / pow(mpl, 2);
        case 'r':
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            if (a > ZEROVALUE) {
                f = (1.0 / 6.0) *
                    (a * (pow(tau, 3) - 3 * tau * pow(tp, 2) + 2 * pow(tp, 3)) + 3 * mpr * (pow(tau, 2) - pow(tp, 2))) /
                    mpr;
            } else {
                f = (1.0 / 6.0) * (-a * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3)) +
                                   3 * mpr * (pow(tau, 2) - pow(tp, 2))) / mpr;
            }
            return a * c * f / g(tp, tm) / pow(mpr, 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::integrate_const_ttau_zeroconvection(char i) const {
    double f{}, c{};
    switch (i) {
        case 'l':
            c = std::signbit(a) ? UNITVALUE : g(tp, tau) - UNITVALUE;
            if (a > ZEROVALUE) {
                f = (-1.0 / 2.0 * a * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2)) + mpl * (tau - tm)) / mpl;
            } else {
                f = ((1.0 / 2.0) * a * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2)) + mpl * (tau - tm)) / mpl;
            }
            return a * c * f / g(tp, tm) / pow(mpl, 2);
        case 'r':
            c = std::signbit(a) ? g(tau, tm) - UNITVALUE : UNITVALUE;
            if (a > ZEROVALUE) {
                f = ((1.0 / 2.0) * a * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2)) + mpr * (tau - tp)) / mpr;
            } else {
                f = (-1.0 / 2.0 * a * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2)) + mpr * (tau - tp)) / mpr;
            }
            return a * c * f / g(tp, tm) / pow(mpr, 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_function(double t) const {
    if (ispositive(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return g(t, tm) * g(tp, tau) / g(tp, tm) / a;
        } else {
            return (UNITVALUE - g(t, tau)) * g(tau, tm) * g(tp, t) / g(tp, tm) / a;
        }
    } else if (isnegative(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return -(UNITVALUE - g(tau, t)) * g(t, tm) * g(tp, tau) / g(tp, tm) / a;
        } else {
            return -g(tau, tm) * g(tp, t) / g(tp, tm) / a;
        }
    } else {
        return Greenfunction::green_function(t);
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_function_t(double t) const {
    if (ispositive(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return (UNITVALUE - g(t, tm)) * g(tp, tau) / g(tp, tm) / mpl;
        } else {
            return (g(t, tau) - UNITVALUE) * g(tau, tm) / g(tp, tm) / mpr;
        }
    } else if (isnegative(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return (UNITVALUE - g(tau, t)) * g(tp, tau) / g(tp, tm) / mpl;
        } else {
            return (g(tp, t) - UNITVALUE) * g(tau, tm) / g(tp, tm) / mpr;
        }
    } else {
        return Greenfunction::green_function_t(t);
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_function_tau(double t) const {
    if (ispositive(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return (g(tp, tau) - UNITVALUE) * g(t, tm) / g(tp, tm) / mpl;
        } else {
            return (UNITVALUE - g(t, tau)) * g(tp, t) / g(tp, tm) / mpr;
        }
    } else if (isnegative(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return (g(tau, t) - UNITVALUE) * g(t, tm) / g(tp, tm) / mpl;
        } else {
            return (UNITVALUE - g(tau, tm)) * g(tp, t) / g(tp, tm) / mpr;
        }
    } else {
        return Greenfunction::green_function_tau(t);
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_function_ttau(double t) const {
    if (ispositive(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return -a * (UNITVALUE - g(tp, tau)) * (UNITVALUE - g(t, tm)) / g(tp, tm) / pow(mpl, 2);
        } else {
            return -a * (UNITVALUE - g(t, tau)) / g(tp, tm) / pow(mpr, 2);
        }
    } else if (isnegative(a / mpl)) {
        if (t < tau || isclose(tm, t)) {
            return a * (UNITVALUE - g(tau, t)) / g(tp, tm) / pow(mpl, 2);
        } else {
            return a * (UNITVALUE - g(tp, t)) * (UNITVALUE - g(tau, tm)) / g(tp, tm) / pow(mpr, 2);
        }
    } else {
        return Greenfunction::green_function_ttau(t);
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_integral(char pos) const {
    double cc{}, cl{}, cr{};
    if (fabs(a / mpl) < NEARZERO) {
        return Greenfunction(tm, tau, tp, mpl, mpr).green_integral(pos);
    } else {
        switch (pos) {
            case 'l':
                cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                return cc * integrate_square('l') - (tau + tp) * cc * integrate_linear('l') +
                       tau * tp * cc * integrate_const('l') +
                       cc * integrate_square('r') - (tau + tp) * cc * integrate_linear('r') +
                       tau * tp * cc * integrate_const('r');
            case 'r':
                cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                return cc * integrate_square('l') - (tau + tm) * cc * integrate_linear('l') +
                       tau * tm * cc * integrate_const('l') +
                       cc * integrate_square('r') - (tau + tm) * cc * integrate_linear('r') +
                       tau * tm * cc * integrate_const('r');
            case 'c':
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cc * integrate_square('l') - (tm + tp) * cc * integrate_linear('l') +
                       tm * tp * cc * integrate_const('l') +
                       cc * integrate_square('r') - (tm + tp) * cc * integrate_linear('r') +
                       tm * tp * cc * integrate_const('r');
            case 'L':
                cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cl * integrate_square('l') - (tau + tp) * cl * integrate_linear('l') +
                       tau * tp * cl * integrate_const('l') + cl * integrate_square('r') -
                       (tau + tp) * cl * integrate_linear('r') + tau * tp * cl * integrate_const('r') +
                       cc * integrate_square('l') - (tm + tp) * cc * integrate_linear('l') +
                       tm * tp * cc * integrate_const('l');
            case 'R':
                cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cr * integrate_square('l') - (tau + tm) * cr * integrate_linear('l') +
                       tau * tm * cr * integrate_const('l') + cr * integrate_square('r') -
                       (tau + tm) * cr * integrate_linear('r') + tau * tm * cr * integrate_const('r') +
                       cc * integrate_square('r') - (tm + tp) * cc * integrate_linear('r') +
                       tm * tp * cc * integrate_const('r');
            default:
                return ZEROVALUE;
        }
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_integral_t(char pos) const {
    double cc{}, cl{}, cr{};
    if (fabs(a / mpl) < NEARZERO) {
        return Greenfunction(tm, tau, tp, mpl, mpr).green_integral_t(pos);
    } else {
        switch (pos) {
            case 'l':
                cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                return cc * integrate_square_t('l') - (tau + tp) * cc * integrate_linear_t('l') +
                       tau * tp * cc * integrate_const_t('l') +
                       cc * integrate_square_t('r') - (tau + tp) * cc * integrate_linear_t('r') +
                       tau * tp * cc * integrate_const_t('r');
            case 'r':
                cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                return cc * integrate_square_t('l') - (tau + tm) * cc * integrate_linear_t('l') +
                       tau * tm * cc * integrate_const_t('l') +
                       cc * integrate_square_t('r') - (tau + tm) * cc * integrate_linear_t('r') +
                       tau * tm * cc * integrate_const_t('r');
            case 'c':
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cc * integrate_square_t('l') - (tm + tp) * cc * integrate_linear_t('l') +
                       tm * tp * cc * integrate_const_t('l') +
                       cc * integrate_square_t('r') - (tm + tp) * cc * integrate_linear_t('r') +
                       tm * tp * cc * integrate_const_t('r');
            case 'L':
                cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cl * integrate_square_t('l') - (tau + tp) * cl * integrate_linear_t('l') +
                       tau * tp * cl * integrate_const_t('l') + cl * integrate_square_t('r') -
                       (tau + tp) * cl * integrate_linear_t('r') + tau * tp * cl * integrate_const_t('r') +
                       cc * integrate_square_t('l') - (tm + tp) * cc * integrate_linear_t('l') +
                       tm * tp * cc * integrate_const_t('l');
            case 'R':
                cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cr * integrate_square_t('l') - (tau + tm) * cr * integrate_linear_t('l') +
                       tau * tm * cr * integrate_const_t('l') + cr * integrate_square_t('r') -
                       (tau + tm) * cr * integrate_linear_t('r') + tau * tm * cr * integrate_const_t('r') +
                       cc * integrate_square_t('r') - (tm + tp) * cc * integrate_linear_t('r') +
                       tm * tp * cc * integrate_const_t('r');
            default:
                return ZEROVALUE;
        }
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_integral_tau(char pos) const {
    double cc{}, cl{}, cr{};
    if (fabs(a / mpl) < NEARZERO) {
        return Greenfunction(tm, tau, tp, mpl, mpr).green_integral_tau(pos);
    } else {
        switch (pos) {
            case 'l':
                cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                return cc * integrate_square_tau('l') - (tau + tp) * cc * integrate_linear_tau('l') +
                       tau * tp * cc * integrate_const_tau('l') +
                       cc * integrate_square_tau('r') - (tau + tp) * cc * integrate_linear_tau('r') +
                       tau * tp * cc * integrate_const_tau('r');
            case 'r':
                cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                return cc * integrate_square_tau('l') - (tau + tm) * cc * integrate_linear_tau('l') +
                       tau * tm * cc * integrate_const_tau('l') +
                       cc * integrate_square_tau('r') - (tau + tm) * cc * integrate_linear_tau('r') +
                       tau * tm * cc * integrate_const_tau('r');
            case 'c':
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cc * integrate_square_tau('l') - (tm + tp) * cc * integrate_linear_tau('l') +
                       tm * tp * cc * integrate_const_tau('l') +
                       cc * integrate_square_tau('r') - (tm + tp) * cc * integrate_linear_tau('r') +
                       tm * tp * cc * integrate_const_tau('r');
            case 'L':
                cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cl * integrate_square_tau('l') - (tau + tp) * cl * integrate_linear_tau('l') +
                       tau * tp * cl * integrate_const_tau('l') + cl * integrate_square_tau('r') -
                       (tau + tp) * cl * integrate_linear_tau('r') + tau * tp * cl * integrate_const_tau('r') +
                       cc * integrate_square_tau('l') - (tm + tp) * cc * integrate_linear_tau('l') +
                       tm * tp * cc * integrate_const_tau('l');
            case 'R':
                cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cr * integrate_const_tau('l') - (tau + tm) * cr * integrate_linear_tau('l') +
                       tau * tm * cr * integrate_const_tau('l') + cr * integrate_square_tau('r') -
                       (tau + tm) * cr * integrate_linear_tau('r') + tau * tm * cr * integrate_const_tau('r') +
                       cc * integrate_square_tau('r') - (tm + tp) * cc * integrate_linear_tau('r') +
                       tm * tp * cc * integrate_const_tau('r');
            default:
                return ZEROVALUE;
        }
    }
}

double AGM::GreenfunctionConvectionDiffusion::green_integral_ttau(char pos) const {
    double cc{}, cl{}, cr{};
    if (fabs(a / mpl) < NEARZERO) {
        return Greenfunction(tm, tau, tp, mpl, mpr).green_integral_ttau(pos);
    } else {
        switch (pos) {
            case 'l':
                cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                return cc * integrate_square_ttau('l') - (tau + tp) * cc * integrate_linear_ttau('l') +
                       tau * tp * cc * integrate_const_ttau('l') +
                       cc * integrate_square_ttau('r') - (tau + tp) * cc * integrate_linear_ttau('r') +
                       tau * tp * cc * integrate_const_ttau('r');
            case 'r':
                cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                return cc * integrate_square_ttau('l') - (tau + tm) * cc * integrate_linear_ttau('l') +
                       tau * tm * cc * integrate_const_ttau('l') +
                       cc * integrate_square_ttau('r') - (tau + tm) * cc * integrate_linear_ttau('r') +
                       tau * tm * cc * integrate_const_ttau('r');
            case 'c':
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cc * integrate_square_ttau('l') - (tm + tp) * cc * integrate_linear_ttau('l') +
                       tm * tp * cc * integrate_const_ttau('l') +
                       cc * integrate_square_ttau('r') - (tm + tp) * cc * integrate_linear_ttau('r') +
                       tm * tp * cc * integrate_const_ttau('r');
            case 'L':
                cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cl * integrate_square_ttau('l') - (tau + tp) * cl * integrate_linear_ttau('l') +
                       tau * tp * cl * integrate_const_ttau('l') + cl * integrate_square_ttau('r') -
                       (tau + tp) * cl * integrate_linear_ttau('r') + tau * tp * cl * integrate_const_ttau('r') +
                       cc * integrate_square_ttau('l') - (tm + tp) * cc * integrate_linear_ttau('l') +
                       tm * tp * cc * integrate_const_ttau('l');
            case 'R':
                cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
                cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
                return cr * integrate_const_ttau('l') - (tau + tm) * cr * integrate_linear_ttau('l') +
                       tau * tm * cr * integrate_const_ttau('l') + cr * integrate_square_ttau('r') -
                       (tau + tm) * cr * integrate_linear_ttau('r') + tau * tm * cr * integrate_const_ttau('r') +
                       cc * integrate_square_ttau('r') - (tm + tp) * cc * integrate_linear_ttau('r') +
                       tm * tp * cc * integrate_const_ttau('r');
            default:
                return ZEROVALUE;
        }
    }
}



