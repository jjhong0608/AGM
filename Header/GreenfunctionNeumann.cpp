//
// Created by NIMS-JUNHONG on 2021/06/16.
//

#include "GreenfunctionNeumann.h"

AGM::GreenfunctionNeumann::GreenfunctionNeumann(double tm, double tau, double tp, double mpl, double mpr)
        : Greenfunction(tm, tau, tp, mpl, mpr) {}

double AGM::GreenfunctionNeumann::integrate_square_ND(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 3.0) * (pow(tau, 3) * (-tau + tp) + pow(tm, 3) * (tau - tp)) / mpr;
        case 'r':
            return (1.0 / 12.0) * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tp + pow(tp, 4)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_linear_ND(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 2.0) * (pow(tau, 2) * (-tau + tp) + pow(tm, 2) * (tau - tp)) / mpr;
        case 'r':
            return (1.0 / 6.0) * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tp + pow(tp, 3)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_const_ND(char i) const {
    switch (i) {
        case 'l':
            return (-tau + tm) * (tau - tp) / mpr;
        case 'r':
            return (1.0 / 2.0) * (pow(tau, 2) - 2 * tau * tp + pow(tp, 2)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_square_t_ND(char i) const {
    switch (i) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (1.0 / 3.0) * (pow(tau, 3) - pow(tp, 3)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_linear_t_ND(char i) const {
    switch (i) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (1.0 / 2.0) * (pow(tau, 2) - pow(tp, 2)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_const_t_ND(char i) const {
    switch (i) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (tau - tp) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_square_tau_ND(char i) const {
    switch (i) {
        case 'r':
            return ZEROVALUE;
        case 'l':
            return (1.0 / 3.0) * (pow(tm, 3) - pow(tp, 3)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_linear_tau_ND(char i) const {
    switch (i) {
        case 'r':
            return ZEROVALUE;
        case 'l':
            return (1.0 / 2.0) * (pow(tm, 2) - pow(tp, 2)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_const_tau_ND(char i) const {
    switch (i) {
        case 'r':
            return ZEROVALUE;
        case 'l':
            return (tm - tp) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_square_ttau_ND(char i) {
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::integrate_linear_ttau_ND(char i) {
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::integrate_const_ttau_ND(char i) {
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::green_function_ND(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (-tau + tp) / mpr;
    } else {
        return (-t + tp) / mpr;
    }
}

double AGM::GreenfunctionNeumann::green_function_t_ND(double t) const {
    if (t < tau || isclose(tm, t)) {
        return ZEROVALUE;
    } else {
        return -1 / mpr;
    }
}

double AGM::GreenfunctionNeumann::green_function_tau_ND(double t) const {
    if (t < tau || isclose(tm, t)) {
        return -1 / mpr;
    } else {
        return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_function_ttau_ND(double t) const {
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::green_integral_ND(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_linear_ND(pos);
    } else {
        return green_integral_square_ND(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_t_ND(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_t_linear_ND(pos);
    } else {
        return green_integral_t_square_ND(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_tau_ND(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_tau_linear_ND(pos);
    } else {
        return green_integral_tau_square_ND(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_ttau_ND(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_ttau_linear_ND(pos);
    } else {
        return green_integral_ttau_square_ND(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_square_ND(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_ND('l') - (tau + tp) * integrate_linear_ND('l') +
                         tau * tp * integrate_const_ND('l') +
                         integrate_square_ND('r') - (tau + tp) * integrate_linear_ND('r') +
                         tau * tp * integrate_const_ND('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_ND('l') - (tau + tm) * integrate_linear_ND('l') +
                         tau * tm * integrate_const_ND('l') +
                         integrate_square_ND('r') - (tau + tm) * integrate_linear_ND('r') +
                         tau * tm * integrate_const_ND('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_ND('l') - (tm + tp) * integrate_linear_ND('l') +
                         tm * tp * integrate_const_ND('l') +
                         integrate_square_ND('r') - (tm + tp) * integrate_linear_ND('r') +
                         tm * tp * integrate_const_ND('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_ND('l') - (tau + tp) * integrate_linear_ND('l') +
                         tau * tp * integrate_const_ND('l') +
                         integrate_square_ND('r') - (tau + tp) * integrate_linear_ND('r') +
                         tau * tp * integrate_const_ND('r')) +
                   cc * (integrate_square_ND('l') - (tm + tp) * integrate_linear_ND('l') +
                         tm * tp * integrate_const_ND('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_ND('l') - (tau + tm) * integrate_linear_ND('l') +
                         tau * tm * integrate_const_ND('l') +
                         integrate_square_ND('r') - (tau + tm) * integrate_linear_ND('r') +
                         tau * tm * integrate_const_ND('r')) +
                   cc * (integrate_square_ND('r') - (tm + tp) * integrate_linear_ND('r') +
                         tm * tp * integrate_const_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_t_square_ND(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc *
                   (integrate_square_t_ND('l') - (tau + tp) * integrate_linear_t_ND('l') +
                    tau * tp * integrate_const_t_ND('l') +
                    integrate_square_t_ND('r') - (tau + tp) * integrate_linear_t_ND('r') +
                    tau * tp * integrate_const_t_ND('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc *
                   (integrate_square_t_ND('l') - (tau + tm) * integrate_linear_t_ND('l') +
                    tau * tm * integrate_const_t_ND('l') +
                    integrate_square_t_ND('r') - (tau + tm) * integrate_linear_t_ND('r') +
                    tau * tm * integrate_const_t_ND('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc *
                   (integrate_square_t_ND('l') - (tm + tp) * integrate_linear_t_ND('l') +
                    tm * tp * integrate_const_t_ND('l') +
                    integrate_square_t_ND('r') - (tm + tp) * integrate_linear_t_ND('r') +
                    tm * tp * integrate_const_t_ND('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl *
                   (integrate_square_t_ND('l') - (tau + tp) * integrate_linear_t_ND('l') +
                    tau * tp * integrate_const_t_ND('l') +
                    integrate_square_t_ND('r') - (tau + tp) * integrate_linear_t_ND('r') +
                    tau * tp * integrate_const_t_ND('r')) +
                   cc *
                   (integrate_square_t_ND('l') - (tm + tp) * integrate_linear_t_ND('l') +
                    tm * tp * integrate_const_t_ND('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr *
                   (integrate_square_t_ND('l') - (tau + tm) * integrate_linear_t_ND('l') +
                    tau * tm * integrate_const_t_ND('l') +
                    integrate_square_t_ND('r') - (tau + tm) * integrate_linear_t_ND('r') +
                    tau * tm * integrate_const_t_ND('r')) +
                   cc *
                   (integrate_square_t_ND('r') - (tm + tp) * integrate_linear_t_ND('r') +
                    tm * tp * integrate_const_t_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_tau_square_ND(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_tau_ND('l') - (tau + tp) * integrate_linear_tau_ND('l') +
                         tau * tp * integrate_const_tau_ND('l') +
                         integrate_square_tau_ND('r') - (tau + tp) * integrate_linear_tau_ND('r') +
                         tau * tp * integrate_const_tau_ND('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_tau_ND('l') - (tau + tm) * integrate_linear_tau_ND('l') +
                         tau * tm * integrate_const_tau_ND('l') +
                         integrate_square_tau_ND('r') - (tau + tm) * integrate_linear_tau_ND('r') +
                         tau * tm * integrate_const_tau_ND('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_tau_ND('l') - (tm + tp) * integrate_linear_tau_ND('l') +
                         tm * tp * integrate_const_tau_ND('l') +
                         integrate_square_tau_ND('r') - (tm + tp) * integrate_linear_tau_ND('r') +
                         tm * tp * integrate_const_tau_ND('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_tau_ND('l') - (tau + tp) * integrate_linear_tau_ND('l') +
                         tau * tp * integrate_const_tau_ND('l') +
                         integrate_square_tau_ND('r') - (tau + tp) * integrate_linear_tau_ND('r') +
                         tau * tp * integrate_const_tau_ND('r')) +
                   cc * (integrate_square_tau_ND('l') - (tm + tp) * integrate_linear_tau_ND('l') +
                         tm * tp * integrate_const_tau_ND('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_tau_ND('l') - (tau + tm) * integrate_linear_tau_ND('l') +
                         tau * tm * integrate_const_tau_ND('l') +
                         integrate_square_tau_ND('r') - (tau + tm) * integrate_linear_tau_ND('r') +
                         tau * tm * integrate_const_tau_ND('r')) +
                   cc * (integrate_square_tau_ND('r') - (tm + tp) * integrate_linear_tau_ND('r') +
                         tm * tp * integrate_const_tau_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_ttau_square_ND(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_ttau_ND('l') - (tau + tp) * integrate_linear_ttau_ND('l') +
                         tau * tp * integrate_const_ttau_ND('l') +
                         integrate_square_ttau_ND('r') - (tau + tp) * integrate_linear_ttau_ND('r') +
                         tau * tp * integrate_const_ttau_ND('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_ttau_ND('l') - (tau + tm) * integrate_linear_ttau_ND('l') +
                         tau * tm * integrate_const_ttau_ND('l') +
                         integrate_square_ttau_ND('r') - (tau + tm) * integrate_linear_ttau_ND('r') +
                         tau * tm * integrate_const_ttau_ND('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_ttau_ND('l') - (tm + tp) * integrate_linear_ttau_ND('l') +
                         tm * tp * integrate_const_ttau_ND('l') +
                         integrate_square_ttau_ND('r') - (tm + tp) * integrate_linear_ttau_ND('r') +
                         tm * tp * integrate_const_ttau_ND('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_ttau_ND('l') - (tau + tp) * integrate_linear_ttau_ND('l') +
                         tau * tp * integrate_const_ttau_ND('l') +
                         integrate_square_ttau_ND('r') - (tau + tp) * integrate_linear_ttau_ND('r') +
                         tau * tp * integrate_const_ttau_ND('r')) +
                   cc * (integrate_square_ttau_ND('l') - (tm + tp) * integrate_linear_ttau_ND('l') +
                         tm * tp * integrate_const_ttau_ND('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_ttau_ND('l') - (tau + tm) * integrate_linear_ttau_ND('l') +
                         tau * tm * integrate_const_ttau_ND('l') +
                         integrate_square_ttau_ND('r') - (tau + tm) * integrate_linear_ttau_ND('r') +
                         tau * tm * integrate_const_ttau_ND('r')) +
                   cc * (integrate_square_ttau_ND('r') - (tm + tp) * integrate_linear_ttau_ND('r') +
                         tm * tp * integrate_const_ttau_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_linear_ND(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_ND('l') - tau * integrate_const_ND('l'));
        case 'r':
            return -cr * (integrate_linear_ND('r') - tau * integrate_const_ND('r'));
        case 'c':
            return cl * (integrate_linear_ND('l') - tm * integrate_const_ND('l')) +
                   cr * (integrate_linear_ND('r') - tp * integrate_const_ND('r'));
        case 'L':
            return cl * (-integrate_linear_ND('l') + tau * integrate_const_ND('l') + integrate_linear_ND('l') -
                         tm * integrate_const_ND('l'));
        case 'R':
            return cr * (-integrate_linear_ND('r') + tau * integrate_const_ND('r') + integrate_linear_ND('r') -
                         tp * integrate_const_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_t_linear_ND(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_t_ND('l') - tau * integrate_const_t_ND('l'));
        case 'r':
            return -cr * (integrate_linear_t_ND('r') - tau * integrate_const_t_ND('r'));
        case 'c':
            return cl * (integrate_linear_t_ND('l') - tm * integrate_const_t_ND('l')) +
                   cr * (integrate_linear_t_ND('r') - tp * integrate_const_t_ND('r'));
        case 'L':
            return cl * (-integrate_linear_t_ND('l') + tau * integrate_const_t_ND('l') + integrate_linear_t_ND('l') -
                         tm * integrate_const_t_ND('l'));
        case 'R':
            return cr * (-integrate_linear_t_ND('r') + tau * integrate_const_t_ND('r') + integrate_linear_t_ND('r') -
                         tp * integrate_const_t_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_tau_linear_ND(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_tau_ND('l') - tau * integrate_const_tau_ND('l'));
        case 'r':
            return -cr * (integrate_linear_tau_ND('r') - tau * integrate_const_tau_ND('r'));
        case 'c':
            return cl * (integrate_linear_tau_ND('l') - tm * integrate_const_tau_ND('l')) +
                   cr * (integrate_linear_tau_ND('r') - tp * integrate_const_tau_ND('r'));
        case 'L':
            return cl *
                   (-integrate_linear_tau_ND('l') + tau * integrate_const_tau_ND('l') + integrate_linear_tau_ND('l') -
                    tm * integrate_const_tau_ND('l'));
        case 'R':
            return cr *
                   (-integrate_linear_tau_ND('r') + tau * integrate_const_tau_ND('r') + integrate_linear_tau_ND('r') -
                    tp * integrate_const_tau_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_ttau_linear_ND(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_ttau_ND('l') - tau * integrate_const_ttau_ND('l'));
        case 'r':
            return -cr * (integrate_linear_ttau_ND('r') - tau * integrate_const_ttau_ND('r'));
        case 'c':
            return cl * (integrate_linear_ttau_ND('l') - tm * integrate_const_ttau_ND('l')) +
                   cr * (integrate_linear_ttau_ND('r') - tp * integrate_const_ttau_ND('r'));
        case 'L':
            return cl * (-integrate_linear_ttau_ND('l') + tau * integrate_const_ttau_ND('l') +
                         integrate_linear_ttau_ND('l') - tm * integrate_const_ttau_ND('l'));
        case 'R':
            return cr * (-integrate_linear_ttau_ND('r') + tau * integrate_const_ttau_ND('r') +
                         integrate_linear_ttau_ND('r') - tp * integrate_const_ttau_ND('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_square_DN(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 12.0) * (3 * pow(tau, 4) - 4 * pow(tau, 3) * tm + pow(tm, 4)) / mpr;
        case 'r':
            return (1.0 / 3.0) * (pow(tau, 3) * (-tau + tm) + pow(tp, 3) * (tau - tm)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_linear_DN(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 6.0) * (2 * pow(tau, 3) - 3 * pow(tau, 2) * tm + pow(tm, 3)) / mpr;
        case 'r':
            return (1.0 / 2.0) * (pow(tau, 2) * (-tau + tm) + pow(tp, 2) * (tau - tm)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_const_DN(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 2.0) * (pow(tau, 2) - 2 * tau * tm + pow(tm, 2)) / mpr;
        case 'r':
            return (-tau + tp) * (tau - tm) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_square_t_DN(char i) const {
    switch (i) {
        case 'r':
            return ZEROVALUE;
        case 'l':
            return (1.0 / 3.0) * (pow(tau, 3) - pow(tm, 3)) / mpl;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_linear_t_DN(char i) const {
    switch (i) {
        case 'r':
            return ZEROVALUE;
        case 'l':
            return (1.0 / 2.0) * (pow(tau, 2) - pow(tm, 2)) / mpl;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_const_t_DN(char i) const {
    switch (i) {
        case 'r':
            return ZEROVALUE;
        case 'l':
            return (tau - tm) / mpl;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_square_tau_DN(char i) const {
    switch (i) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (1.0 / 3.0) * (-pow(tm, 3) + pow(tp, 3)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_linear_tau_DN(char i) const {
    switch (i) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (1.0 / 2.0) * (-pow(tm, 2) + pow(tp, 2)) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_const_tau_DN(char i) const {
    switch (i) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (-tm + tp) / mpr;
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::integrate_square_ttau_DN(char i) {
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::integrate_linear_ttau_DN(char i) {
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::integrate_const_ttau_DN(char i) {
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::green_function_DN(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (t - tm) / mpr;
    } else {
        return (tau - tm) / mpr;
    }
}

double AGM::GreenfunctionNeumann::green_function_t_DN(double t) const {
    if (t < tau || isclose(tm, t)) {
        return 1.0 / mpr;
    } else {
        return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_function_tau_DN(double t) const {
    if (t < tau || isclose(tm, t)) {
        return ZEROVALUE;
    } else {
        return 1.0 / mpr;
    }
}

double AGM::GreenfunctionNeumann::green_function_ttau_DN(double t) const{
    return ZEROVALUE;
}

double AGM::GreenfunctionNeumann::green_integral_DN(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_linear_DN(pos);
    } else {
        return green_integral_square_DN(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_t_DN(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_t_linear_DN(pos);
    } else {
        return green_integral_t_square_DN(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_tau_DN(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_tau_linear_DN(pos);
    } else {
        return green_integral_tau_square_DN(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_ttau_DN(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_ttau_linear_DN(pos);
    } else {
        return green_integral_ttau_square_DN(pos);
    }
}

double AGM::GreenfunctionNeumann::green_integral_square_DN(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_DN('l') - (tau + tp) * integrate_linear_DN('l') +
                         tau * tp * integrate_const_DN('l') +
                         integrate_square_DN('r') - (tau + tp) * integrate_linear_DN('r') +
                         tau * tp * integrate_const_DN('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_DN('l') - (tau + tm) * integrate_linear_DN('l') +
                         tau * tm * integrate_const_DN('l') +
                         integrate_square_DN('r') - (tau + tm) * integrate_linear_DN('r') +
                         tau * tm * integrate_const_DN('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_DN('l') - (tm + tp) * integrate_linear_DN('l') +
                         tm * tp * integrate_const_DN('l') +
                         integrate_square_DN('r') - (tm + tp) * integrate_linear_DN('r') +
                         tm * tp * integrate_const_DN('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_DN('l') - (tau + tp) * integrate_linear_DN('l') +
                         tau * tp * integrate_const_DN('l') +
                         integrate_square_DN('r') - (tau + tp) * integrate_linear_DN('r') +
                         tau * tp * integrate_const_DN('r')) +
                   cc * (integrate_square_DN('l') - (tm + tp) * integrate_linear_DN('l') +
                         tm * tp * integrate_const_DN('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_DN('l') - (tau + tm) * integrate_linear_DN('l') +
                         tau * tm * integrate_const_DN('l') +
                         integrate_square_DN('r') - (tau + tm) * integrate_linear_DN('r') +
                         tau * tm * integrate_const_DN('r')) +
                   cc * (integrate_square_DN('r') - (tm + tp) * integrate_linear_DN('r') +
                         tm * tp * integrate_const_DN('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_t_square_DN(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc *
                   (integrate_square_t_DN('l') - (tau + tp) * integrate_linear_t_DN('l') +
                    tau * tp * integrate_const_t_DN('l') +
                    integrate_square_t_DN('r') - (tau + tp) * integrate_linear_t_DN('r') +
                    tau * tp * integrate_const_t_DN('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc *
                   (integrate_square_t_DN('l') - (tau + tm) * integrate_linear_t_DN('l') +
                    tau * tm * integrate_const_t_DN('l') +
                    integrate_square_t_DN('r') - (tau + tm) * integrate_linear_t_DN('r') +
                    tau * tm * integrate_const_t_DN('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc *
                   (integrate_square_t_DN('l') - (tm + tp) * integrate_linear_t_DN('l') +
                    tm * tp * integrate_const_t_DN('l') +
                    integrate_square_t_DN('r') - (tm + tp) * integrate_linear_t_DN('r') +
                    tm * tp * integrate_const_t_DN('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl *
                   (integrate_square_t_DN('l') - (tau + tp) * integrate_linear_t_DN('l') +
                    tau * tp * integrate_const_t_DN('l') +
                    integrate_square_t_DN('r') - (tau + tp) * integrate_linear_t_DN('r') +
                    tau * tp * integrate_const_t_DN('r')) +
                   cc *
                   (integrate_square_t_DN('l') - (tm + tp) * integrate_linear_t_DN('l') +
                    tm * tp * integrate_const_t_DN('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr *
                   (integrate_square_t_DN('l') - (tau + tm) * integrate_linear_t_DN('l') +
                    tau * tm * integrate_const_t_DN('l') +
                    integrate_square_t_DN('r') - (tau + tm) * integrate_linear_t_DN('r') +
                    tau * tm * integrate_const_t_DN('r')) +
                   cc *
                   (integrate_square_t_DN('r') - (tm + tp) * integrate_linear_t_DN('r') +
                    tm * tp * integrate_const_t_DN('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_tau_square_DN(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_tau_DN('l') - (tau + tp) * integrate_linear_tau_DN('l') +
                         tau * tp * integrate_const_tau_DN('l') +
                         integrate_square_tau_DN('r') - (tau + tp) * integrate_linear_tau_DN('r') +
                         tau * tp * integrate_const_tau_DN('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_tau_DN('l') - (tau + tm) * integrate_linear_tau_DN('l') +
                         tau * tm * integrate_const_tau_DN('l') +
                         integrate_square_tau_DN('r') - (tau + tm) * integrate_linear_tau_DN('r') +
                         tau * tm * integrate_const_tau_DN('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_tau_DN('l') - (tm + tp) * integrate_linear_tau_DN('l') +
                         tm * tp * integrate_const_tau_DN('l') +
                         integrate_square_tau_DN('r') - (tm + tp) * integrate_linear_tau_DN('r') +
                         tm * tp * integrate_const_tau_DN('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_tau_DN('l') - (tau + tp) * integrate_linear_tau_DN('l') +
                         tau * tp * integrate_const_tau_DN('l') +
                         integrate_square_tau_DN('r') - (tau + tp) * integrate_linear_tau_DN('r') +
                         tau * tp * integrate_const_tau_DN('r')) +
                   cc * (integrate_square_tau_DN('l') - (tm + tp) * integrate_linear_tau_DN('l') +
                         tm * tp * integrate_const_tau_DN('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_tau_DN('l') - (tau + tm) * integrate_linear_tau_DN('l') +
                         tau * tm * integrate_const_tau_DN('l') +
                         integrate_square_tau_DN('r') - (tau + tm) * integrate_linear_tau_DN('r') +
                         tau * tm * integrate_const_tau_DN('r')) +
                   cc * (integrate_square_tau_DN('r') - (tm + tp) * integrate_linear_tau_DN('r') +
                         tm * tp * integrate_const_tau_DN('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_ttau_square_DN(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_ttau_DN('l') - (tau + tp) * integrate_linear_ttau_DN('l') +
                         tau * tp * integrate_const_ttau_DN('l') +
                         integrate_square_ttau_DN('r') - (tau + tp) * integrate_linear_ttau_DN('r') +
                         tau * tp * integrate_const_ttau_DN('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_ttau_DN('l') - (tau + tm) * integrate_linear_ttau_DN('l') +
                         tau * tm * integrate_const_ttau_DN('l') +
                         integrate_square_ttau_DN('r') - (tau + tm) * integrate_linear_ttau_DN('r') +
                         tau * tm * integrate_const_ttau_DN('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_ttau_DN('l') - (tm + tp) * integrate_linear_ttau_DN('l') +
                         tm * tp * integrate_const_ttau_DN('l') +
                         integrate_square_ttau_DN('r') - (tm + tp) * integrate_linear_ttau_DN('r') +
                         tm * tp * integrate_const_ttau_DN('r'));
        case 'L':
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_ttau_DN('l') - (tau + tp) * integrate_linear_ttau_DN('l') +
                         tau * tp * integrate_const_ttau_DN('l') +
                         integrate_square_ttau_DN('r') - (tau + tp) * integrate_linear_ttau_DN('r') +
                         tau * tp * integrate_const_ttau_DN('r')) +
                   cc * (integrate_square_ttau_DN('l') - (tm + tp) * integrate_linear_ttau_DN('l') +
                         tm * tp * integrate_const_ttau_DN('l'));
        case 'R':
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_ttau_DN('l') - (tau + tm) * integrate_linear_ttau_DN('l') +
                         tau * tm * integrate_const_ttau_DN('l') +
                         integrate_square_ttau_DN('r') - (tau + tm) * integrate_linear_ttau_DN('r') +
                         tau * tm * integrate_const_ttau_DN('r')) +
                   cc * (integrate_square_ttau_DN('r') - (tm + tp) * integrate_linear_ttau_DN('r') +
                         tm * tp * integrate_const_ttau_DN('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_linear_DN(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_DN('l') - tau * integrate_const_DN('l'));
        case 'r':
            return -cr * (integrate_linear_DN('r') - tau * integrate_const_DN('r'));
        case 'c':
            return cl * (integrate_linear_DN('l') - tm * integrate_const_DN('l')) +
                   cr * (integrate_linear_DN('r') - tp * integrate_const_DN('r'));
        case 'L':
            return cl * (-integrate_linear_DN('l') + tau * integrate_const_DN('l') + integrate_linear_DN('l') -
                         tm * integrate_const_DN('l'));
        case 'R':
            return cr * (-integrate_linear_DN('r') + tau * integrate_const_DN('r') + integrate_linear_DN('r') -
                         tp * integrate_const_DN('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_t_linear_DN(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_t_DN('l') - tau * integrate_const_t_DN('l'));
        case 'r':
            return -cr * (integrate_linear_t_DN('r') - tau * integrate_const_t_DN('r'));
        case 'c':
            return cl * (integrate_linear_t_DN('l') - tm * integrate_const_t_DN('l')) +
                   cr * (integrate_linear_t_DN('r') - tp * integrate_const_t_DN('r'));
        case 'L':
            return cl * (-integrate_linear_t_DN('l') + tau * integrate_const_t_DN('l') + integrate_linear_t_DN('l') -
                         tm * integrate_const_t_DN('l'));
        case 'R':
            return cr * (-integrate_linear_t_DN('r') + tau * integrate_const_t_DN('r') + integrate_linear_t_DN('r') -
                         tp * integrate_const_t_DN('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_tau_linear_DN(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_tau_DN('l') - tau * integrate_const_tau_DN('l'));
        case 'r':
            return -cr * (integrate_linear_tau_DN('r') - tau * integrate_const_tau_DN('r'));
        case 'c':
            return cl * (integrate_linear_tau_DN('l') - tm * integrate_const_tau_DN('l')) +
                   cr * (integrate_linear_tau_DN('r') - tp * integrate_const_tau_DN('r'));
        case 'L':
            return cl *
                   (-integrate_linear_tau_DN('l') + tau * integrate_const_tau_DN('l') + integrate_linear_tau_DN('l') -
                    tm * integrate_const_tau_DN('l'));
        case 'R':
            return cr *
                   (-integrate_linear_tau_DN('r') + tau * integrate_const_tau_DN('r') + integrate_linear_tau_DN('r') -
                    tp * integrate_const_tau_DN('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionNeumann::green_integral_ttau_linear_DN(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_ttau_DN('l') - tau * integrate_const_ttau_DN('l'));
        case 'r':
            return -cr * (integrate_linear_ttau_DN('r') - tau * integrate_const_ttau_DN('r'));
        case 'c':
            return cl * (integrate_linear_ttau_DN('l') - tm * integrate_const_ttau_DN('l')) +
                   cr * (integrate_linear_ttau_DN('r') - tp * integrate_const_ttau_DN('r'));
        case 'L':
            return cl * (-integrate_linear_ttau_DN('l') + tau * integrate_const_ttau_DN('l') +
                         integrate_linear_ttau_DN('l') - tm * integrate_const_ttau_DN('l'));
        case 'R':
            return cr * (-integrate_linear_ttau_DN('r') + tau * integrate_const_ttau_DN('r') +
                         integrate_linear_ttau_DN('r') - tp * integrate_const_ttau_DN('r'));
        default:
            return ZEROVALUE;
    }
}

AGM::GreenfunctionNeumann::~GreenfunctionNeumann() = default;
