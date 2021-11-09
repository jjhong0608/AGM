//
// Created by NIMS-JUNHONG on 2020/12/21.
//

#include "Greenfunction.h"

AGM::Greenfunction::Greenfunction(double tm, double tau, double tp, double mpl, double mpr) : tm(tm), tau(tau), tp(tp),
                                                                                              mpl(mpl), mpr(mpr) {}

AGM::Greenfunction::~Greenfunction() = default;

double AGM::Greenfunction::green_function(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (t - tm) * (tau - tp) / (mpl * (tau - tp) + mpr * (-tau + tm));
    } else {
        return (t - tp) * (tau - tm) / (mpl * (tau - tp) + mpr * (-tau + tm));
    }
}

double AGM::Greenfunction::green_function_t(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (tau - tp) / (mpl * (tau - tp) + mpr * (-tau + tm));
    } else {
        return (tau - tm) / (mpl * (tau - tp) + mpr * (-tau + tm));
    }
}

double AGM::Greenfunction::green_function_tau(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (t - tm) * (mpl * (tau - tp) + mpr * (-tau + tm) - (mpl - mpr) * (tau - tp)) /
               pow(mpl * (tau - tp) + mpr * (-tau + tm), 2);
    } else {
        return (t - tp) * (mpl * (tau - tp) + mpr * (-tau + tm) - (mpl - mpr) * (tau - tm)) /
               pow(mpl * (tau - tp) + mpr * (-tau + tm), 2);
    }
}

double AGM::Greenfunction::green_function_ttau(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (mpl * (tau - tp) + mpr * (-tau + tm) - (mpl - mpr) * (tau - tp)) /
               pow(mpl * (tau - tp) + mpr * (-tau + tm), 2);
    } else {
        return (mpl * (tau - tp) + mpr * (-tau + tm) - (mpl - mpr) * (tau - tm)) /
               pow(mpl * (tau - tp) + mpr * (-tau + tm), 2);
    }
}

double AGM::Greenfunction::integrate_square(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 12.0 * pow(tau - tm, 2) * (tau - tp) * (3 * pow(tau, 2) + 2 * tau * tm + pow(tm, 2)) /
                   (mpl * (tau - tp) - mpr * (tau - tm));
        case 'r':
            return -1.0 / 12.0 * (tau - tm) * pow(tau - tp, 2) * (3 * pow(tau, 2) + 2 * tau * tp + pow(tp, 2)) /
                   (mpl * (tau - tp) - mpr * (tau - tm));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_linear(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 6.0 * pow(tau - tm, 2) * (tau - tp) * (2 * tau + tm) / (mpl * (tau - tp) - mpr * (tau - tm));
        case 'r':
            return -1.0 / 6.0 * (tau - tm) * pow(tau - tp, 2) * (2 * tau + tp) / (mpl * (tau - tp) - mpr * (tau - tm));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_const(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 2.0 * pow(tau - tm, 2) * (tau - tp) / (mpl * (tau - tp) - mpr * (tau - tm));
        case 'r':
            return -1.0 / 2.0 * (tau - tm) * pow(tau - tp, 2) / (mpl * (tau - tp) - mpr * (tau - tm));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_square_t(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 3.0 * (tau - tm) * (tau - tp) * (pow(tau, 2) + tau * tm + pow(tm, 2)) /
                   (mpl * (tau - tp) - mpr * (tau - tm));
        case 'r':
            return -1.0 / 3.0 * (tau - tm) * (tau - tp) * (pow(tau, 2) + tau * tp + pow(tp, 2)) /
                   (mpl * (tau - tp) - mpr * (tau - tm));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_linear_t(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 2.0 * (tau - tm) * (tau + tm) * (tau - tp) / (mpl * (tau - tp) - mpr * (tau - tm));
        case 'r':
            return -1.0 / 2.0 * (tau - tm) * (tau - tp) * (tau + tp) / (mpl * (tau - tp) - mpr * (tau - tm));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_const_t(char i) const {
    switch (i) {
        case 'l':
            return (tau - tm) * (tau - tp) / (mpl * (tau - tp) - mpr * (tau - tm));
        case 'r':
            return -(tau - tm) * (tau - tp) / (mpl * (tau - tp) - mpr * (tau - tm));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_square_tau(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 12.0 * mpr * pow(tau - tm, 2) * (tm - tp) * (3 * pow(tau, 2) + 2 * tau * tm + pow(tm, 2)) /
                   pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        case 'r':
            return -1.0 / 12.0 * mpl * pow(tau - tp, 2) * (tm - tp) * (3 * pow(tau, 2) + 2 * tau * tp + pow(tp, 2)) /
                   pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_linear_tau(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 6.0 * mpr * pow(tau - tm, 2) * (2 * tau + tm) * (tm - tp) /
                   pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        case 'r':
            return -1.0 / 6.0 * mpl * pow(tau - tp, 2) * (2 * tau + tp) * (tm - tp) /
                   pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_const_tau(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 2.0 * mpr * pow(tau - tm, 2) * (tm - tp) / pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        case 'r':
            return -1.0 / 2.0 * mpl * pow(tau - tp, 2) * (tm - tp) / pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_square_ttau(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 3.0 * mpr * (tau - tm) * (tm - tp) * (pow(tau, 2) + tau * tm + pow(tm, 2)) /
                   pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        case 'r':
            return -1.0 / 3.0 * mpl * (tau - tp) * (tm - tp) * (pow(tau, 2) + tau * tp + pow(tp, 2)) /
                   pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_linear_ttau(char i) const {
    switch (i) {
        case 'l':
            return 1.0 / 2.0 * mpr * (tau - tm) * (tau + tm) * (tm - tp) / pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        case 'r':
            return -1.0 / 2.0 * mpl * (tau - tp) * (tau + tp) * (tm - tp) / pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::integrate_const_ttau(char i) const {
    switch (i) {
        case 'l':
            return mpr * (tau - tm) * (tm - tp) / pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        case 'r':
            return -mpl * (tau - tp) * (tm - tp) / pow(mpl * (tau - tp) - mpr * (tau - tm), 2);
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_linear(pos);
    } else {
        return green_integral_square(pos);
    }
}

double AGM::Greenfunction::green_integral_t(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_t_linear(pos);
    } else {
        return green_integral_t_square(pos);
    }
}

double AGM::Greenfunction::green_integral_tau(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_tau_linear(pos);
    } else {
        return green_integral_tau_square(pos);
    }
}

double AGM::Greenfunction::green_integral_ttau(char pos) const {
    if (min(tau - tm, tp - tau) / max(tau - tm, tp - tau) < 0.3) {
        return green_integral_ttau_linear(pos);
    } else {
        return green_integral_ttau_square(pos);
    }
}

double AGM::Greenfunction::green_integral_square(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square('l') - (tau + tp) * integrate_linear('l') + tau * tp * integrate_const('l') +
                         integrate_square('r') - (tau + tp) * integrate_linear('r') + tau * tp * integrate_const('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square('l') - (tau + tm) * integrate_linear('l') + tau * tm * integrate_const('l') +
                         integrate_square('r') - (tau + tm) * integrate_linear('r') + tau * tm * integrate_const('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square('l') - (tm + tp) * integrate_linear('l') + tm * tp * integrate_const('l') +
                         integrate_square('r') - (tm + tp) * integrate_linear('r') + tm * tp * integrate_const('r'));
        case 'L':
            return green_integral_linear(pos);
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square('l') - (tau + tp) * integrate_linear('l') + tau * tp * integrate_const('l') +
                         integrate_square('r') - (tau + tp) * integrate_linear('r') + tau * tp * integrate_const('r')) +
                   cc * (integrate_square('l') - (tm + tp) * integrate_linear('l') + tm * tp * integrate_const('l'));
        case 'R':
            return green_integral_linear(pos);
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square('l') - (tau + tm) * integrate_linear('l') + tau * tm * integrate_const('l') +
                         integrate_square('r') - (tau + tm) * integrate_linear('r') + tau * tm * integrate_const('r')) +
                   cc * (integrate_square('r') - (tm + tp) * integrate_linear('r') + tm * tp * integrate_const('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral_t_square(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc *
                   (integrate_square_t('l') - (tau + tp) * integrate_linear_t('l') + tau * tp * integrate_const_t('l') +
                    integrate_square_t('r') - (tau + tp) * integrate_linear_t('r') + tau * tp * integrate_const_t('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc *
                   (integrate_square_t('l') - (tau + tm) * integrate_linear_t('l') + tau * tm * integrate_const_t('l') +
                    integrate_square_t('r') - (tau + tm) * integrate_linear_t('r') + tau * tm * integrate_const_t('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc *
                   (integrate_square_t('l') - (tm + tp) * integrate_linear_t('l') + tm * tp * integrate_const_t('l') +
                    integrate_square_t('r') - (tm + tp) * integrate_linear_t('r') + tm * tp * integrate_const_t('r'));
        case 'L':
            return green_integral_t_linear(pos);
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl *
                   (integrate_square_t('l') - (tau + tp) * integrate_linear_t('l') + tau * tp * integrate_const_t('l') +
                    integrate_square_t('r') - (tau + tp) * integrate_linear_t('r') +
                    tau * tp * integrate_const_t('r')) +
                   cc *
                   (integrate_square_t('l') - (tm + tp) * integrate_linear_t('l') + tm * tp * integrate_const_t('l'));
        case 'R':
            return green_integral_t_linear(pos);
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr *
                   (integrate_square_t('l') - (tau + tm) * integrate_linear_t('l') + tau * tm * integrate_const_t('l') +
                    integrate_square_t('r') - (tau + tm) * integrate_linear_t('r') +
                    tau * tm * integrate_const_t('r')) +
                   cc *
                   (integrate_square_t('r') - (tm + tp) * integrate_linear_t('r') + tm * tp * integrate_const_t('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral_tau_square(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_tau('l') - (tau + tp) * integrate_linear_tau('l') +
                         tau * tp * integrate_const_tau('l') +
                         integrate_square_tau('r') - (tau + tp) * integrate_linear_tau('r') +
                         tau * tp * integrate_const_tau('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_tau('l') - (tau + tm) * integrate_linear_tau('l') +
                         tau * tm * integrate_const_tau('l') +
                         integrate_square_tau('r') - (tau + tm) * integrate_linear_tau('r') +
                         tau * tm * integrate_const_tau('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_tau('l') - (tm + tp) * integrate_linear_tau('l') +
                         tm * tp * integrate_const_tau('l') +
                         integrate_square_tau('r') - (tm + tp) * integrate_linear_tau('r') +
                         tm * tp * integrate_const_tau('r'));
        case 'L':
            return green_integral_tau_linear(pos);
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_tau('l') - (tau + tp) * integrate_linear_tau('l') +
                         tau * tp * integrate_const_tau('l') +
                         integrate_square_tau('r') - (tau + tp) * integrate_linear_tau('r') +
                         tau * tp * integrate_const_tau('r')) +
                   cc * (integrate_square_tau('l') - (tm + tp) * integrate_linear_tau('l') +
                         tm * tp * integrate_const_tau('l'));
        case 'R':
            return green_integral_tau_linear(pos);
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_tau('l') - (tau + tm) * integrate_linear_tau('l') +
                         tau * tm * integrate_const_tau('l') +
                         integrate_square_tau('r') - (tau + tm) * integrate_linear_tau('r') +
                         tau * tm * integrate_const_tau('r')) +
                   cc * (integrate_square_tau('r') - (tm + tp) * integrate_linear_tau('r') +
                         tm * tp * integrate_const_tau('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral_ttau_square(char pos) const {
    double cc{}, cl{}, cr{};
    switch (pos) {
        case 'l':
            cc = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            return cc * (integrate_square_ttau('l') - (tau + tp) * integrate_linear_ttau('l') +
                         tau * tp * integrate_const_ttau('l') +
                         integrate_square_ttau('r') - (tau + tp) * integrate_linear_ttau('r') +
                         tau * tp * integrate_const_ttau('r'));
        case 'r':
            cc = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            return cc * (integrate_square_ttau('l') - (tau + tm) * integrate_linear_ttau('l') +
                         tau * tm * integrate_const_ttau('l') +
                         integrate_square_ttau('r') - (tau + tm) * integrate_linear_ttau('r') +
                         tau * tm * integrate_const_ttau('r'));
        case 'c':
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cc * (integrate_square_ttau('l') - (tm + tp) * integrate_linear_ttau('l') +
                         tm * tp * integrate_const_ttau('l') +
                         integrate_square_ttau('r') - (tm + tp) * integrate_linear_ttau('r') +
                         tm * tp * integrate_const_ttau('r'));
        case 'L':
            return green_integral_ttau_linear(pos);
            cl = -UNITVALUE / (tau * tm - tau * tp - tm * tm + tm * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cl * (integrate_square_ttau('l') - (tau + tp) * integrate_linear_ttau('l') +
                         tau * tp * integrate_const_ttau('l') +
                         integrate_square_ttau('r') - (tau + tp) * integrate_linear_ttau('r') +
                         tau * tp * integrate_const_ttau('r')) +
                   cc * (integrate_square_ttau('l') - (tm + tp) * integrate_linear_ttau('l') +
                         tm * tp * integrate_const_ttau('l'));
        case 'R':
            return green_integral_ttau_linear(pos);
            cr = UNITVALUE / (tau * tm - tau * tp - tm * tp + tp * tp);
            cc = UNITVALUE / (tau * tau - tau * tm - tau * tp + tm * tp);
            return cr * (integrate_square_ttau('l') - (tau + tm) * integrate_linear_ttau('l') +
                         tau * tm * integrate_const_ttau('l') +
                         integrate_square_ttau('r') - (tau + tm) * integrate_linear_ttau('r') +
                         tau * tm * integrate_const_ttau('r')) +
                   cc * (integrate_square_ttau('r') - (tm + tp) * integrate_linear_ttau('r') +
                         tm * tp * integrate_const_ttau('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral_linear(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear('l') - tau * integrate_const('l'));
        case 'r':
            return -cr * (integrate_linear('r') - tau * integrate_const('r'));
        case 'c':
            return cl * (integrate_linear('l') - tm * integrate_const('l')) +
                   cr * (integrate_linear('r') - tp * integrate_const('r'));
        case 'L':
            return -((integrate_linear('l') + integrate_linear('r')) -
                     tp * (integrate_const('l') + integrate_const('r'))) / (tp - tm);
            return cl * (-integrate_linear('l') + tau * integrate_const('l') + integrate_linear('l') -
                         tm * integrate_const('l'));
        case 'R':
            return -((integrate_linear('l') + integrate_linear('r')) -
                     tm * (integrate_const('l') + integrate_const('r'))) / (tm - tp);
            return cr * (-integrate_linear('r') + tau * integrate_const('r') + integrate_linear('r') -
                         tp * integrate_const('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral_t_linear(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_t('l') - tau * integrate_const_t('l'));
        case 'r':
            return -cr * (integrate_linear_t('r') - tau * integrate_const_t('r'));
        case 'c':
            return cl * (integrate_linear_t('l') - tm * integrate_const_t('l')) +
                   cr * (integrate_linear_t('r') - tp * integrate_const_t('r'));
        case 'L':
            return -((integrate_linear_t('l') + integrate_linear_t('r')) -
                     tp * (integrate_const_t('l') + integrate_const_t('r'))) / (tp - tm);
            return cl * (-integrate_linear_t('l') + tau * integrate_const_t('l') + integrate_linear_t('l') -
                         tm * integrate_const_t('l'));
        case 'R':
            return -((integrate_linear_t('l') + integrate_linear_t('r')) -
                     tm * (integrate_const_t('l') + integrate_const_t('r'))) / (tm - tp);
            return cr * (-integrate_linear_t('r') + tau * integrate_const_t('r') + integrate_linear_t('r') -
                         tp * integrate_const_t('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral_tau_linear(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_tau('l') - tau * integrate_const_tau('l'));
        case 'r':
            return -cr * (integrate_linear_tau('r') - tau * integrate_const_tau('r'));
        case 'c':
            return cl * (integrate_linear_tau('l') - tm * integrate_const_tau('l')) +
                   cr * (integrate_linear_tau('r') - tp * integrate_const_tau('r'));
        case 'L':
            return -((integrate_linear_tau('l') + integrate_linear_tau('r')) -
                     tp * (integrate_const_tau('l') + integrate_const_tau('r'))) / (tp - tm);
            return cl * (-integrate_linear_tau('l') + tau * integrate_const_tau('l') + integrate_linear_tau('l') -
                         tm * integrate_const_tau('l'));
        case 'R':
            return -((integrate_linear_tau('l') + integrate_linear_tau('r')) -
                     tm * (integrate_const_tau('l') + integrate_const_tau('r'))) / (tm - tp);
            return cr * (-integrate_linear_tau('r') + tau * integrate_const_tau('r') + integrate_linear_tau('r') -
                         tp * integrate_const_tau('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::Greenfunction::green_integral_ttau_linear(char pos) const {
    double cl = isclose(tau, tm) ? ZEROVALUE : UNITVALUE / (tau - tm);
    double cr = isclose(tau, tp) ? ZEROVALUE : UNITVALUE / (tau - tp);

    switch (pos) {
        case 'l':
            return -cl * (integrate_linear_ttau('l') - tau * integrate_const_ttau('l'));
        case 'r':
            return -cr * (integrate_linear_ttau('r') - tau * integrate_const_ttau('r'));
        case 'c':
            return cl * (integrate_linear_ttau('l') - tm * integrate_const_ttau('l')) +
                   cr * (integrate_linear_ttau('r') - tp * integrate_const_ttau('r'));
        case 'L':
            return -((integrate_linear_ttau('l') + integrate_linear_ttau('r')) -
                     tp * (integrate_const_ttau('l') + integrate_const_ttau('r'))) / (tp - tm);
            return cl * (-integrate_linear_ttau('l') + tau * integrate_const_ttau('l') + integrate_linear_ttau('l') -
                         tm * integrate_const_ttau('l'));
        case 'R':
            return -((integrate_linear_ttau('l') + integrate_linear_ttau('r')) -
                     tm * (integrate_const_ttau('l') + integrate_const_ttau('r'))) / (tm - tp);
            return cr * (-integrate_linear_ttau('r') + tau * integrate_const_ttau('r') + integrate_linear_ttau('r') -
                         tp * integrate_const_ttau('r'));
        default:
            return ZEROVALUE;
    }
}


