//
// Created by NIMS-JUNHONG on 2021/01/22.
//

#include "GreenfunctionLinear.h"

AGM::GreenfunctionLinear::GreenfunctionLinear(double tm, double tau, double tp, double mpl, double mpr) : Greenfunction(
        tm, tau, tp, mpl, mpr) {}

double AGM::GreenfunctionLinear::green_integral(char pos) const {
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
            return cl * (-integrate_linear('l') + tau * integrate_const('l') + integrate_linear('l') -
                         tm * integrate_const('l'));
        case 'R':
            return cr * (-integrate_linear('r') + tau * integrate_const('r') + integrate_linear('r') -
                         tp * integrate_const('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionLinear::green_integral_t(char pos) const {
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
            return cl * (-integrate_linear_t('l') + tau * integrate_const_t('l') + integrate_linear_t('l') -
                         tm * integrate_const_t('l'));
        case 'R':
            return cr * (-integrate_linear_t('r') + tau * integrate_const_t('r') + integrate_linear_t('r') -
                         tp * integrate_const_t('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionLinear::green_integral_tau(char pos) const {
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
            return cl * (-integrate_linear_tau('l') + tau * integrate_const_tau('l') + integrate_linear_tau('l') -
                         tm * integrate_const_tau('l'));
        case 'R':
            return cr * (-integrate_linear_tau('r') + tau * integrate_const_tau('r') + integrate_linear_tau('r') -
                         tp * integrate_const_tau('r'));
        default:
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionLinear::green_integral_ttau(char pos) const {
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
            return cl * (-integrate_linear_ttau('l') + tau * integrate_const_ttau('l') + integrate_linear_ttau('l') -
                         tm * integrate_const_ttau('l'));
        case 'R':
            return cr * (-integrate_linear_ttau('r') + tau * integrate_const_ttau('r') + integrate_linear_ttau('r') -
                         tp * integrate_const_ttau('r'));
        default:
            return ZEROVALUE;
    }
}
