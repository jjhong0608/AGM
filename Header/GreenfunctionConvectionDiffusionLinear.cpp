//
// Created by JUNHONG-NIMS on 2021/02/10.
//

#include "GreenfunctionConvectionDiffusionLinear.h"

AGM::GreenfunctionConvectionDiffusionLinear::GreenfunctionConvectionDiffusionLinear(double tm, double tau, double tp,
                                                                                    double mpl, double mpr)
        : GreenfunctionConvectionDiffusion(tm, tau, tp, mpl, mpr) {

}

AGM::GreenfunctionConvectionDiffusionLinear::GreenfunctionConvectionDiffusionLinear(double tm, double tau, double tp,
                                                                                    double mpl, double mpr, double a)
        : GreenfunctionConvectionDiffusion(tm, tau, tp, mpl, mpr, a) {

}

double AGM::GreenfunctionConvectionDiffusionLinear::green_integral(char pos) const {
    double rv{}, cc;
    if (fabs(a / mpl) < 1e-7) {
        return GreenfunctionLinear(tm, tau, tp, mpl, mpr).green_integral(pos);
    } else {
        switch (pos) {
            case 'l':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tm - tau);
                    return cc * integrate_linear('l') - tau * cc * integrate_const('l');
                } else {
                    return ZEROVALUE;
                }
            case 'r':
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tp - tau);
                    return cc * integrate_linear('r') - tau * cc * integrate_const('r');
                } else {
                    return ZEROVALUE;
                }
            case 'c':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tau - tm);
                    rv += cc * integrate_linear('l') - tm * cc * integrate_const('l');
                }
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tau - tp);
                    rv += cc * integrate_linear('r') - tp * cc * integrate_const('r');
                }
                return rv;
            default:
                return ZEROVALUE;
        }
    }
}

double AGM::GreenfunctionConvectionDiffusionLinear::green_integral_t(char pos) const {
    double rv{}, cc;
    if (fabs(a / mpl) < 1e-7) {
        return GreenfunctionLinear(tm, tau, tp, mpl, mpr).green_integral_t(pos);
    } else {
        switch (pos) {
            case 'l':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tm - tau);
                    return cc * integrate_linear_t('l') - tau * cc * integrate_const_t('l');
                } else {
                    return ZEROVALUE;
                }
            case 'r':
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tp - tau);
                    return cc * integrate_linear_t('r') - tau * cc * integrate_const_t('r');
                } else {
                    return ZEROVALUE;
                }
            case 'c':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tau - tm);
                    rv += cc * integrate_linear_t('l') - tm * cc * integrate_const_t('l');
                }
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tau - tp);
                    rv += cc * integrate_linear_t('r') - tp * cc * integrate_const_t('r');
                }
                return rv;
            default:
                return ZEROVALUE;
        }
    }
}

double AGM::GreenfunctionConvectionDiffusionLinear::green_integral_tau(char pos) const {
    double rv{}, cc;
    if (fabs(a / mpl) < 1e-7) {
        return GreenfunctionLinear(tm, tau, tp, mpl, mpr).green_integral_tau(pos);
    } else {
        switch (pos) {
            case 'l':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tm - tau);
                    return cc * integrate_linear_tau('l') - tau * cc * integrate_const_tau('l');
                } else {
                    return ZEROVALUE;
                }
            case 'r':
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tp - tau);
                    return cc * integrate_linear_tau('r') - tau * cc * integrate_const_tau('r');
                } else {
                    return ZEROVALUE;
                }
            case 'c':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tau - tm);
                    rv += cc * integrate_linear_tau('l') - tm * cc * integrate_const_tau('l');
                }
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tau - tp);
                    rv += cc * integrate_linear_tau('r') - tp * cc * integrate_const_tau('r');
                }
                return rv;
            default:
                return ZEROVALUE;
        }
    }
}

double AGM::GreenfunctionConvectionDiffusionLinear::green_integral_ttau(char pos) const {
    double rv{}, cc;
    if (fabs(a / mpl) < 1e-7) {
        return GreenfunctionLinear(tm, tau, tp, mpl, mpr).green_integral_ttau(pos);
    } else {
        switch (pos) {
            case 'l':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tm - tau);
                    return cc * integrate_linear_ttau('l') - tau * cc * integrate_const_ttau('l');
                } else {
                    return ZEROVALUE;
                }
            case 'r':
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tp - tau);
                    return cc * integrate_linear_ttau('r') - tau * cc * integrate_const_ttau('r');
                } else {
                    return ZEROVALUE;
                }
            case 'c':
                if (!isclose(tm, tau)) {
                    cc = UNITVALUE / (tau - tm);
                    rv += cc * integrate_linear_ttau('l') - tm * cc * integrate_const_ttau('l');
                }
                if (!isclose(tp, tau)) {
                    cc = UNITVALUE / (tau - tp);
                    rv += cc * integrate_linear_ttau('r') - tp * cc * integrate_const_ttau('r');
                }
                return rv;
            default:
                return ZEROVALUE;
        }
    }
}

