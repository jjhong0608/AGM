//
// Created by NIMS-JUNHONG on 2020/12/21.
//

#ifndef AGM_GREENFUNCTION_H
#define AGM_GREENFUNCTION_H

#include "util.h"

namespace AGM {
    class Greenfunction {
    protected:
        double tm{}, tau{}, tp{}, mpl{}, mpr{};

    public:
        Greenfunction(double tm, double tau, double tp, double mpl, double mpr);

        virtual ~Greenfunction();

        [[nodiscard]] virtual double integrate_square(char i) const;
        [[nodiscard]] virtual double integrate_linear(char i) const;
        [[nodiscard]] virtual double integrate_const(char i) const;

        [[nodiscard]] virtual double integrate_square_t(char i) const;
        [[nodiscard]] virtual double integrate_linear_t(char i) const;
        [[nodiscard]] virtual double integrate_const_t(char i) const;

        [[nodiscard]] virtual double integrate_square_tau(char i) const;
        [[nodiscard]] virtual double integrate_linear_tau(char i) const;
        [[nodiscard]] virtual double integrate_const_tau(char i) const;

        [[nodiscard]] virtual double integrate_square_ttau(char i) const;
        [[nodiscard]] virtual double integrate_linear_ttau(char i) const;
        [[nodiscard]] virtual double integrate_const_ttau(char i) const;

        [[nodiscard]] virtual double green_function(double t) const;
        [[nodiscard]] virtual double green_function_t(double t) const;
        [[nodiscard]] virtual double green_function_tau(double t) const;
        [[nodiscard]] virtual double green_function_ttau(double t) const;

        [[nodiscard]] virtual double green_integral(char pos) const;
        [[nodiscard]] virtual double green_integral_t(char pos) const;
        [[nodiscard]] virtual double green_integral_tau(char pos) const;
        [[nodiscard]] virtual double green_integral_ttau(char pos) const;

        [[nodiscard]] double green_integral_square(char pos) const;
        [[nodiscard]] double green_integral_t_square(char pos) const;
        [[nodiscard]] double green_integral_tau_square(char pos) const;
        [[nodiscard]] double green_integral_ttau_square(char pos) const;

        [[nodiscard]] double green_integral_linear(char pos) const;
        [[nodiscard]] double green_integral_t_linear(char pos) const;
        [[nodiscard]] double green_integral_tau_linear(char pos) const;
        [[nodiscard]] double green_integral_ttau_linear(char pos) const;
    };
}

#endif //AGM_GREENFUNCTION_H
