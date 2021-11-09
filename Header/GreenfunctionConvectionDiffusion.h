//
// Created by NIMS-JUNHONG on 2021/02/08.
//

#ifndef AGM_GREENFUNCTIONCONVECTIONDIFFUSION_H
#define AGM_GREENFUNCTIONCONVECTIONDIFFUSION_H

#include "GreenfunctionReactionDiffusionLinear.h"

namespace AGM {
    class GreenfunctionConvectionDiffusion : public Greenfunction {
    protected:
        double a{};
    public:
        GreenfunctionConvectionDiffusion(double tm, double tau, double tp, double mpl, double mpr);

        GreenfunctionConvectionDiffusion(double tm, double tau, double tp, double mpl, double mpr, double a);

        [[nodiscard]] double g(double u, double v) const;

        [[nodiscard]] double integrate_square(char i) const override;

        [[nodiscard]] double integrate_linear(char i) const override;

        [[nodiscard]] double integrate_const(char i) const override;

        [[nodiscard]] double integrate_square_t(char i) const override;

        [[nodiscard]] double integrate_linear_t(char i) const override;

        [[nodiscard]] double integrate_const_t(char i) const override;

        [[nodiscard]] double integrate_square_tau(char i) const override;

        [[nodiscard]] double integrate_linear_tau(char i) const override;

        [[nodiscard]] double integrate_const_tau(char i) const override;

        [[nodiscard]] double integrate_square_ttau(char i) const override;

        [[nodiscard]] double integrate_linear_ttau(char i) const override;

        [[nodiscard]] double integrate_const_ttau(char i) const override;

        [[nodiscard]] double integrate_square_zeroconvection(char i) const;

        [[nodiscard]] double integrate_linear_zeroconvection(char i) const;

        [[nodiscard]] double integrate_const_zeroconvection(char i) const;

        [[nodiscard]] double integrate_square_t_zeroconvection(char i) const;

        [[nodiscard]] double integrate_linear_t_zeroconvection(char i) const;

        [[nodiscard]] double integrate_const_t_zeroconvection(char i) const;

        [[nodiscard]] double integrate_square_tau_zeroconvection(char i) const;

        [[nodiscard]] double integrate_linear_tau_zeroconvection(char i) const;

        [[nodiscard]] double integrate_const_tau_zeroconvection(char i) const;

        [[nodiscard]] double integrate_square_ttau_zeroconvection(char i) const;

        [[nodiscard]] double integrate_linear_ttau_zeroconvection(char i) const;

        [[nodiscard]] double integrate_const_ttau_zeroconvection(char i) const;

        [[nodiscard]] double green_function(double t) const override;

        [[nodiscard]] double green_function_t(double t) const override;

        [[nodiscard]] double green_function_tau(double t) const override;

        [[nodiscard]] double green_function_ttau(double t) const override;

        [[nodiscard]] double green_integral(char pos) const override;

        [[nodiscard]] double green_integral_t(char pos) const override;

        [[nodiscard]] double green_integral_tau(char pos) const override;

        [[nodiscard]] double green_integral_ttau(char pos) const override;
    };

}

#endif //AGM_GREENFUNCTIONCONVECTIONDIFFUSION_H
