//
// Created by NIMS-JUNHONG on 2021/06/16.
//

#ifndef AGM_GREENFUNCTIONNEUMANN_H
#define AGM_GREENFUNCTIONNEUMANN_H

#include "GreenfunctionLinear.h"


namespace AGM {
    class GreenfunctionNeumann : public Greenfunction {
    public:
        GreenfunctionNeumann(double tm, double tau, double tp, double mpl, double mpr);

        ~GreenfunctionNeumann() override;

        [[nodiscard]] virtual double integrate_square_ND(char i) const;
        [[nodiscard]] virtual double integrate_linear_ND(char i) const;
        [[nodiscard]] virtual double integrate_const_ND(char i) const;

        [[nodiscard]] virtual double integrate_square_t_ND(char i) const;
        [[nodiscard]] virtual double integrate_linear_t_ND(char i) const;
        [[nodiscard]] virtual double integrate_const_t_ND(char i) const;

        [[nodiscard]] virtual double integrate_square_tau_ND(char i) const;
        [[nodiscard]] virtual double integrate_linear_tau_ND(char i) const;
        [[nodiscard]] virtual double integrate_const_tau_ND(char i) const;

        [[nodiscard]] static double integrate_square_ttau_ND(char i) ;
        [[nodiscard]] static double integrate_linear_ttau_ND(char i) ;
        [[nodiscard]] static double integrate_const_ttau_ND(char i) ;

        [[nodiscard]] virtual double green_function_ND(double t) const;
        [[nodiscard]] virtual double green_function_t_ND(double t) const;
        [[nodiscard]] virtual double green_function_tau_ND(double t) const;
        [[nodiscard]] virtual double green_function_ttau_ND(double t) const;

        [[nodiscard]] virtual double green_integral_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_t_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_tau_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_ttau_ND(char pos) const;

        [[nodiscard]] virtual double green_integral_square_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_t_square_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_tau_square_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_ttau_square_ND(char pos) const;

        [[nodiscard]] virtual double green_integral_linear_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_t_linear_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_tau_linear_ND(char pos) const;
        [[nodiscard]] virtual double green_integral_ttau_linear_ND(char pos) const;

        [[nodiscard]] virtual double integrate_square_DN(char i) const;
        [[nodiscard]] virtual double integrate_linear_DN(char i) const;
        [[nodiscard]] virtual double integrate_const_DN(char i) const;

        [[nodiscard]] virtual double integrate_square_t_DN(char i) const;
        [[nodiscard]] virtual double integrate_linear_t_DN(char i) const;
        [[nodiscard]] virtual double integrate_const_t_DN(char i) const;

        [[nodiscard]] virtual double integrate_square_tau_DN(char i) const;
        [[nodiscard]] virtual double integrate_linear_tau_DN(char i) const;
        [[nodiscard]] virtual double integrate_const_tau_DN(char i) const;

        [[nodiscard]] static double integrate_square_ttau_DN(char i) ;
        [[nodiscard]] static double integrate_linear_ttau_DN(char i) ;
        [[nodiscard]] static double integrate_const_ttau_DN(char i) ;

        [[nodiscard]] virtual double green_function_DN(double t) const;
        [[nodiscard]] virtual double green_function_t_DN(double t) const;
        [[nodiscard]] virtual double green_function_tau_DN(double t) const;
        [[nodiscard]] virtual double green_function_ttau_DN(double t) const;

        [[nodiscard]] virtual double green_integral_DN(char pos) const;
        [[nodiscard]] virtual double green_integral_t_DN(char pos) const;
        [[nodiscard]] virtual double green_integral_tau_DN(char pos) const;
        [[nodiscard]] double green_integral_ttau_DN(char pos) const;

        [[nodiscard]] virtual double green_integral_square_DN(char pos) const;
        [[nodiscard]] virtual double green_integral_t_square_DN(char pos) const;
        [[nodiscard]] virtual double green_integral_tau_square_DN(char pos) const;
        [[nodiscard]] double green_integral_ttau_square_DN(char pos) const;

        [[nodiscard]] virtual double green_integral_linear_DN(char pos) const;
        [[nodiscard]] virtual double green_integral_t_linear_DN(char pos) const;
        [[nodiscard]] virtual double green_integral_tau_linear_DN(char pos) const;
        [[nodiscard]] double green_integral_ttau_linear_DN(char pos) const;
    };
}


#endif //AGM_GREENFUNCTIONNEUMANN_H
