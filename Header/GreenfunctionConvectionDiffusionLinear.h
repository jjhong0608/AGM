//
// Created by JUNHONG-NIMS on 2021/02/10.
//

#ifndef AGM_GREENFUNCTIONCONVECTIONDIFFUSIONLINEAR_H
#define AGM_GREENFUNCTIONCONVECTIONDIFFUSIONLINEAR_H

#include "GreenfunctionConvectionDiffusion.h"

namespace AGM {
    class GreenfunctionConvectionDiffusionLinear : public GreenfunctionConvectionDiffusion {
    public:
        GreenfunctionConvectionDiffusionLinear(double tm, double tau, double tp, double mpl, double mpr);

        GreenfunctionConvectionDiffusionLinear(double tm, double tau, double tp, double mpl, double mpr, double a);

        [[nodiscard]] double green_integral(char pos) const override;

        [[nodiscard]] double green_integral_t(char pos) const override;

        [[nodiscard]] double green_integral_tau(char pos) const override;

        [[nodiscard]] double green_integral_ttau(char pos) const override;
    };

}


#endif //AGM_GREENFUNCTIONCONVECTIONDIFFUSIONLINEAR_H
