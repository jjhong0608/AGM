//
// Created by NIMS-JUNHONG on 2021/07/12.
//

#ifndef AGM_GREENFUNCTIONREACTIONDIFFUSIONLINEAR_H
#define AGM_GREENFUNCTIONREACTIONDIFFUSIONLINEAR_H

#include "GreenfunctionReactionDiffusion.h"

namespace AGM {
    class GreenfunctionReactionDiffusionLinear : public GreenfunctionReactionDiffusion {

    public:
        GreenfunctionReactionDiffusionLinear(double tm, double tau, double tp, double mpl, double mpr, double c);

        [[nodiscard]] double green_integral(char pos) const override;

        [[nodiscard]] double green_integral_t(char pos) const override;

        [[nodiscard]] double green_integral_tau(char pos) const override;

        [[nodiscard]] double green_integral_ttau(char pos) const override;
    };
}

#endif //AGM_GREENFUNCTIONREACTIONDIFFUSIONLINEAR_H
