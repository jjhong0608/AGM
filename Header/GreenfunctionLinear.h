//
// Created by NIMS-JUNHONG on 2021/01/22.
//

#ifndef AGM_GREENFUNCTIONLINEAR_H
#define AGM_GREENFUNCTIONLINEAR_H

#include "Greenfunction.h"

namespace AGM {
    class GreenfunctionLinear : public Greenfunction {
    public:
        GreenfunctionLinear(double tm, double tau, double tp, double mpl, double mpr);

        [[nodiscard]] double green_integral(char pos) const override;

        [[nodiscard]] double green_integral_t(char pos) const override;

        [[nodiscard]] double green_integral_tau(char pos) const override;

        [[nodiscard]] double green_integral_ttau(char pos) const override;
    };

}

#endif //AGM_GREENFUNCTIONLINEAR_H
