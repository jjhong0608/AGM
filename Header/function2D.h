//
// Created by 조준홍 on 2021/03/05.
//

#ifndef AGM_FUNCTION2D_H
#define AGM_FUNCTION2D_H

#include "function.h"

namespace AGM {
    class function2D : public function {
    public:
        [[nodiscard]] double u() const override;

        [[nodiscard]] double v() const;

        [[nodiscard]] double p() const;

        [[nodiscard]] double ut() const;

        [[nodiscard]] double vt() const;

        [[nodiscard]] double ux() const;

        [[nodiscard]] double uy() const;

        [[nodiscard]] double vx() const;

        [[nodiscard]] double vy() const;

        [[nodiscard]] double px() const;

        [[nodiscard]] double py() const;

        [[nodiscard]] double f1() const;

        [[nodiscard]] double f2() const;

        [[nodiscard]] double uxx() const;

        [[nodiscard]] double vyy() const;

        [[nodiscard]] double phi() const override;

        [[nodiscard]] double psi() const;
    };

}


#endif //AGM_FUNCTION2D_H
