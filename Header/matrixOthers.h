//
// Created by NIMS-JUNHONG on 2021/02/01.
//

#ifndef AGM_MATRIXOTHERS_H
#define AGM_MATRIXOTHERS_H

#include "matrix.h"

namespace AGM {
    template<typename pt>
    class matrixTwoDimension : public matrix<pt> {
    private:
        std::vector<pt> *pts0{};
    public:
        matrixTwoDimension();

        matrixTwoDimension(std::vector<pt> *pts, std::vector<pt> *pts0);

        ~matrixTwoDimension() override;

        [[nodiscard]] std::vector<pt> *getPts0() const;

        void setPts0(std::vector<pt> *pVector);

        void calculateMatrix() override;

        void calculateMatrixWithPressure(std::vector<Point> *pressure);
    };

}


#endif //AGM_MATRIXOTHERS_H
