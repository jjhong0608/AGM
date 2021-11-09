//
// Created by NIMS-JUNHONG on 2021/04/07.
//

#ifndef AGM_MATRIXEIGEN_H
#define AGM_MATRIXEIGEN_H

#include "matrixOthers.h"
#include <Eigen/Sparse>

namespace AGM {
    template<typename pt>
    class matrixEIGEN : public matrix<pt> {
    private:
        Eigen::SparseMatrix<double, Eigen::RowMajor> A{};
        int fixedPtrIdx{};
        int *ia_t{}, *ja_t{};
        double *ent_t{};
        int n_entry{};
    public:
        matrixEIGEN();

        explicit matrixEIGEN(std::vector<pt> *pts);

        matrixEIGEN(std::vector<pt> *pts, int fixedPtrIdx);

        ~matrixEIGEN() override;

        [[nodiscard]] int getFixedPtrIdx() const;

        void setFixedPtrIdx(int i);

        void makeMatrix() override;

        void calculateMatrix() override;

        void exportTranseposedMatrix();

        void exportNormatMatrix();
    };

}

#endif //AGM_MATRIXEIGEN_H
