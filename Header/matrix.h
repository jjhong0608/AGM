//
// Created by NIMS-JUNHONG on 2020/12/31.
//

#ifndef AGM_MATRIX_H
#define AGM_MATRIX_H

#include "function2D.h"
#include "mkl_pardiso.h"

namespace AGM {
    template<typename pt>
    class matrix {
    protected:
        int *ia{}, *ja{};
        double *ent{};
        std::vector<pt> *pts{};
    public:
        matrix();

        explicit matrix(std::vector<pt> *pts);

        virtual ~matrix();

        [[nodiscard]] std::vector<pt> *getPts() const;

        void setPts(std::vector<pt> *pVector);

        [[nodiscard]] int *getIa() const;

        [[nodiscard]] int *getJa() const;

        [[nodiscard]] double *getEnt() const;

        virtual void makeMatrix();

        virtual void calculateMatrix();

        void exportMatrix();

        void deleteMatrix();
    };

}


#endif //AGM_MATRIX_H
