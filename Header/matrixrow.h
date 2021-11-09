//
// Created by NIMS-JUNHONG on 2020/12/24.
//

#ifndef AGM_MATRIXROW_H
#define AGM_MATRIXROW_H

#include "AxialLine.h"

namespace AGM {
    struct matrix_element {
        int idx{};
        double value{};
    };

    class matrix_row : public std::vector<matrix_element> {
    public:
        void remove(int i);
        double &operator[](int i);
        matrix_row operator+(const matrix_row &src) const;
        matrix_row operator-(const matrix_row &src) const;
        matrix_row operator*(double d) const;
        matrix_row operator+=(const matrix_row &src);
        matrix_row operator-=(const matrix_row &src);
    };
}


#endif //AGM_MATRIXROW_H
