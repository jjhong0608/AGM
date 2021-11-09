//
// Created by NIMS-JUNHONG on 2020/12/24.
//

#include "matrixrow.h"

void AGM::matrix_row::remove(int i) {
    for (int j = 0; j < size(); ++j) {
        if (at(j).idx == i) {
            erase(begin() + j);
        }
    }
}

double &AGM::matrix_row::operator[](int i) {
    if (empty()) {
        matrix_element args;
        args.idx = i;
        args.value = ZEROVALUE;
        emplace_back(args);
        return front().value;
    } else {
        for (int j = 0; j < size(); ++j) {
            if (at(j).idx == i) {
                return at(j).value;
            } else if (at(j).idx > i) {
                matrix_element args;
                args.idx = i;
                args.value = ZEROVALUE;
                emplace(begin() + j, args);
                return at(j).value;
            }
        }
    }
    matrix_element args;
    args.idx = i;
    args.value = ZEROVALUE;
    emplace_back(args);
    return back().value;
}

AGM::matrix_row AGM::matrix_row::operator+(const AGM::matrix_row &src) const {
    auto row = matrix_row();
    row = *this;
    for (const auto &i : src) {
        row[i.idx] += i.value;
    }
    return row;
}

AGM::matrix_row AGM::matrix_row::operator-(const AGM::matrix_row &src) const {
    auto row = matrix_row();
    row = *this;
    for (const auto &i : src) {
        row[i.idx] -= i.value;
    }
    return row;
}

AGM::matrix_row AGM::matrix_row::operator*(double d) const {
    auto row = matrix_row();
    row = *this;
    for (auto &i : row) {
        row[i.idx] *= d;
    }
    return row;
}

AGM::matrix_row AGM::matrix_row::operator+=(const AGM::matrix_row &src) {
    for (const auto &i : src) {
        (*this)[i.idx] += i.value;
    }
    return *this;
}

AGM::matrix_row AGM::matrix_row::operator-=(const AGM::matrix_row &src) {
    for (const auto &i : src) {
        (*this)[i.idx] -= i.value;
    }
    return *this;
}
