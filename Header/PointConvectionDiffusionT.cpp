//
// Created by NIMS-JUNHONG on 2021/03/04.
//

#include "PointConvectionDiffusionT.h"

std::vector<AGM::Value> *AGM::PointConvectionDiffusionT::two_levels_up_values;

std::vector<AGM::PointConvectionDiffusionT> *AGM::PointConvectionDiffusionT::ptsCDT;

std::vector<AGM::PointConvectionDiffusionT> *AGM::PointConvectionDiffusionT::getPtsCdt() {
    return ptsCDT;
}

void AGM::PointConvectionDiffusionT::setPtsCdt(std::vector<PointConvectionDiffusionT> *ptsCdt) {
    ptsCDT = ptsCdt;
}

std::vector<AGM::Coordinate> *AGM::PointConvectionDiffusionT::convections;

AGM::PointConvectionDiffusionT::PointConvectionDiffusionT() = default;

const AGM::Coordinate &AGM::PointConvectionDiffusionT::getConvection() const {
    return convection;
}

void AGM::PointConvectionDiffusionT::setConvection(const AGM::Coordinate &coefficient) {
    PointConvectionDiffusionT::convection = coefficient;
}

std::vector<AGM::Coordinate> *AGM::PointConvectionDiffusionT::getConvections() {
    return convections;
}

std::vector<AGM::Value> *AGM::PointConvectionDiffusionT::getTwoLevelsUpValues() {
    return two_levels_up_values;
}

void AGM::PointConvectionDiffusionT::setTwoLevelsUpValues(std::vector<Value> *twoLevelsUpValues) {
    two_levels_up_values = twoLevelsUpValues;
}

void AGM::PointConvectionDiffusionT::setConvections(std::vector<Coordinate> *vector) {
    PointConvectionDiffusionT::convections = vector;
}

void AGM::PointConvectionDiffusionT::findElement(std::vector<Point> *src) {
    int index{};
    auto ele = std::array<Point *, 12>{nullptr,};
    for (auto &item : src->at(idx).getElement().getElement()) {
        if (item) {
            ele.at(index) = &(ptsCDT->at(item->getIdx()));
        }
        ++index;
    }
    element.setElement(ele);
}

void AGM::PointConvectionDiffusionT::calcRepresentationFormula_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, CN * mp, CN * mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, CN * mp, CN * mp, convection[1]);
    std::array<matrix_row, 2> row;

    auto eraseInterface = [&](Point *pt, int i) -> void {
        auto checkInterface = [&](Point *pt) -> bool {
            auto getEachMp = [&](Point *pt, Point *ptl, Point *ptr) -> double {
                if (pt) {
                    return pt->getMp();
                } else {
                    if (ptl->getCondition() == 'C') return ptl->getMp();
                    else if (ptr->getCondition() == 'C') return ptr->getMp();
                    else printError("getEachMp");
                }
                return pt->getMp();
            };

            return (pt->getCondition() == 'I') &&
                   !(isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                             getEachMp((*pt)[W], (*pt)[WN], (*pt)[WS])) &&
                     isclose(getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW]),
                             getEachMp((*pt)[S], (*pt)[SE], (*pt)[SW])) &&
                     isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                             getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW])));
        };

        if (checkInterface(pt)) {
            row[i][idx + ptsnum] += row[i][pt->getIdx() + ptsnum];
            row[i].remove(pt->getIdx() + ptsnum);
        }
    };

    row[0][idx] = -UNITVALUE - HALFVALUE * UNITVALUE * gfunc_x.green_integral('c') / delta;
    row[0][element[W]->getIdx()] =
            CN * mp * gfunc_x.green_function_t(xm) - HALFVALUE * UNITVALUE * gfunc_x.green_integral('l') / delta;
    row[0][element[E]->getIdx()] =
            -CN * mp * gfunc_x.green_function_t(xp) - HALFVALUE * UNITVALUE * gfunc_x.green_integral('r') / delta;

    row[0][element[W]->getIdx() + ptsnum] = CN * gfunc_x.green_integral('l');
    row[0][idx + ptsnum] = CN * gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + ptsnum] = CN * gfunc_x.green_integral('r');

    rhsMatrixRow[0][element[W]->getIdx()] = HALFVALUE * gfunc_x.green_integral('l');
    rhsMatrixRow[0][idx] = HALFVALUE * gfunc_x.green_integral('c');
    rhsMatrixRow[0][element[E]->getIdx()] = HALFVALUE * gfunc_x.green_integral('r');

    row[1][idx] = -UNITVALUE - HALFVALUE * UNITVALUE * gfunc_y.green_integral('c') / delta;
    row[1][element[S]->getIdx()] =
            CN * mp * gfunc_y.green_function_t(ym) - HALFVALUE * UNITVALUE * gfunc_y.green_integral('l') / delta;
    row[1][element[N]->getIdx()] =
            -CN * mp * gfunc_y.green_function_t(yp) - HALFVALUE * UNITVALUE * gfunc_y.green_integral('r') / delta;

    row[1][element[S]->getIdx() + ptsnum] = -CN * gfunc_y.green_integral('l');
    row[1][idx + ptsnum] = -CN * gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + ptsnum] = -CN * gfunc_y.green_integral('r');

    rhsMatrixRow[1][element[S]->getIdx()] = HALFVALUE * gfunc_y.green_integral('l');
    rhsMatrixRow[1][idx] = HALFVALUE * gfunc_y.green_integral('c');
    rhsMatrixRow[1][element[N]->getIdx()] = HALFVALUE * gfunc_y.green_integral('r');

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);

    matrixRow[0] = row[0] + row[1];
    matrixRow[1] = row[0] - row[1];

    for (const auto &j : rhsMatrixRow[0]) {
        rb[0] -= j.value * rhsParallelToX(j.idx);
        rb[1] -= j.value * rhsParallelToX(j.idx);
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * rhsParallelToY(j.idx);
        rb[1] += j.value * rhsParallelToY(j.idx);
    }
}

void AGM::PointConvectionDiffusionT::calcRepresentationFormula_neumann() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, CN * mp, CN * mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, CN * mp, CN * mp, convection[1]);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusion(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * UNITVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * UNITVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = CN * sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = CN * sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = HALFVALUE * func.green_integral('R');

        return mRow * coefficient;
    };

    auto approximatePhi = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto approximateRhs = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 2 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 2 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionConvectionDiffusionLinear *func,
                            double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] +=
                    mp0 * func->green_function_ttau(d) - HALFVALUE * UNITVALUE * func->green_integral_tau(c) / delta;
            row[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral_tau(c);
            row[i][pt->getIdx() + 2 * ptsnum] += HALFVALUE * func->green_integral_tau(c);
        } else {
            row[i] += approximateSol(ptl, ptr,
                                     mp0 * func->green_function_ttau(d) -
                                     HALFVALUE * UNITVALUE * func->green_integral_tau(c) / delta,
                                     i0, std::abs(mp0));
            row[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral_tau(c), i0);
            row[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral_tau(c), i0);
        }
    };

    row[0][idx] = -HALFVALUE * UNITVALUE * gfunc_x.green_integral_tau('c') / delta;
    row[0][idx + ptsnum] = CN * gfunc_x.green_integral_tau('c');
    row[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -CN * mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], CN * mp, &gfunc_x, xm, 0, 1, 'l');

    row[1][idx] = -HALFVALUE * UNITVALUE * gfunc_y.green_integral_tau('c') / delta;
    row[1][idx + ptsnum] = -CN * gfunc_y.green_integral_tau('c');
    row[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -CN * mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], CN * mp, &gfunc_y, ym, 1, 0, 'l');

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] += row[i].back().value * normal[i];
            row[i].pop_back();
        }
    }

    matrixRow[0] = row[0] * normal[0] + row[1] * normal[1];
    rb[0] = value["bdv"];
    for (const auto &j : rhsMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToX(j.idx);
        } else {
            rb[0] -= j.value * rhsParallelToY(j.idx - ptsnum);
        }
    }
    for (const auto &j : rhsMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToY(j.idx);
        } else {
            rb[0] -= j.value * rhsParallelToX(j.idx - ptsnum);
        }
    }

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

void AGM::PointConvectionDiffusionT::calcRepresentationFormula_interface() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto getEachMp = [&](Point *pt, Point *ptl, Point *ptr) -> double {
        if (pt) {
            return pt->getMp();
        } else {
            if (ptl->getCondition() == 'C') return ptl->getMp();
            else if (ptr->getCondition() == 'C') return ptr->getMp();
            else printError("getEachMp");
        }
        return pt->getMp();
    };

    double mpw{getEachMp(element[W], element[WN], element[WS])};
    double mpe{getEachMp(element[E], element[EN], element[ES])};
    double mps{getEachMp(element[S], element[SE], element[SW])};
    double mpn{getEachMp(element[N], element[NE], element[NW])};

    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, CN * mpw, CN * mpe, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, CN * mps, CN * mpn, convection[1]);

    auto checkInterface = [&](Point *pt) -> bool {
        return (pt->getCondition() == 'I') &&
               !(isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                         getEachMp((*pt)[W], (*pt)[WN], (*pt)[WS])) &&
                 isclose(getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW]),
                         getEachMp((*pt)[S], (*pt)[SE], (*pt)[SW])) &&
                 isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                         getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW])));
    };

    auto checkMatrixRow = [&](matrix_row *row, Point *ptl, Point *ptr) -> void {
        if (checkInterface(ptl)) {
            (*row)[ptr->getIdx() + ptsnum] += (*row)[ptl->getIdx() + ptsnum];
            row->remove(ptl->getIdx() + ptsnum);
        }
        if (checkInterface(ptr)) {
            (*row)[ptl->getIdx() + ptsnum] += (*row)[ptr->getIdx() + ptsnum];
            row->remove(ptr->getIdx() + ptsnum);
        }
    };

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusion(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * UNITVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * UNITVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = CN * sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = CN * sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 3 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = HALFVALUE * func.green_integral('R');

        return mRow * coefficient;
    };

    auto approximatePhi = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + ptsnum] = (tb - tm) / (tp - tm);

        checkMatrixRow(&mRow, ptl, ptr);

        return mRow * coefficient;
    };

    auto approximateRhs = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 2 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 2 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionConvectionDiffusion *func,
                            double d, int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] +=
                    mp0 * func->green_function_t(d) - HALFVALUE * UNITVALUE * func->green_integral(c) / delta;
            if (!checkInterface(this)) row[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral(c);
            else row[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral(C);
            row[i][pt->getIdx() + 2 * ptsnum] += HALFVALUE * func->green_integral(c);
        } else {
            row[i] += approximateSol(ptl, ptr,
                                     mp0 * func->green_function_t(d) -
                                     HALFVALUE * UNITVALUE * func->green_integral(c) / delta, i0,
                                     std::abs(mp0));
            if (!checkInterface(this)) row[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral(c), i0);
            else row[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral(C), i0);
            row[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral(c), i0);
        }
    };

    row[0][idx] = -UNITVALUE - HALFVALUE * UNITVALUE * gfunc_x.green_integral('c') / delta;
    if (!checkInterface(this)) row[0][idx + ptsnum] = CN * gfunc_x.green_integral('c');
    row[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral('c');

    assignMatrix(element[E], element[ES], element[EN], -CN * mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], CN * mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    row[1][idx] = -UNITVALUE - HALFVALUE * UNITVALUE * gfunc_y.green_integral('c') / delta;
    if (!checkInterface(this)) row[1][idx + ptsnum] = -CN * gfunc_y.green_integral('c');
    row[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral('c');

    assignMatrix(element[N], element[NW], element[NE], -CN * mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], CN * mps, &gfunc_y, ym, 1, 0, 'l', 'L');

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] += row[i].back().value;
            row[i].pop_back();
        }
    }

    matrixRow[0] = row[0] + row[1];
    for (const auto &j : rhsMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToX(j.idx);
        } else {
            rb[0] -= j.value * rhsParallelToY(j.idx - ptsnum);
        }
    }
    for (const auto &j : rhsMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToY(j.idx);
        } else {
            rb[0] -= j.value * rhsParallelToX(j.idx - ptsnum);
        }
    }
    if (isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mpn)) {
        matrixRow[1] = row[0] - row[1];
        for (const auto &j : rhsMatrixRow[0]) {
            if (j.idx < ptsnum) {
                rb[1] -= j.value * rhsParallelToX(j.idx);
            } else {
                rb[1] -= j.value * rhsParallelToY(j.idx - ptsnum);
            }
        }
        for (const auto &j : rhsMatrixRow[1]) {
            if (j.idx < ptsnum) {
                rb[1] += j.value * rhsParallelToY(j.idx);
            } else {
                rb[1] += j.value * rhsParallelToX(j.idx - ptsnum);
            }
        }
    } else {
        matrixRow[1][idx + ptsnum] = UNITVALUE;
    }
}

void AGM::PointConvectionDiffusionT::makeDifferentiation_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, CN * mp, CN * mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, CN * mp, CN * mp, convection[1]);

    auto eraseInterface = [&](Point *pt, int i) -> void {
        auto checkInterface = [&](Point *pt) -> bool {
            auto getEachMp = [&](Point *pt, Point *ptl, Point *ptr) -> double {
                if (pt) {
                    return pt->getMp();
                } else {
                    if (ptl->getCondition() == 'C') return ptl->getMp();
                    else if (ptr->getCondition() == 'C') return ptr->getMp();
                    else printError("getEachMp");
                }
                return pt->getMp();
            };

            return (pt->getCondition() == 'I') &&
                   !(isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                             getEachMp((*pt)[W], (*pt)[WN], (*pt)[WS])) &&
                     isclose(getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW]),
                             getEachMp((*pt)[S], (*pt)[SE], (*pt)[SW])) &&
                     isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                             getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW])));
        };

        if (checkInterface(pt)) {
            dMatrixRow[i][idx + ptsnum] += dMatrixRow[i][pt->getIdx() + ptsnum];
            dMatrixRow[i].remove(pt->getIdx() + ptsnum);
        }
    };

    dMatrixRow[0][idx] = -HALFVALUE * UNITVALUE * gfunc_x.green_integral_tau('c') / delta;
    dMatrixRow[0][element[W]->getIdx()] = CN * mp * gfunc_x.green_function_ttau(xm)
                                          - HALFVALUE * UNITVALUE * gfunc_x.green_integral_tau('l') / delta;
    dMatrixRow[0][element[E]->getIdx()] = -CN * mp * gfunc_x.green_function_ttau(xp)
                                          - HALFVALUE * UNITVALUE * gfunc_x.green_integral_tau('r') / delta;

    dMatrixRow[0][element[W]->getIdx() + ptsnum] = CN * gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + ptsnum] = CN * gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->getIdx() + ptsnum] = CN * gfunc_x.green_integral_tau('r');

    dMatrixRow[0][element[W]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('r');

    dMatrixRow[1][idx] = -HALFVALUE * UNITVALUE * gfunc_y.green_integral_tau('c') / delta;
    dMatrixRow[1][element[S]->getIdx()] = CN * mp * gfunc_y.green_function_ttau(ym)
                                          - HALFVALUE * UNITVALUE * gfunc_y.green_integral_tau('l') / delta;
    dMatrixRow[1][element[N]->getIdx()] = -CN * mp * gfunc_y.green_function_ttau(yp)
                                          - HALFVALUE * UNITVALUE * gfunc_y.green_integral_tau('r') / delta;

    dMatrixRow[1][element[S]->getIdx() + ptsnum] = -CN * gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + ptsnum] = -CN * gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + ptsnum] = -CN * gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('r');

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);
}

void AGM::PointConvectionDiffusionT::makeDifferentiation_boundary() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, CN * mp, CN * mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, CN * mp, CN * mp, convection[1]);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusion(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * UNITVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * UNITVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = CN * sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = CN * sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = HALFVALUE * func.green_integral('R');

        return mRow * coefficient;
    };

    auto approximatePhi = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto approximateRhs = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 2 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 2 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionConvectionDiffusionLinear *func,
                            double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] += mp0 * func->green_function_ttau(d)
                                           - HALFVALUE * UNITVALUE * func->green_integral_tau(c) / delta;
            dMatrixRow[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral_tau(c);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] += HALFVALUE * func->green_integral_tau(c);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr,
                                            mp0 * func->green_function_ttau(d) -
                                            HALFVALUE * UNITVALUE * func->green_integral_tau(c) / delta,
                                            i0, std::abs(mp0));
            dMatrixRow[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral_tau(c), i0);
        }
    };

    dMatrixRow[0][idx] = -HALFVALUE * UNITVALUE * gfunc_x.green_integral_tau('c') / delta;
    dMatrixRow[0][idx + ptsnum] = CN * gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -CN * mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], CN * mp, &gfunc_x, xm, 0, 1, 'l');

    dMatrixRow[1][idx] = -HALFVALUE * UNITVALUE * gfunc_y.green_integral_tau('c') / delta;
    dMatrixRow[1][idx + ptsnum] = -CN * gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -CN * mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], CN * mp, &gfunc_y, ym, 1, 0, 'l');
}

void AGM::PointConvectionDiffusionT::makeDifferentiation_interface() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto getEachMp = [&](Point *pt, Point *ptl, Point *ptr) -> double {
        if (pt) {
            return pt->getMp();
        } else {
            if (ptl->getCondition() == 'C') return ptl->getMp();
            else if (ptr->getCondition() == 'C') return ptr->getMp();
            else printError("getEachMp");
        }
        return pt->getMp();
    };

    double mpw{getEachMp(element[W], element[WN], element[WS])};
    double mpe{getEachMp(element[E], element[EN], element[ES])};
    double mps{getEachMp(element[S], element[SE], element[SW])};
    double mpn{getEachMp(element[N], element[NE], element[NW])};

    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, CN * mpw, CN * mpe, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, CN * mps, CN * mpn, convection[1]);

    auto checkInterface = [&](Point *pt) -> bool {
        return (pt->getCondition() == 'I') &&
               !(isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                         getEachMp((*pt)[W], (*pt)[WN], (*pt)[WS])) &&
                 isclose(getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW]),
                         getEachMp((*pt)[S], (*pt)[SE], (*pt)[SW])) &&
                 isclose(getEachMp((*pt)[E], (*pt)[EN], (*pt)[ES]),
                         getEachMp((*pt)[N], (*pt)[NE], (*pt)[NW])));
    };

    auto checkMatrixRow = [&](matrix_row *row, Point *ptl, Point *ptr) -> void {
        if (checkInterface(ptl)) {
            (*row)[ptr->getIdx() + ptsnum] += (*row)[ptl->getIdx() + ptsnum];
            row->remove(ptl->getIdx() + ptsnum);
        }
        if (checkInterface(ptr)) {
            (*row)[ptl->getIdx() + ptsnum] += (*row)[ptr->getIdx() + ptsnum];
            row->remove(ptr->getIdx() + ptsnum);
        }
    };

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusion(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * UNITVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * UNITVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = CN * sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = CN * sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        dMatrixRow[(i + 1) % 2][ptl->getIdx() + 3 * ptsnum] += HALFVALUE * func.green_integral('L') * coefficient;
        dMatrixRow[(i + 1) % 2][ptr->getIdx() + 3 * ptsnum] += HALFVALUE * func.green_integral('R') * coefficient;

        return mRow * coefficient;
    };

    auto approximatePhi = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + ptsnum] = (tb - tm) / (tp - tm);

        checkMatrixRow(&mRow, ptl, ptr);

        return mRow * coefficient;
    };

    auto approximateRhs = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 2 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 2 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionConvectionDiffusion *func,
                            double d, int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] +=
                    mp0 * func->green_function_ttau(d) - HALFVALUE * UNITVALUE * func->green_integral_tau(c) / delta;
            if (!checkInterface(this)) dMatrixRow[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral_tau(c);
            else dMatrixRow[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral_tau(C);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] += HALFVALUE * func->green_integral_tau(c);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr,
                                            mp0 * func->green_function_ttau(d) -
                                            HALFVALUE * UNITVALUE * func->green_integral_tau(c) / delta, i0, std::abs(mp0));
            if (!checkInterface(this))
                dMatrixRow[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral_tau(c), i0);
            else dMatrixRow[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral_tau(C), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral_tau(c), i0);
        }
    };

    dMatrixRow[0][idx] += -HALFVALUE * UNITVALUE * gfunc_x.green_integral_tau('c') / delta;
    if (!checkInterface(this)) dMatrixRow[0][idx + ptsnum] += CN * gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] += HALFVALUE * gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -CN * mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], CN * mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    dMatrixRow[1][idx] += -HALFVALUE * UNITVALUE * gfunc_y.green_integral_tau('c') / delta;
    if (!checkInterface(this)) dMatrixRow[1][idx + ptsnum] += -CN * gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] += HALFVALUE * gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -CN * mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], CN * mps, &gfunc_y, ym, 1, 0, 'l', 'L');
}

void AGM::PointConvectionDiffusionT::calcDifferentiation() {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * ptsCDT->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * ptsCDT->at(item.idx - ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dx += item.value * rhsParallelToX(item.idx - 2 * ptsnum);
        } else {
            dx += item.value * rhsParallelToY(item.idx - 3 * ptsnum);
        }
    });
    value["dx"] = dx;
    value["dxx"] = (convection[0] * dx -
                    (CN * value["phi"] - HALFVALUE * UNITVALUE * value["sol"] / delta + HALFVALUE * rhsParallelToX(idx))) /
                   (CN * mp);

    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * ptsCDT->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * ptsCDT->at(item.idx % ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dy += item.value * rhsParallelToY(item.idx - 2 * ptsnum);
        } else {
            dy += item.value * rhsParallelToX(item.idx - 3 * ptsnum);
        }
    });
    value["dy"] = dy;
    value["dyy"] = (convection[1] * dy -
                    (-CN * value["phi"] - HALFVALUE * UNITVALUE * value["sol"] / delta + HALFVALUE * rhsParallelToY(idx))) /
                   (CN * mp);

//    if (getCondition() == 'D') {
//        if (axialLine[0]) {
//            if (!iszero((*this) - *(element[E]))) {
//                value["dx"] = element[E]->getValue()["dx"];
//                value["dy"] = element[E]->getValue()["dy"];
//            } else if (!iszero((*this) - *(element[W]))) {
//                value["dx"] = element[W]->getValue()["dx"];
//                value["dy"] = element[W]->getValue()["dy"];
//            } else {
//                printError("AGM::Point::calcContinuity");
//            }
//        }
//        if (axialLine[1]) {
//            if (!iszero((*this) - *(element[N]))) {
//                value["dx"] = element[N]->getValue()["dx"];
//                value["dy"] = element[N]->getValue()["dy"];
//            } else if (!iszero((*this) - *(element[S]))) {
//                value["dx"] = element[S]->getValue()["dx"];
//                value["dy"] = element[S]->getValue()["dy"];
//            } else {
//                printError("AGM::Point::calcContinuity");
//            }
//        }
//    }

}

