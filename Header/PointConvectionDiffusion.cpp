//
// Created by NIMS-JUNHONG on 2021/01/26.
//

#include "PointConvectionDiffusion.h"

std::vector<AGM::PointConvectionDiffusion> *AGM::PointConvectionDiffusion::ptsCD;

std::vector<AGM::Coordinate> *AGM::PointConvectionDiffusion::convections;

AGM::PointConvectionDiffusion::PointConvectionDiffusion() = default;

std::vector<AGM::Coordinate> *AGM::PointConvectionDiffusion::getConvections() {
    return convections;
}

void AGM::PointConvectionDiffusion::setConvections(std::vector<Coordinate> *vector) {
    PointConvectionDiffusion::convections = vector;
}

const AGM::Coordinate &AGM::PointConvectionDiffusion::getConvection() const {
    return convection;
}

void AGM::PointConvectionDiffusion::setConvection(const AGM::Coordinate &coefficient) {
    PointConvectionDiffusion::convection = coefficient;
}

std::vector<AGM::PointConvectionDiffusion> *AGM::PointConvectionDiffusion::getPtsCd() {
    return ptsCD;
}

void AGM::PointConvectionDiffusion::setPtsCd(std::vector<PointConvectionDiffusion> *ptsCd) {
    ptsCD = ptsCd;
}

void AGM::PointConvectionDiffusion::findElement(std::vector<Point> *src) {
    int index{};
    auto ele = std::array<Point *, 12> {nullptr, };
    for (auto &item : src->at(idx).getElement().getElement()) {
        if (item) {
            ele.at(index) = &(ptsCD->at(item->getIdx()));
        }
        ++index;
    }
    element.setElement(ele);
}

void AGM::PointConvectionDiffusion::calcRepresentationFormula_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, mp, mp, convection[1]);
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

    row[0][idx] = -UNITVALUE;
    row[0][element[W]->getIdx()] = mp * gfunc_x.green_function_t(xm);
    row[0][element[E]->getIdx()] = -mp * gfunc_x.green_function_t(xp);

    row[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral('r');

    rhsMatrixRow[0][element[W]->getIdx()] = HALFVALUE * gfunc_x.green_integral('l');
    rhsMatrixRow[0][idx] = HALFVALUE * gfunc_x.green_integral('c');
    rhsMatrixRow[0][element[E]->getIdx()] = HALFVALUE * gfunc_x.green_integral('r');

    row[1][idx] = -UNITVALUE;
    row[1][element[S]->getIdx()] = mp * gfunc_y.green_function_t(ym);
    row[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_t(yp);

    row[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral('l');
    row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral('r');

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
        rb[0] -= j.value * ptsCD->at(j.idx)["rhs"];
        rb[1] -= j.value * ptsCD->at(j.idx)["rhs"];
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * ptsCD->at(j.idx)["rhs"];
        rb[1] += j.value * ptsCD->at(j.idx)["rhs"];
    }
}

void AGM::PointConvectionDiffusion::calcRepresentationFormula_neumann() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convection[1]);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusionLinear(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('R');

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
            row[i][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->getIdx() + 2 * ptsnum] += HALFVALUE * func->green_integral_tau(c);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            row[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral_tau(c), i0);
        }
    };

    row[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    row[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    row[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx % ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
    }

    matrixRow[0] = row[0] * normal[0] + row[1] * normal[1];

    rb[0] = value["bdv"];
    for (const auto &i : rhsMatrixRow) {
        for (const auto &j : i) {
            rb[0] -= j.value * ptsCD->at(j.idx).getValue()["rhs"];
        }
    }

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

void AGM::PointConvectionDiffusion::calcRepresentationFormula_interface() {
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

    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, mpw, mpe, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, mps, mpn, convection[1]);

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
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('R');

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
            row[i][pt->getIdx()] += mp0 * func->green_function_t(d);
            if (!checkInterface(this)) row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(c);
            else row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(C);
            row[i][pt->getIdx() + 2 * ptsnum] += HALFVALUE * func->green_integral(c);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_t(d), i0, std::abs(mp0));
            if (!checkInterface(this)) row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(c), i0);
            else row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(C), i0);
            row[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral(c), i0);
        }
    };

    row[0][idx] = -UNITVALUE;
    if (!checkInterface(this)) row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral('c');

    assignMatrix(element[E], element[ES], element[EN], -mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    row[1][idx] = -UNITVALUE;
    if (!checkInterface(this)) row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral('c');

    assignMatrix(element[N], element[NW], element[NE], -mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], mps, &gfunc_y, ym, 1, 0, 'l', 'L');

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx % ptsnum] = row[i].back().value;
            row[i].pop_back();
        }
    }

    matrixRow[0] = row[0] + row[1];
    for (const auto &j : rhsMatrixRow[0]) {
        rb[0] -= j.value * ptsCD->at(j.idx)["rhs"];
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * ptsCD->at(j.idx)["rhs"];
    }
    if (isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mpn)) {
        matrixRow[1] = row[0] - row[1];
        for (const auto &j : rhsMatrixRow[0]) {
            rb[1] -= j.value * ptsCD->at(j.idx)["rhs"];
        }
        for (const auto &j : rhsMatrixRow[1]) {
            rb[1] += j.value * ptsCD->at(j.idx)["rhs"];
        }
    } else {
        matrixRow[1][idx + ptsnum] = UNITVALUE;
    }
}

void AGM::PointConvectionDiffusion::makeDifferentiation_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, mp, mp, convection[1]);

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

    dMatrixRow[0][element[W]->getIdx()] = mp * gfunc_x.green_function_ttau(xm);
    dMatrixRow[0][element[E]->getIdx()] = -mp * gfunc_x.green_function_ttau(xp);

    dMatrixRow[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral_tau('r');

    dMatrixRow[0][element[W]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('r');

    dMatrixRow[1][element[S]->getIdx()] = mp * gfunc_y.green_function_ttau(ym);
    dMatrixRow[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_ttau(yp);

    dMatrixRow[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('r');

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);
}

void AGM::PointConvectionDiffusion::makeDifferentiation_boundary() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convection[1]);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -1.0E0 : 1.0E0;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusion(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('R');

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
        double sign = i ? -1.0E0 : 1.0E0;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            dMatrixRow[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] += HALFVALUE * func->green_integral_tau(c);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral_tau(c), i0);
        }
    };

    dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');
}

void AGM::PointConvectionDiffusion::makeDifferentiation_interface() {
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

    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, mpw, mpe, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, mps, mpn, convection[1]);

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
        double sign = i ? -1.0E0 : 1.0E0;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusion(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('R');

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
        double sign = i ? -1.0E0 : 1.0E0;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] = mp0 * func->green_function_ttau(d);
            if (!checkInterface(this)) dMatrixRow[i][pt->getIdx() + ptsnum] = sign * func->green_integral_tau(c);
            else dMatrixRow[i][pt->getIdx() + ptsnum] = sign * func->green_integral_tau(C);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] = HALFVALUE * func->green_integral_tau(c);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            if (!checkInterface(this))
                dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            else dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(C), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, HALFVALUE * func->green_integral_tau(c), i0);
        }
    };

    if (!checkInterface(this)) dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] = HALFVALUE * gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    if (!checkInterface(this)) dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] = HALFVALUE * gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], mps, &gfunc_y, ym, 1, 0, 'l', 'L');
}

void AGM::PointConvectionDiffusion::calcDifferentiation() {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * ptsCD->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * ptsCD->at(item.idx % ptsnum)["phi"];
        } else {
            dx += item.value * ptsCD->at(item.idx % ptsnum)["rhs"];
        }
    });
    value["dx"] = dx;
    value["dxx"] = (convection[0] * dx - (value["phi"] + HALFVALUE * value["rhs"])) / mp;

    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * ptsCD->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * ptsCD->at(item.idx % ptsnum)["phi"];
        } else {
            dy += item.value * ptsCD->at(item.idx % ptsnum)["rhs"];
        }
    });
    value["dy"] = dy;
    value["dyy"] = (convection[1] * dy - (-value["phi"] + HALFVALUE * value["rhs"])) / mp;
}

void AGM::PointConvectionDiffusion::calcContinuity_cross(std::vector<Point> *uvel, std::vector<Point> *vvel) {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convection[1]);

    auto u = std::array<double, 4>(), v = std::array<double, 4>(), p = std::array<double, 4>();
    auto phi = std::array<double, 4>(), psi = std::array<double, 4>();
    auto f1 = std::array<double, 4>(), f2 = std::array<double, 4>();

    auto assignVal = [&](EWNS e) -> void {
        u[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["sol"];
        v[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["sol"];
        phi[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["phi"];
        psi[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["phi"];
        f1[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["rhs"];
        f2[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["rhs"];
        p[e] = element[e]->getValue()["sol"];
    };

    assignVal(E);
    assignVal(W);
    assignVal(N);
    assignVal(S);

    value["sol"] = -HALFVALUE * mp *
                   (mp * gfunc_x.green_function_ttau(xm) * u[W] -
                    mp * gfunc_x.green_function_ttau(xp) * u[E] +
                    mp * gfunc_y.green_function_ttau(ym) * v[S] -
                    mp * gfunc_y.green_function_ttau(yp) * v[N] +
                    gfunc_x.green_integral_tau('l') * phi[W] +
                    gfunc_x.green_integral_tau('c') * uvel->at(idx)["phi"] +
                    gfunc_x.green_integral_tau('r') * phi[E] -
                    gfunc_y.green_integral_tau('l') * psi[S] -
                    gfunc_y.green_integral_tau('c') * vvel->at(idx)["phi"] -
                    gfunc_y.green_integral_tau('r') * psi[N] +
                    HALFVALUE * gfunc_x.green_integral_tau('l') * f1[W] +
                    HALFVALUE * gfunc_x.green_integral_tau('c') * uvel->at(idx)["rhs"] +
                    HALFVALUE * gfunc_x.green_integral_tau('r') * f1[E] +
                    HALFVALUE * gfunc_y.green_integral_tau('l') * f2[S] +
                    HALFVALUE * gfunc_y.green_integral_tau('c') * vvel->at(idx)["rhs"] +
                    HALFVALUE * gfunc_y.green_integral_tau('r') * f2[N] +
                    gfunc_x.green_integral_ttau('l') * p[W] +
                    gfunc_x.green_integral_ttau('c') * value["sol"] +
                    gfunc_x.green_integral_ttau('r') * p[E] +
                    gfunc_y.green_integral_ttau('l') * p[S] +
                    gfunc_y.green_integral_ttau('c') * value["sol"] +
                    gfunc_y.green_integral_ttau('r') * p[N]);
}

void AGM::PointConvectionDiffusion::calcContinuity_boundary(std::vector<Point> *uvel, std::vector<Point> *vvel) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convection[1]);

    auto approximateSol = [&](Point *ptl, Point *ptr, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusionLinear(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('R');

        mRow[ptl->getIdx() + 4 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 4 * ptsnum] = func.green_integral_t('R');

        return mRow;
    };

    auto approximate = [&](Point *ptl, Point *ptr, const std::string &string, int i) -> double {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        return ptl->getValue()[string] * (tp - tb) / (tp - tm) + ptr->getValue()[string] * (tb - tm) / (tp - tm);
    };

    auto assignSol = [&](Point *ptl, Point *ptr, int i0, std::vector<Point> *vector, int i) -> double {
        auto row = approximateSol(ptl, ptr, i0, mp);
        double sol{};
        for (const auto &item : row) {
            if (item.idx < ptsnum) {
                sol += item.value * vector->at(item.idx).getValue()["sol"];
            } else if (item.idx < 2 * ptsnum) {
                sol += item.value * vector->at(item.idx - ptsnum).getValue()["phi"];
            } else if (item.idx < 3 * ptsnum) {
                sol += item.value * vector->at(item.idx - 2 * ptsnum).getValue()["rhs"];
            }
            if (i == i0) {
                sol += item.value * value["sol"];
            }
        }
        return sol;
    };

    auto u = std::array<double, 4>(), v = std::array<double, 4>(), p = std::array<double, 4>();
    auto phi = std::array<double, 4>(), psi = std::array<double, 4>();
    auto f1 = std::array<double, 4>(), f2 = std::array<double, 4>();

    auto assignVal = [&](EWNS e, EWNS e0, EWNS e1, int i, int j) -> void {
        if (element[e]) {
            u[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["sol"];
            v[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["sol"];
            phi[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["phi"];
            psi[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["phi"];
            f1[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["rhs"];
            f2[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["rhs"];
            p[e] = element[e]->getValue()["sol"];
        } else {
            u[e] = assignSol(uvel->at(idx).getElement().getElement()[e0], uvel->at(idx).getElement().getElement()[e1],
                             j, uvel, 0);
            v[e] = assignSol(vvel->at(idx).getElement().getElement()[e0], vvel->at(idx).getElement().getElement()[e1],
                             j, vvel, 1);
            phi[e] = approximate(uvel->at(idx).getElement().getElement()[e0],
                                 uvel->at(idx).getElement().getElement()[e1], "phi", i);
            psi[e] = approximate(vvel->at(idx).getElement().getElement()[e0],
                                 vvel->at(idx).getElement().getElement()[e1], "phi", i);
            f1[e] = approximate(uvel->at(idx).getElement().getElement()[e0],
                                uvel->at(idx).getElement().getElement()[e1], "rhs", i);
            f2[e] = approximate(vvel->at(idx).getElement().getElement()[e0],
                                vvel->at(idx).getElement().getElement()[e1], "rhs", i);
            p[e] = approximate(element[e0], element[e1], "sol", i);
        }
    };

    assignVal(E, EN, ES, 0, 1);
    assignVal(W, WN, WS, 0, 1);
    assignVal(N, NE, NW, 1, 0);
    assignVal(S, SE, SW, 1, 0);

    value["sol"] = -HALFVALUE * mp *
                   (mp * gfunc_x.green_function_ttau(xm) * u[W] -
                    mp * gfunc_x.green_function_ttau(xp) * u[E] +
                    mp * gfunc_y.green_function_ttau(ym) * v[S] -
                    mp * gfunc_y.green_function_ttau(yp) * v[N] +
                    gfunc_x.green_integral_tau('l') * phi[W] +
                    gfunc_x.green_integral_tau('c') * uvel->at(idx)["phi"] +
                    gfunc_x.green_integral_tau('r') * phi[E] -
                    gfunc_y.green_integral_tau('l') * psi[S] -
                    gfunc_y.green_integral_tau('c') * vvel->at(idx)["phi"] -
                    gfunc_y.green_integral_tau('r') * psi[N] +
                    HALFVALUE * gfunc_x.green_integral_tau('l') * f1[W] +
                    HALFVALUE * gfunc_x.green_integral_tau('c') * uvel->at(idx)["rhs"] +
                    HALFVALUE * gfunc_x.green_integral_tau('r') * f1[E] +
                    HALFVALUE * gfunc_y.green_integral_tau('l') * f2[S] +
                    HALFVALUE * gfunc_y.green_integral_tau('c') * vvel->at(idx)["rhs"] +
                    HALFVALUE * gfunc_y.green_integral_tau('r') * f2[N] +
                    gfunc_x.green_integral_ttau('l') * p[W] +
                    gfunc_x.green_integral_ttau('c') * value["sol"] +
                    gfunc_x.green_integral_ttau('r') * p[E] +
                    gfunc_y.green_integral_ttau('l') * p[S] +
                    gfunc_y.green_integral_ttau('c') * value["sol"] +
                    gfunc_y.green_integral_ttau('r') * p[N]);
}

void AGM::PointConvectionDiffusion::calcContinuity_interface(std::vector<Point> *uvel, std::vector<Point> *vvel) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, mp, mp, convection[1]);

    auto approximateSol = [&](Point *ptl, Point *ptr, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        double cv = (convections->at(ptl->getIdx())[i] * (tp - tb) + convections->at(ptr->getIdx())[i] * (tb - tm)) /
                    (tp - tm);
        auto func = GreenfunctionConvectionDiffusion(tm, tb, tp, d, d, cv);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('L');
        mRow[ptr->getIdx() + 2 * ptsnum] = HALFVALUE * func.green_integral('R');

        mRow[ptl->getIdx() + 4 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 4 * ptsnum] = func.green_integral_t('R');

        return mRow;
    };

    auto approximate = [&](Point *ptl, Point *ptr, const std::string &string, int i) -> double {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        return ptl->getValue()[string] * (tp - tb) / (tp - tm) + ptr->getValue()[string] * (tb - tm) / (tp - tm);
    };

    auto assignSol = [&](Point *ptl, Point *ptr, int i0, std::vector<Point> *vector, int i) -> double {
        auto row = approximateSol(ptl, ptr, i0, mp);
        double sol{};
        for (const auto &item : row) {
            if (item.idx < ptsnum) {
                sol += item.value * vector->at(item.idx).getValue()["sol"];
            } else if (item.idx < 2 * ptsnum) {
                sol += item.value * vector->at(item.idx - ptsnum).getValue()["phi"];
            } else if (item.idx < 3 * ptsnum) {
                sol += item.value * vector->at(item.idx - 2 * ptsnum).getValue()["rhs"];
            }
            if (i == i0) {
                sol += item.value * value["sol"];
            }
        }
        return sol;
    };

    auto u = std::array<double, 4>(), v = std::array<double, 4>(), p = std::array<double, 4>();
    auto phi = std::array<double, 4>(), psi = std::array<double, 4>();
    auto f1 = std::array<double, 4>(), f2 = std::array<double, 4>();

    auto assignVal = [&](EWNS e, EWNS e0, EWNS e1, int i, int j) -> void {
        if (element[e]) {
            u[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["sol"];
            v[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["sol"];
            phi[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["phi"];
            psi[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["phi"];
            f1[e] = uvel->at(idx).getElement().getElement()[e]->getValue()["rhs"];
            f2[e] = vvel->at(idx).getElement().getElement()[e]->getValue()["rhs"];
            p[e] = element[e]->getValue()["sol"];
        } else {
            u[e] = assignSol(uvel->at(idx).getElement().getElement()[e0], uvel->at(idx).getElement().getElement()[e1],
                             j, uvel, 0);
            v[e] = assignSol(vvel->at(idx).getElement().getElement()[e0], vvel->at(idx).getElement().getElement()[e1],
                             j, vvel, 1);
            phi[e] = approximate(uvel->at(idx).getElement().getElement()[e0],
                                 uvel->at(idx).getElement().getElement()[e1], "phi", i);
            psi[e] = approximate(vvel->at(idx).getElement().getElement()[e0],
                                 vvel->at(idx).getElement().getElement()[e1], "phi", i);
            f1[e] = approximate(uvel->at(idx).getElement().getElement()[e0],
                                uvel->at(idx).getElement().getElement()[e1], "rhs", i);
            f2[e] = approximate(vvel->at(idx).getElement().getElement()[e0],
                                vvel->at(idx).getElement().getElement()[e1], "rhs", i);
            p[e] = approximate(element[e0], element[e1], "sol", i);
        }
    };

    assignVal(E, EN, ES, 0, 1);
    assignVal(W, WN, WS, 0, 1);
    assignVal(N, NE, NW, 1, 0);
    assignVal(S, SE, SW, 1, 0);

    value["sol"] = -HALFVALUE * mp *
                   (mp * gfunc_x.green_function_ttau(xm) * u[W] -
                    mp * gfunc_x.green_function_ttau(xp) * u[E] +
                    mp * gfunc_y.green_function_ttau(ym) * v[S] -
                    mp * gfunc_y.green_function_ttau(yp) * v[N] +
                    gfunc_x.green_integral_tau('l') * phi[W] +
                    gfunc_x.green_integral_tau('c') * uvel->at(idx)["phi"] +
                    gfunc_x.green_integral_tau('r') * phi[E] -
                    gfunc_y.green_integral_tau('l') * psi[S] -
                    gfunc_y.green_integral_tau('c') * vvel->at(idx)["phi"] -
                    gfunc_y.green_integral_tau('r') * psi[N] +
                    HALFVALUE * gfunc_x.green_integral_tau('l') * f1[W] +
                    HALFVALUE * gfunc_x.green_integral_tau('c') * uvel->at(idx)["rhs"] +
                    HALFVALUE * gfunc_x.green_integral_tau('r') * f1[E] +
                    HALFVALUE * gfunc_y.green_integral_tau('l') * f2[S] +
                    HALFVALUE * gfunc_y.green_integral_tau('c') * vvel->at(idx)["rhs"] +
                    HALFVALUE * gfunc_y.green_integral_tau('r') * f2[N] +
                    gfunc_x.green_integral_ttau('l') * p[W] +
                    gfunc_x.green_integral_ttau('c') * value["sol"] +
                    gfunc_x.green_integral_ttau('r') * p[E] +
                    gfunc_y.green_integral_ttau('l') * p[S] +
                    gfunc_y.green_integral_ttau('c') * value["sol"] +
                    gfunc_y.green_integral_ttau('r') * p[N]);
}

void AGM::PointConvectionDiffusion::makePressureTerm_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convection[1]);

    pMatrixRow[0][element[W]->getIdx()] = gfunc_x.green_integral_t('l');
    pMatrixRow[0][getIdx()] = gfunc_x.green_integral_t('c');
    pMatrixRow[0][element[E]->getIdx()] = gfunc_x.green_integral_t('r');

    pMatrixRow[1][element[S]->getIdx()] = gfunc_y.green_integral_t('l');
    pMatrixRow[1][getIdx()] = gfunc_y.green_integral_t('c');
    pMatrixRow[1][element[N]->getIdx()] = gfunc_y.green_integral_t('r');
}

void AGM::PointConvectionDiffusion::makePressureTerm_neumann() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convection[1]);

    auto approximateP = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, GreenfunctionConvectionDiffusionLinear *func, int i,
                            int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            pMatrixRow[i][pt->getIdx()] += func->green_integral_ttau(c);
        } else {
            pMatrixRow[i] += approximateP(ptl, ptr, sign * func->green_integral_ttau(c), i0);
        }
    };

    pMatrixRow[0][idx] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[E], element[ES], element[EN], &gfunc_x, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], &gfunc_x, 0, 1, 'l');

    pMatrixRow[1][idx] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[N], element[NW], element[NE], &gfunc_y, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], &gfunc_y, 1, 0, 'l');}

void AGM::PointConvectionDiffusion::makePressureTerm_interface() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, mp, mp, convection[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, mp, mp, convection[1]);

    auto approximateP = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, GreenfunctionConvectionDiffusion *func, int i, int i0,
                            char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            pMatrixRow[i][pt->getIdx()] += func->green_integral_ttau(c);
        } else {
            pMatrixRow[i] += approximateP(ptl, ptr, sign * func->green_integral_ttau(c), i0);
        }
    };

    pMatrixRow[0][idx] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[E], element[ES], element[EN], &gfunc_x, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], &gfunc_x, 0, 1, 'l');

    pMatrixRow[1][idx] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[N], element[NW], element[NE], &gfunc_y, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], &gfunc_y, 1, 0, 'l');
}

void AGM::PointConvectionDiffusion::updateRb_cross() {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        rb[0] -= j.value * ptsCD->at(j.idx)["rhs"];
        rb[1] -= j.value * ptsCD->at(j.idx)["rhs"];
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * ptsCD->at(j.idx)["rhs"];
        rb[1] += j.value * ptsCD->at(j.idx)["rhs"];
    }
}

void AGM::PointConvectionDiffusion::updateRb_neumann() {
    rb[0] = value["bdv"];
    for (const auto &i : rhsMatrixRow) {
        for (const auto &j : i) {
            rb[0] -= j.value * ptsCD->at(j.idx).getValue()["rhs"];
        }
    }
}

void AGM::PointConvectionDiffusion::updateRb_interface() {
    updateRb_cross();
}
