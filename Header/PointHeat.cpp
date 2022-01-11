//
// Created by NIMS-JUNHONG on 2021/01/25.
//

#include "PointHeat.h"

double AGM::PointHeat::time;
double AGM::PointHeat::delta;
double AGM::PointHeat::CN;
double AGM::PointHeat::BDF;
std::vector<AGM::Value> *AGM::PointHeat::previous_values;
std::function<double(int)> AGM::PointHeat::rhsParallelToX;
std::function<double(int)> AGM::PointHeat::rhsParallelToY;
std::vector<AGM::PointHeat> *AGM::PointHeat::ptsH;

double AGM::PointHeat::getTime() {
    return time;
}

void AGM::PointHeat::setTime(double d) {
    PointHeat::time = d;
}

double AGM::PointHeat::getDelta() {
    return delta;
}

void AGM::PointHeat::setDelta(double d) {
    PointHeat::delta = d;
}

double AGM::PointHeat::getCn() {
    return CN;
}

void AGM::PointHeat::setCn(double cn) {
    CN = cn;
}

double AGM::PointHeat::getBdf() {
    return BDF;
}

void AGM::PointHeat::setBdf(double bfd) {
    BDF = bfd;
}

std::vector<AGM::Value> *AGM::PointHeat::getPreviousValues() {
    return previous_values;
}

void AGM::PointHeat::setPreviousValues(std::vector<AGM::Value> *previousValues) {
    previous_values = previousValues;
}

const std::function<double(int)> &AGM::PointHeat::getRhsParallelToX() {
    return rhsParallelToX;
}

void AGM::PointHeat::setRhsParallelToX(const std::function<double(int)> &function) {
    PointHeat::rhsParallelToX = function;
}

const std::function<double(int)> &AGM::PointHeat::getRhsParallelToY() {
    return rhsParallelToY;
}

void AGM::PointHeat::setRhsParallelToY(const std::function<double(int)> &function) {
    PointHeat::rhsParallelToY = function;
}

std::vector<AGM::PointHeat> *AGM::PointHeat::getPtsH() {
    return ptsH;
}

void AGM::PointHeat::setPtsH(std::vector<PointHeat> *pVector) {
    PointHeat::ptsH = pVector;
}

void AGM::PointHeat::findElement(std::vector<AGM::Point> *src) {
    int index{};
    auto ele = std::array<Point *, 12> {nullptr, };
    for (auto &item : src->at(idx).getElement().getElement()) {
        if (item) {
            ele.at(index) = &(ptsH->at(item->getIdx()));
        }
        ++index;
    }
    element.setElement(ele);
}

void AGM::PointHeat::findElement1(std::vector<AGM::Point> *src) {
    int index{};
    auto ele = std::array<Point *, 12> {nullptr, };
    for (auto &item : src->at(idx).getElement1().getElement()) {
        if (item) {
            ele.at(index) = &(ptsH->at(item->getIdx()));
        }
        ++index;
    }
    element1.setElement(ele);
}

void AGM::PointHeat::calcRepresentationFormula_cross() {
    std::array<matrix_row, 2> row;
    calcRepresentationFormula_cross_RD();
    row[0] = matrixRow[0];
    row[1] = matrixRow[1];
//    std::string string{};
//    if (std::any_of(matrixRow[0].begin(), matrixRow[0].end(), [](matrix_element &matrixElement) -> bool {
//        if (matrixElement.idx < ptsnum) {
//            return std::fabs(matrixElement.value) < 1.0E-10;
//        } else {
//            return false;
//        }
//    })) {
//        string += "x";
//    }
//    if (std::any_of(matrixRow[1].begin(), matrixRow[1].end(), [](matrix_element &matrixElement) -> bool {
//        if (matrixElement.idx < ptsnum) {
//            return std::fabs(matrixElement.value) < 1.0E-10;
//        } else {
//            return false;
//        }
//    })) {
//        string += "y";
//    }
//    if (string.empty()) {
//        row[0] = matrixRow[0];
//        row[1] = matrixRow[1];
//    } else {
//        if (string == "x") {
//            row[0] = calcRepresentationFormula_cross_E(string);
//            row[1] = matrixRow[1];
//        } else if (string == "y") {
//            row[0] = matrixRow[0];
//            row[1] = calcRepresentationFormula_cross_E(string);
//        } else {
//            matrixRow[0] = matrix_row{};
//            matrixRow[1] = matrix_row{};
//            calcRepresentationFormula_cross_E();
//            row[0] = matrixRow[0];
//            row[1] = matrixRow[1];
//        }
//    }

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 4 * ptsnum) {
            pMatrixRow[i][row[i].back().idx - 4 * ptsnum] = row[i].back().value;
            row[i].pop_back();
        }
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value;
            row[i].pop_back();
        }
    }

    matrixRow[0] = row[0] + row[1];
    matrixRow[1] = row[0] - row[1];
}

void AGM::PointHeat::calcRepresentationFormula_cross_E() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);
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
                return ZEROVALUE;
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

    row[0][idx] = -UNITVALUE - gfunc_x.green_integral('c') / delta;
    row[0][element[W]->getIdx()] = mp * gfunc_x.green_function_t(xm) - gfunc_x.green_integral('l') / delta;
    row[0][element[E]->getIdx()] = -mp * gfunc_x.green_function_t(xp) - gfunc_x.green_integral('r') / delta;

    row[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('l');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');
    row[0][element[E]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('r');

    row[1][idx] = -UNITVALUE - gfunc_y.green_integral('c') / delta;
    row[1][element[S]->getIdx()] = mp * gfunc_y.green_function_t(ym) - gfunc_y.green_integral('l') / delta;
    row[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_t(yp) - gfunc_y.green_integral('r') / delta;

    row[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral('l');
    row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('l');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('l');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');
    row[1][element[N]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('r');

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);

    matrixRow[0] = row[0];
    matrixRow[1] = row[1];
}

void AGM::PointHeat::calcRepresentationFormula_cross_RD() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);
    std::array<matrix_row, 2> row;

    row[0][idx] = -UNITVALUE;
    row[0][element[W]->getIdx()] = mp * gfunc_x.green_function_t(xm);
    row[0][element[E]->getIdx()] = -mp * gfunc_x.green_function_t(xp);

    row[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('l');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');
    row[0][element[E]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('r');

    row[1][idx] = -UNITVALUE;
    row[1][element[S]->getIdx()] = mp * gfunc_y.green_function_t(ym);
    row[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_t(yp);

    row[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral('l');
    row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('l');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('l');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');
    row[1][element[N]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('r');

    matrixRow[0] = row[0];
    matrixRow[1] = row[1];
}

AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_cross_E(std::string &string) {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);
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
                return ZEROVALUE;
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

    row[0][idx] = -UNITVALUE - gfunc_x.green_integral('c') / delta;
    row[0][element[W]->getIdx()] = mp * gfunc_x.green_function_t(xm) - gfunc_x.green_integral('l') / delta;
    row[0][element[E]->getIdx()] = -mp * gfunc_x.green_function_t(xp) - gfunc_x.green_integral('r') / delta;

    row[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('l');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');
    row[0][element[E]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('r');

    row[1][idx] = -UNITVALUE - gfunc_y.green_integral('c') / delta;
    row[1][element[S]->getIdx()] = mp * gfunc_y.green_function_t(ym) - gfunc_y.green_integral('l') / delta;
    row[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_t(yp) - gfunc_y.green_integral('r') / delta;

    row[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral('l');
    row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('l');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('l');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');
    row[1][element[N]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('r');

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);

    if (string == "x") {
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_cross_E(std::string &string)");
    }
    return row[0];
}

AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_cross_RD(std::string &string) {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);
    std::array<matrix_row, 2> row;

    row[0][idx] = -UNITVALUE;
    row[0][element[W]->getIdx()] = mp * gfunc_x.green_function_t(xm);
    row[0][element[E]->getIdx()] = -mp * gfunc_x.green_function_t(xp);

    row[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral('r');

    row[0][element[W]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('l');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');
    row[0][element[E]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t('r');

    row[1][idx] = -UNITVALUE;
    row[1][element[S]->getIdx()] = mp * gfunc_y.green_function_t(ym);
    row[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_t(yp);

    row[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral('l');
    row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('l');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][element[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral('r');

    row[1][element[S]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('l');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');
    row[1][element[N]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t('r');

    if (string == "x") {
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_cross_RD(std::string &string)");
    }
    return row[0];
}

void AGM::PointHeat::calcRepresentationFormula_neumann() {
//    calcRepresentationFormula_neumann_NN_RD();
    std::string string{};
    if (element[E] == this || element[W] == this) {
        string += "x";
    }
    if (element[N] == this || element[S] == this) {
        string += "y";
    }
    makeDifferentiation_boundary();
    if (!string.empty()) {
        std::array<matrix_row, 2> row{};
        if (string == "x") {
            if (this->getAxialLine('x')) {
                row[0] = calcRepresentationFormula_neumann_NN_RD(string);
            } else {
                row[0] = calcRepresentationFormula_neumann_D_RD(string);
            }
        } else if (string == "y") {
            if (this->getAxialLine('y')) {
                row[1] = calcRepresentationFormula_neumann_NN_RD(string);
            } else {
                row[1] = calcRepresentationFormula_neumann_D_RD(string);
            }
        } else {
            string = "x";
            if (this->getAxialLine('x')) {
                row[0] = calcRepresentationFormula_neumann_NN_RD(string);
            } else {
                row[0] = calcRepresentationFormula_neumann_D_RD(string);
            }
            string = "y";
            if (this->getAxialLine('y')) {
                row[1] = calcRepresentationFormula_neumann_NN_RD(string);
            } else {
                row[1] = calcRepresentationFormula_neumann_D_RD(string);
            }
        }
        for (int i = 0; i < 2; ++i) {
            if (!row[i].empty() && !iszero(normal[i])) {
                while (row[i].back().idx >= 4 * ptsnum) {
                    pMatrixRow[i][row[i].back().idx - 4 * ptsnum] = row[i].back().value * normal[i];
                    row[i].pop_back();
                }
                while (row[i].back().idx >= 2 * ptsnum) {
                    rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value * normal[i];
                    row[i].pop_back();
                }
            }
        }

        if (iszero(normal[0])) {
            matrixRow[0] = row[1] * normal[1];
        } else if (iszero(normal[1])) {
            matrixRow[0] = row[0] * normal[0];
        } else {
            matrixRow[0] = row[0] * normal[0] + row[1] * normal[1];
        }

        AGM::Point *pt{};
        EWNS ewns = findPoint_dirichlet_and_Neumann(pt);
        calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
    } else {
        calcRepresentationFormula_neumann_D();
    }
}

void AGM::PointHeat::calcRepresentationFormula_neumann_D_E() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d, int i, int i0,
                            char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] +=
                    mp0 * func->green_function_ttau(d) - HALFVALUE * func->green_integral_tau(c) / delta;
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d) -
                                               HALFVALUE * func->green_integral_tau(c) / delta, i0, std::abs(mp0));
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_ttau(c) + sign1 * func->green_function_tau(d),
                                       i0);
        }
    };

    row[0][idx] = -gfunc_x.green_integral_tau('c') / delta;
    row[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c');

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], CN * mp, &gfunc_x, xm, 0, 1, 'l');

    row[1][idx] = -gfunc_y.green_integral_tau('c') / delta;
    row[1][idx + ptsnum] = -CN * gfunc_y.green_integral_tau('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -CN * mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], CN * mp, &gfunc_y, ym, 1, 0, 'l');

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] += row[i].back().value * normal[i];
            row[i].pop_back();
        }
    }

    matrixRow[0] = row[0] * normal[0] + row[1] * normal[1];

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

void AGM::PointHeat::calcRepresentationFormula_neumann_D_RD() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
            int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_ttau(c) + sign1 * func->green_function_tau(d),
                                       i0);
        }
    };

    row[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    row[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 4 * ptsnum) {
            if (!iszero(normal[i])) pMatrixRow[i][row[i].back().idx - 4 * ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
        while (row[i].back().idx >= 2 * ptsnum) {
            if (!iszero(normal[i])) rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
    }

    if (iszero(normal[0])) {
        matrixRow[0] = row[1] * normal[1];
    } else if (iszero(normal[1])) {
        matrixRow[0] = row[0] * normal[0];
    } else {
        matrixRow[0] = row[0] * normal[0] + row[1] * normal[1];
    }

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);

    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

void AGM::PointHeat::calcRepresentationFormula_neumann_NN_E() {
    double xm{}, xb{}, xp{}, ym{}, yb{}, yp{};
    EWNS nx{NW}, ny{WS};
    double xx{}, yy{};

    auto findPts = [&](EWNS ewns, EWNS ewns1, EWNS ewns2, int i, double *tm, double *tb, double *tp, EWNS *nn) -> void {
        if (element1[ewns] == this || element1[ewns1] == this) {
            if (element1[ewns] == this) {
                *tm = element1[ewns1] ? element1[ewns1]->getCoordinate()[i] : element1[ewns2]->getCoordinate()[i];
                *tb = element[ewns1] ? element[ewns1]->getCoordinate()[i] : element[ewns2]->getCoordinate()[i];
                *tp = coordinate[i];
                *nn = ewns;
            }
        } else {
            if (i) {
                ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
                yb = getCoordinate()[1];
                yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
                ny = N;
            } else {
                xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
                xb = getCoordinate()[0];
                xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
                nx = E;
            }
        }
    };

    findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
    findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
    findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);

    if (iszero(normal[0])) {
        xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
        xb = getCoordinate()[0];
        xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
        nx = E;
    }
    if (iszero(normal[1])) {
        ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
        yb = getCoordinate()[1];
        yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
        ny = N;
    }

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };


    std::array<matrix_row, 2> row{};

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE + HALFVALUE * func->green_integral_ND(c) / delta;
            } else {
                row[i][pt->getIdx()] += UNITVALUE - HALFVALUE * func->green_integral_ND(c) / delta;
            }

            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_ND(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE - HALFVALUE * func->green_integral_ND(c) / delta, i0,
                                         std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE - HALFVALUE * func->green_integral_ND(c) / delta, i0,
                                         std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_ND(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_ND(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_ND(c), i0);
        }
    };

    auto assignMatrix_DN = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE + HALFVALUE * func->green_integral_DN(c) / delta;
            } else {
                row[i][pt->getIdx()] += UNITVALUE - HALFVALUE * func->green_integral_DN(c) / delta;
            }
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_DN(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE - HALFVALUE * func->green_integral_DN(c) / delta, i0,
                                         std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE - HALFVALUE * func->green_integral_DN(c) / delta, i0,
                                         std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_DN(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_DN(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_DN(c), i0);
        }
    };

    if (nx == E) {
        xx = mp * gfunc_x.green_function_DN(xp);
        row[0][idx] = -HALFVALUE * gfunc_x.green_integral_DN('r') / delta;
        row[0][idx + ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_DN('r') - gfunc_x.green_function_DN(xp);

        assignMatrix_DN(element[W], element[WS], element[WN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_DN(element1[W], element1[WS], element1[WN], mp, &gfunc_x, xm, 0, 1, 'l');
    } else if (nx == W) {
        xx = -mp * gfunc_x.green_function_ND(xm);
        row[0][idx] = -HALFVALUE * gfunc_x.green_integral_ND('l') / delta;
        row[0][idx + ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_ND('l') + gfunc_x.green_function_ND(xm);

        assignMatrix_ND(element[E], element[ES], element[EN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_ND(element1[E], element1[ES], element1[EN], mp, &gfunc_x, xp, 0, 1, 'r');
    } else {
        printError("AGM::PointHeat::calcRepresentationFormula_neumann1()");
    }

    if (ny == N) {
        yy = mp * gfunc_y.green_function_DN(yp);
        row[1][idx] = -HALFVALUE * gfunc_y.green_integral_DN('r') / delta;
        row[1][idx + ptsnum] = -gfunc_y.green_integral_DN('r');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_DN('r');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_DN('r') - gfunc_y.green_function_DN(yp);

        assignMatrix_DN(element[S], element[SW], element[SE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_DN(element1[S], element1[SW], element1[SE], mp, &gfunc_y, ym, 1, 0, 'l');
    } else if (ny == S) {
        yy = -mp * gfunc_y.green_function_ND(ym);
        row[1][idx] = -HALFVALUE * gfunc_y.green_integral_ND('r') / delta;
        row[1][idx + ptsnum] = -gfunc_y.green_integral_ND('l');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_ND('l');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_ND('l') + gfunc_y.green_function_ND(ym);

        assignMatrix_ND(element[N], element[NW], element[NE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_ND(element1[N], element1[NW], element1[NE], mp, &gfunc_y, yp, 1, 0, 'r');
    } else {
        printError("AGM::PointHeat::calcRepresentationFormula_neumann1()");
    }

    double c{};
    for (int i = 0; i < 2; ++i) {
        c = i ? -yy : -xx;
        for (auto &item : row[i]) {
            item.value /= c;
        }
    }


    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 4 * ptsnum) {
            if (!iszero(normal[i])) pMatrixRow[i][row[i].back().idx - 4 * ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
        while (row[i].back().idx >= 2 * ptsnum) {
            if (!iszero(normal[i])) rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
    }

    if (iszero(normal[0])) {
        matrixRow[0] = row[1] * normal[1];
    } else if (iszero(normal[1])) {
        matrixRow[0] = row[0] * normal[0];
    } else {
        matrixRow[0] = row[0] * normal[0] + row[1] * normal[1];
    }

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);
//    if (pt->getCondition() == 'C') {
//        calcRepresentationFormula_boundary();
//    } else {
//        calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
//        printError("calcRepresentationFormula_neumann1()");
//    }
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

void AGM::PointHeat::calcRepresentationFormula_neumann_NN_RD() {
    double xm{}, xb{}, xp{}, ym{}, yb{}, yp{};
    EWNS nx{NW}, ny{WS};

    auto findPts = [&](EWNS ewns, EWNS ewns1, EWNS ewns2, int i, double *tm, double *tb, double *tp, EWNS *nn) -> void {
        if (element1[ewns] == this || element1[ewns1] == this) {
            if (element1[ewns] == this) {
                *tm = element1[ewns1] ? element1[ewns1]->getCoordinate()[i] : element1[ewns2]->getCoordinate()[i];
                *tb = element[ewns1] ? element[ewns1]->getCoordinate()[i] : element[ewns2]->getCoordinate()[i];
                *tp = coordinate[i];
                *nn = ewns;
            }
        } else {
            if (i) {
                ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
                yb = getCoordinate()[1];
                yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
                ny = N;
            } else {
                xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
                xb = getCoordinate()[0];
                xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
                nx = E;
            }
        }
    };

    findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
    findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
    findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);

    if (iszero(normal[0])) {
        xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
        xb = getCoordinate()[0];
        xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
        nx = E;
    }
    if (iszero(normal[1])) {
        ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
        yb = getCoordinate()[1];
        yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
        ny = N;
    }

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                               double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE;
            } else {
                row[i][pt->getIdx()] += -mp * func->green_function_t_ND(d);
            }

            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_ND(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, -mp * func->green_function_t_ND(d), i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_ND(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_ND(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_ND(c), i0);
        }
    };

    auto assignMatrix_DN = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                               double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE;
            } else {
                row[i][pt->getIdx()] += mp * func->green_function_t_DN(d);
            }
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_DN(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, mp * func->green_function_t_DN(d), i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_DN(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_DN(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_DN(c), i0);
        }
    };

    if (nx == E) {
        row[0][idx] = mp * gfunc_x.green_function_DN(xp);
        row[0][idx + ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_DN('r') - gfunc_x.green_function_DN(xp);

        assignMatrix_DN(element[W], element[WS], element[WN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_DN(element1[W], element1[WS], element1[WN], mp, &gfunc_x, xm, 0, 1, 'l');
    } else if (nx == W) {
        row[0][idx] = -mp * gfunc_x.green_function_ND(xm);
        row[0][idx + ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_ND('l') + gfunc_x.green_function_ND(xm);

        assignMatrix_ND(element[E], element[ES], element[EN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_ND(element1[E], element1[ES], element1[EN], mp, &gfunc_x, xp, 0, 1, 'r');
    } else {
        printError("AGM::Point::calcRepresentationFormula_neumann1()");
    }

    if (ny == N) {
        row[1][idx] = mp * gfunc_y.green_function_DN(yp);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_DN('r');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_DN('r');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_DN('r') - gfunc_y.green_function_DN(yp);

        assignMatrix_DN(element[S], element[SW], element[SE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_DN(element1[S], element1[SW], element1[SE], mp, &gfunc_y, ym, 1, 0, 'l');
    } else if (ny == S) {
        row[1][idx] = -mp * gfunc_y.green_function_ND(ym);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_ND('l');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_ND('l');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_ND('l') + gfunc_y.green_function_ND(ym);

        assignMatrix_ND(element[N], element[NW], element[NE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_ND(element1[N], element1[NW], element1[NE], mp, &gfunc_y, yp, 1, 0, 'r');
    } else {
        printError("AGM::Point::calcRepresentationFormula_neumann1()");
    }

    double c{};
    for (int i = 0; i < 2; ++i) {
        c = -row[i][idx];
        for (auto &item : row[i]) {
            item.value /= c;
        }
        row[i][idx] = ZEROVALUE;
    }

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 4 * ptsnum) {
            if (!iszero(normal[i])) pMatrixRow[i][row[i].back().idx - 4 * ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
        while (row[i].back().idx >= 2 * ptsnum) {
            if (!iszero(normal[i])) rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
    }

    if (iszero(normal[0])) {
        matrixRow[0] = row[1] * normal[1];
    } else if (iszero(normal[1])) {
        matrixRow[0] = row[0] * normal[0];
    } else {
        matrixRow[0] = row[0] * normal[0] + row[1] * normal[1];
    }

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);

//    if (pt->getCondition() == 'C') {
//        calcRepresentationFormula_boundary();
//    } else {
//        calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
//        printError("calcRepresentationFormula_neumann1()");
//    }
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_neumann_D_RD(std::string &string) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                            double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_ttau(c) + sign1 * func->green_function_tau(d),
                                       i0);
        }
    };

    row[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    row[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');

    if (string == "x") {
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_neumann_D_RD(std::string &string)");
    }
    return row[0];
}

AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_neumann_NN_E(std::string &string) {
    double xm{}, xb{}, xp{}, ym{}, yb{}, yp{};
    EWNS nx{NW}, ny{WS};
    double xx{}, yy{};

    auto findPts = [&](EWNS ewns, EWNS ewns1, EWNS ewns2, int i, double *tm, double *tb, double *tp, EWNS *nn) -> void {
        if (element1[ewns] == this || element1[ewns1] == this) {
            if (element1[ewns] == this) {
                *tm = element1[ewns1] ? element1[ewns1]->getCoordinate()[i] : element1[ewns2]->getCoordinate()[i];
                *tb = element[ewns1] ? element[ewns1]->getCoordinate()[i] : element[ewns2]->getCoordinate()[i];
                *tp = coordinate[i];
                *nn = ewns;
            }
        } else {
            if (i) {
                ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
                yb = getCoordinate()[1];
                yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
                ny = N;
            } else {
                xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
                xb = getCoordinate()[0];
                xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
                nx = E;
            }
        }
    };

    findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
    findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
    findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);

    if (iszero(normal[0])) {
        xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
        xb = getCoordinate()[0];
        xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
        nx = E;
    }
    if (iszero(normal[1])) {
        ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
        yb = getCoordinate()[1];
        yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
        ny = N;
    }

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE + func->green_integral_ND(c) / delta;
            } else {
                row[i][pt->getIdx()] += UNITVALUE - func->green_integral_ND(c) / delta;
            }

            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_ND(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE - func->green_integral_ND(c) / delta, i0,
                                         std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE - func->green_integral_ND(c) / delta, i0,
                                         std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_ND(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_ND(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_ND(c), i0);
        }
    };

    auto assignMatrix_DN = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE + func->green_integral_DN(c) / delta;
            } else {
                row[i][pt->getIdx()] += UNITVALUE - func->green_integral_DN(c) / delta;
            }
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_DN(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE - func->green_integral_DN(c) / delta, i0,
                                         std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE - func->green_integral_DN(c) / delta, i0,
                                         std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_DN(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_DN(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_DN(c), i0);
        }
    };

    if (nx == E) {
        xx = mp * gfunc_x.green_function_DN(xp);
        row[0][idx] = -gfunc_x.green_integral_DN('r') / delta;
        row[0][idx + ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_DN('r') - gfunc_x.green_function_DN(xp);

        assignMatrix_DN(element[W], element[WS], element[WN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_DN(element1[W], element1[WS], element1[WN], mp, &gfunc_x, xm, 0, 1, 'l');
    } else if (nx == W) {
        xx = -mp * gfunc_x.green_function_ND(xm);
        row[0][idx] = -gfunc_x.green_integral_ND('l') / delta;
        row[0][idx + ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_ND('l') + gfunc_x.green_function_ND(xm);

        assignMatrix_ND(element[E], element[ES], element[EN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_ND(element1[E], element1[ES], element1[EN], mp, &gfunc_x, xp, 0, 1, 'r');
    } else {
        printError("AGM::PointHeat::calcRepresentationFormula_neumann1()");
    }

    if (ny == N) {
        yy = mp * gfunc_y.green_function_DN(yp);
        row[1][idx] = -gfunc_y.green_integral_DN('r') / delta;
        row[1][idx + ptsnum] = -gfunc_y.green_integral_DN('r');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_DN('r');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_DN('r') - gfunc_y.green_function_DN(yp);

        assignMatrix_DN(element[S], element[SW], element[SE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_DN(element1[S], element1[SW], element1[SE], mp, &gfunc_y, ym, 1, 0, 'l');
    } else if (ny == S) {
        yy = -mp * gfunc_y.green_function_ND(ym);
        row[1][idx] = -gfunc_y.green_integral_ND('r') / delta;
        row[1][idx + ptsnum] = -gfunc_y.green_integral_ND('l');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_ND('l');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_ND('l') + gfunc_y.green_function_ND(ym);

        assignMatrix_ND(element[N], element[NW], element[NE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_ND(element1[N], element1[NW], element1[NE], mp, &gfunc_y, yp, 1, 0, 'r');
    } else {
        printError("AGM::PointHeat::calcRepresentationFormula_neumann1()");
    }

    double c{};
    for (int i = 0; i < 2; ++i) {
        c = i ? -yy : -xx;
        for (auto &item : row[i]) {
            item.value /= c;
        }
    }

    if (string == "x") {
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_neumann_E_RD(std::string &string)");
    }
    return row[0];
}

AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_neumann_NN_RD(std::string &string) {
    double xm{}, xb{}, xp{}, ym{}, yb{}, yp{};
    EWNS nx{NW}, ny{WS};

    auto findPts = [&](EWNS ewns, EWNS ewns1, EWNS ewns2, int i, double *tm, double *tb, double *tp, EWNS *nn) -> void {
        if (element1[ewns] == this || element1[ewns1] == this) {
            if (element1[ewns] == this) {
                *tm = element1[ewns1] ? element1[ewns1]->getCoordinate()[i] : element1[ewns2]->getCoordinate()[i];
                *tb = element[ewns1] ? element[ewns1]->getCoordinate()[i] : element[ewns2]->getCoordinate()[i];
                *tp = coordinate[i];
                *nn = ewns;
            }
        } else {
            if (i) {
                ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
                yb = getCoordinate()[1];
                yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
                ny = N;
            } else {
                xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
                xb = getCoordinate()[0];
                xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
                nx = E;
            }
        }
    };

    findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
    findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
    findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);

    if (iszero(normal[0])) {
        xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
        xb = getCoordinate()[0];
        xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
        nx = E;
    }
    if (iszero(normal[1])) {
        ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
        yb = getCoordinate()[1];
        yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
        ny = N;
    }

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                               double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE;
            } else {
                row[i][pt->getIdx()] += -mp * func->green_function_t_ND(d);
            }

            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_ND(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, -mp * func->green_function_t_ND(d), i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_ND(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_ND(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_ND(c), i0);
        }
    };

    auto assignMatrix_DN = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                               double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE;
            } else {
                row[i][pt->getIdx()] += mp * func->green_function_t_DN(d);
            }
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_DN(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, mp * func->green_function_t_DN(d), i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_DN(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_DN(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_DN(c), i0);
        }
    };

    if (nx == E) {
        row[0][idx] = mp * gfunc_x.green_function_DN(xp);
        row[0][idx + ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_DN('r') - gfunc_x.green_function_DN(xp);

        assignMatrix_DN(element[W], element[WS], element[WN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_DN(element1[W], element1[WS], element1[WN], mp, &gfunc_x, xm, 0, 1, 'l');
    } else if (nx == W) {
        row[0][idx] = -mp * gfunc_x.green_function_ND(xm);
        row[0][idx + ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_ND('l') + gfunc_x.green_function_ND(xm);

        assignMatrix_ND(element[E], element[ES], element[EN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_ND(element1[E], element1[ES], element1[EN], mp, &gfunc_x, xp, 0, 1, 'r');
    } else {
        printError("AGM::Point::calcRepresentationFormula_neumann1()");
    }

    if (ny == N) {
        row[1][idx] = mp * gfunc_y.green_function_DN(yp);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_DN('r');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_DN('r');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_DN('r') - gfunc_y.green_function_DN(yp);

        assignMatrix_DN(element[S], element[SW], element[SE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_DN(element1[S], element1[SW], element1[SE], mp, &gfunc_y, ym, 1, 0, 'l');
    } else if (ny == S) {
        row[1][idx] = -mp * gfunc_y.green_function_ND(ym);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_ND('l');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_ND('l');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_ND('l') + gfunc_y.green_function_ND(ym);

        assignMatrix_ND(element[N], element[NW], element[NE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_ND(element1[N], element1[NW], element1[NE], mp, &gfunc_y, yp, 1, 0, 'r');
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_neumann_NN_RD(std::string &string)");
    }

    std::array<double, 2> c{};
    for (int i = 0; i < 2; ++i) {
        c[i] = -row[i][idx];
        for (auto &item : row[i]) {
            item.value /= c[i];
        }
        row[i][idx] = ZEROVALUE;
    }

    if (string == "x") {
        if (fabs(c[0]) < 1e-10) {
            return calcRepresentationFormula_neumann_NN_E(string);
        }
        return row[0];
    } else if (string == "y") {
        if (fabs(c[1]) < 1e-10) {
            return calcRepresentationFormula_neumann_NN_E(string);
        }
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_neumann_NN_RD(std::string &string)");
    }
    return row[0];
}

void AGM::PointHeat::calcRepresentationFormula_interface() {
    auto getEachMp = [&](Point *pt, Point *ptl, Point *ptr) -> double {
        if (pt) {
            return pt->getMp();
        } else {
            if (ptl->getCondition() == 'C') return ptl->getMp();
            else if (ptr->getCondition() == 'C') return ptr->getMp();
            else printError("getEachMp");
        }
        return ZEROVALUE;
    };

    double mpw{getEachMp(element[W], element[WN], element[WS])};
    double mpe{getEachMp(element[E], element[EN], element[ES])};
    double mps{getEachMp(element[S], element[SE], element[SW])};
    double mpn{getEachMp(element[N], element[NE], element[NW])};

    std::array<matrix_row, 2> row;

    calcRepresentationFormula_interface_RD();
    row[0] = matrixRow[0];
    row[1] = matrixRow[1];
//    std::string string{};
//    if (std::any_of(matrixRow[0].begin(), matrixRow[0].end(), [](matrix_element &matrixElement) -> bool {
//        if (matrixElement.idx < ptsnum) {
//            return std::fabs(matrixElement.value) < 1.0E-10;
//        } else {
//            return false;
//        }
//    })) {
//        string += "x";
//    }
//    if (std::any_of(matrixRow[1].begin(), matrixRow[1].end(), [](matrix_element &matrixElement) -> bool {
//        if (matrixElement.idx < ptsnum) {
//            return std::fabs(matrixElement.value) < 1.0E-10;
//        } else {
//            return false;
//        }
//    })) {
//        string += "y";
//    }
//    if (string.empty()) {
//        row[0] = matrixRow[0];
//        row[1] = matrixRow[1];
//    } else {
//        if (string == "x") {
//            row[0] = calcRepresentationFormula_interface_E(string);
//            row[1] = matrixRow[1];
//        } else if (string == "y") {
//            row[0] = matrixRow[0];
//            row[1] = calcRepresentationFormula_interface_E(string);
//        } else {
//            matrixRow[0] = matrix_row{};
//            matrixRow[1] = matrix_row{};
//            calcRepresentationFormula_interface_E();
//            row[0] = matrixRow[0];
//            row[1] = matrixRow[1];
//        }
//    }

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 4 * ptsnum) {
            pMatrixRow[i][row[i].back().idx - 4 * ptsnum] = row[i].back().value;
            row[i].pop_back();
        }
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value;
            row[i].pop_back();
        }
    }

    matrixRow[0] = row[0] + row[1];
    if (isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mpn)) {
        matrixRow[1] = row[0] - row[1];
    } else {
        exit(2);
        matrixRow[1][idx + ptsnum] = UNITVALUE;
    }
}

void AGM::PointHeat::calcRepresentationFormula_interface_E() {
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
        return ZEROVALUE;
    };

    double mpw{getEachMp(element[W], element[WN], element[WS])};
    double mpe{getEachMp(element[E], element[EN], element[ES])};
    double mps{getEachMp(element[S], element[SE], element[SW])};
    double mpn{getEachMp(element[N], element[NE], element[NW])};

    auto gfunc_x = Greenfunction(xm, xb, xp, mpw, mpe);
    auto gfunc_y = Greenfunction(ym, yb, yp, mps, mpn);

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

//    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
//        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
//        double sign = i ? -UNITVALUE : UNITVALUE;
//        auto func = Greenfunction(tm, tb, tp, d, d);
//        auto mRow = matrix_row();
//        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - func.green_integral('L') / delta;
//        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - func.green_integral('R') / delta;
//
//        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
//        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');
//
//        checkMatrixRow(&mRow, ptl, ptr);
//
//        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
//        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');
//
//        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
//        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');
//
//        return mRow * coefficient;
//    };

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
                            int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] += mp0 * func->green_function_t(d) - func->green_integral(c) / delta;
            if (!checkInterface(this)) row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(c);
            else row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(C);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t(c);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_t(d) -
                                               func->green_integral(c) / delta, i0, std::abs(mp0));
            if (!checkInterface(this)) row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(c), i0);
            else row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(C), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t(c), i0);
        }
    };

    row[0][idx] = -UNITVALUE - gfunc_x.green_integral('c') / delta;
    if (!checkInterface(this)) row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');

    assignMatrix(element[E], element[ES], element[EN], -mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    row[1][idx] = -UNITVALUE - gfunc_y.green_integral('c') / delta;
    if (!checkInterface(this)) row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');

    assignMatrix(element[N], element[NW], element[NE], -mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], mps, &gfunc_y, ym, 1, 0, 'l', 'L');

    matrixRow[0] = row[0];
    matrixRow[1] = row[1];
}

void AGM::PointHeat::calcRepresentationFormula_interface_RD() {
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

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mpw, mpe, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mps, mpn, UNITVALUE / delta);

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
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                            double d, int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] += mp0 * func->green_function_t(d);
            if (!checkInterface(this)) row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(c);
            else row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(C);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t(c);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_t(d), i0, std::abs(mp0));
            if (!checkInterface(this)) row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(c), i0);
            else row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(C), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t(c), i0);
        }
    };

    row[0][idx] = -UNITVALUE;
    if (!checkInterface(this)) row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');

    assignMatrix(element[E], element[ES], element[EN], -mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    row[1][idx] = -UNITVALUE;
    if (!checkInterface(this)) row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');

    assignMatrix(element[N], element[NW], element[NE], -mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], mps, &gfunc_y, ym, 1, 0, 'l', 'L');

    matrixRow[0] = row[0];
    matrixRow[1] = row[1];
}

AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_interface_E(std::string &string) {
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
        return ZEROVALUE;
    };

    double mpw{getEachMp(element[W], element[WN], element[WS])};
    double mpe{getEachMp(element[E], element[EN], element[ES])};
    double mps{getEachMp(element[S], element[SE], element[SW])};
    double mpn{getEachMp(element[N], element[NE], element[NW])};

    auto gfunc_x = Greenfunction(xm, xb, xp, mpw, mpe);
    auto gfunc_y = Greenfunction(ym, yb, yp, mps, mpn);

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

//    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
//        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
//        double sign = i ? -UNITVALUE : UNITVALUE;
//        auto func = Greenfunction(tm, tb, tp, d, d);
//        auto mRow = matrix_row();
//        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - func.green_integral('L') / delta;
//        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - func.green_integral('R') / delta;
//
//        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
//        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');
//
//        checkMatrixRow(&mRow, ptl, ptr);
//
//        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
//        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');
//
//        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
//        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');
//
//        return mRow * coefficient;
//    };

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
                            int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] += mp0 * func->green_function_t(d) - func->green_integral(c) / delta;
            if (!checkInterface(this)) row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(c);
            else row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(C);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t(c);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_t(d) -
                                               func->green_integral(c) / delta, i0, std::abs(mp0));
            if (!checkInterface(this)) row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(c), i0);
            else row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(C), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t(c), i0);
        }
    };

    row[0][idx] = -UNITVALUE - gfunc_x.green_integral('c') / delta;
    if (!checkInterface(this)) row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');

    assignMatrix(element[E], element[ES], element[EN], -mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    row[1][idx] = -UNITVALUE - gfunc_y.green_integral('c') / delta;
    if (!checkInterface(this)) row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');

    assignMatrix(element[N], element[NW], element[NE], -mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], mps, &gfunc_y, ym, 1, 0, 'l', 'L');

    if (string == "x") {
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_interface_E(std::string &string)");
    }
    return row[0];
}

AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_interface_RD(std::string &string) {
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

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mpw, mpe, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mps, mpn, UNITVALUE / delta);

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
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                            double d, int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] += mp0 * func->green_function_t(d);
            if (!checkInterface(this)) row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(c);
            else row[i][pt->getIdx() + ptsnum] += sign * func->green_integral(C);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t(c);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_t(d), i0, std::abs(mp0));
            if (!checkInterface(this)) row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(c), i0);
            else row[i] += approximatePhi(ptl, ptr, sign * func->green_integral(C), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t(c), i0);
        }
    };

    row[0][idx] = -UNITVALUE;
    if (!checkInterface(this)) row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t('c');

    assignMatrix(element[E], element[ES], element[EN], -mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    row[1][idx] = -UNITVALUE;
    if (!checkInterface(this)) row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral('c');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t('c');

    assignMatrix(element[N], element[NW], element[NE], -mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], mps, &gfunc_y, ym, 1, 0, 'l', 'L');

    if (string == "x") {
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::PointHeat::calcRepresentationFormula_interface_RD(std::string &string)");
    }
    return row[0];
}

void AGM::PointHeat::calcRepresentationFormula_boundary() {
    matrixRow[1][idx + ptsnum] = UNITVALUE;
    return;


    double xm{}, xb{}, xp{}, ym{}, yb{}, yp{};
    EWNS nx{NW}, ny{WS};

    auto findPts = [&](EWNS ewns, EWNS ewns1, EWNS ewns2, int i, double *tm, double *tb, double *tp, EWNS *nn) -> void {
        if (element1[ewns] == this || element1[ewns1] == this) {
            if (element1[ewns] == this) {
                *tm = element1[ewns1] ? element1[ewns1]->getCoordinate()[i] : element1[ewns2]->getCoordinate()[i];
                *tb = element[ewns1] ? element[ewns1]->getCoordinate()[i] : element[ewns2]->getCoordinate()[i];
                *tp = coordinate[i];
                *nn = ewns;
            }
        }
    };

    if (axialLine[0] && (element1[E] == this || element1[W] == this)) {
        findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
        findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    } else if (axialLine[1] && (element1[N] == this || element1[S] == this)) {
        findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
        findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);
    } else {
        printError("AGM::Point::calcRepresentationFormula_boundary1()");
    }

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    std::array<matrix_row, 2> row{};

    if (nx == E) {
        row[0][idx] = -mp * gfunc_x.green_function_t_ND(xp);
        row[0][idx + ptsnum] = gfunc_x.green_integral_ND('r');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_ND('r');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_ND('r');

        row[0][element[W]->getIdx()] = -UNITVALUE;
        row[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral_ND('c');
        row[0][element[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_ND('c');
        row[0][element[W]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t_ND('c');

        row[0][element1[W]->getIdx() + ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][element1[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][element1[W]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t_ND('l') + gfunc_x.green_function_ND(xm);

        row[0] += element1[W]->getDMatrixRow()[0] * -mp * gfunc_x.green_function_ND(xm);
    } else if (nx == W) {
        row[0][idx] = mp * gfunc_x.green_function_t_DN(xm);
        row[0][idx + ptsnum] = gfunc_x.green_integral_DN('l');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_DN('l');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_DN('l');

        row[0][element[E]->getIdx()] = -UNITVALUE;
        row[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral_DN('c');
        row[0][element[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_DN('c');
        row[0][element[E]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t_DN('c');

        row[0][element1[E]->getIdx() + ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][element1[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][element1[E]->getIdx() + 4 * ptsnum] = gfunc_x.green_integral_t_DN('r') - gfunc_x.green_function_DN(xp);

        row[0] += element1[E]->getDMatrixRow()[0] * mp * gfunc_x.green_function_DN(xp);
    } else if (ny == N) {
        row[1][idx] = -mp * gfunc_y.green_function_t_ND(yp);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_ND('r');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_ND('r');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_ND('r');

        row[1][element[S]->getIdx()] = -UNITVALUE;
        row[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral_ND('c');
        row[1][element[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_ND('c');
        row[1][element[S]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t_ND('c');

        row[1][element1[S]->getIdx() + ptsnum] = -gfunc_y.green_integral_ND('l');
        row[1][element1[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_ND('l');
        row[1][element1[S]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t_ND('l') + gfunc_y.green_function_ND(ym);

        row[1] += element1[S]->getDMatrixRow()[1] * -mp * gfunc_y.green_function_ND(ym);
    } else if (ny == S) {
        row[1][idx] = mp * gfunc_y.green_function_t_DN(ym);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_DN('l');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_DN('l');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_DN('l');

        row[1][element[N]->getIdx()] = -UNITVALUE;
        row[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral_DN('c');
        row[1][element[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_DN('c');
        row[1][element[N]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t_DN('c');

        row[1][element1[N]->getIdx() + ptsnum] = -gfunc_y.green_integral_DN('r');
        row[1][element1[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_DN('r');
        row[1][element1[N]->getIdx() + 4 * ptsnum] = gfunc_y.green_integral_t_DN('r') - gfunc_y.green_function_DN(yp);

        row[1] += element1[N]->getDMatrixRow()[1] * mp * gfunc_y.green_function_DN(yp);
    } else {
        printError("AGM::PointHeat::calcRepresentationFormula_boundary()");
    }

    for (int i = 0; i < 2; ++i) {
        if (!row[i].empty()) {
            while (row[i].back().idx >= 2 * ptsnum) {
                phiMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value;
                row[i].pop_back();
            }
        }
    }

    if (axialLine[0] && (element1[E] == this || element1[W] == this)) {
        matrixRow[1] = row[0];
    } else if (axialLine[1] && (element1[N] == this || element1[S] == this)) {
        matrixRow[1] = row[1];
    } else {
        printError("AGM::PointHeat::calcRepresentationFormula_boundary1()");
    }
}

void AGM::PointHeat::updateRb_cross() {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        rb[0] -= j.value * rhsParallelToX(j.idx);
        rb[1] -= j.value * rhsParallelToX(j.idx);
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * rhsParallelToY(j.idx);
        rb[1] += j.value * rhsParallelToY(j.idx);
    }
}

void AGM::PointHeat::updateRb_dirichlet() {
    rb[0] = value["bdv"];
//    for (const auto &i : phiMatrixRow[0]) {
//        rb[1] -= i.value * rhsParallelToX(i.idx);
//    }
//    for (const auto &i : phiMatrixRow[1]) {
//        rb[1] -= i.value * rhsParallelToY(i.idx);
//    }
}

void AGM::PointHeat::updateRb_neumann() {
    rb[0] = value["bdv"];
    rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToX(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * rhsParallelToY(j.idx - ptsnum);
        } else {
            printError("AGM::PointHeat::updateRb_neumann()");
        }
    }
    for (const auto &j : rhsMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToY(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * rhsParallelToX(j.idx - ptsnum);
        } else {
            printError("AGM::PointHeat::updateRb_neumann()");
        }
    }
    for (const auto &item : phiMatrixRow[0]) {
        rb[1] -= item.value * rhsParallelToX(item.idx);
    }
    for (const auto &item : phiMatrixRow[1]) {
        rb[1] -= item.value * rhsParallelToY(item.idx);
    }

}

void AGM::PointHeat::updateRb_interface() {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToX(j.idx);
            rb[1] -= j.value * rhsParallelToX(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * rhsParallelToY(j.idx - ptsnum);
            rb[1] -= j.value * rhsParallelToY(j.idx - ptsnum);
        }
    }
    for (const auto &j : rhsMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * rhsParallelToY(j.idx);
            rb[1] += j.value * rhsParallelToY(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * rhsParallelToX(j.idx - ptsnum);
            rb[1] += j.value * rhsParallelToX(j.idx - ptsnum);
        }
    }
}

void AGM::PointHeat::makeDifferentiation_cross() {
    makeDifferentiation_cross_RD();
}

void AGM::PointHeat::makeDifferentiation_cross_RD() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

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

    dMatrixRow[0][element[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_tau('r');

    dMatrixRow[0][element[W]->getIdx() + 4 * ptsnum] =
            gfunc_x.green_integral_ttau('l') + gfunc_x.green_function_tau(xm);
    dMatrixRow[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    dMatrixRow[0][element[E]->getIdx() + 4 * ptsnum] =
            gfunc_x.green_integral_ttau('r') - gfunc_x.green_function_tau(xp);

    dMatrixRow[1][element[S]->getIdx()] = mp * gfunc_y.green_function_ttau(ym);
    dMatrixRow[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_ttau(yp);

    dMatrixRow[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->getIdx() + 4 * ptsnum] =
            gfunc_y.green_integral_ttau('l') + gfunc_y.green_function_tau(ym);
    dMatrixRow[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    dMatrixRow[1][element[N]->getIdx() + 4 * ptsnum] =
            gfunc_y.green_integral_ttau('r') - gfunc_y.green_function_tau(yp);

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);
}

void AGM::PointHeat::makeDifferentiation_cross_E() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

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

    dMatrixRow[0][idx] = -HALFVALUE * gfunc_x.green_integral_tau('c') / delta;
    dMatrixRow[0][element[W]->getIdx()] = mp * gfunc_x.green_function_ttau(xm)
                                          - HALFVALUE * gfunc_x.green_integral_tau('l') / delta;
    dMatrixRow[0][element[E]->getIdx()] = -mp * gfunc_x.green_function_ttau(xp)
                                          - HALFVALUE * gfunc_x.green_integral_tau('r') / delta;

    dMatrixRow[0][element[W]->getIdx() + ptsnum] = gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->getIdx() + ptsnum] = gfunc_x.green_integral_tau('r');

    dMatrixRow[0][element[W]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->getIdx() + 2 * ptsnum] = gfunc_x.green_integral_tau('r');

    dMatrixRow[0][element[W]->getIdx() + 4 * ptsnum] =
            gfunc_x.green_integral_ttau('l') + gfunc_x.green_function_tau(xm);
    dMatrixRow[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    dMatrixRow[0][element[E]->getIdx() + 4 * ptsnum] =
            gfunc_x.green_integral_ttau('r') - gfunc_x.green_function_tau(xp);

    dMatrixRow[1][idx] = -HALFVALUE * gfunc_y.green_integral_tau('c') / delta;
    dMatrixRow[1][element[S]->getIdx()] = mp * gfunc_y.green_function_ttau(ym)
                                          - HALFVALUE * gfunc_y.green_integral_tau('l') / delta;
    dMatrixRow[1][element[N]->getIdx()] = -mp * gfunc_y.green_function_ttau(yp)
                                          - HALFVALUE * gfunc_y.green_integral_tau('r') / delta;

    dMatrixRow[1][element[S]->getIdx() + ptsnum] = -gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + ptsnum] = -gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->getIdx() + 2 * ptsnum] = gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->getIdx() + 4 * ptsnum] =
            gfunc_y.green_integral_ttau('l') + gfunc_y.green_function_tau(ym);
    dMatrixRow[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    dMatrixRow[1][element[N]->getIdx() + 4 * ptsnum] =
            gfunc_y.green_integral_ttau('r') - gfunc_y.green_function_tau(yp);

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);
}

void AGM::PointHeat::makeDifferentiation_boundary() {
    makeDifferentiation_boundary_RD();
}

void AGM::PointHeat::makeDifferentiation_boundary_E() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionLinear(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionLinear(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = CN * sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = CN * sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
                            int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] +=
                    mp0 * func->green_function_ttau(d) - HALFVALUE * func->green_integral_tau(c) / delta;
            dMatrixRow[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            dMatrixRow[i][pt->getIdx() + 4 * ptsnum] +=
                    func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d) -
                                                      HALFVALUE * func->green_integral_tau(c) / delta, i0,
                                            std::abs(mp0));
            dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs_t(ptl, ptr,
                                              func->green_integral_ttau(c) + sign1 * func->green_function_tau(d), i0);
        }
    };

    dMatrixRow[0][idx] = -HALFVALUE * gfunc_x.green_integral_tau('c') / delta;
    dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    dMatrixRow[1][idx] = -HALFVALUE * gfunc_y.green_integral_tau('c') / delta;
    dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');
}

void AGM::PointHeat::makeDifferentiation_boundary_RD() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                            double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            dMatrixRow[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            dMatrixRow[i][pt->getIdx() + 4 * ptsnum] +=
                    func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs_t(ptl, ptr,
                                              func->green_integral_ttau(c) + sign1 * func->green_function_tau(d), i0);
        }
    };

    dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');
}

void AGM::PointHeat::makeDifferentiation_boundary(const std::string &string) {
    if (string == "x") {
        if (this->getAxialLine('x')) makeDifferentiation_boundary_RD_N(string);
        else makeDifferentiation_boundary_RD(string);
    } else if (string == "y") {
        if (this->getAxialLine('y')) makeDifferentiation_boundary_RD_N(string);
        else makeDifferentiation_boundary_RD(string);
    } else {
        if (this->getAxialLine('x')) makeDifferentiation_boundary_RD_N(std::string("x"));
        else makeDifferentiation_boundary_RD(std::string("x"));

        if (this->getAxialLine('y')) makeDifferentiation_boundary_RD_N(std::string("y"));
        else makeDifferentiation_boundary_RD(std::string("y"));
    }
}

void AGM::PointHeat::makeDifferentiation_boundary_E(const std::string &string) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};
    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
                            int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] +=
                    mp0 * func->green_function_ttau(d) - HALFVALUE * func->green_integral_tau(c) / delta;
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->getIdx() + 4 * ptsnum] +=
                    func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d) -
                                               HALFVALUE * func->green_integral_tau(c) / delta, i0, std::abs(mp0));
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_ttau(c) + sign1 * func->green_function_tau(d),
                                       i0);
        }
    };

    row[0][idx] = -HALFVALUE * gfunc_x.green_integral_tau('c') / delta;
    row[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    row[1][idx] = -HALFVALUE * gfunc_y.green_integral_tau('c') / delta;
    row[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');

    if (string == "x") dMatrixRow[0] = row[0];
    else if (string == "y") dMatrixRow[1] = row[1];
    else if (string == "xy") {
        dMatrixRow[0] = row[0];
        dMatrixRow[1] = row[1];
    }
}

void AGM::PointHeat::makeDifferentiation_boundary_E_N(const std::string &string) {
    double xm{}, xb{}, xp{}, ym{}, yb{}, yp{};
    EWNS nx{NW}, ny{WS};
    double xx{}, yy{};

    auto findPts = [&](EWNS ewns, EWNS ewns1, EWNS ewns2, int i, double *tm, double *tb, double *tp, EWNS *nn) -> void {
        if (element1[ewns] == this || element1[ewns1] == this) {
            if (element1[ewns] == this) {
                *tm = element1[ewns1] ? element1[ewns1]->getCoordinate()[i] : element1[ewns2]->getCoordinate()[i];
                *tb = element[ewns1] ? element[ewns1]->getCoordinate()[i] : element[ewns2]->getCoordinate()[i];
                *tp = coordinate[i];
                *nn = ewns;
            }
        } else {
            if (i) {
                ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
                yb = getCoordinate()[1];
                yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
                ny = N;
            } else {
                xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
                xb = getCoordinate()[0];
                xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
                nx = E;
            }
        }
    };

    findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
    findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
    findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - HALFVALUE * func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - HALFVALUE * func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE + HALFVALUE * func->green_integral_ND(c) / delta;
            } else {
                row[i][pt->getIdx()] += UNITVALUE - HALFVALUE * func->green_integral_ND(c) / delta;
            }

            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_ND(c);

        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr,
                                         -UNITVALUE - HALFVALUE * func->green_integral_ND(c) / delta, i0,
                                         std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE - HALFVALUE * func->green_integral_ND(c) / delta,
                                         i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_ND(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_ND(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_ND(c), i0);
        }
    };

    auto assignMatrix_DN = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE + HALFVALUE * func->green_integral_DN(c) / delta;
            } else {
                row[i][pt->getIdx()] += UNITVALUE - HALFVALUE * func->green_integral_DN(c) / delta;
            }
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_DN(c);

        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr,
                                         -UNITVALUE - HALFVALUE * func->green_integral_DN('c') / delta, i0,
                                         std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE - HALFVALUE * func->green_integral_DN(c) / delta,
                                         i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_DN(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_DN(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_DN(c), i0);
        }
    };

    if (nx == E) {
        xx = mp * gfunc_x.green_function_DN(xp);
        row[0][idx] = -HALFVALUE * gfunc_x.green_integral_DN('r') / delta;
        row[0][idx + ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_DN('r') - gfunc_x.green_function_DN(xp);

        assignMatrix_DN(element[W], element[WS], element[WN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_DN(element1[W], element1[WS], element1[WN], mp, &gfunc_x, xm, 0, 1, 'l');
    } else if (nx == W) {
        xx = -mp * gfunc_x.green_function_ND(xm);
        row[0][idx] = -HALFVALUE * gfunc_x.green_integral_ND('l') / delta;
        row[0][idx + ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_ND('l') + gfunc_x.green_function_ND(xm);

        assignMatrix_ND(element[E], element[ES], element[EN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_ND(element1[E], element1[ES], element1[EN], mp, &gfunc_x, xp, 0, 1, 'r');
    } else {
        printError("AGM::PointHeat::makeDifferentiation_boundary()");
    }

    if (ny == N) {
        yy = mp * gfunc_y.green_function_DN(yp);
        row[1][idx] = -HALFVALUE * gfunc_y.green_integral_DN('r') / delta;
        row[1][idx + ptsnum] = -gfunc_y.green_integral_DN('r');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_DN('r');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_DN('r') - gfunc_y.green_function_DN(yp);

        assignMatrix_DN(element[S], element[SW], element[SE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_DN(element1[S], element1[SW], element1[SE], mp, &gfunc_y, ym, 1, 0, 'l');
    } else if (ny == S) {
        yy = -mp * gfunc_y.green_function_ND(ym);
        row[1][idx] = -HALFVALUE * gfunc_y.green_integral_ND('l') / delta;
        row[1][idx + ptsnum] = -gfunc_y.green_integral_ND('l');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_ND('l');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_ND('l') + gfunc_y.green_function_ND(ym);

        assignMatrix_ND(element[N], element[NW], element[NE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_ND(element1[N], element1[NW], element1[NE], mp, &gfunc_y, yp, 1, 0, 'r');
    } else {
        printError("AGM::PointHeat::makeDifferentiation_boundary()");
    }

    double c{};
    for (int i = 0; i < 2; ++i) {
        c = i ? -yy : -xx;
        for (auto &item : row[i]) {
            item.value /= c;
        }
    }

    if (string == "x") dMatrixRow[0] = row[0];
    else if (string == "y") dMatrixRow[1] = row[1];
    else if (string == "xy") {
        dMatrixRow[0] = row[0];
        dMatrixRow[1] = row[1];
    }
}

void AGM::PointHeat::makeDifferentiation_boundary_RD(const std::string &string) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                            double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            row[i][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->getIdx() + 4 * ptsnum] +=
                    func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            row[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            row[i] += approximateRhs_t(ptl, ptr,
                                       func->green_integral_ttau(c) + sign1 * func->green_function_tau(d), i0);
        }
    };

    row[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    row[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');

    if (string == "x") {
        dMatrixRow[0] = row[0];
    } else if (string == "y") {
        dMatrixRow[1] = row[1];
    } else if (string == "xy") {
        dMatrixRow[0] = row[0];
        dMatrixRow[1] = row[1];
    }
}

void AGM::PointHeat::makeDifferentiation_boundary_RD_N(const std::string &string) {
    double xm{}, xb{}, xp{}, ym{}, yb{}, yp{};
    EWNS nx{NW}, ny{WS};

    auto findPts = [&](EWNS ewns, EWNS ewns1, EWNS ewns2, int i, double *tm, double *tb, double *tp, EWNS *nn) -> void {
        if (element1[ewns] == this || element1[ewns1] == this) {
            if (element1[ewns] == this) {
                *tm = element1[ewns1] ? element1[ewns1]->getCoordinate()[i] : element1[ewns2]->getCoordinate()[i];
                *tb = element[ewns1] ? element[ewns1]->getCoordinate()[i] : element[ewns2]->getCoordinate()[i];
                *tp = coordinate[i];
                *nn = ewns;
            }
        } else {
            if (i) {
                ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
                yb = getCoordinate()[1];
                yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
                ny = N;
            } else {
                xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
                xb = getCoordinate()[0];
                xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
                nx = E;
            }
        }
    };

    findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
    findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
    findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    std::array<matrix_row, 2> row{};

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                               double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE;
            } else {
                row[i][pt->getIdx()] += -mp * func->green_function_t_ND(d);
            }
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_ND(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, -mp * func->green_function_t_ND(d), i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_ND(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_ND(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_ND(c), i0);
        }
    };

    auto assignMatrix_DN = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                               double d, int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->getIdx()] -= UNITVALUE;
            } else {
                row[i][pt->getIdx()] += mp * func->green_function_t_DN(d);
            }
            row[i][pt->getIdx() + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->getIdx() + 4 * ptsnum] += func->green_integral_t_DN(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, mp * func->green_function_t_DN(d), i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_DN(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_DN(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_DN(c), i0);
        }
    };

    if (nx == E) {
        row[0][idx] = mp * gfunc_x.green_function_DN(xp);
        row[0][idx + ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_DN('r');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_DN('r') - gfunc_x.green_function_DN(xp);

        assignMatrix_DN(element[W], element[WS], element[WN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_DN(element1[W], element1[WS], element1[WN], mp, &gfunc_x, xm, 0, 1, 'l');
    } else if (nx == W) {
        row[0][idx] = -mp * gfunc_x.green_function_ND(xm);
        row[0][idx + ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 2 * ptsnum] = gfunc_x.green_integral_ND('l');
        row[0][idx + 4 * ptsnum] = gfunc_x.green_integral_t_ND('l') + gfunc_x.green_function_ND(xm);

        assignMatrix_ND(element[E], element[ES], element[EN], mp, &gfunc_x, xb, 0, 1, 'c');
        assignMatrix_ND(element1[E], element1[ES], element1[EN], mp, &gfunc_x, xp, 0, 1, 'r');
    } else {
        printError("AGM::PointHeat::makeDifferentiation_boundary1()");
    }

    if (ny == N) {
        row[1][idx] = mp * gfunc_y.green_function_DN(yp);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_DN('r');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_DN('r');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_DN('r') - gfunc_y.green_function_DN(yp);

        assignMatrix_DN(element[S], element[SW], element[SE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_DN(element1[S], element1[SW], element1[SE], mp, &gfunc_y, ym, 1, 0, 'l');
    } else if (ny == S) {
        row[1][idx] = -mp * gfunc_y.green_function_ND(ym);
        row[1][idx + ptsnum] = -gfunc_y.green_integral_ND('l');
        row[1][idx + 2 * ptsnum] = gfunc_y.green_integral_ND('l');
        row[1][idx + 4 * ptsnum] = gfunc_y.green_integral_t_ND('l') + gfunc_y.green_function_ND(ym);

        assignMatrix_ND(element[N], element[NW], element[NE], mp, &gfunc_y, yb, 1, 0, 'c');
        assignMatrix_ND(element1[N], element1[NW], element1[NE], mp, &gfunc_y, yp, 1, 0, 'r');
    } else {
        printError("AGM::Point::makeDifferentiation_boundary1()");
    }

    for (int i = 0; i < 2; ++i) {
        dv[i] = -row[i][idx];
        row[i][idx] = ZEROVALUE;
    }

    if (string == "x") {
        dv[1] = ZEROVALUE;
        dMatrixRow[0] = row[0];
    }
    else if (string == "y") {
        dv[0] = ZEROVALUE;
        dMatrixRow[1] = row[1];
    }
    else if (string == "xy") {
        dMatrixRow[0] = row[0];
        dMatrixRow[1] = row[1];
    }
}

void AGM::PointHeat::makeDifferentiation_interface() {
    makeDifferentiation_interface_RD();
}

void AGM::PointHeat::makeDifferentiation_interface_E() {
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

    auto gfunc_x = Greenfunction(xm, xb, xp, CN * mpw, CN * mpe);
    auto gfunc_y = Greenfunction(ym, yb, yp, CN * mps, CN * mpn);

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
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm) - func.green_integral('L') / delta;
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp) - func.green_integral('R') / delta;

        mRow[ptl->getIdx() + ptsnum] = CN * sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = CN * sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        dMatrixRow[(i + 1) % 2][ptl->getIdx() + 3 * ptsnum] += func.green_integral('L') * coefficient;
        dMatrixRow[(i + 1) % 2][ptr->getIdx() + 3 * ptsnum] += func.green_integral('R') * coefficient;

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

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
                            int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] +=
                    mp0 * func->green_function_ttau(d) - func->green_integral_tau(c) / delta;
            if (!checkInterface(this)) dMatrixRow[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral_tau(c);
            else dMatrixRow[i][pt->getIdx() + ptsnum] += CN * sign * func->green_integral_tau(C);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr,
                                            mp0 * func->green_function_ttau(d) -
                                            func->green_integral_tau(c) / delta, i0, std::abs(mp0));
            if (!checkInterface(this))
                dMatrixRow[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral_tau(c), i0);
            else dMatrixRow[i] += approximatePhi(ptl, ptr, CN * sign * func->green_integral_tau(C), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
        }
    };

    dMatrixRow[0][idx] += -gfunc_x.green_integral_tau('c') / delta;
    if (!checkInterface(this)) dMatrixRow[0][idx + ptsnum] += CN * gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] += gfunc_x.green_integral_tau('c');

    assignMatrix(element[E], element[ES], element[EN], -CN * mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], CN * mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    dMatrixRow[1][idx] += -gfunc_y.green_integral_tau('c') / delta;
    if (!checkInterface(this)) dMatrixRow[1][idx + ptsnum] += -CN * gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] += gfunc_y.green_integral_tau('c');

    assignMatrix(element[N], element[NW], element[NE], -CN * mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], CN * mps, &gfunc_y, ym, 1, 0, 'l', 'L');
}

void AGM::PointHeat::makeDifferentiation_interface_RD() {
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

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mpw, mpe, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mps, mpn, UNITVALUE / delta);

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
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, UNITVALUE / delta);
        auto mRow = matrix_row();
        mRow[ptl->getIdx()] = d * func.green_function_t(tm);
        mRow[ptr->getIdx()] = -d * func.green_function_t(tp);

        mRow[ptl->getIdx() + ptsnum] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + ptsnum] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptl, ptr);

        mRow[ptl->getIdx() + 3 * ptsnum] = func.green_integral('L');
        mRow[ptr->getIdx() + 3 * ptsnum] = func.green_integral('R');

        mRow[ptl->getIdx() + 5 * ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + 5 * ptsnum] = func.green_integral_t('R');

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

    auto approximateRhs_t = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx() + 4 * ptsnum] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + 4 * ptsnum] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionReactionDiffusion *func,
                            double d, int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;
        if (pt) {
            dMatrixRow[i][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            if (!checkInterface(this)) dMatrixRow[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(c);
            else dMatrixRow[i][pt->getIdx() + ptsnum] += sign * func->green_integral_tau(C);
            dMatrixRow[i][pt->getIdx() + 2 * ptsnum] += func->green_integral_tau(c);
            dMatrixRow[i][pt->getIdx() + 4 * ptsnum] +=
                    func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            if (!checkInterface(this))
                dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            else dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(C), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs_t(ptl, ptr,
                                              func->green_integral_ttau(c) + sign1 * func->green_function_tau(d), i0);
        }
    };

    if (!checkInterface(this)) dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[E], element[ES], element[EN], -mpe, &gfunc_x, xp, 0, 1, 'r', 'R');
    assignMatrix(element[W], element[WS], element[WN], mpw, &gfunc_x, xm, 0, 1, 'l', 'L');

    if (!checkInterface(this)) dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    dMatrixRow[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;

    assignMatrix(element[N], element[NW], element[NE], -mpn, &gfunc_y, yp, 1, 0, 'r', 'R');
    assignMatrix(element[S], element[SW], element[SE], mps, &gfunc_y, ym, 1, 0, 'l', 'L');
}

void AGM::PointHeat::calcDifferentiation() {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * ptsH->at(item.idx)["sol"];
            std::cout << "[" << item.idx << "] = " << item.value * ptsH->at(item.idx)["sol"] << "\n";
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * ptsH->at(item.idx - ptsnum)["phi"];
            std::cout << "[" << item.idx << "] = " << item.value * ptsH->at(item.idx - ptsnum)["phi"] << "\n";
        } else if (item.idx < 3 * ptsnum) {
            dx += item.value * rhsParallelToX(item.idx - 2 * ptsnum);
            std::cout << "[" << item.idx << "] = " << item.value * rhsParallelToX(item.idx - 2 * ptsnum) << "\n";
        } else if (item.idx < 4 * ptsnum) {
            dx += item.value * rhsParallelToY(item.idx - 3 * ptsnum);
            std::cout << "[" << item.idx << "] = " << item.value * rhsParallelToY(item.idx - 3 * ptsnum) << "\n";
        }
    });
    value["dx"] = dx;
    exit(1);
    value["dxx"] =
            -(CN * value["phi"] - HALFVALUE * value["sol"] / delta + HALFVALUE * rhsParallelToX(idx)) / (CN * mp);
    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * ptsH->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * ptsH->at(item.idx % ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dy += item.value * rhsParallelToY(item.idx - 2 * ptsnum);
        } else if (item.idx < 4 * ptsnum) {
            dy += item.value * rhsParallelToX(item.idx - 3 * ptsnum);
        }
    });
    value["dy"] = dy;
    value["dyy"] =
            -(-CN * value["phi"] - HALFVALUE * value["sol"] / delta + HALFVALUE * rhsParallelToY(idx)) / (CN * mp);
}

void AGM::PointHeat::calcDifferentiation(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * ptsH->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * ptsH->at(item.idx % ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dx += item.value * f(item.idx % ptsnum);
        } else if (item.idx < 4 * ptsnum) {
            dx += item.value * g(item.idx % ptsnum);
        }
    });
    value["dx"] = dx;
    value["dxx"] = -(value["phi"] + HALFVALUE * f(idx)) / mp;
    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * ptsH->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * ptsH->at(item.idx % ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dy += item.value * g(item.idx % ptsnum);
        } else if (item.idx < 4 * ptsnum) {
            dy += item.value * f(item.idx % ptsnum);
        }
    });
    value["dy"] = dy;
    value["dyy"] = -(-value["phi"] + HALFVALUE * g(idx)) / mp;
}

void AGM::PointHeat::calcDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                         const std::function<double(int)> &gx, const std::function<double(int)> &gy) {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * ptsH->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * ptsH->at(item.idx % ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dx += item.value * fx(item.idx % ptsnum);
        } else if (item.idx < 4 * ptsnum) {
            dx += item.value * fy(item.idx % ptsnum);
        } else if (item.idx < 5 * ptsnum) {
            dx += item.value * gx(item.idx % ptsnum);
        } else {
            dx += item.value * gy(item.idx % ptsnum);
        }
    });
    value["dx"] = dx;
    if (!iszero(dv[0])) value["dx"] /= dv[0];
    value["dxx"] = -(value["phi"] + HALFVALUE * fx(idx)) / mp;
    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * ptsH->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * ptsH->at(item.idx % ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dy += item.value * fy(item.idx % ptsnum);
        } else if (item.idx < 4 * ptsnum) {
            dy += item.value * fx(item.idx % ptsnum);
        } else if (item.idx < 5 * ptsnum) {
            dy += item.value * gy(item.idx % ptsnum);
        } else {
            dy += item.value * gx(item.idx % ptsnum);
        }
    });
    value["dy"] = dy;
    if (!iszero(dv[1])) value["dy"] /= dv[1];
    value["dyy"] = -(-value["phi"] + HALFVALUE * fy(idx)) / mp;
}

void AGM::PointHeat::calcDifferentiationWithPressure(std::vector<Point> *pressure, char i) {
    double dx{}, dy{};
    if (i == 'u') {
        std::for_each(pdMatrixRow[0].begin(), pdMatrixRow[0].end(), [&](matrix_element &item) {
            dx += item.value * pressure->at(item.idx)["sol"];
        });
        value["dx"] += dx;
        value["dxx"] += pressure->at(idx)["dx"] / (CN * mp);
    } else if (i == 'v') {
        std::for_each(pdMatrixRow[1].begin(), pdMatrixRow[1].end(), [&](matrix_element &item) {
            dy += item.value * pressure->at(item.idx)["sol"];
        });
        value["dy"] += dy;
        value["dyy"] += pressure->at(idx)["dy"] / (CN * mp);
    }
}

void AGM::PointHeat::makePressureTerm_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    pMatrixRow[0][element[W]->getIdx()] = gfunc_x.green_integral_t('l');
    pMatrixRow[0][idx] = gfunc_x.green_integral_t('c');
    pMatrixRow[0][element[E]->getIdx()] = gfunc_x.green_integral_t('r');

    pMatrixRow[1][element[S]->getIdx()] = gfunc_y.green_integral_t('l');
    pMatrixRow[1][idx] = gfunc_y.green_integral_t('c');
    pMatrixRow[1][element[N]->getIdx()] = gfunc_y.green_integral_t('r');
}

void AGM::PointHeat::makePressureTerm_neumann() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionLinear(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionLinear(ym, yb, yp, mp, mp);

    auto approximateP = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx()] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx()] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, GreenfunctionLinear *func, int i, int i0,
                            char c) -> void {
        if (pt) {
            pMatrixRow[i][pt->getIdx()] += func->green_integral_ttau(c);
        } else {
            pMatrixRow[i] += approximateP(ptl, ptr, func->green_integral_ttau(c), i0);
        }
    };

    pMatrixRow[0][idx] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[E], element[ES], element[EN], &gfunc_x, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], &gfunc_x, 0, 1, 'l');

    pMatrixRow[1][idx] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[N], element[NW], element[NE], &gfunc_y, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], &gfunc_y, 1, 0, 'l');

    for (int i = 0; i < 2; ++i) {
        for (auto &item : pMatrixRow[i]) {
            item.value *= normal[i];
        }
    }
}

void AGM::PointHeat::makePressureTerm_interface() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta);
    auto gfunc_y = GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto func = GreenfunctionReactionDiffusion(tm, tb, tp, d, d, HALFVALUE * BDF / delta);
        auto mRow = matrix_row();
//        mRow[ptl->getIdx() + ptsnum] = func.green_integral_t('L');
//        mRow[ptr->getIdx() + ptsnum] = func.green_integral_t('R');

        mRow[ptl->getIdx() + ptsnum] = func.green_integral_t('l') + func.green_integral_t('c') * (tp - tb) / (tp - tm);
        mRow[ptr->getIdx() + ptsnum] = func.green_integral_t('r') + func.green_integral_t('c') * (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto approximatePhi = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx()] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx()] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d, int i, int i0,
                            char c) -> void {
        if (pt) {
            pMatrixRow[i][pt->getIdx()] += func->green_integral_t(c);
        } else {
            pMatrixRow[i0] += approximateSol(ptl, ptr, mp0 * func->green_function_t(d), i0, std::abs(mp0));
//            pMatrixRow[i0] += approximateSol(ptl, ptr, UNITVALUE, i0, std::abs(mp0));
            pMatrixRow[i] += approximatePhi(ptl, ptr, func->green_integral_t(c), i0);
        }
    };

    pMatrixRow[0][idx] += gfunc_x.green_integral_t('c');
    assignMatrix(element[E], element[ES], element[EN], -mp, &gfunc_x, xp, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], mp, &gfunc_x, xm, 0, 1, 'l');

    pMatrixRow[1][idx] += gfunc_y.green_integral_t('c');
    assignMatrix(element[N], element[NW], element[NE], -mp, &gfunc_y, yp, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], mp, &gfunc_y, ym, 1, 0, 'l');
}
