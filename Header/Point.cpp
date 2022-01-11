//
// Created by NIMS-JUNHONG on 2020/12/23.
//

#include "Point.h"
#include "function2D.h"

AGM::Element::Element() = default;

AGM::Element::~Element() {
    for (int i = 0; i < 12; ++i) {
        element[i] = nullptr;
    }
}

const std::array<AGM::Point *, 12> &AGM::Element::getElement() const {
    return element;
}

void AGM::Element::setElement(const std::array<Point *, 12> &array) {
    Element::element = array;
}

AGM::Point *&AGM::Element::operator[](int i) {
    return element[i];
}

AGM::Point *&AGM::Element::operator[](AGM::EWNS ewns) {
    return element[ewns];
}

AGM::Element &AGM::Element::operator=(const AGM::Element &rhs) {
    if (this != &rhs) {
        element = rhs.element;
    }
    return *this;
}

int AGM::Point::ptsnum;
std::vector<AGM::Point> *AGM::Point::pts;

AGM::Point::Point() = default;

AGM::Point::Point(const AGM::Point &src) {
    idx = src.idx;
    coordinate = src.coordinate;
    normal = src.normal;
    mp = src.mp;
    condition = src.condition;
    element = src.element;
    value = src.value;
}

AGM::Point::Point(const AGM::Coordinate &coordinate) : coordinate(coordinate) {}

AGM::Point::Point(const AGM::Coordinate &coordinate, double mp) : coordinate(coordinate), mp(mp) {}

AGM::Point::Point(int idx) : idx(idx) {}

AGM::Point::Point(int idx, const AGM::Coordinate &coordinate) : idx(idx), coordinate(coordinate) {}

AGM::Point::Point(int idx, const AGM::Coordinate &coordinate, double mp) : idx(idx), coordinate(coordinate), mp(mp) {}

AGM::Point::~Point() = default;

int AGM::Point::getIdx() const {
    return idx;
}

void AGM::Point::setIdx(int i) {
    Point::idx = i;
}

int AGM::Point::getPtsnum() {
    return ptsnum;
}

void AGM::Point::setPtsnum(int i) {
    Point::ptsnum = i;
}

const AGM::Coordinate &AGM::Point::getCoordinate() const {
    return coordinate;
}

void AGM::Point::setCoordinate(const AGM::Coordinate &src) {
    Point::coordinate = src;
}

const AGM::Element &AGM::Point::getElement() const {
    return element;
}

void AGM::Point::setElement(const AGM::Element &src) {
    Point::element = src;
}

const AGM::Element &AGM::Point::getElement1() const {
    return element1;
}

void AGM::Point::setElement1(const AGM::Element &src) {
    Point::element1 = src;
}

double AGM::Point::getMp() const {
    return mp;
}

void AGM::Point::setMp(double d) {
    Point::mp = d;
}

char AGM::Point::getCondition() const {
    return condition;
}

void AGM::Point::setCondition(char i) {
    Point::condition = i;
}

void AGM::Point::setNormal(const AGM::Coordinate &n) {
    Point::normal = n;
}

const AGM::Coordinate &AGM::Point::getNormal() const {
    return normal;
}

const AGM::Value &AGM::Point::getValue() const {
    return value;
}

AGM::Value &AGM::Point::getValue() {
    return value;
}

void AGM::Point::setValue(const AGM::Value &v) {
    Point::value = v;
}

const std::array<AGM::matrix_row, 2> &AGM::Point::getMatrixRow() const {
    return matrixRow;
}

void AGM::Point::setMatrixRow(const std::array<matrix_row, 2> &row) {
    Point::matrixRow = row;
}

const std::array<AGM::matrix_row, 2> &AGM::Point::getRhsMatrixRow() const {
    return rhsMatrixRow;
}

void AGM::Point::setRhsMatrixRow(const std::array<matrix_row, 2> &row) {
    Point::rhsMatrixRow = row;
}

const std::array<AGM::matrix_row, 2> &AGM::Point::getPhiMatrixRow() const {
    return phiMatrixRow;
}

void AGM::Point::setPhiMatrixRow(const std::array<matrix_row, 2> &row) {
    Point::phiMatrixRow = row;
}

const std::array<AGM::matrix_row, 2> &AGM::Point::getDMatrixRow() const {
    return dMatrixRow;
}

void AGM::Point::setDMatrixRow(const std::array<matrix_row, 2> &row) {
    Point::dMatrixRow = row;
}

const std::array<AGM::matrix_row, 2> &AGM::Point::getPMatrixRow() const {
    return pMatrixRow;
}

void AGM::Point::setPMatrixRow(const std::array<matrix_row, 2> &row) {
    Point::pMatrixRow = row;
}

const std::array<AGM::matrix_row, 2> &AGM::Point::getPdMatrixRow() const {
    return pdMatrixRow;
}

void AGM::Point::setPdMatrixRow(const std::array<matrix_row, 2> &Row) {
    Point::pdMatrixRow = Row;
}

const std::array<double, 2> &AGM::Point::getRb() const {
    return rb;
}

void AGM::Point::setRb(const std::array<double, 2> &array) {
    Point::rb = array;
}

const std::array<AGM::AxialLine *, 2> &AGM::Point::getAxialLine() const {
    return axialLine;
}

AGM::AxialLine *&AGM::Point::getAxialLine(char i) {
    if (i == 'x') return axialLine[0];
    if (i == 'y') return axialLine[1];
    printError("AGM::AxialLine *&AGM::Point::getAxialLine", "index (which = %c) is wrong", i);
    return axialLine[0];
}

void AGM::Point::setAxialLine(const std::array<AxialLine *, 2> &array) {
    Point::axialLine = array;
}

void AGM::Point::setAxialLine(AGM::AxialLine *line, char i) {
    if (i == 'x') axialLine[0] = line;
    if (i == 'y') axialLine[1] = line;
}

std::vector<AGM::Point> *AGM::Point::getPts() {
    return pts;
}

void AGM::Point::setPts(std::vector<Point> *pVector) {
    Point::pts = pVector;
}

void AGM::Point::findElement(std::vector<AxialLine> *xline, std::vector<AxialLine> *yline) {
    if (condition == 'C') return;

    AxialLine *alin{};
    Point *ptl{}, *ptr{};

    auto isleft = [](AxialLine *line, double pt, double qt) -> bool { return (qt < (*line)[0]) && ((*line)[0] < pt); };
    auto isright = [](AxialLine *line, double pt, double qt) -> bool { return (pt < (*line)[0]) && ((*line)[0] < qt); };
    auto iscontain = [](AxialLine *line, double pt) -> bool {
        return ((*line)[1] - 1.0E-10 < pt) && (pt < ((*line)[2] + 1.0E-10));
    };
    auto isleftpt = [](double x, double y, double z) -> bool { return (z < x) && (x < (y + 1.0E-10)); };
    auto isrightpt = [](double x, double y, double z) -> bool { return ((y - 1.0E-10) < x) && (x < z); };

    auto assignElement = [&](const EWNS &ewns, double qt, std::vector<AxialLine> *line, auto func, double qtl,
                             double qtr, int idx0, int idx1, const EWNS &ewns0,
                             const EWNS &ewns1, double n) -> void {
        if (element[ewns] == this && n > 1.0E-10) return;
        alin = nullptr;
//        double lineDist{std::numeric_limits<double>::max()};
//        double sign = ewns == E || ewns == N ? UNITVALUE : -UNITVALUE;
        if ((!element[ewns]) || element[ewns] == this) {
            for (auto &i : *line) {
                if (axialLine[idx1] == &i) {
                    continue;
                }
//                if ((sign * (i[0] - (*axialLine[idx1])[0]) > NEARZERO) && (sign * (i[0] - (*axialLine[idx1])[0]) < lineDist)) {
//                    lineDist = fabs(i[0] - (*axialLine[idx1])[0]);
//                }
                if (func(&i, coordinate[idx0], qt) && iscontain(&i, coordinate[idx1])) {
                    alin = &i;
                    qt = (*alin)[0];
                }
            }
//            if ((condition == 'D' || condition == 'N') && fabs(lineDist - sign * (qt - (*axialLine[idx1])[0])) > NEARZERO) {
//                alin = nullptr;
//            }
            if (alin) {
                for (auto &i : *alin) {
                    if (isleftpt((*i)[idx1], coordinate[idx1], qtl)) {
                        ptl = i;
                        qtl = (*ptl)[idx1];
                    }
                    if (isrightpt((*i)[idx1], coordinate[idx1], qtr)) {
                        ptr = i;
                        qtr = (*ptr)[idx1];
                    }
                }
            }
            if (ptr && ptl && ptr == ptl) {
                element[ewns] = ptr;
                ptr = ptl = nullptr;
            }
            if (ptr) {
                element[ewns0] = ptr;
                element[ewns] = nullptr;
            }
            if (ptl) {
                element[ewns1] = ptl;
                element[ewns] = nullptr;
            }
        }
    };

    assignElement(E, std::numeric_limits<double>::max(), yline, isright, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 0, 1, EN, ES, normal[0]);
    assignElement(W, -std::numeric_limits<double>::max(), yline, isleft, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 0, 1, WN, WS, -normal[0]);

    if (axialLine[0] && (!element[E] || !element[W])) {
        axialLine[0] = nullptr;
        element[E] = nullptr;
        element[EN] = nullptr;
        element[ES] = nullptr;
        element[W] = nullptr;
        element[WN] = nullptr;
        element[WS] = nullptr;
        element1[E] = nullptr;
        element1[EN] = nullptr;
        element1[ES] = nullptr;
        element1[W] = nullptr;
        element1[WN] = nullptr;
        element1[WS] = nullptr;
        assignElement(E, std::numeric_limits<double>::max(), yline, isright, -std::numeric_limits<double>::max(),
                      std::numeric_limits<double>::max(), 0, 1, EN, ES, normal[0]);
        assignElement(W, -std::numeric_limits<double>::max(), yline, isleft, -std::numeric_limits<double>::max(),
                      std::numeric_limits<double>::max(), 0, 1, WN, WS, -normal[0]);
    }

    assignElement(N, std::numeric_limits<double>::max(), xline, isright, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 1, 0, NE, NW, normal[1]);
    assignElement(S, -std::numeric_limits<double>::max(), xline, isleft, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 1, 0, SE, SW, -normal[1]);

    if (axialLine[1] && (!element[N] || !element[S])) {
        axialLine[1] = nullptr;
        element[N] = nullptr;
        element[NE] = nullptr;
        element[NW] = nullptr;
        element[S] = nullptr;
        element[SE] = nullptr;
        element[SW] = nullptr;
        element1[N] = nullptr;
        element1[NE] = nullptr;
        element1[NW] = nullptr;
        element1[S] = nullptr;
        element1[SE] = nullptr;
        element1[SW] = nullptr;
        assignElement(N, std::numeric_limits<double>::max(), xline, isright, -std::numeric_limits<double>::max(),
                      std::numeric_limits<double>::max(), 1, 0, NE, NW, normal[1]);
        assignElement(S, -std::numeric_limits<double>::max(), xline, isleft, -std::numeric_limits<double>::max(),
                      std::numeric_limits<double>::max(), 1, 0, SE, SW, -normal[1]);
    }

    auto assignNoElement = [&](int i, EWNS ewns, EWNS ewns0, EWNS ewns1, EWNS ewns2) -> void {
        if (axialLine[i] && !(element[ewns] || element[ewns0] || element[ewns1]) && (element[ewns2] == this)) {
            element[ewns] = this;
            if (i) {
                setNormal(Coordinate(ZEROVALUE, normal[1] / std::fabs(normal[1])));
            } else {
                setNormal(Coordinate(normal[0] / std::fabs(normal[0]), ZEROVALUE));
            }
        }
    };
    assignNoElement(0, N, NE, NW, S);
    assignNoElement(0, S, SE, SW, N);
    assignNoElement(1, E, EN, ES, W);
    assignNoElement(1, W, WN, WS, E);
}

void AGM::Point::findElement1(std::vector<AxialLine> *xline, std::vector<AxialLine> *yline) {
    if (condition == 'C') return;

    AxialLine *alin{}, *alin1;
    Point *ptl{}, *ptr{};

    auto isleft = [](AxialLine *line, double pt, double qt) -> bool { return (qt < (*line)[0]) && ((*line)[0] < pt); };
    auto isright = [](AxialLine *line, double pt, double qt) -> bool { return (pt < (*line)[0]) && ((*line)[0] < qt); };
    auto iscontain = [](AxialLine *line, double pt) -> bool {
        return ((*line)[1] - 1.0E-10 < pt) && (pt < ((*line)[2] + 1.0E-10));
    };
    auto isleftpt = [](double x, double y, double z) -> bool { return (z < x) && (x < (y + 1.0E-10)); };
    auto isrightpt = [](double x, double y, double z) -> bool { return ((y - 1.0E-10) < x) && (x < z); };

    auto assignElement = [&](const EWNS &ewns, double qt, std::vector<AxialLine> *line, auto func, double qtl,
                             double qtr, int idx0, int idx1, const EWNS &ewns0,
                             const EWNS &ewns1, double n) -> void {
        if (element1[ewns] == this && n > 1.0E-10) {
            return;
        }
        double qt1 = qt;
        alin = nullptr;
        alin1 = nullptr;

        if ((!element1[ewns]) || element1[ewns] == this) {
            for (auto &i : *line) {
                if (axialLine[idx1] == &i) continue;
                if (func(&i, coordinate[idx0], qt) && iscontain(&i, coordinate[idx1])) {
                    alin = &i;
                    qt = (*alin)[0];
                }
            }
            if (alin) {
                for (auto &i : *line) {
                    if (axialLine[idx1] == &i) continue;
                    if (func(&i, (*alin)[0], qt1) && iscontain(&i, coordinate[idx1])) {
                        alin1 = &i;
                        qt1 = (*alin1)[0];
                    }
                }
            }
            if (alin1) {
                for (auto &i : *alin1) {
                    if (isleftpt((*i)[idx1], coordinate[idx1], qtl)) {
                        ptl = i;
                        qtl = (*ptl)[idx1];
                    }
                    if (isrightpt((*i)[idx1], coordinate[idx1], qtr)) {
                        ptr = i;
                        qtr = (*ptr)[idx1];
                    }
                }
            }

            if (ptr && ptl && ptr == ptl) {
                element1[ewns] = ptr;
                ptr = ptl = nullptr;
            }
            if (ptr) {
                element1[ewns0] = ptr;
                element1[ewns] = nullptr;
            }
            if (ptl) {
                element1[ewns1] = ptl;
                element1[ewns] = nullptr;
            }
        }
    };

    assignElement(E, std::numeric_limits<double>::max(), yline, isright, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 0, 1, EN, ES, normal[0]);
    assignElement(W, -std::numeric_limits<double>::max(), yline, isleft, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 0, 1, WN, WS, -normal[0]);
    assignElement(N, std::numeric_limits<double>::max(), xline, isright, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 1, 0, NE, NW, normal[1]);
    assignElement(S, -std::numeric_limits<double>::max(), xline, isleft, -std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(), 1, 0, SE, SW, -normal[1]);

    auto assignNoElement = [&](const EWNS &ewns, const EWNS &ewns0, const EWNS &ewns1) -> void {
        if (!(element1[ewns] || element1[ewns0] || element1[ewns1])) {
            element1[ewns] = element[ewns];
            element1[ewns0] = element[ewns0];
            element1[ewns1] = element[ewns1];
        }
    };

    assignNoElement(E, EN, ES);
    assignNoElement(W, WN, WS);
    assignNoElement(N, NE, NW);
    assignNoElement(S, SE, SW);

}

void AGM::Point::findElement(std::vector<Point> *src) {
    int index{};
    auto ele = std::array<Point *, 12> {nullptr, };
    for (auto &item : src->at(idx).getElement().getElement()) {
        if (item) {
            ele.at(index) = &(pts->at(item->getIdx()));
        }
        ++index;
    }
    element.setElement(ele);
}

void AGM::Point::findElement1(std::vector<Point> *src) {
    int index{};
    auto ele = std::array<Point *, 12> {nullptr, };
    for (auto &item : src->at(idx).getElement1().getElement()) {
        if (item) {
            ele.at(index) = &(pts->at(item->getIdx()));
        }
        ++index;
    }
    element1.setElement(ele);
}

void AGM::Point::calcRepresentationFormula_cross() {
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
            row[i][idx + ptsnum] += row[i][pt->idx + ptsnum];
            row[i].remove(pt->idx + ptsnum);
        }
    };

    row[0][idx] = -UNITVALUE;
    row[0][element[W]->idx] = mp * gfunc_x.green_function_t(xm);
    row[0][element[E]->idx] = -mp * gfunc_x.green_function_t(xp);

    row[0][element[W]->idx + ptsnum] = gfunc_x.green_integral('l');
    row[0][idx + ptsnum] = gfunc_x.green_integral('c');
    row[0][element[E]->idx + ptsnum] = gfunc_x.green_integral('r');

    rhsMatrixRow[0][element[W]->idx] = gfunc_x.green_integral('l');
    rhsMatrixRow[0][idx] = gfunc_x.green_integral('c');
    rhsMatrixRow[0][element[E]->idx] = gfunc_x.green_integral('r');

    pMatrixRow[0][element[W]->idx] = gfunc_x.green_integral_t('l');
    pMatrixRow[0][idx] = gfunc_x.green_integral_t('c');
    pMatrixRow[0][element[E]->idx] = gfunc_x.green_integral_t('r');

    row[1][idx] = -UNITVALUE;
    row[1][element[S]->idx] = mp * gfunc_y.green_function_t(ym);
    row[1][element[N]->idx] = -mp * gfunc_y.green_function_t(yp);

    row[1][element[S]->idx + ptsnum] = -gfunc_y.green_integral('l');
    row[1][idx + ptsnum] = -gfunc_y.green_integral('c');
    row[1][element[N]->idx + ptsnum] = -gfunc_y.green_integral('r');

    rhsMatrixRow[1][element[S]->idx] = gfunc_y.green_integral('l');
    rhsMatrixRow[1][idx] = gfunc_y.green_integral('c');
    rhsMatrixRow[1][element[N]->idx] = gfunc_y.green_integral('r');

    pMatrixRow[1][element[S]->idx] = gfunc_y.green_integral_t('l');
    pMatrixRow[1][idx] = gfunc_y.green_integral_t('c');
    pMatrixRow[1][element[N]->idx] = gfunc_y.green_integral_t('r');

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);

    matrixRow[0] = row[0] + row[1];
    matrixRow[1] = row[0] - row[1];
}

void AGM::Point::calcRepresentationFormula_dirichlet() {
    matrixRow[0][idx] = UNITVALUE;

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

void AGM::Point::calcRepresentationFormula_neumann() {
//    calcRepresentationFormula_neumann_D();
//    calcRepresentationFormula_neumann_N();
    std::string string{};
    if (element[E] == this || element[W] == this) {
        string += "x";
    }
    if (element[N] == this || element[S] == this) {
        string += "y";
    }
    if (!string.empty()) {
        std::array<matrix_row, 2> row{};
        if (string == "x") {
            if (this->getAxialLine('x')) {
                row[0] = calcRepresentationFormula_neumann_N(string);
            } else {
                row[0] = calcRepresentationFormula_neumann_D(string);
            }
        } else if (string == "y") {
            if (this->getAxialLine('y')) {
                row[1] = calcRepresentationFormula_neumann_N(string);
            } else {
                row[1] = calcRepresentationFormula_neumann_D(string);
            }
        } else {
            string = "x";
            if (this->getAxialLine('x')) {
                row[0] = calcRepresentationFormula_neumann_N(string);
            } else {
                row[0] = calcRepresentationFormula_neumann_D(string);
            }
            string = "y";
            if (this->getAxialLine('y')) {
                row[1] = calcRepresentationFormula_neumann_N(string);
            } else {
                row[1] = calcRepresentationFormula_neumann_D(string);
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

void AGM::Point::calcRepresentationFormula_neumann_D() {
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
            row[i][pt->idx] += mp0 * func->green_function_ttau(d);
            row[i][pt->idx + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
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

    for (int i = 0; i < 2; ++i) {
        while (row[i].back().idx >= 4 * ptsnum) {
            pMatrixRow[i][row[i].back().idx - 4 * ptsnum] = row[i].back().value * normal[i];
            row[i].pop_back();
        }
        while (row[i].back().idx >= 2 * ptsnum) {
            rhsMatrixRow[i][row[i].back().idx - 2 * ptsnum] = row[i].back().value * normal[i];
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

//    calcRepresentationFormula_boundary();

    AGM::Point *pt{};
    std::string string{};
    EWNS ewns = findPoint_dirichlet_and_Neumann(pt);
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

void AGM::Point::calcRepresentationFormula_neumann_N() {
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

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
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

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->idx] -= UNITVALUE;
            } else {
                row[i][pt->idx] += UNITVALUE;
            }
            row[i][pt->idx + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_t_ND(c);

        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE, i0, std::abs(mp0));
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
                row[i][pt->idx] -= UNITVALUE;
            } else {
                row[i][pt->idx] += UNITVALUE;
            }
            row[i][pt->idx + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_t_DN(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE, i0, std::abs(mp0));
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
    calcRepresentationFormula_dirichlet_and_Neumann(pt, ewns, 1);
}

AGM::matrix_row AGM::Point::calcRepresentationFormula_neumann_D(std::string &string) {
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
            row[i][pt->idx] += mp0 * func->green_function_ttau(d);
            row[i][pt->idx + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
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
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::Point::calcRepresentationFormula_neumann_D(std::string &string)");
    }
    return row[0];
}

AGM::matrix_row AGM::Point::calcRepresentationFormula_neumann_N(std::string &string) {
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

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
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

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->idx] -= UNITVALUE;
            } else {
                row[i][pt->idx] += UNITVALUE;
            }
            row[i][pt->idx + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_t_ND(c);

        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE, i0, std::abs(mp0));
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
                row[i][pt->idx] -= UNITVALUE;
            } else {
                row[i][pt->idx] += UNITVALUE;
            }
            row[i][pt->idx + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_t_DN(c);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE, i0, std::abs(mp0));
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

    if (string == "x") {
        return row[0];
    } else if (string == "y") {
        return row[1];
    } else {
        printError("AGM::matrix_row AGM::Point::calcRepresentationFormula_neumann_N(std::string &string)");
    }
    return row[0];
}

void AGM::Point::calcRepresentationFormula_dirichlet_and_Neumann(AGM::Point *pt, const EWNS &ewns, int order) {
    Point point = *pt;
    double xm = ewns == E || ewns == W ? point[W]->getCoordinate()[0] : ZEROVALUE;
    double xb = point.getCoordinate()[0];
    double xp = ewns == E || ewns == W ? point[E]->getCoordinate()[0] : ZEROVALUE;
    double ym = ewns == N || ewns == S ? point[S]->getCoordinate()[1] : ZEROVALUE;
    double yb = point.getCoordinate()[1];
    double yp = ewns == N || ewns == S ? point[N]->getCoordinate()[1] : ZEROVALUE;

    auto secondOrderExtrapolation = [&](int i, double tm, double tb, double tp, EWNS ewns1, EWNS ewns2) -> void {
        double d = (tp - tm) * (tp - tb) * (tb - tm);
        double t0 = coordinate[i];
        auto firstTerm = [t0, d](double w0, double w1) -> double { return t0 * t0 * (w0 - w1) / d; };
        auto secondTerm = [t0, d](double w0, double w1) -> double { return t0 * (w0 * w0 - w1 * w1) / d; };
        auto thirdTerm = [d](double w0, double w1) -> double { return w0 * w1 * (w0 - w1) / d; };
        auto secondOrder = [&](double w0, double w1) -> double {
            return firstTerm(w0, w1) - secondTerm(w0, w1) + thirdTerm(w0, w1);
        };
        matrixRow[1][point[ewns1]->getIdx() + ptsnum] = secondOrder(tp, tb);
        matrixRow[1][point.getIdx() + ptsnum] = -secondOrder(tp, tm);
        matrixRow[1][point[ewns2]->getIdx() + ptsnum] = secondOrder(tb, tm);
        matrixRow[1][idx + ptsnum] = -UNITVALUE;
    };
    auto firstOrderExtrapolation = [&](int i, double t1, double t2, EWNS ewns0) -> void {
        double d = t2 - t1;
        double t0 = coordinate[i];
        auto firstTerm = [t0, d]() -> double { return t0 / d; };
        auto secondTerm = [d](double t) -> double { return t / d; };
        auto firstOrder = [&](double t) -> double { return -firstTerm() + secondTerm(t); };
        matrixRow[1][point[ewns0]->getIdx() + ptsnum] = -firstOrder(t1);
        matrixRow[1][point.getIdx() + ptsnum] = firstOrder(t2);
        matrixRow[1][idx + ptsnum] = -UNITVALUE;
    };
    auto zeroOrderExtrapolation = [&](EWNS ewns1) -> void {
        matrixRow[1][point[ewns1]->getIdx() + ptsnum] = UNITVALUE;
        matrixRow[1][idx + ptsnum] = -UNITVALUE;
    };
    if (ewns == E || ewns == W) {
        if (order == 2) {
            secondOrderExtrapolation(0, xm, xb, xp, W, E);
        } else if (order == 1) {
            if (ewns == E) {
                firstOrderExtrapolation(0, xb, xp, E);
            } else {
                firstOrderExtrapolation(0, xb, xm, W);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(ewns);
        } else {
            printError("AGM::Point::calcRepresentationFormula_dirichlet_and_Neumann", "order (which is %d) error",
                       order);
        }
    } else if (ewns == N || ewns == S) {
        if (order == 2) {
            secondOrderExtrapolation(1, ym, yb, yp, S, N);
        } else if (order == 1) {
            if (ewns == N) {
                firstOrderExtrapolation(1, yb, yp, N);
            } else {
                firstOrderExtrapolation(1, yb, ym, S);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(ewns);
        } else {
            printError("AGM::Point::calcRepresentationFormula_dirichlet_and_Neumann", "order (which is %d) erorr",
                       order);
        }
    } else {
        printError("AGM::Point::calcRepresentationFormula_dirichlet_and_Neumann", "ewns (which is %d) erorr", ewns);
    }
}

void AGM::Point::calcRepresentationFormula_phi() {
    Point *pt0{}, *pt1{}, *pt2{};
    auto assignPoints = [&](int i0, int i1, EWNS e0, EWNS e1, EWNS e2, EWNS e3, EWNS e20, EWNS e21, EWNS e30, EWNS e31) {
        if (axialLine[i0]) {
            if (iszero((*this) - *(getElement().getElement()[e0]))) {
                pt0 = getElement().getElement()[e1];
                pt1 = getElement1().getElement()[e1];
            } else if (iszero((*this) - *(getElement().getElement()[e1]))) {
                pt0 = getElement().getElement()[e0];
                pt1 = getElement1().getElement()[e0];
            }

            if (axialLine[i1]) {
                if (iszero((*this) - *(getElement().getElement()[e2]))) {
                    pt2 = getElement().getElement()[e3];
                } else if (iszero((*this) - *(getElement().getElement()[e3]))) {
                    pt2 = getElement().getElement()[e2];
                }
            } else {
                if (getElement().getElement()[e2] && getElement().getElement()[e3]) {
                    if (!iszero((*this) - *(getElement().getElement()[e2])) &&
                    !iszero((*this) - *(getElement().getElement()[e3]))) {
                        pt2 = getElement().getElement()[e3];
                    } else if (iszero((*this) - *(getElement().getElement()[e2]))) {
                        pt2 = getElement().getElement()[e3];
                    } else if (iszero((*this) - *(getElement().getElement()[e3]))) {
                        pt2 = getElement().getElement()[e2];
                    }
                } else if (!getElement().getElement()[e2]) {
                    pt2 = getElement().getElement()[e21];
                } else if (!getElement().getElement()[e3]) {
                    pt2 = getElement().getElement()[e31];
                }
            }
        }
    };

    if (axialLine[0]) {
        assignPoints(0, 1, E, W, N, S, NE, NW, SE, SW);
    } else if (axialLine[1]) {
        assignPoints(1, 0, N, S, E, W, EN, ES, WN, WS);
    }
    calcRepresentationFormula_phi_plane(pt0, pt1, pt2);
}

void AGM::Point::calcRepresentationFormula_phi_plane(Point *pt0, Point *pt1, Point *pt2) {
    double x{coordinate[0]}, y{coordinate[1]};
    double x0{pt0->getCoordinate()[0]}, x1{pt1->getCoordinate()[0]}, x2{pt2->getCoordinate()[0]};
    double y0{pt0->getCoordinate()[1]}, y1{pt1->getCoordinate()[1]}, y2{pt2->getCoordinate()[1]};
    double cz{x0 * y1 - x0 * y2 - x1 * y0 + x1 * y2 + x2 * y0 - x2 * y1};
    double cz0{-x * y1 + x * y2 + x1 * y - x1 * y2 - x2 * y + x2 * y1};
    double cz1{x * y0 - x * y2 - x0 * y + x0 * y2 + x2 * y - x2 * y0};
    double cz2{-x * y0 + x * y1 + x0 * y - x0 * y1 - x1 * y + x1 * y0};

//    matrixRow[0][idx + ptsnum] = -cz;
//    matrixRow[0][pt0->getIdx() + ptsnum] = -cz0;
//    matrixRow[0][pt1->getIdx() + ptsnum] = -cz1;
//    matrixRow[0][pt2->getIdx() + ptsnum] = -cz2;

    matrixRow[1][idx + ptsnum] = cz;
    matrixRow[1][pt0->getIdx() + ptsnum] = cz0;
    matrixRow[1][pt1->getIdx() + ptsnum] = cz1;
    matrixRow[1][pt2->getIdx() + ptsnum] = cz2;
}

void AGM::Point::calcRepresentationFormula_interface() {
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

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
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
            row[i][pt->idx] += mp0 * func->green_function_t(d);
            if (!checkInterface(this)) row[i][pt->idx + ptsnum] += sign * func->green_integral(c);
            else row[i][pt->idx + ptsnum] += sign * func->green_integral(C);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_t(c);
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
        matrixRow[1][idx + ptsnum] = UNITVALUE;
    }
}

void AGM::Point::calcRepresentationFormula_boundary() {
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

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    std::array<matrix_row, 2> row{};

    if (nx == E) {
        row[0][idx] = UNITVALUE;
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
        row[0][idx] = UNITVALUE;
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
        row[1][idx] = UNITVALUE;
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
        row[1][idx] = UNITVALUE;
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
        printError("AGM::Point::calcRepresentationFormula_boundary()");
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
        printError("AGM::Point::calcRepresentationFormula_boundary1()");
    }
}

void AGM::Point::calcRepresentationFormula() {
    switch (condition) {
        case 'C':
            calcRepresentationFormula_cross();
            break;
        case 'D':
            calcRepresentationFormula_dirichlet();
            break;
        case 'N':
            calcRepresentationFormula_neumann();
            break;
        case 'I':
            calcRepresentationFormula_interface();
            break;
        case 'F':
            std::cout << "F" << std::endl;
            break;
        default:
            printError("AGM::Point::calcRepresentationFormula", "boundary condition (which is %c) is wrong", condition);
    }
}

double AGM::Point::checkRepresentationFormula() {
    double sum_value{};
    auto func = function();
    for (const auto &i : matrixRow[0]) {
        if (i.idx < ptsnum) {
            func.setPoint(pts->at(i.idx));
            func.setEps(this->getMp());
            func.setTime(0.1);
            sum_value += i.value * func.u();
        } else if (i.idx < 2 * ptsnum) {
            func.setPoint(pts->at(i.idx % ptsnum));
            func.setEps(this->getMp());
            func.setTime(0.1);
            sum_value += i.value * func.phi();
        }
    }
    func.setX(0.02);
    func.setY(0.98);
    func.setEps(this->getMp());
    func.setTime(0.1);
    return sum_value - rb[0];
}

double AGM::Point::checkDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                        const std::function<double(int)> &gx, const std::function<double(int)> &gy) {
    double sum_value{};
    double test{};
    auto func = function();
    for (const auto &item : dMatrixRow[0]) {
        if (item.idx < ptsnum) {
            func.setPoint(pts->at(item.idx));
            func.setEps(this->getMp());
            func.setTime(0.1);
            sum_value += item.value * func.u();
        } else if (item.idx < 2 * ptsnum) {
            func.setPoint(pts->at(item.idx % ptsnum));
            func.setEps(this->getMp());
            func.setTime(0.1);
            sum_value += item.value * func.phi();
        } else if (item.idx < 3 * ptsnum) {
            sum_value += item.value * fx(item.idx % ptsnum);
        } else if (item.idx < 4 * ptsnum) {
            sum_value += item.value * fy(item.idx % ptsnum);
        } else if (item.idx < 5 * ptsnum) {
            sum_value += item.value * gx(item.idx % ptsnum);
        } else if (item.idx < 6 * ptsnum) {
            sum_value += item.value * gy(item.idx % ptsnum);
        }
    }
    return sum_value;
}

void AGM::Point::makeDifferentiation() {
    switch (condition) {
        case 'C':
            makeDifferentiation_cross();
            break;
        case 'I':
            makeDifferentiation_interface();
            break;
        default:
            std::string string{};
            if (element[E] == this || element[W] == this) {
                string += "x";
            }
            if (element[N] == this || element[S] == this) {
                string += "y";
            }
            makeDifferentiation_boundary();
            if (!string.empty()) makeDifferentiation_boundary(string);
            break;
    }
}

void AGM::Point::makeDifferentiation_cross() {
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
            dMatrixRow[i][idx + ptsnum] += dMatrixRow[i][pt->idx + ptsnum];
            dMatrixRow[i].remove(pt->idx + ptsnum);
        }
    };

    dMatrixRow[0][element[W]->idx] = mp * gfunc_x.green_function_ttau(xm);
    dMatrixRow[0][element[E]->idx] = -mp * gfunc_x.green_function_ttau(xp);

    dMatrixRow[0][element[W]->idx + ptsnum] = gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->idx + ptsnum] = gfunc_x.green_integral_tau('r');

    dMatrixRow[0][element[W]->idx + 2 * ptsnum] = gfunc_x.green_integral_tau('l');
    dMatrixRow[0][idx + 2 * ptsnum] = gfunc_x.green_integral_tau('c');
    dMatrixRow[0][element[E]->idx + 2 * ptsnum] = gfunc_x.green_integral_tau('r');

    dMatrixRow[0][element[W]->idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('l') + gfunc_x.green_function_tau(xm);
    dMatrixRow[0][idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    dMatrixRow[0][element[E]->idx + 4 * ptsnum] = gfunc_x.green_integral_ttau('r') - gfunc_x.green_function_tau(xp);

    dMatrixRow[1][element[S]->idx] = mp * gfunc_y.green_function_ttau(ym);
    dMatrixRow[1][element[N]->idx] = -mp * gfunc_y.green_function_ttau(yp);

    dMatrixRow[1][element[S]->idx + ptsnum] = -gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + ptsnum] = -gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->idx + ptsnum] = -gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->idx + 2 * ptsnum] = gfunc_y.green_integral_tau('l');
    dMatrixRow[1][idx + 2 * ptsnum] = gfunc_y.green_integral_tau('c');
    dMatrixRow[1][element[N]->idx + 2 * ptsnum] = gfunc_y.green_integral_tau('r');

    dMatrixRow[1][element[S]->idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('l') + gfunc_y.green_function_tau(ym);
    dMatrixRow[1][idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    dMatrixRow[1][element[N]->idx + 4 * ptsnum] = gfunc_y.green_integral_ttau('r') - gfunc_y.green_function_tau(yp);

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);
}

void AGM::Point::makeDifferentiation_boundary() {
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

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
                            int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'l' ? UNITVALUE : -UNITVALUE;

        if (pt) {
            dMatrixRow[i][pt->idx] += mp0 * func->green_function_ttau(d);
            dMatrixRow[i][pt->idx + ptsnum] += sign * func->green_integral_tau(c);
            dMatrixRow[i][pt->idx + 2 * ptsnum] += func->green_integral_tau(c);
            dMatrixRow[i][pt->idx + 4 * ptsnum] += func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);
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

void AGM::Point::makeDifferentiation_boundary(const std::string &string) {
    if (string == "x") {
        if (this->getAxialLine('x')) makeDifferentiation_boundary_N(string);
        else makeDifferentiation_boundary_D(string);
    } else if (string == "y") {
        if (this->getAxialLine('y')) makeDifferentiation_boundary_N(string);
        else makeDifferentiation_boundary_D(string);
    } else {
        if (this->getAxialLine('x')) makeDifferentiation_boundary_N("x");
        else makeDifferentiation_boundary_D("x");

        if (this->getAxialLine('y')) makeDifferentiation_boundary_N("y");
        else makeDifferentiation_boundary_D("y");
    }
}

void AGM::Point::makeDifferentiation_boundary_D(const std::string &string) {
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
            row[i][pt->idx] += mp0 * func->green_function_ttau(d);
            row[i][pt->idx + ptsnum] += sign * func->green_integral_tau(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_tau(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_ttau(c) + sign1 * func->green_function_tau(d);

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

void AGM::Point::makeDifferentiation_boundary_N(const std::string &string) {
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

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
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

    auto assignMatrix_ND = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'c' ? ZEROVALUE : -UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->idx] -= UNITVALUE;
            } else {
                row[i][pt->idx] += UNITVALUE;
            }

            row[i][pt->idx + ptsnum] += sign * func->green_integral_ND(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_ND(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_t_ND(c) + sign1 * func->green_function_ND(d);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE, i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_ND(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_ND(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_ND(c) + sign1 * func->green_function_ND(d), i0);
        }
    };

    auto assignMatrix_DN = [&](Point *pt, Point *ptl, Point *ptr, double mp0, GreenfunctionNeumann *func, double d,
                               int i, int i0, char c) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        double sign1 = c == 'c' ? ZEROVALUE : UNITVALUE;
        if (pt) {
            if (c == 'c') {
                row[i][pt->idx] -= UNITVALUE;
            } else {
                row[i][pt->idx] += UNITVALUE;
            }
            row[i][pt->idx + ptsnum] += sign * func->green_integral_DN(c);
            row[i][pt->idx + 2 * ptsnum] += func->green_integral_DN(c);
            row[i][pt->idx + 4 * ptsnum] += func->green_integral_t_DN(c) + sign1 * func->green_function_DN(d);
        } else {
            if (c == 'c') {
                row[i] += approximateSol(ptl, ptr, -UNITVALUE, i0, std::abs(mp0));
            } else {
                row[i] += approximateSol(ptl, ptr, UNITVALUE, i0, std::abs(mp0));
            }
            row[i] += approximatePhi(ptl, ptr, sign * func->green_integral_DN(c), i0);
            row[i] += approximateRhs(ptl, ptr, func->green_integral_DN(c), i0);
            row[i] += approximateRhs_t(ptl, ptr, func->green_integral_t_DN(c) + sign1 * func->green_function_DN(d), i0);
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
        printError("AGM::Point::makeDifferentiation_boundary1()");
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
    } else if (string == "y") {
        dv[0] = ZEROVALUE;
        dMatrixRow[1] = row[1];
    } else if (string == "xy") {
        dMatrixRow[0] = row[0];
        dMatrixRow[1] = row[1];
    }
}

void AGM::Point::makeDifferentiation_interface() {
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

    auto approximateSol = [&](Point *ptl, Point *ptr, double coefficient, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
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

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, double mp0, Greenfunction *func, double d,
                            int i, int i0, char c, char C) -> void {
        double sign = i ? -UNITVALUE : UNITVALUE;
        if (pt) {
            dMatrixRow[i][pt->idx] += mp0 * func->green_function_ttau(d);
            if (!checkInterface(this)) dMatrixRow[i][pt->idx + ptsnum] += sign * func->green_integral_tau(c);
            else dMatrixRow[i][pt->idx + ptsnum] += sign * func->green_integral_tau(C);
            dMatrixRow[i][pt->idx + 2 * ptsnum] += func->green_integral_tau(c);
            dMatrixRow[i][pt->idx + 4 * ptsnum] += func->green_integral_ttau(c);
        } else {
            dMatrixRow[i] += approximateSol(ptl, ptr, mp0 * func->green_function_ttau(d), i0, std::abs(mp0));
            if (!checkInterface(this))
                dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(c), i0);
            else dMatrixRow[i] += approximatePhi(ptl, ptr, sign * func->green_integral_tau(C), i0);
            dMatrixRow[i] += approximateRhs(ptl, ptr, func->green_integral_tau(c), i0);
            dMatrixRow[i] += approximateRhs_t(ptl, ptr, func->green_integral_ttau(c), i0);
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

void AGM::Point::calcDifferentiation() {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * pts->at(item.idx % ptsnum)["phi"];
        } else {
            dx += item.value * HALFVALUE * pts->at(item.idx % ptsnum)["rhs"];
        }
    });
    value["dx"] = dx;
    value["dxx"] = -(value["phi"] + HALFVALUE * value["rhs"]) / mp;
    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * pts->at(item.idx % ptsnum)["phi"];
        } else {
            dy += item.value * HALFVALUE * pts->at(item.idx % ptsnum)["rhs"];
        }
    });
    value["dy"] = dy;
    value["dyy"] = -(-value["phi"] + HALFVALUE * value["rhs"]) / mp;
}

void AGM::Point::calcDifferentiation(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * pts->at(item.idx % ptsnum)["phi"];
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
            dy += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * pts->at(item.idx % ptsnum)["phi"];
        } else if (item.idx < 3 * ptsnum) {
            dy += item.value * g(item.idx % ptsnum);
        } else if (item.idx < 4 * ptsnum) {
            dy += item.value * f(item.idx % ptsnum);
        }
    });
    value["dy"] = dy;
    value["dyy"] = -(-value["phi"] + HALFVALUE * g(idx)) / mp;
}

void AGM::Point::calcDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                     const std::function<double(int)> &gx, const std::function<double(int)> &gy) {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * pts->at(item.idx % ptsnum)["phi"];
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
    value["dxx"] = -(value["phi"] + fx(idx)) / mp;
    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * pts->at(item.idx % ptsnum)["phi"];
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
    value["dyy"] = -(-value["phi"] + fy(idx)) / mp;
}

void AGM::Point::calcDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                     const std::function<double(int)> &gx, const std::function<double(int)> &gy,
                                     const std::function<double(int)> &hx, const std::function<double(int)> &hy) {
    double dx{}, dy{};
    std::for_each(dMatrixRow[0].begin(), dMatrixRow[0].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dx += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dx += item.value * pts->at(item.idx % ptsnum)["phi"];
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
    value["dxx"] = -(value["phi"] + hx(idx)) / mp;
    std::for_each(dMatrixRow[1].begin(), dMatrixRow[1].end(), [&](matrix_element &item) {
        if (item.idx < ptsnum) {
            dy += item.value * pts->at(item.idx)["sol"];
        } else if (item.idx < 2 * ptsnum) {
            dy += item.value * pts->at(item.idx % ptsnum)["phi"];
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
    value["dyy"] = -(-value["phi"] + hy(idx)) / mp;
}

void AGM::Point::calcDifferentiation_dirichlet_and_Neumann(AGM::Point *pt, const EWNS &ewns, int order,
                                                           const std::string &string) {
    Point point = *pt;
    double xm = point[W]->getCoordinate()[0];
    double xb = point.getCoordinate()[0];
    double xp = point[E]->getCoordinate()[0];
    double ym = point[S]->getCoordinate()[1];
    double yb = point.getCoordinate()[1];
    double yp = point[N]->getCoordinate()[1];

    auto secondOrderExtrapolation = [&](int i, double tm, double tb, double tp, EWNS ewns1, EWNS ewns2) -> void {
        double d = (tp - tm) * (tp - tb) * (tb - tm);
        double t0 = coordinate[i];
        auto firstTerm = [t0, d](double w0, double w1) -> double { return t0 * t0 * (w0 - w1) / d; };
        auto secondTerm = [t0, d](double w0, double w1) -> double { return t0 * (w0 * w0 - w1 * w1) / d; };
        auto thirdTerm = [d](double w0, double w1) -> double { return w0 * w1 * (w0 - w1) / d; };
        auto secondOrder = [&](double w0, double w1) -> double {
            return firstTerm(w0, w1) - secondTerm(w0, w1) + thirdTerm(w0, w1);
        };
        this->getValue()[string] =
                point[ewns1]->getValue()[string] * secondOrder(tp, tb) - point[string] * secondOrder(tp, tm) +
                point[ewns2]->getValue()[string] * secondOrder(tb, tm);
    };
    auto firstOrderExtrapolation = [&](int i, double t1, double t2, EWNS ewns0) -> void {
        double d = t2 - t1;
        double t0 = coordinate[i];
        auto firstTerm = [t0, d]() -> double { return t0 / d; };
        auto secondTerm = [d](double t) -> double { return t / d; };
        auto firstOrder = [&](double t) -> double { return -firstTerm() + secondTerm(t); };
        this->getValue()[string] = -point[ewns0]->getValue()[string] * firstOrder(t1) + point[string] * firstOrder(t2);
    };
    auto zeroOrderExtrapolation = [&](EWNS ewns1) -> void {
        this->getValue()[string] = point[ewns1]->getValue()[string];
    };
    if (ewns == E || ewns == W) {
        if (order == 2) {
            secondOrderExtrapolation(0, xm, xb, xp, W, E);
        } else if (order == 1) {
            if (ewns == E) {
                firstOrderExtrapolation(0, xb, xp, E);
            } else {
                firstOrderExtrapolation(0, xb, xm, W);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(ewns);
        } else {
            printError("AGM::Point::calcDifferentiation_dirichlet_and_Neumann", "order (which is %d) error",
                       order);
        }
    } else if (ewns == N || ewns == S) {
        if (order == 2) {
            secondOrderExtrapolation(1, ym, yb, yp, S, N);
        } else if (order == 1) {
            if (ewns == N) {
                firstOrderExtrapolation(1, yb, yp, N);
            } else {
                firstOrderExtrapolation(1, yb, ym, S);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(ewns);
        } else {
            printError("AGM::Point::calcDifferentiation_dirichlet_and_Neumann", "order (which is %d) erorr",
                       order);
        }
    } else {
        printError("AGM::Point::calcDifferentiation_dirichlet_and_Neumann", "ewns (which is %d) erorr", ewns);
    }
}

void AGM::Point::calcDifferentiationWithPressure(std::vector<Point> *pressure, char i) {
    double dx{}, dy{};
    if (i == 'u') {
        std::for_each(pdMatrixRow[0].begin(), pdMatrixRow[0].end(), [&](matrix_element &item) {
            dx += item.value * pressure->at(item.idx)["sol"];
        });
        value["dx"] += dx;
        value["dxx"] = -(value["phi"] + HALFVALUE * value["rhs"] + dx) / mp;
    } else if (i == 'v') {
        std::for_each(pdMatrixRow[1].begin(), pdMatrixRow[1].end(), [&](matrix_element &item) {
            dy += item.value * pressure->at(item.idx)["sol"];
        });
        value["dy"] += dy;
        value["dyy"] = -(-value["phi"] + HALFVALUE * value["rhs"] + dy) / mp;
    }
}

AGM::EWNS AGM::Point::findPoint_dirichlet_and_Neumann(AGM::Point *&pt) {
    auto func = [&](const EWNS &ewns) -> bool {
        if (*this - *element[ewns] > 1.0E-10) {
            pt = (*(element[ewns]))[ewns];
            if (*this - *element[ewns] < 1.0E-10) {
                pt = (*pt)[ewns];
            }
            return true;
        } else {
            return false;
        }
    };
    for (const auto &i : axialLine) {
        if (i) {
            if (!(i->empty())) {
                if (i->getMark() == 'x') {
                    if (func(E)) return W;
                    if (func(W)) return E;
                }
                if (i->getMark() == 'y') {
                    if (func(N)) return S;
                    if (func(S)) return N;
                }
            }
        }
    }
//    /* ------ */
//    if (!(axialLine[0] || axialLine[1])) {
//        pt = this;
//        return SW;
//    }
//    /* ------ */

    printInformation();
    printError("AGM::Point::findPoint_dirichlet_and_Neumann");
    return E;
}

void AGM::Point::findPoint_non_axial_lines(Point *&pt0, Point *&pt1, Point *&pt2) {
    auto find_point = [&](EWNS ewns0, EWNS ewns1, EWNS ewns2, EWNS ewns3, EWNS ewns4, EWNS ewns5) -> void {
        if (element.getElement()[ewns0] && element.getElement()[ewns0] != this && (*element.getElement()[ewns0] - *this) > NEARZERO) {
            if (pt0) {
                if (pt1) pt2 = element.getElement()[ewns0];
                else pt1 = element.getElement()[ewns0];
            } else {
                pt0 = element.getElement()[ewns0];
            }
        } else if (element.getElement()[ewns1] && element.getElement()[ewns2]) {
            if (pt0) {
                pt1 = element.getElement()[ewns1];
                pt2 = element.getElement()[ewns2];
            } else {
                pt0 = element.getElement()[ewns1];
                pt1 = element.getElement()[ewns2];
            }
        } else if (element.getElement()[ewns3] && element.getElement()[ewns3] != this && (*element.getElement()[ewns3] - *this) > NEARZERO) {
            if (pt0) {
                if (pt1) pt2 = element.getElement()[ewns3];
                else pt1 = element.getElement()[ewns3];
            } else {
                pt0 = element.getElement()[ewns3];
            }
        } else if (element.getElement()[ewns4] && element.getElement()[ewns5]) {
            if (pt0) {
                pt1 = element.getElement()[ewns4];
                pt2 = element.getElement()[ewns5];
            } else {
                pt0 = element.getElement()[ewns4];
                pt1 = element.getElement()[ewns5];
            }
        }
    };
    auto find_min_point = [&]() -> Point* {
        double min_value{1e3};
        Point *return_point{};

        for (auto &item: getElement().getElement()) {
            if (!item) {
                continue;
            }
            if (item->getIdx() == this->getIdx()) {
                continue;
            }
            if (*item - *this < min_value) {
                min_value = *item - *this;
                return_point = item;
            }
        }
        return return_point;
    };

    pt0 = find_min_point();
//    find_point(E, EN, ES, W, WN, WS);
//    find_point(N, NE, NW, S, SE, SW);
}

AGM::Point &AGM::Point::operator=(const AGM::Point &rhs) {
    if (this != &rhs) {
        idx = rhs.idx;
        coordinate = rhs.coordinate;
        normal = rhs.normal;
        value = rhs.value;
        element = rhs.element;
        condition = rhs.condition;
        mp = rhs.mp;
        matrixRow = rhs.matrixRow;
        rhsMatrixRow = rhs.rhsMatrixRow;
        dMatrixRow = rhs.dMatrixRow;
        axialLine = rhs.axialLine;
        rb = rhs.rb;
    }
    return *this;
}

double &AGM::Point::operator[](int i) {
    return coordinate[i];
}

const double &AGM::Point::operator[](int i) const {
    return coordinate[i];
}

AGM::Point *&AGM::Point::operator[](AGM::EWNS ewns) {
    return element[ewns];
}

double &AGM::Point::operator[](const std::string &string) {
    return value[string];
}

const double &AGM::Point::operator[](const std::string &string) const {
    return value[string];
}

double AGM::Point::operator-(const AGM::Point &rhs) {
    return (coordinate - rhs.coordinate).norm();
}

void AGM::Point::updateRb() {
    switch (condition) {
        case 'C':
            updateRb_cross();
            break;
        case 'D':
            updateRb_dirichlet();
            break;
        case 'N':
            updateRb_neumann();
            break;
        case 'I':
            updateRb_interface();
            break;
        case 'F':
            std::cout << "F" << std::endl;
            break;
    }
}

void AGM::Point::updateRb_cross() {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        rb[0] -= j.value * HALFVALUE * pts->at(j.idx)["rhs"];
        rb[1] -= j.value * HALFVALUE * pts->at(j.idx)["rhs"];
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * HALFVALUE * pts->at(j.idx)["rhs"];
        rb[1] += j.value * HALFVALUE * pts->at(j.idx)["rhs"];
    }
}

void AGM::Point::updateRb_dirichlet() {
    rb[0] = value["bdv"];
    for (const auto &row : phiMatrixRow) {
        for (const auto &item : row) {
            rb[1] -= item.value * HALFVALUE * pts->at(item.idx)["rhs"];
        }
    }
}

void AGM::Point::updateRb_neumann() {
    rb[0] = value["bdv"];
    for (const auto &row : rhsMatrixRow) {
        for (const auto &item : row) {
            if (item.idx < ptsnum) {
                rb[0] -= item.value * HALFVALUE * pts->at(item.idx)["rhs"];
            }
        }
    }
    rb[1] = ZEROVALUE;
    for (const auto &row : phiMatrixRow) {
        for (const auto &item : row) {
            rb[1] -= item.value * HALFVALUE * pts->at(item.idx)["rhs"];
        }
    }
}

void AGM::Point::updateRb_interface() {
    updateRb_cross();
}

void AGM::Point::updateRb(const std::function<double(int)>& f) {
    switch (condition) {
        case 'C':
            updateRb_cross(f);
            break;
        case 'D':
            updateRb_dirichlet(f);
            break;
        case 'N':
            updateRb_neumann(f);
            break;
        case 'I':
            updateRb_interface(f);
            break;
        case 'F':
            std::cout << "F" << std::endl;
            break;
    }
}

void AGM::Point::updateRb_cross(const std::function<double(int)>& f) {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        rb[0] -= j.value * f(j.idx);
        rb[1] -= j.value * f(j.idx);
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * f(j.idx);
        rb[1] += j.value * f(j.idx);
    }
}

void AGM::Point::updateRb_dirichlet(const std::function<double(int)>& f) {
    rb[0] = value["bdv"];
    for (const auto &row : phiMatrixRow) {
        for (const auto &item : row) {
            rb[1] -= item.value * f(item.idx);
        }
    }
}

void AGM::Point::updateRb_neumann(const std::function<double(int)>& f) {
    rb[0] = value["bdv"];
    for (const auto &row : rhsMatrixRow) {
        for (const auto &item : row) {
            if (item.idx < ptsnum) {
                rb[0] -= item.value * f(item.idx);
            }
        }
    }
    rb[1] = ZEROVALUE;
    for (const auto &row : phiMatrixRow) {
        for (const auto &item : row) {
            rb[1] -= item.value * f(item.idx);
        }
    }
}

void AGM::Point::updateRb_interface(const std::function<double(int)>& f) {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
            rb[1] -= j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
            rb[1] -= j.value * f(j.idx - ptsnum);
        }
    }
    for (const auto &j : rhsMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
            rb[1] += j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
            rb[1] += j.value * f(j.idx - ptsnum);
        }
    }
}

void AGM::Point::updateRb_t(const std::function<double(int)>& f) {
    switch (condition) {
        case 'C':
            updateRb_t_cross(f);
            break;
        case 'D':
            updateRb_t_dirichlet(f);
            break;
        case 'N':
            updateRb_t_neumann(f);
            break;
        case 'I':
            updateRb_t_interface(f);
            break;
        case 'F':
            std::cout << "F" << std::endl;
            break;
    }
}

void AGM::Point::updateRb_t_cross(const std::function<double(int)>& f) {
    for (const auto &j : pMatrixRow[0]) {
        rb[0] -= j.value * f(j.idx);
        rb[1] -= j.value * f(j.idx);
    }
    for (const auto &j : pMatrixRow[1]) {
        rb[0] -= j.value * f(j.idx);
        rb[1] += j.value * f(j.idx);
    }
}

void AGM::Point::updateRb_t_dirichlet(const std::function<double(int)>& f) {

}

void AGM::Point::updateRb_t_neumann(const std::function<double(int)>& f) {
    for (const auto &j : pMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
        }
    }
    for (const auto &j : pMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
        }
    }
}

void AGM::Point::updateRb_t_interface(const std::function<double(int)>& f) {
    for (const auto &j : pMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
            rb[1] -= j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
            rb[1] -= j.value * f(j.idx - ptsnum);
        }
    }
    for (const auto &j : pMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
            rb[1] += j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
            rb[1] += j.value * f(j.idx - ptsnum);
        }
    }
}

void AGM::Point::updateRb(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    switch (condition) {
        case 'C':
            updateRb_cross(f, g);
            break;
        case 'D':
            updateRb_dirichlet(f, g);
            break;
        case 'N':
            updateRb_neumann(f, g);
            break;
        case 'I':
            updateRb_interface(f, g);
            break;
        case 'F':
            std::cout << "F" << std::endl;
            break;
    }
}

void AGM::Point::updateRb_cross(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        rb[0] -= j.value * f(j.idx);
        rb[1] -= j.value * f(j.idx);
    }
    for (const auto &j : rhsMatrixRow[1]) {
        rb[0] -= j.value * g(j.idx);
        rb[1] += j.value * g(j.idx);
    }
}

void AGM::Point::updateRb_dirichlet(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    rb[0] = value["bdv"];
    rb[1] = ZEROVALUE;
}

void AGM::Point::updateRb_neumann(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    rb[0] = value["bdv"];
    rb[1] = ZEROVALUE;
    for (const auto &item : rhsMatrixRow[0]) {
        if (item.idx < ptsnum) {
            rb[0] -= item.value * f(item.idx);
        } else if (item.idx < 2 * ptsnum) {
            rb[0] -= item.value * g(item.idx - ptsnum);
        }
    }
    for (const auto &item : rhsMatrixRow[1]) {
        if (item.idx < ptsnum) {
            rb[0] -= item.value * g(item.idx);
        } else if (item.idx < 2 * ptsnum) {
            rb[0] -= item.value * f(item.idx - ptsnum);
        }
    }
}

void AGM::Point::updateRb_interface(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &j : rhsMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
            rb[1] -= j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * g(j.idx - ptsnum);
            rb[1] -= j.value * g(j.idx - ptsnum);
        } else {
            printError("AGM::Point::updateRb_interface");
        }
    }
    for (const auto &j : rhsMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * g(j.idx);
            rb[1] += j.value * g(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
            rb[1] += j.value * f(j.idx - ptsnum);
        } else {
            printError("AGM::Point::updateRb_interface");
        }
    }
}

void AGM::Point::updateRb_t(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    switch (condition) {
        case 'C':
            updateRb_t_cross(f, g);
            break;
        case 'D':
            updateRb_t_dirichlet(f, g);
            break;
        case 'N':
            updateRb_t_neumann(f, g);
            break;
        case 'I':
            updateRb_t_interface(f, g);
            break;
        case 'F':
            std::cout << "F" << std::endl;
            break;
    }
}

void AGM::Point::updateRb_t_cross(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    for (const auto &j : pMatrixRow[0]) {
        rb[0] -= j.value * f(j.idx);
        rb[1] -= j.value * f(j.idx);
    }
    for (const auto &j : pMatrixRow[1]) {
        rb[0] -= j.value * g(j.idx);
        rb[1] += j.value * g(j.idx);
    }
}

void AGM::Point::updateRb_t_dirichlet(const std::function<double(int)>& f, const std::function<double(int)>& g) {
}

void AGM::Point::updateRb_t_neumann(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    for (const auto &item : pMatrixRow[0]) {
        if (item.idx < ptsnum) {
            rb[0] -= item.value * f(item.idx);
        } else if (item.idx < 2 * ptsnum) {
            rb[0] -= item.value * g(item.idx - ptsnum);
        }
    }
    for (const auto &item : pMatrixRow[1]) {
        if (item.idx < ptsnum) {
            rb[0] -= item.value * g(item.idx);
        } else if (item.idx < 2 * ptsnum) {
            rb[0] -= item.value * f(item.idx - ptsnum);
        }
    }
}

void AGM::Point::updateRb_t_interface(const std::function<double(int)>& f, const std::function<double(int)>& g) {
    for (const auto &j : pMatrixRow[0]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * f(j.idx);
            rb[1] -= j.value * f(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * g(j.idx - ptsnum);
            rb[1] -= j.value * g(j.idx - ptsnum);
        } else {
            printError("AGM::Point::updateRb_t_interface");
        }
    }
    for (const auto &j : pMatrixRow[1]) {
        if (j.idx < ptsnum) {
            rb[0] -= j.value * g(j.idx);
            rb[1] += j.value * g(j.idx);
        } else if (j.idx < 2 * ptsnum) {
            rb[0] -= j.value * f(j.idx - ptsnum);
            rb[1] += j.value * f(j.idx - ptsnum);
        } else {
            printError("AGM::Point::updateRb_t_interface");
        }
    }
}

void AGM::Point::calculate_finite_difference(char i) {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];

    if (i == 'x') {
        value["dx"] =
                (pts->at(element[E]->getIdx()).getValue()["sol"] - pts->at(element[W]->getIdx()).getValue()["sol"]) /
                (xp - xm);
    } else if (i == 'y') {
        value["dy"] =
                (pts->at(element[N]->getIdx()).getValue()["sol"] - pts->at(element[S]->getIdx()).getValue()["sol"]) /
                (yp - ym);
    } else {
        printError("AGM::Point::calculate_finite_difference", "i = %c", i);
    }
}

void AGM::Point::calculate_cross_difference() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];

    double value0 = (pts->at(element[E]->getIdx()).getValue()["dy"] - pts->at(element[W]->getIdx()).getValue()["dy"]) /
                    (xp - xm);
    double value1 = (pts->at(element[N]->getIdx()).getValue()["dx"] - pts->at(element[S]->getIdx()).getValue()["dx"]) /
                    (yp - ym);

    value["dxy"] = HALFVALUE * (value0 + value1);
}

double AGM::Point::integrate() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    return value["rhs"] * (xp - xm) * (yp - ym) * HALFVALUE * HALFVALUE;
//    return value["rhs"];
}

void AGM::Point::translate_rhs(double d) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    value["rhs"] -= d;
}

void AGM::Point::calcContinuity(std::vector<Point> *uvel, std::vector<Point> *vvel) {
    switch (condition) {
        case 'C':
            calcContinuity_cross(uvel, vvel);
            break;
        case 'D':
            if (axialLine[0]) {
                if (!iszero((*this) - *(element[E]))) {
                    value["sol"] = element[E]->getValue()["sol"];
                } else if (!iszero((*this) - *(element[W]))) {
                    value["sol"] = element[W]->getValue()["sol"];
                } else {
                    printError("AGM::Point::calcContinuity");
                }
            }
            if (axialLine[1]) {
                if (!iszero((*this) - *(element[N]))) {
                    value["sol"] = element[N]->getValue()["sol"];
                } else if (!iszero((*this) - *(element[S]))) {
                    value["sol"] = element[S]->getValue()["sol"];
                } else {
                    printError("AGM::Point::calcContinuity");
                }
            }
//            calcContinuity_boundary(uvel, vvel);
            break;
        case 'I':
            calcContinuity_interface(uvel, vvel);
            break;
        case 'N':
            calcContinuity_boundary(uvel, vvel);
            break;
        default:
            printError("AGM::Point::calcContinuity", "condition (which is %c) is wrong.", condition);
    }
}

void AGM::Point::calcContinuity_cross(std::vector<Point> *uvel, std::vector<Point> *vvel) {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

    auto u = std::array<double, 4>(), v = std::array<double, 4>(), p = std::array<double, 4>();
    auto phi = std::array<double, 4>(), psi = std::array<double, 4>();
    auto f1 = std::array<double, 4>(), f2 = std::array<double, 4>();

    auto assignVal = [&](EWNS e) -> void {
        u[e] = uvel->at(idx).element[e]->getValue()["sol"];
        v[e] = vvel->at(idx).element[e]->getValue()["sol"];
        phi[e] = uvel->at(idx).element[e]->getValue()["phi"];
        psi[e] = vvel->at(idx).element[e]->getValue()["phi"];
        f1[e] = uvel->at(idx).element[e]->getValue()["rhs"];
        f2[e] = vvel->at(idx).element[e]->getValue()["rhs"];
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

void AGM::Point::calcContinuity_boundary(std::vector<Point> *uvel, std::vector<Point> *vvel) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionLinear(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionLinear(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = GreenfunctionLinear(tm, tb, tp, d, d);
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
            u[e] = uvel->at(idx).element[e]->getValue()["sol"];
            v[e] = vvel->at(idx).element[e]->getValue()["sol"];
            phi[e] = uvel->at(idx).element[e]->getValue()["phi"];
            psi[e] = vvel->at(idx).element[e]->getValue()["phi"];
            f1[e] = uvel->at(idx).element[e]->getValue()["rhs"];
            f2[e] = vvel->at(idx).element[e]->getValue()["rhs"];
            p[e] = element[e]->getValue()["sol"];
        } else {
            u[e] = assignSol(uvel->at(idx).element[e0], uvel->at(idx).element[e1], j, uvel, 0);
            v[e] = assignSol(vvel->at(idx).element[e0], vvel->at(idx).element[e1], j, vvel, 1);
            phi[e] = approximate(uvel->at(idx).element[e0], uvel->at(idx).element[e1], "phi", i);
            psi[e] = approximate(vvel->at(idx).element[e0], vvel->at(idx).element[e1], "phi", i);
            f1[e] = approximate(uvel->at(idx).element[e0], uvel->at(idx).element[e1], "rhs", i);
            f2[e] = approximate(vvel->at(idx).element[e0], vvel->at(idx).element[e1], "rhs", i);
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

void AGM::Point::calcContinuity_interface(std::vector<Point> *uvel, std::vector<Point> *vvel) {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

    auto approximateSol = [&](Point *ptl, Point *ptr, int i, double d) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        double sign = i ? -UNITVALUE : UNITVALUE;
        auto func = Greenfunction(tm, tb, tp, d, d);
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
            u[e] = uvel->at(idx).element[e]->getValue()["sol"];
            v[e] = vvel->at(idx).element[e]->getValue()["sol"];
            phi[e] = uvel->at(idx).element[e]->getValue()["phi"];
            psi[e] = vvel->at(idx).element[e]->getValue()["phi"];
            f1[e] = uvel->at(idx).element[e]->getValue()["rhs"];
            f2[e] = vvel->at(idx).element[e]->getValue()["rhs"];
            p[e] = element[e]->getValue()["sol"];
        } else {
            u[e] = assignSol(uvel->at(idx).element[e0], uvel->at(idx).element[e1], j, uvel, 0);
            v[e] = assignSol(vvel->at(idx).element[e0], vvel->at(idx).element[e1], j, vvel, 1);
            phi[e] = approximate(uvel->at(idx).element[e0], uvel->at(idx).element[e1], "phi", i);
            psi[e] = approximate(vvel->at(idx).element[e0], vvel->at(idx).element[e1], "phi", i);
            f1[e] = approximate(uvel->at(idx).element[e0], uvel->at(idx).element[e1], "rhs", i);
            f2[e] = approximate(vvel->at(idx).element[e0], vvel->at(idx).element[e1], "rhs", i);
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

void AGM::Point::makePressureTerm() {
    switch (condition) {
        case 'C':
            makePressureTerm_cross();
            break;
        case 'D':
            break;
        case 'N':
            makePressureTerm_neumann();
            break;
        case 'I':
            makePressureTerm_interface();
            break;
        default:
            printError("AGM::Point::makePressureTerm", "condition (which is %c) is wrong.", condition);
    }
}

void AGM::Point::makePressureTerm_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

    pMatrixRow[0][element[W]->idx] = gfunc_x.green_integral_t('l');
    pMatrixRow[0][idx] = gfunc_x.green_integral_t('c');
    pMatrixRow[0][element[E]->idx] = gfunc_x.green_integral_t('r');

    pMatrixRow[1][element[S]->idx] = gfunc_y.green_integral_t('l');
    pMatrixRow[1][idx] = gfunc_y.green_integral_t('c');
    pMatrixRow[1][element[N]->idx] = gfunc_y.green_integral_t('r');
}

void AGM::Point::makePressureTerm_neumann() {
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
            pMatrixRow[i][pt->idx] += func->green_integral_ttau(c);
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

void AGM::Point::makePressureTerm_interface() {
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
        auto func = Greenfunction(tm, tb, tp, d, d);
        auto mRow = matrix_row();
        mRow[ptl->getIdx() + ptsnum] = func.green_integral_t('L');
        mRow[ptr->getIdx() + ptsnum] = func.green_integral_t('R');

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
            pMatrixRow[i][pt->idx] += func->green_integral_t(c);
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

void AGM::Point::calcPressureTerm(std::vector<Point> *pressure, char i) {
    switch (condition) {
        case 'C':
            calcPressureTerm_cross(pressure, i);
            break;
        case 'D':
            break;
        case 'N':
            calcPressureTerm_neumman(pressure, i);
            break;
        case 'I':
            calcPressureTerm_interface(pressure, i);
            break;
        default:
            printError("AGM::Point::calcPressureTerm", "condition (whidch is %d) is wrong.", condition);
    }
}

void AGM::Point::calcPressureTerm_cross(std::vector<Point> *pressure, char i) {
    if (i == 'u') {
        for (const auto &item : pMatrixRow[0]) {
            rb[0] -= item.value * pressure->at(item.idx)["sol"];
            rb[1] -= item.value * pressure->at(item.idx)["sol"];
        }
    } else if (i == 'v') {
        for (const auto &item : pMatrixRow[1]) {
            rb[0] -= item.value * pressure->at(item.idx)["sol"];
            rb[1] += item.value * pressure->at(item.idx)["sol"];
        }
    } else {
        printError("AGM::Point::calcPressureTerm_cross", "i (which is %c) is wrong.", i);
    }
}

void AGM::Point::calcPressureTerm_neumman(std::vector<Point> *pressure, char i) {
    if (i == 'u') {
        for (const auto &item : pMatrixRow[0]) {
            rb[0] -= item.value * pressure->at(item.idx)["sol"] * normal[0];
        }
    } else if (i == 'v') {
        for (const auto &item : pMatrixRow[1]) {
            rb[0] -= item.value * pressure->at(item.idx)["sol"] * normal[1];
        }
    } else {
        printError("AGM::Point::calcPressureTerm_neumman", "i (which is %c) is wrong.", i);
    }
}

void AGM::Point::calcPressureTerm_interface(std::vector<Point> *pressure, char i) {
    if (i == 'u') {
        for (const auto &item : pMatrixRow[0]) {
            rb[0] -= item.value * pressure->at(item.idx)["sol"];
            rb[1] -= item.value * pressure->at(item.idx)["sol"];
        }
    } else if (i == 'v') {
        for (const auto &item : pMatrixRow[1]) {
            rb[0] -= item.value * pressure->at(item.idx)["sol"];
            rb[1] += item.value * pressure->at(item.idx)["sol"];
        }
    } else {
        printError("AGM::Point::calcPressureTerm_interface", "i (which is %c) is wrong.", i);
    }
}

void AGM::Point::makePressureDifferenceTerm() {
    switch (condition) {
        case 'C':
            makePressrueDifferenceTerm_cross();
            break;
        case 'D':
            makePressrueDifferenceTerm_boundary();
            break;
        case 'I':
            makePressrueDifferenceTerm_interface();
            break;
        case 'N':
            makePressrueDifferenceTerm_boundary();
            break;
        default:
            printError("AGM::Point::makePressureDifferenceTerm", "condition (which is %c) is wrong.", condition);
    }

}

void AGM::Point::makePressrueDifferenceTerm_cross() {
    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

    pdMatrixRow[0][element[W]->idx] = gfunc_x.green_integral_ttau('l');
    pdMatrixRow[0][idx] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    pdMatrixRow[0][element[E]->idx] = gfunc_x.green_integral_ttau('r');

    pdMatrixRow[1][element[S]->idx] = gfunc_y.green_integral_ttau('l');
    pdMatrixRow[1][idx] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    pdMatrixRow[1][element[N]->idx] = gfunc_y.green_integral_ttau('r');

}

void AGM::Point::makePressrueDifferenceTerm_boundary() {
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

    if (axialLine[0]) {
        findPts(E, W, WN, 0, &xm, &xb, &xp, &nx);
        findPts(W, E, EN, 0, &xp, &xb, &xm, &nx);
    } else if (axialLine[1]) {
        findPts(N, S, SE, 1, &ym, &yb, &yp, &ny);
        findPts(S, N, NE, 1, &yp, &yb, &ym, &ny);
    } else {
        printError("AGM::Point::makePressrueDifferenceTerm_boundary()");
    }

    auto gfunc_x = GreenfunctionNeumann(xm, xb, xp, mp, mp);
    auto gfunc_y = GreenfunctionNeumann(ym, yb, yp, mp, mp);

    std::array<matrix_row, 2> row{};

    if (nx == E) {
        row[0][idx] = gfunc_x.green_integral_t_ND('r');
        row[0][element[W]->getIdx()] = gfunc_x.green_integral_t_ND('c');
        row[0][element1[W]->getIdx()] = gfunc_x.green_integral_t_ND('l');

        if (element1[W]->getIdx() > idx) {
            printError("AGM::Point::makePressrueDifferenceTerm_boundary()",
                       "element1[W] (which is %d) is greater than this (which is %d)", element1[W]->getIdx(), idx);
        }
        row[0] += element1[W]->getPdMatrixRow()[0] * -mp * gfunc_x.green_function_ND(xm);
    } else if (nx == W) {
        row[0][idx] = gfunc_x.green_integral_t_DN('l');
        row[0][element[E]->getIdx()] = gfunc_x.green_integral_t_DN('c');
        row[0][element1[E]->getIdx()] = gfunc_x.green_integral_t_DN('r');

        if (element1[E]->getIdx() > idx) {
            printError("AGM::Point::makePressrueDifferenceTerm_boundary()",
                       "element1[E] (which is %d) is greater than this (which is %d)", element1[E]->getIdx(), idx);
        }
        row[0] += element1[E]->getPdMatrixRow()[0] * mp * gfunc_x.green_function_DN(xp);
    } else if (ny == N) {
        row[1][idx] = gfunc_y.green_integral_t_ND('r');
        row[1][element[S]->getIdx()] = gfunc_y.green_integral_t_ND('c');
        row[1][element1[S]->getIdx()] = gfunc_y.green_integral_t_ND('l');

        if (element1[S]->getIdx() > idx) {
            printError("AGM::Point::makePressrueDifferenceTerm_boundary()",
                       "element1[S] (which is %d) is greater than this (which is %d)", element1[S]->getIdx(), idx);
        }
        row[1] += element1[S]->getPdMatrixRow()[1] * -mp * gfunc_y.green_function_ND(ym);
    } else if (ny == S) {
        row[1][idx] = gfunc_y.green_integral_t_DN('l');
        row[1][element[N]->getIdx()] = gfunc_y.green_integral_t_DN('c');
        row[1][element1[N]->getIdx()] = gfunc_y.green_integral_t_DN('r');

        if (element1[N]->getIdx() > idx) {
            printError("AGM::Point::makePressrueDifferenceTerm_boundary()",
                       "element1[N] (which is %d) is greater than this (which is %d)", element1[N]->getIdx(), idx);
        }
        row[1] += element1[N]->getPdMatrixRow()[1] * mp * gfunc_y.green_function_DN(yp);
    } else {
        printError("AGM::Point::makePressrueDifferenceTerm_boundary()");
    }

    pdMatrixRow[0] = row[0];
    pdMatrixRow[1] = row[1];
}

void AGM::Point::makePressrueDifferenceTerm_interface() {
    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];

    auto gfunc_x = Greenfunction(xm, xb, xp, mp, mp);
    auto gfunc_y = Greenfunction(ym, yb, yp, mp, mp);

    auto approximateP = [&](Point *ptl, Point *ptr, double coefficient, int i) -> matrix_row {
        double tm{ptl->getCoordinate()[i]}, tb{coordinate[i]}, tp{ptr->getCoordinate()[i]};
        auto mRow = matrix_row();

        mRow[ptl->getIdx()] = (tp - tb) / (tp - tm);
        mRow[ptr->getIdx()] = (tb - tm) / (tp - tm);

        return mRow * coefficient;
    };

    auto assignMatrix = [&](Point *pt, Point *ptl, Point *ptr, Greenfunction *func, int i, int i0, char c) -> void {
        if (pt) {
            pdMatrixRow[i][pt->idx] += func->green_integral_ttau(c);
        } else {
            pdMatrixRow[i] += approximateP(ptl, ptr, func->green_integral_ttau(c), i0);
        }
    };

    pdMatrixRow[0][idx] = gfunc_x.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[E], element[ES], element[EN], &gfunc_x, 0, 1, 'r');
    assignMatrix(element[W], element[WS], element[WN], &gfunc_x, 0, 1, 'l');

    pdMatrixRow[1][idx] = gfunc_y.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(element[N], element[NW], element[NE], &gfunc_y, 1, 0, 'r');
    assignMatrix(element[S], element[SW], element[SE], &gfunc_y, 1, 0, 'l');
}

void AGM::Point::calcContinuity(std::vector<Point *> *uvel, std::vector<Point *> *vvel,
                                std::vector<Coordinate> *convections) {

    switch (condition) {
        case 'C':
            calcContinuity_cross(uvel, vvel, convections);
            break;
        case 'D':
            if (axialLine[0]) {
                if (!iszero((*this) - *(element[E]))) {
                    value["sol"] = element[E]->getValue()["sol"];
                } else if (!iszero((*this) - *(element[W]))) {
                    value["sol"] = element[W]->getValue()["sol"];
                } else {
                    printError("AGM::Point::calcContinuity");
                }
            }
            if (axialLine[1]) {
                if (!iszero((*this) - *(element[N]))) {
                    value["sol"] = element[N]->getValue()["sol"];
                } else if (!iszero((*this) - *(element[S]))) {
                    value["sol"] = element[S]->getValue()["sol"];
                } else {
                    printError("AGM::Point::calcContinuity");
                }
            }
//            calcContinuity_boundary(uvel, vvel, convections);
            break;
        case 'I':
            calcContinuity_interface(uvel, vvel, convections);
            break;
        case 'N':
            calcContinuity_boundary(uvel, vvel, convections);
            break;
        default:
            printError("AGM::Point::calcContinuity", "condition (which is %c) is wrong.", condition);
    }
}

void AGM::Point::calcContinuity_cross(std::vector<Point *> *uvel, std::vector<Point *> *vvel,
                                      std::vector<Coordinate> *convections) {

    double xm = element[W]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E]->getCoordinate()[0];
    double ym = element[S]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convections->at(idx)[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convections->at(idx)[1]);

    auto u = std::array<double, 4>(), v = std::array<double, 4>(), p = std::array<double, 4>();
    auto phi = std::array<double, 4>(), psi = std::array<double, 4>();
    auto f1 = std::array<double, 4>(), f2 = std::array<double, 4>();

    auto assignVal = [&](EWNS e) -> void {
        u[e] = uvel->at(idx)->element[e]->getValue()["sol"];
        v[e] = vvel->at(idx)->element[e]->getValue()["sol"];
        phi[e] = uvel->at(idx)->element[e]->getValue()["phi"];
        psi[e] = vvel->at(idx)->element[e]->getValue()["phi"];
        f1[e] = uvel->at(idx)->element[e]->getValue()["rhs"];
        f2[e] = vvel->at(idx)->element[e]->getValue()["rhs"];
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
                    gfunc_x.green_integral_tau('c') * uvel->at(idx)->getValue()["phi"] +
                    gfunc_x.green_integral_tau('r') * phi[E] -
                    gfunc_y.green_integral_tau('l') * psi[S] -
                    gfunc_y.green_integral_tau('c') * vvel->at(idx)->getValue()["phi"] -
                    gfunc_y.green_integral_tau('r') * psi[N] +
                    HALFVALUE * gfunc_x.green_integral_tau('l') * f1[W] +
                    HALFVALUE * gfunc_x.green_integral_tau('c') * uvel->at(idx)->getValue()["rhs"] +
                    HALFVALUE * gfunc_x.green_integral_tau('r') * f1[E] +
                    HALFVALUE * gfunc_y.green_integral_tau('l') * f2[S] +
                    HALFVALUE * gfunc_y.green_integral_tau('c') * vvel->at(idx)->getValue()["rhs"] +
                    HALFVALUE * gfunc_y.green_integral_tau('r') * f2[N] +
                    gfunc_x.green_integral_ttau('l') * p[W] +
                    gfunc_x.green_integral_ttau('c') * value["sol"] +
                    gfunc_x.green_integral_ttau('r') * p[E] +
                    gfunc_y.green_integral_ttau('l') * p[S] +
                    gfunc_y.green_integral_ttau('c') * value["sol"] +
                    gfunc_y.green_integral_ttau('r') * p[N]);
}

void AGM::Point::calcContinuity_boundary(std::vector<Point *> *uvel, std::vector<Point *> *vvel,
                                         std::vector<Coordinate> *convections) {

    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusionLinear(xm, xb, xp, mp, mp, convections->at(idx)[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusionLinear(ym, yb, yp, mp, mp, convections->at(idx)[1]);

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

    auto assignSol = [&](Point *ptl, Point *ptr, int i0, std::vector<Point *> *vector, int i) -> double {
        auto row = approximateSol(ptl, ptr, i0, mp);
        double sol{};
        for (const auto &item : row) {
            if (item.idx < ptsnum) {
                sol += item.value * vector->at(item.idx)->getValue()["sol"];
            } else if (item.idx < 2 * ptsnum) {
                sol += item.value * vector->at(item.idx - ptsnum)->getValue()["phi"];
            } else if (item.idx < 3 * ptsnum) {
                sol += item.value * vector->at(item.idx - 2 * ptsnum)->getValue()["rhs"];
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
            u[e] = uvel->at(idx)->element[e]->getValue()["sol"];
            v[e] = vvel->at(idx)->element[e]->getValue()["sol"];
            phi[e] = uvel->at(idx)->element[e]->getValue()["phi"];
            psi[e] = vvel->at(idx)->element[e]->getValue()["phi"];
            f1[e] = uvel->at(idx)->element[e]->getValue()["rhs"];
            f2[e] = vvel->at(idx)->element[e]->getValue()["rhs"];
            p[e] = element[e]->getValue()["sol"];
        } else {
            u[e] = assignSol(uvel->at(idx)->element[e0], uvel->at(idx)->element[e1], j, uvel, 0);
            v[e] = assignSol(vvel->at(idx)->element[e0], vvel->at(idx)->element[e1], j, vvel, 1);
            phi[e] = approximate(uvel->at(idx)->element[e0], uvel->at(idx)->element[e1], "phi", i);
            psi[e] = approximate(vvel->at(idx)->element[e0], vvel->at(idx)->element[e1], "phi", i);
            f1[e] = approximate(uvel->at(idx)->element[e0], uvel->at(idx)->element[e1], "rhs", i);
            f2[e] = approximate(vvel->at(idx)->element[e0], vvel->at(idx)->element[e1], "rhs", i);
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
                    gfunc_x.green_integral_tau('c') * uvel->at(idx)->getValue()["phi"] +
                    gfunc_x.green_integral_tau('r') * phi[E] -
                    gfunc_y.green_integral_tau('l') * psi[S] -
                    gfunc_y.green_integral_tau('c') * vvel->at(idx)->getValue()["phi"] -
                    gfunc_y.green_integral_tau('r') * psi[N] +
                    HALFVALUE * gfunc_x.green_integral_tau('l') * f1[W] +
                    HALFVALUE * gfunc_x.green_integral_tau('c') * uvel->at(idx)->getValue()["rhs"] +
                    HALFVALUE * gfunc_x.green_integral_tau('r') * f1[E] +
                    HALFVALUE * gfunc_y.green_integral_tau('l') * f2[S] +
                    HALFVALUE * gfunc_y.green_integral_tau('c') * vvel->at(idx)->getValue()["rhs"] +
                    HALFVALUE * gfunc_y.green_integral_tau('r') * f2[N] +
                    gfunc_x.green_integral_ttau('l') * p[W] +
                    gfunc_x.green_integral_ttau('c') * value["sol"] +
                    gfunc_x.green_integral_ttau('r') * p[E] +
                    gfunc_y.green_integral_ttau('l') * p[S] +
                    gfunc_y.green_integral_ttau('c') * value["sol"] +
                    gfunc_y.green_integral_ttau('r') * p[N]);
}

void AGM::Point::calcContinuity_interface(std::vector<Point *> *uvel, std::vector<Point *> *vvel,
                                          std::vector<Coordinate> *convections) {

    double xm = element[W] ? element[W]->getCoordinate()[0] : element[WN]->getCoordinate()[0];
    double xb = getCoordinate()[0];
    double xp = element[E] ? element[E]->getCoordinate()[0] : element[EN]->getCoordinate()[0];
    double ym = element[S] ? element[S]->getCoordinate()[1] : element[SE]->getCoordinate()[1];
    double yb = getCoordinate()[1];
    double yp = element[N] ? element[N]->getCoordinate()[1] : element[NE]->getCoordinate()[1];
    auto gfunc_x = GreenfunctionConvectionDiffusion(xm, xb, xp, mp, mp, convections->at(idx)[0]);
    auto gfunc_y = GreenfunctionConvectionDiffusion(ym, yb, yp, mp, mp, convections->at(idx)[1]);

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

    auto assignSol = [&](Point *ptl, Point *ptr, int i0, std::vector<Point *> *vector, int i) -> double {
        auto row = approximateSol(ptl, ptr, i0, mp);
        double sol{};
        for (const auto &item : row) {
            if (item.idx < ptsnum) {
                sol += item.value * vector->at(item.idx)->getValue()["sol"];
            } else if (item.idx < 2 * ptsnum) {
                sol += item.value * vector->at(item.idx - ptsnum)->getValue()["phi"];
            } else if (item.idx < 3 * ptsnum) {
                sol += item.value * vector->at(item.idx - 2 * ptsnum)->getValue()["rhs"];
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
            u[e] = uvel->at(idx)->element[e]->getValue()["sol"];
            v[e] = vvel->at(idx)->element[e]->getValue()["sol"];
            phi[e] = uvel->at(idx)->element[e]->getValue()["phi"];
            psi[e] = vvel->at(idx)->element[e]->getValue()["phi"];
            f1[e] = uvel->at(idx)->element[e]->getValue()["rhs"];
            f2[e] = vvel->at(idx)->element[e]->getValue()["rhs"];
            p[e] = element[e]->getValue()["sol"];
        } else {
            u[e] = assignSol(uvel->at(idx)->element[e0], uvel->at(idx)->element[e1], j, uvel, 0);
            v[e] = assignSol(vvel->at(idx)->element[e0], vvel->at(idx)->element[e1], j, vvel, 1);
            phi[e] = approximate(uvel->at(idx)->element[e0], uvel->at(idx)->element[e1], "phi", i);
            psi[e] = approximate(vvel->at(idx)->element[e0], vvel->at(idx)->element[e1], "phi", i);
            f1[e] = approximate(uvel->at(idx)->element[e0], uvel->at(idx)->element[e1], "rhs", i);
            f2[e] = approximate(vvel->at(idx)->element[e0], vvel->at(idx)->element[e1], "rhs", i);
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
                    gfunc_x.green_integral_tau('c') * uvel->at(idx)->getValue()["phi"] +
                    gfunc_x.green_integral_tau('r') * phi[E] -
                    gfunc_y.green_integral_tau('l') * psi[S] -
                    gfunc_y.green_integral_tau('c') * vvel->at(idx)->getValue()["phi"] -
                    gfunc_y.green_integral_tau('r') * psi[N] +
                    HALFVALUE * gfunc_x.green_integral_tau('l') * f1[W] +
                    HALFVALUE * gfunc_x.green_integral_tau('c') * uvel->at(idx)->getValue()["rhs"] +
                    HALFVALUE * gfunc_x.green_integral_tau('r') * f1[E] +
                    HALFVALUE * gfunc_y.green_integral_tau('l') * f2[S] +
                    HALFVALUE * gfunc_y.green_integral_tau('c') * vvel->at(idx)->getValue()["rhs"] +
                    HALFVALUE * gfunc_y.green_integral_tau('r') * f2[N] +
                    gfunc_x.green_integral_ttau('l') * p[W] +
                    gfunc_x.green_integral_ttau('c') * value["sol"] +
                    gfunc_x.green_integral_ttau('r') * p[E] +
                    gfunc_y.green_integral_ttau('l') * p[S] +
                    gfunc_y.green_integral_ttau('c') * value["sol"] +
                    gfunc_y.green_integral_ttau('r') * p[N]);
}

void AGM::Point::printInformation() {
    std::string ordinal{}, bd{};
    int ordinalIdx = idx % 10;
    if (ordinalIdx == 1) ordinal = "st";
    else if (ordinalIdx == 2) ordinal = "nd";
    else if (ordinalIdx == 3) ordinal = "rd";
    else ordinal = "th";

    if (condition == 'C') bd = "Cross point";
    else if (condition == 'D') bd = "Dirichlet boundary";
    else if (condition == 'N') bd = "Neumann boundary";
    else if (condition == 'I') bd = "Interface point";

    std::cout.precision(16);
    std::cout << std::scientific;
    std::cout << idx << ordinal << " point information:" << '\n';
    std::cout << '\t' << "(x, y) = (" << coordinate[0] << ", " << coordinate[1] << ")\n";
    std::cout << '\t' << "boundary condition = " << bd << '\n';
    if (condition == 'D' || condition == 'N') {
        std::cout << '\t' << "boundary value = " << value["bdv"] << '\n';
    }
    std::cout << '\t' << "Right-hand side = " << value["rhs"] << '\n';
    if (condition == 'N') {
        std::cout << '\t' << "normal vector = (" << normal[0] << ", " << normal[1] << ")" << '\n';
    }
    std::array<std::string, 12> elementName{"East", "West", "North", "South", "East-North", "East-South", "West-North",
                                            "West-South", "North-East", "North-West", "South-East", "South-West"};
    for (int i = 0; i < 12; ++i) {
        std::cout << '\t' << elementName[i] << ":\n";
        if (element[i]) {
            std::cout << "\t\t" << "(x, y) = (" << element[i]->getCoordinate()[0] << ", "
                      << element[i]->getCoordinate()[1] << "), condition = " << element[i]->getCondition()
                      << ", distance = " << (*element[i]) - *this << '\n';
        } else {
            std::cout << "\t\tIt does not exist\n";
        }
    }
}

void AGM::Point::roundMatrix() {

    auto roundValue = [&](double d) -> double {
        char char0[1024], char1[1024];
        sprintf(char0, "%23.16e", d);
        auto string0 = std::string(char0);
        double double0 = std::round(std::stod(string0.substr(1, 18)) * 1e11);
        sprintf(char1, "%23.16e", double0);
        auto string1 = std::string(char1);
        return std::stod(string0[0] + string1.substr(1, 18) + "E" + string0.substr(20, 23));
    };

    for (auto &item : matrixRow) {
        for (auto it = item.begin(); it != item.end();) {
            if (iszero(it->value)) {
                it = item.erase(it);
            } else {
                ++it;
            }
        }
        for (auto &item1 : item) {
            item1.value = roundValue(item1.value);
        }
    }

    for (auto &item : rhsMatrixRow) {
        for (auto it = item.begin(); it != item.end();) {
            if (iszero(it->value)) {
                it = item.erase(it);
            } else {
                ++it;
            }
        }
        for (auto &item1 : item) {
            item1.value = roundValue(item1.value);
        }
    }

    for (auto &item : dMatrixRow) {
        for (auto it = item.begin(); it != item.end();) {
            if (iszero(it->value)) {
                it = item.erase(it);
            } else {
                ++it;
            }
        }
        for (auto &item1 : item) {
            item1.value = roundValue(item1.value);
        }
    }

    for (auto &item : pMatrixRow) {
        for (auto it = item.begin(); it != item.end();) {
            if (iszero(it->value)) {
                it = item.erase(it);
            } else {
                ++it;
            }
        }
        for (auto &item1 : item) {
            item1.value = roundValue(item1.value);
        }
    }

    if (condition == 'N' && isclose(coordinate[0], 0) && isclose(coordinate[1], 0.5)) {
        std::cout << '\n';
    }
}
