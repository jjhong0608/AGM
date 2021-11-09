//
// Created by NIMS-JUNHONG on 2021/01/04.
//

#include "WriteFile.h"

template<typename T>
AGM::WriteFile<T>::WriteFile() = default;

template<typename T>
AGM::WriteFile<T>::WriteFile(std::vector<T> *pts) : pts(pts) {}

template<typename T>
T *AGM::WriteFile<T>::getPt() const {
    return pt;
}

template<typename T>
std::vector<T> *AGM::WriteFile<T>::getPts() const {
    return pts;
}

template<typename T>
void AGM::WriteFile<T>::setPoint(T *pPoint) {
    pt = pPoint;
}

template<typename T>
void AGM::WriteFile<T>::setPoints(std::vector<T> *points) {
    pts = points;
}

template<typename T>
auto AGM::WriteFile<T>::calculateErrorAtPoint() {
    double xm = (*pt)[W] ? (*pt)[W]->getCoordinate()[0] : (*pt)[WN]->getCoordinate()[0];
    double xp = (*pt)[E] ? (*pt)[E]->getCoordinate()[0] : (*pt)[EN]->getCoordinate()[0];
    double ym = (*pt)[S] ? (*pt)[S]->getCoordinate()[1] : (*pt)[SE]->getCoordinate()[1];
    double yp = (*pt)[N] ? (*pt)[N]->getCoordinate()[1] : (*pt)[NE]->getCoordinate()[1];
    auto func = function(*pt);
    func.setTime(AGM::PointHeat::getTime() + AGM::PointHeat::getDelta());
    func.setEps(0.01);
    double numerator{pow(pt->getValue()["sol"] - func.u(), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
    double denominator{pow(func.u(), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
    return std::array<double, 2>{numerator, denominator};
}

template<typename T>
double AGM::WriteFile<T>::calculateError() {
    double numerator{}, denominator{};
    std::array<double, 2> error{};
    for (auto &i : *pts) {
        pt = &i;
        error = calculateErrorAtPoint();
        numerator += error[0];
        denominator += error[1];
    }

    return sqrt(numerator) / sqrt(denominator);
}

template<typename T>
std::array<double, 2> AGM::WriteFile<T>::calculateErrorAtPoint(const std::string &string) {
    double xm = (*pt)[W] ? (*pt)[W]->getCoordinate()[0] : (*pt)[WN]->getCoordinate()[0];
    double xp = (*pt)[E] ? (*pt)[E]->getCoordinate()[0] : (*pt)[EN]->getCoordinate()[0];
    double ym = (*pt)[S] ? (*pt)[S]->getCoordinate()[1] : (*pt)[SE]->getCoordinate()[1];
    double yp = (*pt)[N] ? (*pt)[N]->getCoordinate()[1] : (*pt)[NE]->getCoordinate()[1];
    auto func = function(*pt);
    func.setTime(AGM::PointHeat::getTime() + AGM::PointHeat::getDelta());
    auto f = [&]() -> double {
        if (string == "sol") return func.u();
        if (string == "phi") return func.phi();
        if (string == "dx") return func.dx();
        if (string == "dy") return func.dy();
        if (string == "grad") return func.grad();
        if (string == "dxx") return func.dxx();
        if (string == "dyy") return func.dyy();
        printError("WriteFile::calculateErrorAtPoint");
    };
    double value = string == "grad" ? sqrt(pow(pt->getValue()["dx"], 2) + pow(pt->getValue()["dy"], 2))
                                    : pt->getValue()[string];
    double numerator{pow(value - f(), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
    double denominator{pow(f(), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
    return std::array<double, 2>{numerator, denominator};
}

template<typename T>
double AGM::WriteFile<T>::calculateError(const std::string &string) {
    double numerator{}, denominator{};
    std::array<double, 2> error{};
    for (auto &i : *pts) {
        pt = &i;
        error = calculateErrorAtPoint(string);
        numerator += error[0];
        denominator += error[1];
    }
    return sqrt(numerator) / sqrt(denominator);
}

template<typename T>
auto AGM::WriteFile<T>::calculateIterationErrorAtPoint(AGM::Value *value) {
    double xm = (*pt)[W] ? (*pt)[W]->getCoordinate()[0] : (*pt)[WN]->getCoordinate()[0];
    double xp = (*pt)[E] ? (*pt)[E]->getCoordinate()[0] : (*pt)[EN]->getCoordinate()[0];
    double ym = (*pt)[S] ? (*pt)[S]->getCoordinate()[1] : (*pt)[SE]->getCoordinate()[1];
    double yp = (*pt)[N] ? (*pt)[N]->getCoordinate()[1] : (*pt)[NE]->getCoordinate()[1];
    double numerator{pow(value->getSol() - pt->getValue()["sol"], 2) * (xp - xm) * (yp - ym) * HALFVALUE * HALFVALUE};
    double denominator{pow(pt->getValue()["sol"], 2) * (xp - xm) * (yp - ym) * HALFVALUE * HALFVALUE};
    return std::array<double, 2>{numerator, denominator};
}

template<typename T>
double AGM::WriteFile<T>::calculateIterationError(std::vector<Value> *values) {
    double numerator{}, denominator{};
    std::array<double ,2> error{};
    auto iter = values->begin();
    for (auto &item : *pts) {
        pt = &item;
        error = calculateIterationErrorAtPoint(iter.base());
        numerator += error[0];
        denominator += error[1];
        ++iter;
    }
    return sqrt(numerator) / sqrt(denominator);
}

template<typename T>
void AGM::WriteFile<T>::writeResult(const std::string &string) {
    std::ofstream f(string);
    f.precision(16);

    std::for_each(pts->begin(), pts->end(), [&](Point &point) {
        f << std::scientific;
        f << point.getIdx() << '\t';
        f << point[0] << '\t';
        f << point[1] << '\t';
        f << point["sol"] << '\t';
        f << point["phi"] << '\t';
        f << point["dx"] << '\t';
        f << point["dy"] << '\t';
        f << point["dxx"] << '\t';
        f << point["dyy"] << std::endl;
    });
    f.close();
}

template class AGM::WriteFile<AGM::Point>;
template class AGM::WriteFile<AGM::PointHeat>;
template class AGM::WriteFile<AGM::PointConvectionDiffusion>;
template class AGM::WriteFile<AGM::PointConvectionDiffusionT>;

template<typename T1, typename T2>
AGM::WriteFileNS<T1, T2>::WriteFileNS() = default;

template<typename T1, typename T2>
AGM::WriteFileNS<T1, T2>::WriteFileNS(std::vector<T1> *pts, std::vector<T2> *ptsU, std::vector<T2> *ptsV) :
        WriteFile<T1>(pts), ptsU(ptsU), ptsV(ptsV) {}

template<typename T1, typename T2>
void AGM::WriteFileNS<T1, T2>::writeResult(const std::string &string) {
    std::ofstream f(string);
    auto iterU = ptsU->begin(), iterV = ptsV->begin();
    f.precision(16);
    f << std::scientific;
    std::for_each(WriteFile<T1>::pts->begin(), WriteFile<T1>::pts->end(), [&](Point &point) {
        f << point.getIdx() << '\t';
        f << point[0] << '\t';
        f << point[1] << '\t';
        f << (*iterU)["sol"] << '\t';
        f << (*iterV)["sol"] << '\t';
        f << point["sol"] << '\t';
        f << (*iterU)["dx"] << '\t';
        f << (*iterU)["dy"] << '\t';
        f << (*iterV)["dx"] << '\t';
        f << (*iterV)["dy"] << '\t';
        f << (*iterU)["phi"] << '\t';
        f << (*iterV)["phi"] << '\t';
        f << point["dx"] << '\t';
        f << point["dy"] << '\n';
        ++iterU;
        ++iterV;
    });
    f.close();
}

template<typename T1, typename T2>
std::array<double, 2> AGM::WriteFileNS<T1, T2>::calculateErrorAtPoint(const std::string &string) {
    double xm = (*ptU)[W] ? (*ptU)[W]->getCoordinate()[0] : (*ptU)[WN]->getCoordinate()[0];
    double xp = (*ptU)[E] ? (*ptU)[E]->getCoordinate()[0] : (*ptU)[EN]->getCoordinate()[0];
    double ym = (*ptU)[S] ? (*ptU)[S]->getCoordinate()[1] : (*ptU)[SE]->getCoordinate()[1];
    double yp = (*ptU)[N] ? (*ptU)[N]->getCoordinate()[1] : (*ptU)[NE]->getCoordinate()[1];
    auto func = function2D();
    func.setPoint(*ptU);
    func.setTime(AGM::PointHeat::getTime() + AGM::PointHeat::getDelta());
    auto f = [&]() -> double {
        if (string == "uvel") return func.u();
        if (string == "vvel") return func.v();
        if (string == "pressure") return func.p();
        if (string == "ux") return func.ux();
        if (string == "uy") return func.uy();
        if (string == "vx") return func.vx();
        if (string == "vy") return func.vy();
        printError("WriteFileNS::calculateErrorAtPoint");
    };
    auto g = [&]() -> double {
        if (string == "uvel") return ptU->getValue()["sol"];
        if (string == "vvel") return ptV->getValue()["sol"];
        if (string == "pressure") return this->getPt()->getValue()["sol"];
        if (string == "ux") return ptU->getValue()["dx"];
        if (string == "uy") return ptU->getValue()["dy"];
        if (string == "vx") return ptV->getValue()["dx"];
        if (string == "vy") return ptV->getValue()["dy"];
        printError("WriteFileNS::calculateErrorAtPoint");
    };
    double numerator{pow(g() - f(), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
    double denominator{pow(f(), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
    return std::array<double, 2>{numerator, denominator};
}

template<typename T1, typename T2>
double AGM::WriteFileNS<T1, T2>::calculateError(const std::string &string) {
    double numerator{}, denominator{};
    std::array<double, 2> error{};
    for (int j = 0; j < ptsU->size(); ++j) {
        ptU = &(ptsU->at(j));
        ptV = &(ptsV->at(j));
        this->setPoint(&(this->getPts()->at(j)));
        error = calculateErrorAtPoint(string);
        numerator += error[0];
        denominator += error[1];
    }
    return sqrt(numerator) / sqrt(denominator);
}

template class AGM::WriteFileNS<AGM::Point, AGM::Point>;
template class AGM::WriteFileNS<AGM::Point, AGM::PointHeat>;
template class AGM::WriteFileNS<AGM::Point, AGM::PointConvectionDiffusion>;
template class AGM::WriteFileNS<AGM::Point, AGM::PointConvectionDiffusionT>;
