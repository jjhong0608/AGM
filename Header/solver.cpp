//
// Created by NIMS-JUNHONG on 2021/01/04.
//

#include <vector>
#include "solver.h"
#include "PointConvectionDiffusion.h"

AGM::solver::solver(std::vector<Point> *pts) : pts(pts) {}

AGM::solver::solver(std::vector<AGM::Point> *pts, std::vector<AGM::AxialLine> *pXline,
                    std::vector<AGM::AxialLine> *pYline) : pts(pts), pXline(pXline), pYline(pYline) {}

AGM::solver::~solver() = default;

void AGM::solver::setPoints(std::vector<Point> *points) {
    pts = points;
}

std::vector<AGM::AxialLine> *AGM::solver::getPXline() const {
    return pXline;
}

void AGM::solver::setPXline(std::vector<AGM::AxialLine> *xline) {
    solver::pXline = xline;
}

std::vector<AGM::AxialLine> *AGM::solver::getPYline() const {
    return pYline;
}

void AGM::solver::setPYline(std::vector<AGM::AxialLine> *yline) {
    solver::pYline = yline;
}

void AGM::solver::ellipticSolver() {
    AGM::Point::setPtsnum(int(pts->size()));
    AGM::Point::setPts(pts);

    auto uvel = std::vector<Value>(Point::getPtsnum());
    auto vvel = std::vector<Value>(Point::getPtsnum());

    auto f = function();
    int fixed_idx{};

    auto rhs = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhs_x_t = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhs_y_t = [&](int i) -> double {
        return ZEROVALUE;
    };

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        f.setPoint(pt);
        f.setEps(pt.getMp());
        if (pt.getCondition() == 'D') {
            pt["bdv"] = f.u();
        }
        if (pt.getCondition() == 'N') {
            pt["bdv"] = f.dx() * pt.getNormal()[0] + f.dy() * pt.getNormal()[1];
        }
    });

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        pt.makeDifferentiation();
        pt.calcRepresentationFormula();
        pt.updateRb(rhs, rhs);
        pt.updateRb_t(rhs_x_t, rhs_y_t);
    });

    auto matrix = AGM::matrix<Point>(pts);
//    auto matrix = AGM::matrixEIGEN<Point>(pts, fixed_idx);
    matrix.makeMatrix();
    matrix.calculateMatrix();

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        pt.calcDifferentiation(rhs, rhs, rhs_x_t, rhs_y_t);
    });

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        f.setPoint(pt);
        f.setEps(pt.getMp());
        if (pt.getCondition() == 'D') {
            pt["dx"] = f.dx();
            pt["dy"] = f.dy();
        }
        if (pt.getCondition() == 'N') {
            pt["dx"] = f.dx();
            pt["dy"] = f.dy();
        }
    });
}

void AGM::solver::streamSolver() {
    AGM::Point::setPtsnum(int(pts->size()));
    AGM::Point::setPts(pts);

    auto uvel = std::vector<Value>(Point::getPtsnum());
    auto vvel = std::vector<Value>(Point::getPtsnum());
    auto rightbound = std::vector<double>();

    auto f = function();
    int fixed_idx{};

    auto rhs = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhs_x_t = [&](int i) -> double {
        return -vvel.at(i)["sol"];
    };
    auto rhs_y_t = [&](int i) -> double {
        return uvel.at(i)["sol"];
    };

    auto NS = std::ifstream("/home/jjhong0608/docker1/AGM_test/Navier-Stokes/BFS/Re800/AGM_Result");
    if (NS.fail()) {
        printError("file is not opened");
    }
    int idx{};
    double x{}, y{}, p{}, px{}, py{};

    int bdnum{};
    auto assignBoundary = [&](int i) -> double {
        double sum{};
        for (int j = 0; j < i + 1; ++j) {
            sum += rightbound.at(j);
        }
        return (sum - HALFVALUE * rightbound.at(i)) / 80.0;
    };

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        NS >> idx >> x >> y >> uvel.at(idx)["sol"] >> vvel.at(idx)["sol"] >> p >> uvel.at(idx)["dx"]
           >> uvel.at(idx)["dy"] >> vvel.at(idx)["dx"] >> vvel.at(idx)["dy"] >> uvel.at(idx)["phi"]
           >> vvel.at(idx)["phi"] >> px >> py;

        if (isclose(x, 1.5e1)) {
            rightbound.emplace_back(uvel.at(idx)["sol"]);
        }

        f.setPoint(pt);
        f.setEps(pt.getMp());
        if (pt.getCondition() == 'D') {
            pt["bdv"] = f.u();
        }
        if (pt.getCondition() == 'N') {
            pt.setCondition('D');
//            pt["bdv"] = f.u();
            pt["bdv"] = assignBoundary(bdnum++);
        }
        pt.setMp(UNITVALUE);
    });

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        pt.makeDifferentiation();
        pt.calcRepresentationFormula();
        pt.updateRb(rhs, rhs);
        pt.updateRb_t(rhs_x_t, rhs_y_t);
    });

    auto matrix = AGM::matrix<Point>(pts);
//    auto matrix = AGM::matrixEIGEN<Point>(pts, fixed_idx);
    matrix.makeMatrix();
    matrix.calculateMatrix();

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        pt.calcDifferentiation(rhs, rhs, rhs_x_t, rhs_y_t);
    });
}

void AGM::solver::heatSolver() {
    Point::setPtsnum(int(pts->size()));
    Point::setPts(pts);
    auto ptsH = std::vector<PointHeat>(Point::getPtsnum());
    PointHeat::setPtsH(&ptsH);
    PointHeat::setTime(ZEROVALUE);
    PointHeat::setDelta(0.01);
    auto previousValue = std::vector<Value>(Point::getPtsnum());
    PointHeat::setPreviousValues(&previousValue);
    PointHeat::setCn(UNITVALUE);
    PointHeat::setBdf(UNITVALUE);
    auto rhsHeat_x = [&](int i) -> double {
        return HALFVALUE * (ptsH.at(i)["rhs"] + previousValue.at(i)["rhs"]) +
               previousValue.at(i)["sol"] / PointHeat::getDelta() + previousValue.at(i)["phi"];
    };
    auto rhsHeat_y = [&](int i) -> double {
        return HALFVALUE * (ptsH.at(i)["rhs"] + previousValue.at(i)["rhs"]) +
               previousValue.at(i)["sol"] / PointHeat::getDelta() - previousValue.at(i)["phi"];
    };
    auto rhsHeat_t_x = [&](int i) -> double {
        return -ptsH.at(i).getMp() * previousValue.at(i)["dx"];
    };
    auto rhsHeat_t_y = [&](int i) -> double {
        return -ptsH.at(i).getMp() * previousValue.at(i)["dy"];
    };
    int nIter{};
    PointHeat::setRhsParallelToX(rhsHeat_x);
    PointHeat::setRhsParallelToY(rhsHeat_y);
    PointHeat::setBdf(UNITVALUE);

    auto assignInitial = [&]() -> void {
        auto f = function();
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            ptsH.at(j).Point::operator=(pts->at(j));
            f.setPoint(ptsH.at(j));
            f.setTime(PointHeat::getTime());
            previousValue.at(j)["sol"] = f.u();
            previousValue.at(j)["rhs"] = f.f();
            previousValue.at(j)["phi"] = f.phi();
            previousValue.at(j)["dx"] = f.dx();
            previousValue.at(j)["dy"] = f.dy();
            f.setTime(PointHeat::getTime() + PointHeat::getDelta());
            if (ptsH.at(j).getCondition() == 'D') {
                ptsH.at(j)["bdv"] = f.u();
            } else if (ptsH.at(j).getCondition() == 'N') {
                ptsH.at(j)["bdv"] = f.dx() * ptsH.at(j).getNormal()[0] + f.dy() * ptsH.at(j).getNormal()[1];
            }
            ptsH.at(j)["rhs"] = f.f();
        }
    };

    auto assignBoundaryvalue = [&]() -> void {
        auto f = function();
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            f.setPoint(ptsH.at(j));
            f.setTime(PointHeat::getTime() + PointHeat::getDelta());
            if (ptsH.at(j).getCondition() == 'D') {
                ptsH.at(j)["bdv"] = f.u();
            } else if (ptsH.at(j).getCondition() == 'N') {
                ptsH.at(j)["bdv"] = f.dx() * ptsH.at(j).getNormal()[0] + f.dy() * ptsH.at(j).getNormal()[1];
            }
            ptsH.at(j)["rhs"] = f.f();
        }
    };

    auto findElement = [&]() -> void {
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            ptsH.at(j).findElement(pts);
            ptsH.at(j).findElement1(pts);
        }
    };

    auto makeMatrix = [&]() -> void {
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            ptsH.at(j).calcRepresentationFormula();
            ptsH.at(j).makeDifferentiation();
            ptsH.at(j).updateRb(rhsHeat_x, rhsHeat_y);
            ptsH.at(j).updateRb_t(rhsHeat_t_x, rhsHeat_t_y);
        }
    };

    auto calculateDifferentiation = [&]() -> void {
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            ptsH.at(j).calcDifferentiation(rhsHeat_x, rhsHeat_y, rhsHeat_t_x, rhsHeat_t_y);
        }
    };

    auto updateRb = [&]() -> void {
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            ptsH.at(j).updateRb(rhsHeat_x, rhsHeat_y);
            ptsH.at(j).updateRb_t(rhsHeat_t_x, rhsHeat_t_y);
        }
    };

    auto updateValue = [&]() -> void {
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            previousValue.at(j) = ptsH.at(j).getValue();
        }
    };

    auto checkMatrix = [&]() -> void {
        double v0{}, v1{};
        int i0{};
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            v0 = ptsH.at(j).checkRepresentationFormula();
            std::cout << "idx = " << ptsH.at(j).getIdx() << ", v0 = " << v0  << ", condition = " << ptsH.at(j).getCondition() << "\n";
            if (fabs(v1) < fabs(v0)) {
                v1 = v0;
                i0 = ptsH.at(j).getIdx();
            }
        }
        std::cout << "matrix value = " << v1 << '\n';
        ptsH.at(i0).printInformation();
    };

    auto checkDifferentiation = [&]() -> void {
        auto f = function();
        double v0{}, v1{};
        int i0{};
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            v0 = ptsH.at(j).checkDifferentiation(rhsHeat_x, rhsHeat_y, rhsHeat_t_x, rhsHeat_t_y);
            f.setPoint(ptsH.at(j));
            f.setEps(ptsH.at(j).getMp());
            f.setTime(0.1);
            v0 = fabs(v0 - f.dx());
            if (fabs(v1) < fabs(v0)) {
                v1 = v0;
                i0 = j;
            }
        }
        std::cout << "max difference = " << v1 << '\n';
        ptsH.at(i0).printInformation();
    };

    auto checkCalculation = [&]() -> void {
        double v0{}, v1{};
        int i0{};
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            for (int k = 0; k < 2; ++k) {
                v0 = -ptsH.at(j).getRb()[k];
                for (const auto &item : ptsH.at(j).getMatrixRow()[k]) {
                    if (item.idx < Point::getPtsnum()) {
                        v0 += item.value * ptsH.at(item.idx)["sol"];
                    } else if (item.idx < Point::getPtsnum() * 2) {
                        v0 += item.value * ptsH.at(item.idx - Point::getPtsnum())["phi"];
                    }
                }
                if (fabs(v1) < fabs(v0)) {
                    v1 = v0;
                    i0 = ptsH.at(j).getIdx();
                }
            }
        }
        std::cout << "calculation value = " << v1 << '\n';
        ptsH.at(i0).printInformation();
    };

    auto updateTime = [&]() -> void {
        PointHeat::setTime(PointHeat::getTime() + PointHeat::getDelta());
        std::cout << "current time = " << PointHeat::getTime() << '\n';
        ++nIter;
    };

    assignInitial();
    findElement();
    makeMatrix();
    auto matrix = AGM::matrix<PointHeat>(&ptsH);
    matrix.makeMatrix();
    matrix.calculateMatrix();
    calculateDifferentiation();
//    checkDifferentiation();
//    checkMatrix();
    updateValue();
//    checkCalculation();

    auto wf = WriteFile<PointHeat>(&ptsH);

    double min_value{100.0}, vv{};


    while (PointHeat::getTime() + PointHeat::getDelta() < ZEROVALUE - HALFVALUE * PointHeat::getDelta()) {
        updateTime();
        assignBoundaryvalue();
        updateRb();
        matrix.calculateMatrix();
        calculateDifferentiation();
        updateValue();

        if (nIter % 1 == 0) {
            std::cout.precision(16);
            vv = wf.calculateError("sol");
            std::cout << std::setw(32) << std::left << "Relative error of solution = " << std::scientific
                      << vv << std::endl;
            if (min_value > vv) {
                min_value = vv;
            }
            std::cout << "min value = " << min_value << "\n";
//            wf.writeResult("AGM_Result_" + std::to_string(nIter));
//            system(("cp AGM_Result_" + std::to_string(nIter) + " ../docker/AGM_test/").c_str());
        }
    }
    checkDifferentiation();
    for (int j = 0; j < Point::getPtsnum(); ++j) {
        pts->at(j).setValue(ptsH.at(j).getValue());
    }
}

void AGM::solver::convectiondiffusionSolver() {
    auto ptsCD = std::vector<PointConvectionDiffusion>(pts->size());
    auto convections = std::vector<Coordinate>(pts->size());
    PointConvectionDiffusion::setConvections(&convections);
    PointConvectionDiffusion::setPtsnum(int(pts->size()));
    PointConvectionDiffusion::setPts(pts);
    PointConvectionDiffusion::setPtsCd(&ptsCD);
    auto iter = pts->begin();
    auto iterConvections = convections.begin();

    auto assignInitial = [&]() {
        auto f = function();
        auto f2 = function2D();
        f2.setTime(UNITVALUE);
        std::for_each(ptsCD.begin(), ptsCD.end(), [&](PointConvectionDiffusion &pt) -> void {
            pt.Point::operator=(*iter);
            f.setPoint(pt);
            f2.setPoint(pt);
            *iterConvections = Coordinate(f2.u(), f2.v());
            pt["bdv"] = f.u();
            pt["rhs"] = f.f();
            ++iterConvections, ++iter;
        });
        iterConvections = convections.begin(), iter = pts->begin();
    };

    auto findElement = [&]() {
        std::for_each(ptsCD.begin(), ptsCD.end(), [&](PointConvectionDiffusion &pt) -> void {
            pt.findElement(pts);
        });
    };

    auto makeMatrix = [&]() {
        std::for_each(ptsCD.begin(), ptsCD.end(), [&](PointConvectionDiffusion &pt) -> void {
            pt.setConvection(*iterConvections);
            pt.calcRepresentationFormula();
            pt.makeDifferentiation();
            ++iterConvections;
        });
        iterConvections = convections.begin();
    };

    auto calcDifferentiation = [&]() {
        std::for_each(ptsCD.begin(), ptsCD.end(), [&](PointConvectionDiffusion &pt) -> void {
            pt.calcDifferentiation();
        });
    };

    assignInitial();
    findElement();
    makeMatrix();

    auto matrix = AGM::matrix<PointConvectionDiffusion>(&ptsCD);
    matrix.makeMatrix();
    matrix.calculateMatrix();
    calcDifferentiation();

    for (int j = 0; j < pts->size(); ++j) {
        pts->at(j).setValue(ptsCD.at(j).getValue());
    }
}

void AGM::solver::unsteadyconvectiondiffusionSolver() {
    auto ptsCDT = std::vector<PointConvectionDiffusionT>(pts->size());
    auto previousValues = std::vector<Value>(pts->size());
    auto twoLevelsUpValues = std::vector<Value>(pts->size());
    auto convections = std::vector<Coordinate>(pts->size());
    PointConvectionDiffusionT::setConvections(&convections);
    PointConvectionDiffusionT::setPtsnum(int(pts->size()));
    PointConvectionDiffusionT::setPtsCdt(&ptsCDT);
    PointConvectionDiffusionT::setPreviousValues(&previousValues);
    PointConvectionDiffusionT::setTwoLevelsUpValues(&twoLevelsUpValues);

    PointConvectionDiffusionT::setTime(ZEROVALUE);
    PointConvectionDiffusionT::setDelta(0.01);
    PointConvectionDiffusionT::setCn(UNITVALUE);

    PointHeat::setRhsParallelToX([&](int i) -> double {
        return PointConvectionDiffusionT::getCn() * PointConvectionDiffusionT::getPtsCdt()->at(i)["rhs"] +
               2.0E0 * PointConvectionDiffusionT::getPreviousValues()->at(i)["sol"] /
               PointConvectionDiffusionT::getDelta() -
               HALFVALUE * PointConvectionDiffusionT::getTwoLevelsUpValues()->at(i)["sol"] /
               PointConvectionDiffusionT::getDelta();
    });
    PointConvectionDiffusionT::setRhsParallelToY([&](int i) -> double {
        return PointConvectionDiffusionT::getCn() * PointConvectionDiffusionT::getPtsCdt()->at(i)["rhs"] +
               2.0E0 * PointConvectionDiffusionT::getPreviousValues()->at(i)["sol"] /
               PointConvectionDiffusionT::getDelta() -
               HALFVALUE * PointConvectionDiffusionT::getTwoLevelsUpValues()->at(i)["sol"] /
               PointConvectionDiffusionT::getDelta();
    });

    auto previousIter = previousValues.begin();
    auto twoLevelsIter = twoLevelsUpValues.begin();
    auto iterConvections = convections.begin();
    auto iter = pts->begin();

    auto assignInitial = [&]() {
        auto f = function();
        auto f2 = function2D();
        std::for_each(ptsCDT.begin(), ptsCDT.end(), [&](PointConvectionDiffusionT &pt) -> void {
            pt.Point::operator=(*iter);
            f.setPoint(pt);
            f2.setPoint(pt);
            f.setTime(PointConvectionDiffusionT::getTime());
            f2.setTime(PointConvectionDiffusionT::getTime());
            (*previousIter)["sol"] = f.u();

            f.setTime(PointConvectionDiffusionT::getTime() - PointConvectionDiffusionT::getDelta());
            f2.setTime(PointConvectionDiffusionT::getTime() - PointConvectionDiffusionT::getDelta());
            (*twoLevelsIter)["sol"] = f.u();

            f.setTime(PointConvectionDiffusionT::getTime() + PointConvectionDiffusionT::getDelta());
            *iterConvections = Coordinate(-UNITVALUE, -UNITVALUE);
            pt["bdv"] = f.u();
            pt["rhs"] = f.f();

            ++previousIter, ++twoLevelsIter, ++iterConvections, ++iter;
        });
        previousIter = previousValues.begin(), twoLevelsIter = twoLevelsUpValues.begin();
        iterConvections = convections.begin(), iter = pts->begin();
    };

    auto assignBoundary = [&]() {
        auto f = function();
        auto f2 = function2D();
        std::for_each(ptsCDT.begin(), ptsCDT.end(), [&](PointConvectionDiffusionT &pt) -> void {
            f.setPoint(pt);
            f2.setPoint(pt);

            f.setTime(PointConvectionDiffusionT::getTime() + PointConvectionDiffusionT::getDelta());
            f2.setTime(PointConvectionDiffusionT::getTime() + PointConvectionDiffusionT::getDelta());
            pt["bdv"] = f.u();
            pt["rhs"] = f.f();
            *iterConvections = Coordinate(-UNITVALUE, -UNITVALUE);
            ++iterConvections;
        });
        iterConvections = convections.begin();
    };

    auto findElement = [&]() {
        std::for_each(ptsCDT.begin(), ptsCDT.end(), [&](PointConvectionDiffusionT &pt) -> void {
            pt.findElement(pts);
        });
    };

    auto makeMatrix = [&]() {
        std::for_each(ptsCDT.begin(), ptsCDT.end(), [&](PointConvectionDiffusionT &pt) -> void {
            pt.setConvection(*iterConvections);
            pt.calcRepresentationFormula();
            pt.makeDifferentiation();
            ++iterConvections;
        });
        iterConvections = convections.begin();
    };

    auto eraseMatrix = [&]() {
        std::for_each(ptsCDT.begin(), ptsCDT.end(), [&](PointConvectionDiffusionT &pt) -> void {
            pt.setMatrixRow(std::array<matrix_row, 2>());
            pt.setRhsMatrixRow(std::array<matrix_row, 2>());
            pt.setDMatrixRow(std::array<matrix_row, 2>());
            pt.setRb(std::array<double, 2>());
        });
    };

    auto calcDifferentiation = [&]() {
        std::for_each(ptsCDT.begin(), ptsCDT.end(), [&](PointConvectionDiffusionT &pt) -> void {
            pt.calcDifferentiation();
        });
    };

    auto updateValues = [&]() {
        std::for_each(ptsCDT.begin(), ptsCDT.end(), [&](PointConvectionDiffusionT &pt) -> void {
            *twoLevelsIter = *previousIter;
            *previousIter = pt.getValue();
            ++twoLevelsIter;
            ++previousIter;
        });
        previousIter = previousValues.begin(), twoLevelsIter = twoLevelsUpValues.begin();
    };

    auto updateTime = [&]() {
        PointConvectionDiffusionT::setTime(
                PointConvectionDiffusionT::getTime() + PointConvectionDiffusionT::getDelta());
        std::cout.precision(5);
        std::cout << std::defaultfloat << "current time = " << PointConvectionDiffusionT::getTime() << std::endl;
    };

    assignInitial();
    findElement();
    makeMatrix();

    auto matrix = AGM::matrix<PointConvectionDiffusionT>(&ptsCDT);
    matrix.makeMatrix();
    matrix.calculateMatrix();

    calcDifferentiation();
    updateValues();

    auto wf = WriteFile<PointConvectionDiffusionT>(&ptsCDT);
    auto filename = "errors_" + std::to_string(PointConvectionDiffusionT::getDelta());
    auto output = std::ofstream(filename);
    output.precision(16);
    output << PointConvectionDiffusionT::getTime() << " " << wf.calculateError("sol") << std::endl;

    while (PointConvectionDiffusionT::getTime() + PointConvectionDiffusionT::getDelta() <
           UNITVALUE - HALFVALUE * PointConvectionDiffusionT::getDelta()) {
        updateTime();
        assignBoundary();
        eraseMatrix();
        makeMatrix();

        matrix.deleteMatrix();
        matrix.makeMatrix();
        matrix.calculateMatrix();

        calcDifferentiation();
        updateValues();
        output << PointConvectionDiffusionT::getTime() << " " << wf.calculateError("sol") << std::endl;
    }

    output.close();
    auto mv = "mv errors_" + std::to_string(PointConvectionDiffusionT::getDelta()) + " ../docker/AGM_test/";
    system(mv.c_str());

    for (int j = 0; j < pts->size(); ++j) {
        pts->at(j).setValue(ptsCDT.at(j).getValue());
    }

}

std::pair<std::vector<AGM::PointHeat>, std::vector<AGM::PointHeat>> AGM::solver::NavierStokesSolver() {
    std::cout << "Reynolds number = " << UNITVALUE / pts->at(0).getMp() << "\n";
    std::cout << "The number of threads = " << NT << "\n";
    std::cout << "The number of Points = " << pts->size() << "\n";
    omp_set_dynamic(0);
    omp_set_num_threads(NT);
    Point::setPtsnum(int(pts->size()));
    Point::setPts(pts);
    auto uvel = std::vector<PointHeat>(Point::getPtsnum());
    auto vvel = std::vector<PointHeat>(Point::getPtsnum());
    PointHeat::setTime(0);
    PointHeat::setDelta(0.005);
    double terminal{251};
    auto puvel = std::vector<Value>(Point::getPtsnum());
    auto pvvel = std::vector<Value>(Point::getPtsnum());
    auto ppvel = std::vector<Value>(Point::getPtsnum());
    auto ppuvel = std::vector<Value>(Point::getPtsnum());
    auto ppvvel = std::vector<Value>(Point::getPtsnum());
    auto uvalue = std::vector<Value>(Point::getPtsnum());
    auto vvalue = std::vector<Value>(Point::getPtsnum());
    PointHeat::setCn(UNITVALUE);
    PointHeat::setBdf(UNITVALUE);
    auto rhs_u_x = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / PointHeat::getDelta() +
               puvel.at(i)["phi"];
    };
    auto rhs_u_x1 = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / PointHeat::getDelta() +
               puvel.at(i)["phi"];
    };
    auto rhs_u_y = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / PointHeat::getDelta() -
               puvel.at(i)["phi"] + 2 * puvel.at(i)["sol"] * pvvel.at(i)["dy"];
    };
    auto rhs_u_y1 = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / PointHeat::getDelta() -
               puvel.at(i)["phi"] + uvalue.at(i)["sol"] * vvalue.at(i)["dy"] + puvel.at(i)["sol"] * pvvel.at(i)["dy"];
    };
    auto rhs_v_x = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / PointHeat::getDelta() +
               pvvel.at(i)["phi"] + 2 * pvvel.at(i)["sol"] * puvel.at(i)["dx"];
    };
    auto rhs_v_x1 = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / PointHeat::getDelta() +
               pvvel.at(i)["phi"] + vvalue.at(i)["sol"] * uvalue.at(i)["dx"] + pvvel.at(i)["sol"] * puvel.at(i)["dx"];
    };
    auto rhs_v_y = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / PointHeat::getDelta() -
               pvvel.at(i)["phi"];
    };
    auto rhs_v_y1 = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / PointHeat::getDelta() -
               pvvel.at(i)["phi"];
    };
    auto rhs_p_x = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhs_p_y = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhs_u_x_t = [&](int i) -> double {
        return pow(puvel.at(i)["sol"], 2) - uvel.at(i).getMp() * puvel.at(i)["dx"];
    };
    auto rhs_u_x_t1 = [&](int i) -> double {
        return HALFVALUE * (pow(uvalue.at(i)["sol"], 2) + pow(puvel.at(i)["sol"], 2)) -
               uvel.at(i).getMp() * puvel.at(i)["dx"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
    };
    auto rhs_u_y_t = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
    };
    auto rhs_u_y_t1 = [&](int i) -> double {
        return uvalue.at(i)["sol"] * vvalue.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"] -
               uvel.at(i).getMp() * puvel.at(i)["dy"];
    };
    auto rhs_v_x_t = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
    };
    auto rhs_v_x_t1 = [&](int i) -> double {
        return uvalue.at(i)["sol"] * vvalue.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"] -
               vvel.at(i).getMp() * pvvel.at(i)["dx"];
    };
    auto rhs_v_y_t = [&](int i) -> double {
        return pow(pvvel.at(i)["sol"], 2) - vvel.at(i).getMp() * pvvel.at(i)["dy"];
    };
    auto rhs_v_y_t1 = [&](int i) -> double {
        return HALFVALUE * (pow(vvalue.at(i)["sol"], 2) + pow(pvvel.at(i)["sol"], 2)) -
               vvel.at(i).getMp() * pvvel.at(i)["dy"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
    };
    auto rhs_p_x_t = [&](int i) -> double {
        return uvel.at(i)["sol"] / PointHeat::getDelta();
    };
    auto rhs_p_y_t = [&](int i) -> double {
        return vvel.at(i)["sol"] / PointHeat::getDelta();
    };
    auto rhs_p_x_diff = [&](int i) -> double {
        return -uvel.at(i)["dx"] / PointHeat::getDelta();
    };
    auto rhs_p_y_diff = [&](int i) -> double {
        return -vvel.at(i)["dy"] / PointHeat::getDelta();
    };
    int nIter{}, fixed{}, nDiv{};

//    auto NS = std::ifstream("/home/jjhong0608/docker/AGM_test/ExternalFlow/test/Re100_7/AGM_Result_2");
//    if (NS.fail()) {
//        printError("file is not opened");
//    }
//    int idx{};
//    double x{}, y{}, p{}, px{}, py{};

    auto assignInitial = [&]() -> void {
        auto f = function2D();
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            f.setTime(PointHeat::getTime() - PointHeat::getDelta());
            uvel.at(j).Point::operator=(pts->at(j));
            vvel.at(j).Point::operator=(pts->at(j));
            f.setPoint(uvel.at(j));
            ppuvel.at(j)["sol"] = f.u();
            ppuvel.at(j)["rhs"] = f.f1();
            ppuvel.at(j)["phi"] = f.phi();
            ppuvel.at(j)["dx"] = f.ux();
            ppuvel.at(j)["dy"] = f.uy();
            f.setPoint(vvel.at(j));
            ppvvel.at(j)["sol"] = f.v();
            ppvvel.at(j)["rhs"] = f.f2();
            ppvvel.at(j)["phi"] = f.psi();
            ppvvel.at(j)["dx"] = f.vx();
            ppvvel.at(j)["dy"] = f.vy();
            f.setTime(PointHeat::getTime());
            f.setPoint(uvel.at(j));
            puvel.at(j)["sol"] = UNITVALUE;
            puvel.at(j)["rhs"] = f.f1();
            puvel.at(j)["phi"] = f.phi();
            puvel.at(j)["dx"] = f.ux();
            puvel.at(j)["dy"] = f.uy();
            f.setPoint(vvel.at(j));
            pvvel.at(j)["sol"] = ZEROVALUE;
            pvvel.at(j)["rhs"] = f.f2();
            pvvel.at(j)["phi"] = f.psi();
            pvvel.at(j)["dx"] = f.vx();
            pvvel.at(j)["dy"] = f.vy();

//            NS >> idx >> x >> y >> puvel.at(idx)["sol"] >> pvvel.at(idx)["sol"] >> ppvel.at(idx)["sol"]
//               >> puvel.at(idx)["dx"] >> puvel.at(idx)["dy"] >> pvvel.at(idx)["dx"] >> pvvel.at(idx)["dy"]
//               >> puvel.at(idx)["phi"] >> pvvel.at(idx)["phi"] >> ppvel.at(idx)["dx"] >> ppvel.at(idx)["dy"];


            f.setTime(PointHeat::getTime() + PointHeat::getDelta());
            f.setPoint(uvel.at(j));
            uvel.at(j)["sol"] = f.u();
            uvel.at(j)["bdv"] = f.u();
            uvel.at(j)["rhs"] = f.f1();
            uvel.at(j)["dx"] = f.ux();
            uvel.at(j)["dy"] = f.uy();
            f.setPoint(vvel.at(j));
            vvel.at(j)["sol"] = f.v();
            vvel.at(j)["bdv"] = f.v();
            vvel.at(j)["rhs"] = f.f2();
            vvel.at(j)["dx"] = f.vx();
            vvel.at(j)["dy"] = f.vy();
            pts->at(j)["sol"] = f.p();
            pts->at(j)["bdv"] = ZEROVALUE;
            pts->at(j)["dx"] = f.px();
            pts->at(j)["dy"] = f.py();
            pts->at(j).setMp(UNITVALUE);
            if (pts->at(j).getCondition() == 'D') {
                pts->at(j).setCondition('N');
                pts->at(j)["bdv"] = ZEROVALUE;
            } else if (pts->at(j).getCondition() == 'N') {
                uvel.at(j)["bdv"] = ZEROVALUE;
                vvel.at(j)["bdv"] = ZEROVALUE;
                pts->at(j)["bdv"] = ZEROVALUE;
                pts->at(j).setCondition('D');
            }
            if (isclose(pts->at(j)[0], UNITVALUE) && isclose(pts->at(j)[1], UNITVALUE)) {
                fixed = j;
            }
        }
    };

    auto assignBoundaryvalue = [&]() -> void {
        #pragma omp single
        {
            auto f = function2D();
            for (int j = 0; j < Point::getPtsnum(); ++j) {
                f.setTime(PointHeat::getTime() + PointHeat::getDelta());
                f.setPoint(uvel.at(j));
                uvel.at(j)["bdv"] = f.u();
                uvel.at(j)["rhs"] = f.f1();
                if (uvel.at(j).getCondition() == 'N') {
                    uvel.at(j)["bdv"] = ZEROVALUE;
                }
                f.setPoint(vvel.at(j));
                vvel.at(j)["bdv"] = f.v();
                vvel.at(j)["rhs"] = f.f2();
                if (vvel.at(j).getCondition() == 'N') {
                    vvel.at(j)["bdv"] = ZEROVALUE;
                }
            }
        }
    };

    auto findElement = [&]() -> void {
        PointHeat::setPtsH(&uvel);
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            uvel.at(j).findElement(pts);
            uvel.at(j).findElement1(pts);
        }
        PointHeat::setPtsH(&vvel);
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            vvel.at(j).findElement(pts);
            vvel.at(j).findElement1(pts);
        }
    };

    auto makeMatrix = [&]() -> void {
        #pragma omp single
        PointHeat::setPtsH(&uvel);
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            uvel.at(j).calcRepresentationFormula();
            uvel.at(j).makeDifferentiation();
            uvel.at(j).updateRb(rhs_u_x, rhs_u_y);
            uvel.at(j).updateRb_t(rhs_u_x_t, rhs_u_y_t);

        }
        #pragma omp barrier
        #pragma omp single
        PointHeat::setPtsH(&vvel);
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            vvel.at(j).calcRepresentationFormula();
            vvel.at(j).makeDifferentiation();
            vvel.at(j).updateRb(rhs_v_x, rhs_v_y);
            vvel.at(j).updateRb_t(rhs_v_x_t, rhs_v_y_t);
        }
    };

    auto makeMatrix_p = [&]() -> void {
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            pts->at(j).calcRepresentationFormula();
            pts->at(j).makeDifferentiation();
            pts->at(j).updateRb(rhs_p_x, rhs_p_y);
            pts->at(j).updateRb_t(rhs_p_x_t, rhs_p_y_t);
        }
    };

    auto calculateDifferentiation = [&]() -> void {
        #pragma omp single
        PointHeat::setPtsH(&uvel);
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            uvel.at(j).calcDifferentiation(rhs_u_x, rhs_u_y, rhs_u_x_t, rhs_u_y_t);
        }
        #pragma omp barrier
        #pragma omp single
        PointHeat::setPtsH(&vvel);
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            vvel.at(j).calcDifferentiation(rhs_v_x, rhs_v_y, rhs_v_x_t, rhs_v_y_t);
        }
    };

    auto calculateDifferentiation1 = [&]() -> void {
        #pragma omp single
        PointHeat::setPtsH(&uvel);
        double d{};
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            d = uvel.at(j)["dx"];
            uvel.at(j).calcDifferentiation(rhs_u_x1, rhs_u_y1, rhs_u_x_t1, rhs_u_y_t1);
            uvel.at(j)["dx"] = d;
        }
        #pragma omp barrier
        #pragma omp single
        PointHeat::setPtsH(&vvel);
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            d = vvel.at(j)["dy"];
            vvel.at(j).calcDifferentiation(rhs_v_x1, rhs_v_y1, rhs_v_x_t1, rhs_v_y_t1);
            vvel.at(j)["dy"] = d;
        }
    };

    auto calculateDifferentiation_p = [&]() -> void {
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            pts->at(j).calcDifferentiation(rhs_p_x, rhs_p_y, rhs_p_x_t, rhs_p_y_t, rhs_p_x_diff, rhs_p_y_diff);
        }
    };

    auto updateRb = [&]() -> void {
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            uvel.at(j).updateRb(rhs_u_x, rhs_u_y);
            uvel.at(j).updateRb_t(rhs_u_x_t, rhs_u_y_t);
            vvel.at(j).updateRb(rhs_v_x, rhs_v_y);
            vvel.at(j).updateRb_t(rhs_v_x_t, rhs_v_y_t);
        }
    };

    auto updateRb1 = [&]() -> void {
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            uvel.at(j).updateRb(rhs_u_x1, rhs_u_y1);
            uvel.at(j).updateRb_t(rhs_u_x_t1, rhs_u_y_t1);
            vvel.at(j).updateRb(rhs_v_x1, rhs_v_y1);
            vvel.at(j).updateRb_t(rhs_v_x_t1, rhs_v_y_t1);
        }
    };

    auto updateRb_p = [&]() -> void {
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            pts->at(j).updateRb(rhs_p_x, rhs_p_y);
            pts->at(j).updateRb_t(rhs_p_x_t, rhs_p_y_t);
        }
    };

    auto updateValue = [&]() -> void {
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            if (pts->at(j).getCondition() != 'N') {
                uvel.at(j)["sol"] -= pts->at(j)["dx"] * PointHeat::getDelta();
                vvel.at(j)["sol"] -= pts->at(j)["dy"] * PointHeat::getDelta();
            }
            pts->at(j)["sol"] -= HALFVALUE * uvel.at(j).getMp() * (uvel.at(j)["dx"] + vvel.at(j)["dy"]);
            uvel.at(j)["dx"] -= pts->at(j)["dxx"] * PointHeat::getDelta();
            vvel.at(j)["dy"] -= pts->at(j)["dyy"] * PointHeat::getDelta();
            uvalue.at(j) = uvel.at(j).getValue();
            vvalue.at(j) = vvel.at(j).getValue();
        }
    };

    auto updateValue1 = [&]() -> void {
        #pragma omp for
        for (int j = 0; j < Point::getPtsnum(); ++j) {
            ppuvel.at(j) = puvel.at(j);
            ppvvel.at(j) = pvvel.at(j);
            puvel.at(j) = uvel.at(j).getValue();
            pvvel.at(j) = vvel.at(j).getValue();
            ppvel.at(j) = pts->at(j).getValue();
        }
    };

    auto updateTime = [&]() -> void {
        #pragma omp single
        {
            PointHeat::setTime(PointHeat::getTime() + PointHeat::getDelta());
            std::cout << "current time = " << PointHeat::getTime() << '\n';
            ++nIter;
        }
    };

    double cfl0{}, cfl1{};
    bool cfl_bool{false};
    auto checkDiv = [&]() -> bool {
        auto findMax = [&](std::vector<PointHeat> *vector) -> int {
            double maxVal{};
            int maxIdx{};
            for (int j = 0; j < Point::getPtsnum(); ++j) {
                if (std::fabs(vector->at(j)["sol"]) > maxVal) {
                    maxVal = vector->at(j)["sol"];
                    maxIdx = vector->at(j).getIdx();
                }
            }
            return maxIdx;
        };

        int maxIdxU = findMax(&uvel), maxIdxV = findMax(&vvel);
        double maxH = 0.04;
        double vv = max(std::fabs(uvel.at(maxIdxU)["sol"]), std::fabs(vvel.at(maxIdxV)["sol"])) * PointHeat::getDelta();
        std::cout << "0 = " << vv << "\n";
        std::cout << "1 = " << maxH << "\n";
        if (std::isnan(vv)) {
            return true;
        }
        return vv > maxH;
    };

    assignInitial();
    findElement();

    auto matrix = matrixTwoDimension<PointHeat>(&uvel, &vvel);
    auto matrixP = AGM::matrix<Point>(pts);
//    auto matrixP = matrixEIGEN<Point>(pts, fixed);
    #pragma omp parallel
    {
        makeMatrix();
    }
    matrix.makeMatrix();
    matrix.calculateMatrix();

    #pragma omp parallel
    {
        calculateDifferentiation();
        makeMatrix_p();
    }

    matrixP.makeMatrix();
    matrixP.calculateMatrix();

    #pragma omp parallel
    {
        calculateDifferentiation_p();
        updateValue();
        updateRb1();
    }
    matrix.calculateMatrix();
    #pragma omp parallel
    {
        calculateDifferentiation1();
        updateValue1();
    }

    auto wf = WriteFileNS<Point, PointHeat>(pts, &uvel, &vvel);

    system("mkdir -p /home/jjhong0608/docker/AGM_test/ExternalFlow/test/Re100_10");

    while (PointHeat::getTime() + PointHeat::getDelta() < terminal - HALFVALUE * PointHeat::getDelta()) {
        #pragma omp parallel
        {
            updateTime();
//            assignBoundaryvalue();
            updateRb();
        }
        matrix.calculateMatrix();
        #pragma omp parallel
        {
            calculateDifferentiation();
            updateRb_p();
        }
        matrixP.calculateMatrix();
        #pragma omp parallel
        {
            calculateDifferentiation_p();
            updateValue();
            updateRb1();
        }
        matrix.calculateMatrix();
        #pragma omp parallel
        {
            calculateDifferentiation1();
            updateValue1();
        }

//        if (checkDiv()) {
//            wf.writeResult("/home/jjhong0608/docker/AGM_test/ExternalFlow/test/Re240/AGM_Div_" + std::to_string(nDiv));
//            ++nDiv;
//            if (std::isnan(wf.calculateError("uvel"))) {
//                printError("NaN value0");
//                exit(1);
//            }
//        }
        if (nIter % 20 == 0) {
            wf.writeResult("/home/jjhong0608/docker/AGM_test/ExternalFlow/test/Re100_10/AGM_Temp_" + std::to_string(nIter / 20));
            system("python /home/jjhong0608/docker/AGM_test/ExternalFlow/test/calculateDrag.py");
        }

        if (nIter % 200 == 0) {
            wf.writeResult("/home/jjhong0608/docker/AGM_test/ExternalFlow/test/Re100_10/AGM_Result_" + std::to_string(nIter / 200));
            std::cout << std::setw(32) << std::left << "Relative error of u-velocity = " << std::scientific
            << wf.calculateError("uvel") << std::endl;
            std::cout << std::setw(32) << std::left << "Relative error of v-velocity = " << std::scientific
            << wf.calculateError("vvel") << std::endl;
            std::cout << std::setw(32) << std::left << "Relative error of pressure = " << std::scientific
            << wf.calculateError("pressure") << std::endl;
            if (std::isnan(wf.calculateError("uvel"))) {
                printError("NaN value");
                exit(1);
            }
        }
    }

    return std::make_pair(uvel, vvel);
}

std::pair<std::vector<AGM::Point>, std::vector<AGM::Point>> AGM::solver::StokesSolver() {
    Point::setPtsnum(int(pts->size()));
    Point::setPts(pts);
    auto ptsU = std::vector<Point>(Point::getPtsnum()), ptsV = std::vector<Point>(Point::getPtsnum());
    auto iter = pts->begin(), iterU = ptsU.begin(), iterV = ptsV.begin();
    auto previousValue = std::vector<Value>(Point::getPtsnum());
    auto previousIter = previousValue.begin();
    double pError{UNITVALUE}, tol{2e-2}, fixedValue{};
    int nIter{}, fixedPoint{};
    auto f = function2D();
    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        f.setPoint(pt);
        pt["sol"] = f.p();
    });
    Point::setPts(&ptsU);
    std::for_each(ptsU.begin(), ptsU.end(), [&](Point &pt) -> void {
        pt = *iter;
        f.setPoint(pt);
        pt.findElement(pts);
        pt["sol"] = f.u();
        pt["bdv"] = f.u();
        ++iter;
    });
    iter = pts->begin();
    Point::setPts(&ptsV);
    std::for_each(ptsV.begin(), ptsV.end(), [&](Point &pt) -> void {
        pt = *iter;
        f.setPoint(pt);
        pt.findElement(pts);
        pt["sol"] = f.v();
        pt["bdv"] = f.v();
        ++iter;
    });
    iter = pts->begin();
    std::for_each(ptsU.begin(), ptsU.end(), [&](Point &pt) -> void {
        pt.calcRepresentationFormula();
        pt.makePressureTerm();
        pt.calcPressureTerm(pts, 'u');
        pt.makeDifferentiation();
    });
    std::for_each(ptsV.begin(), ptsV.end(), [&](Point &pt) -> void {
        pt.setMatrixRow(iterU->getMatrixRow());
        pt.setDMatrixRow(iterU->getDMatrixRow());
        pt.setPMatrixRow(iterU->getPMatrixRow());
        pt.setRhsMatrixRow(iterU->getRhsMatrixRow());
        pt.updateRb();
        pt.calcPressureTerm(pts, 'v');
        ++iterU;
    });
    iterU = ptsU.begin();

    auto matrix = AGM::matrixTwoDimension<Point>(&ptsU, &ptsV);
    matrix.makeMatrix();
    matrix.calculateMatrix();

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        pt.calcContinuity(&ptsU, &ptsV);
        if (isclose(pt[0], HALFVALUE) && isclose(pt[1], HALFVALUE)) {
            fixedPoint = pt.getIdx();
        }
    });

    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        pt["sol"] = pt["sol"] * 1.5E0 - previousIter->getSol() * HALFVALUE;
        ++previousIter;
    });
    previousIter = previousValue.begin();

    f.setPoint(pts->at(fixedPoint));
    fixedValue = pts->at(fixedPoint)["sol"] - f.p();
    std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
        pt["sol"] -= fixedValue;
    });

    auto wf = WriteFile<Point>(pts);
    pError = wf.calculateIterationError(&previousValue);

    std::cout.precision(16);
    std::cout << "Relative Pressure error of " << ++nIter << " Iterations: " << std::scientific << pError;

    while (pError > tol) {
        std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
            previousIter->setSol(pt["sol"]);
            ++previousIter;
        });
        previousIter = previousValue.begin();

        std::for_each(ptsU.begin(), ptsU.end(), [&](Point &pt) -> void {
            pt.updateRb();
            pt.calcPressureTerm(pts, 'u');
        });
        std::for_each(ptsV.begin(), ptsV.end(), [&](Point &pt) -> void {
            pt.updateRb();
            pt.calcPressureTerm(pts, 'v');
        });

        matrix.calculateMatrix();

        std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
            pt.calcContinuity(&ptsU, &ptsV);
        });

        std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
            pt["sol"] = pt["sol"] * 1.5E0 - previousIter->getSol() * HALFVALUE;
            ++previousIter;
        });
        previousIter = previousValue.begin();

        fixedValue = pts->at(fixedPoint)["sol"] - f.p();
        std::for_each(pts->begin(), pts->end(), [&](Point &pt) -> void {
            pt["sol"] -= fixedValue;
        });

        pError = wf.calculateIterationError(&previousValue);
        std::cout << "Relative Pressure error of " << ++nIter << " Iterations: " << pError << '\n';
    }

    std::for_each(ptsU.begin(), ptsU.end(), [&](Point &pt) {
        pt.calcDifferentiation();
    });
    std::for_each(ptsV.begin(), ptsV.end(), [&](Point &pt) {
        pt.calcDifferentiation();
    });

    return std::make_pair(ptsU, ptsV);
}
