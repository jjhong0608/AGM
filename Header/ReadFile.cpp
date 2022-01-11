//
// Created by NIMS-JUNHONG on 2020/12/29.
//

#include "ReadFile.h"

void AGM::ReadFile::loadAxialData(const std::string& filename, std::vector<AGM::Point> *pts,
                                  std::vector<AxialLine> *aline_x, std::vector<AxialLine> *aline_y) {
    std::ifstream f(filename);
    std::string line{};
    std::array<std::string, 3> temp_string{};
    int nRegion{}, idx{}, index{}, temp_integer{}, present_idx{};
    int nCross{}, nBound{}, nXaxialPoints{}, nYaxialPoints{};
    int countXaxialPoitns{}, countYaxialPoints{};
    std::array<double, 2> temp_pt{};
    double mp{}, bv{};
    char bc{};
    AGM::Point pt{};
    std::vector<std::vector<int>> xline{}, yline{};
    std::vector<int> aline{};

    if (!f.is_open()) {
        printError("AGM::ReadFile::loadAxialData", "No Axial data file: \"%s\"\nPlease check file name",
                   filename.c_str());
    }

    std::cout << "Axial file: \"" << filename << "\" open\n";

    while (!f.eof()) {
        std::getline(f, line);
        if (line.empty()) std::getline(f, line);
        if (line.find("ENDREGION") != std::string::npos) {
            idx = 0;
            std::getline(f, line);
        }
        if (line.find("REGION") != std::string::npos) {
            nRegion += 1;
            std::cout << "REGION " << nRegion << " reading...\n";
            std::getline(f, line);
            f >> temp_string[0] >> temp_string[1] >> temp_string[2] >> mp;
        }
        if (line.find('=') != std::string::npos) {
            idx += 1;
            switch (idx) {
                case 1:
                    nCross = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    for (int i = 0; i < nCross; ++i) {
                        f >> temp_integer >> pt[0] >> pt[1];
                        pt.setMp(mp);
                        pt.setCondition('C');
                        pt.setIdx(index++);
                        pts->emplace_back(pt);
                    }
                    break;
                case 2:
                    nBound = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    for (int i = 0; i < nBound; ++i) {
                        f >> temp_integer >> pt[0] >> pt[1] >> bc >> bv >> temp_pt[0] >> temp_pt[1];
                        pt.setMp(mp);
                        if (bc != 'D' && bc != 'N' && bc != 'I' && bc != 'F') {
                            printError("AGM::ReadFile::loadAxialData", "boundary condition (which is %s) is wrong", bc);
                        }
                        pt.setCondition(bc);
                        pt.getValue()["bdv"] = bv;
                        pt.setNormal(AGM::Coordinate(temp_pt[0], temp_pt[1]));
                        pt.setIdx(index++);
                        pts->emplace_back(pt);
                    }
                    break;
                case 3:
                    countXaxialPoitns = 0;
                    nXaxialPoints = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    while (countXaxialPoitns < nXaxialPoints) {
                        f >> present_idx;
                        countXaxialPoitns += 1;
                        if (pts->at(present_idx).getCondition() == 'C') {
                            printError("AGM::ReadFile::loadAxialData",
                                       "The first Point of the x-axial line is the cross Point");
                        }
                        aline.emplace_back(present_idx);
                        f >> present_idx;
                        countXaxialPoitns += 1;
                        while (pts->at(present_idx).getCondition() == 'C') {
                            aline.emplace_back(present_idx);
                            f >> present_idx;
                            countXaxialPoitns += 1;
                        }
                        aline.emplace_back(present_idx);
                        xline.emplace_back(aline);
                        aline.clear();
                    }
                    break;
                case 4:
                    countYaxialPoints = 0;
                    nYaxialPoints = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    while (countYaxialPoints < nYaxialPoints) {
                        f >> present_idx;
                        countYaxialPoints += 1;
                        if (pts->at(present_idx).getCondition() == 'C') {
                            printError("AGM::ReadFile::loadAxialData",
                                       "The first Point of the y-axial line is the cross Point");
                        }
                        aline.emplace_back(present_idx);
                        f >> present_idx;
                        countYaxialPoints += 1;
                        while (pts->at(present_idx).getCondition() == 'C') {
                            aline.emplace_back(present_idx);
                            f >> present_idx;
                            countYaxialPoints += 1;
                        }
                        aline.emplace_back(present_idx);
                        yline.emplace_back(aline);
                        aline.clear();
                    }
                    break;
                default:
                    printError("AGM::ReadFile::loadAxialData", "switch index (which is %d) error", idx);
                    break;
            }
        }
    }
    f.close();

    auto xAxialLine = AGM::AxialLine('x');
    auto yAxialLine = AGM::AxialLine('y');

    for (const auto &i : xline) {
        xAxialLine[0] = pts->at(i.front())[1];
        xAxialLine[1] = pts->at(i.front())[0];
        xAxialLine[2] = pts->at(i.back())[0];
        for (const auto &j : i) {
            xAxialLine.emplace_back(&(pts->at(j)));
        }
        aline_x->emplace_back(xAxialLine);
        xAxialLine.clear();
    }

    for (const auto &i : yline) {
        yAxialLine[0] = pts->at(i.front())[0];
        yAxialLine[1] = pts->at(i.front())[1];
        yAxialLine[2] = pts->at(i.back())[1];
        for (const auto &j : i) {
            yAxialLine.emplace_back(&(pts->at(j)));
        }
        aline_y->emplace_back(yAxialLine);
        yAxialLine.clear();
    }

    for (auto &i : *aline_x) {
        for (const auto &j : i) {
            j->setAxialLine(&i, 'x');
        }
    }

    for (auto &i : *aline_y) {
        for (const auto &j : i) {
            j->setAxialLine(&i, 'y');
        }
    }

    AGM::Point *prev{};

    for (const auto &i : *aline_x) {
        for (const auto &j : i) {
            if (j == i.front()) {
                if (j->getCondition() != 'I') {
                    (*j)[W] = j;
                    if (j->getNormal()[1] > ZEROVALUE && !j->getAxialLine('y')) (*j)[N] = j;
                    if (j->getNormal()[1] < ZEROVALUE && !j->getAxialLine('y')) (*j)[S] = j;
                    if (isclose(j->getNormal()[1], ZEROVALUE) && !j->getAxialLine('y')) (*j)[N] = j, (*j)[S] = j;
                }
                prev = j;
            } else if (j == i.back()) {
                (*prev)[E] = j;
                if (j->getCondition() != 'I') {
                    (*j)[E] = j;
                    if (j->getNormal()[1] > ZEROVALUE && !j->getAxialLine('y')) (*j)[N] = j;
                    if (j->getNormal()[1] < ZEROVALUE && !j->getAxialLine('y')) (*j)[S] = j;
                    if (isclose(j->getNormal()[1], ZEROVALUE) && !j->getAxialLine('y')) (*j)[N] = j, (*j)[S] = j;
                }
                (*j)[W] = prev;
            } else {
                (*prev)[E] = j;
                (*j)[W] = prev;
                prev = j;
            }
        }
    }

    for (const auto &i : *aline_y) {
        for (const auto &j : i) {
            if (j == i.front()) {
                if (j->getCondition() != 'I') {
                    (*j)[S] = j;
                    if (j->getNormal()[0] > ZEROVALUE && !j->getAxialLine('x')) (*j)[E] = j;
                    if (j->getNormal()[0] < ZEROVALUE && !j->getAxialLine('x')) (*j)[W] = j;
                    if (isclose(j->getNormal()[0], ZEROVALUE) && !j->getAxialLine('x')) (*j)[E] = j, (*j)[W] = j;
                }
                prev = j;
            } else if (j == i.back()) {
                (*prev)[N] = j;
                if (j->getCondition() != 'I') {
                    (*j)[N] = j;
                    if (j->getNormal()[0] > ZEROVALUE && !j->getAxialLine('x')) (*j)[E] = j;
                    if (j->getNormal()[0] < ZEROVALUE && !j->getAxialLine('x')) (*j)[W] = j;
                    if (isclose(j->getNormal()[0], ZEROVALUE) && !j->getAxialLine('x')) (*j)[E] = j, (*j)[W] = j;
                }
                (*j)[S] = prev;
            } else {
                (*prev)[N] = j;
                (*j)[S] = prev;
                prev = j;
            }
        }
    }

    auto element = Element();
    for (auto &item : *pts) {
        element = item.getElement();
        if (item.getAxialLine('x')) {
            if (element[E]) element[E] = item[E]->getElement().getElement()[E];
            if (element[W]) element[W] = item[W]->getElement().getElement()[W];
        }
        if (item.getAxialLine('y')) {
            if (element[N]) element[N] = item[N]->getElement().getElement()[N];
            if (element[S]) element[S] = item[S]->getElement().getElement()[S];
        }
        item.setElement1(element);
    }
}
