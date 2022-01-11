//
// Created by NIMS-JUNHONG on 2020/12/28.
//

#include "AxialLine.h"

AGM::AxialLine::AxialLine() = default;

AGM::AxialLine::AxialLine(char mark) : mark(mark) {}

char AGM::AxialLine::getMark() const {
    return mark;
}

void AGM::AxialLine::setMark(char i) {
    AxialLine::mark = i;
}

double &AGM::AxialLine::operator[](int i) {
    return coordinate[i];
}

AGM::AxialLine::~AxialLine() = default;
