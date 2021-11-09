//
// Created by JUNHONG-NIMS on 2021/01/09.
//

#include "value.h"

AGM::Value::Value() = default;

AGM::Value::~Value() = default;

double AGM::Value::getSol() const {
    return sol;
}

void AGM::Value::setSol(double d) {
    Value::sol = d;
}

double AGM::Value::getPhi() const {
    return phi;
}

void AGM::Value::setPhi(double d) {
    Value::phi = d;
}

double AGM::Value::getRhs() const {
    return rhs;
}

void AGM::Value::setRhs(double d) {
    Value::rhs = d;
}

double AGM::Value::getBdv() const {
    return bdv;
}

void AGM::Value::setBdv(double d) {
    Value::bdv = d;
}

double AGM::Value::getDx() const {
    return dx;
}

void AGM::Value::setDx(double d) {
    Value::dx = d;
}

double AGM::Value::getDy() const {
    return dy;
}

void AGM::Value::setDy(double d) {
    Value::dy = d;
}

double AGM::Value::getDxx() const {
    return dxx;
}

void AGM::Value::setDxx(double d) {
    Value::dxx = d;
}

double AGM::Value::getDyy() const {
    return dyy;
}

void AGM::Value::setDyy(double d) {
    Value::dyy = d;
}

double AGM::Value::getDxy() const {
    return dxy;
}

void AGM::Value::setDxy(double d) {
    Value::dxy = d;
}

double &AGM::Value::operator[](const std::string &string) {
    if (string == "sol") return sol;
    if (string == "phi") return phi;
    if (string == "rhs") return rhs;
    if (string == "bdv") return bdv;
    if (string == "dx") return dx;
    if (string == "dy") return dy;
    if (string == "dxx") return dxx;
    if (string == "dyy") return dyy;
    if (string == "dxy") return dxy;
    printError("Value::operator[]", "input = %s", string.c_str());
}

const double &AGM::Value::operator[](const std::string &string) const {
    if (string == "sol") return sol;
    if (string == "phi") return phi;
    if (string == "rhs") return rhs;
    if (string == "bdv") return bdv;
    if (string == "dx") return dx;
    if (string == "dy") return dy;
    if (string == "dxx") return dxx;
    if (string == "dyy") return dyy;
    if (string == "dxy") return dxy;
    printError("Value::operator[]", "input = %s", string.c_str());
}

AGM::Value &AGM::Value::operator=(const AGM::Value &value) {
    if (this != &value) {
        sol = value.sol;
        phi = value.phi;
        rhs = value.rhs;
        bdv = value.bdv;
        dx = value.dx;
        dy = value.dy;
        dxx = value.dxx;
        dyy = value.dyy;
        dxy = value.dxy;
    }
    return *this;
}
