//
// Created by NIMS-JUNHONG on 2020/12/22.
//
#ifndef AGM_UTIL_H
#define AGM_UTIL_H

#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <limits>
#include <sstream>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <numeric>
#include <functional>
#include <chrono>
#include <ctime>
#include <omp.h>

#ifndef UNITVALUE
#define UNITVALUE 1.0000000000000000E0
#endif
#ifndef HALFVALUE
#define HALFVALUE 5.0000000000000000E-1
#endif
#ifndef ZEROVALUE
#define ZEROVALUE 0.0000000000000000E0
#endif
#ifndef NEARZERO
#define NEARZERO 1.0000000000000000E-10
#endif
#ifndef NT
#define NT 10
#endif

namespace AGM {
    bool isclose(double x, double y, double eps = NEARZERO);

    void printError(const std::string &function_name);

    void printError(const char *function_name, const char *fmt, ...);

    bool iszero(double x, double eps = NEARZERO);

    double sgn(double d);

    bool ispositive(double d);

    bool isnegative(double d);

    [[nodiscard]] double min(double d0, double d1);

    [[nodiscard]] double max(double d0, double d1);

    enum EWNS {
        E, W, N, S, EN, ES, WN, WS, NE, NW, SE, SW
    };
}

extern std::ofstream ElapsedTime;

#endif //AGM_UTIL_H
