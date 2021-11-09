//
// Created by NIMS-JUNHONG on 2020/12/28.
//

#ifndef AGM_AXIALLINE_H
#define AGM_AXIALLINE_H

#include "GreenfunctionConvectionDiffusionLinear.h"

namespace AGM {
    class Point;

    class AxialLine : public std::vector<Point *> {
    private:
        char mark{};
        std::array<double, 3> coordinate{};
    public:
        AxialLine();

        explicit AxialLine(char mark);

        virtual ~AxialLine();

        [[nodiscard]] char getMark() const;

        void setMark(char i);

        double &operator[](int i);
    };
}


#endif //AGM_AXIALLINE_H
