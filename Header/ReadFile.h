//
// Created by NIMS-JUNHONG on 2020/12/29.
//

#ifndef AGM_READFILE_H
#define AGM_READFILE_H

#include "solver.h"

namespace AGM {
    class ReadFile {
    public:
        static void
        loadAxialData(const std::string &filename, std::vector<AGM::Point> *pts, std::vector<AxialLine> *aline_x,
                      std::vector<AxialLine> *aline_y);
    };

}


#endif //AGM_READFILE_H
