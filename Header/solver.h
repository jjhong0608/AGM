//
// Created by NIMS-JUNHONG on 2021/01/04.
//

#ifndef AGM_SOLVER_H
#define AGM_SOLVER_H

#include <vector>
#include "WriteFile.h"

namespace AGM {
    class solver {
    private:
        std::vector<Point> *pts{};
        std::vector<AxialLine> *pXline{}, *pYline{};
    public:
        explicit solver(std::vector<Point> *pts);

        solver(std::vector<Point> *pts, std::vector<AxialLine> *pXline, std::vector<AxialLine> *pYline);

        virtual ~solver();

        void setPoints(std::vector<Point> *points);

        [[nodiscard]] std::vector<AxialLine> *getPXline() const;

        void setPXline(std::vector<AxialLine> *xline);

        [[nodiscard]] std::vector<AxialLine> *getPYline() const;

        void setPYline(std::vector<AxialLine> *yline);

        void ellipticSolver();

        void streamSolver();

        void heatSolver();

        void convectiondiffusionSolver();

        void unsteadyconvectiondiffusionSolver();

        std::pair<std::vector<AGM::PointHeat>, std::vector<AGM::PointHeat>> NavierStokesSolver();

        std::pair<std::vector<AGM::Point>, std::vector<AGM::Point>> StokesSolver();
    };

}


#endif //AGM_SOLVER_H
