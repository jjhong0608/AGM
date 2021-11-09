//
// Created by NIMS-JUNHONG on 2021/03/04.
//

#ifndef AGM_POINTCONVECTIONDIFFUSIONT_H
#define AGM_POINTCONVECTIONDIFFUSIONT_H

#include "PointConvectionDiffusion.h"

namespace AGM {
    class PointConvectionDiffusionT : public PointHeat {
    private:
        static std::vector<Value> *two_levels_up_values;
        static std::vector<Coordinate> *convections;
        Coordinate convection{};
        static std::vector<PointConvectionDiffusionT> *ptsCDT;

    public:
        PointConvectionDiffusionT();

        static std::vector<PointConvectionDiffusionT> *getPtsCdt();

        static void setPtsCdt(std::vector<PointConvectionDiffusionT> *ptsCdt);

        [[nodiscard]] const Coordinate &getConvection() const;

        void setConvection(const Coordinate &coefficient);

        static std::vector<Coordinate> *getConvections();

        static void setConvections(std::vector<Coordinate> *vector);

        static std::vector<Value> *getTwoLevelsUpValues();

        static void setTwoLevelsUpValues(std::vector<Value> *twoLevelsUpValues);

        void findElement(std::vector<Point> *src) override;

        void calcRepresentationFormula_cross() override;

        void calcRepresentationFormula_neumann() override;

        void calcRepresentationFormula_interface() override;

        void makeDifferentiation_cross() override;

        void makeDifferentiation_boundary() override;

        void makeDifferentiation_interface() override;

        void calcDifferentiation() override;
    };

}


#endif //AGM_POINTCONVECTIONDIFFUSIONT_H
