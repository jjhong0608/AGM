//
// Created by NIMS-JUNHONG on 2021/01/26.
//

#ifndef AGM_POINTCONVECTIONDIFFUSION_H
#define AGM_POINTCONVECTIONDIFFUSION_H

#include "PointHeat.h"

namespace AGM {
    class PointConvectionDiffusion : public Point {
    private:
        static std::vector<Coordinate> *convections;
        static std::vector<PointConvectionDiffusion> *ptsCD;
        Coordinate convection{};
    public:
        PointConvectionDiffusion();

        static std::vector<Coordinate> *getConvections();

        static void setConvections(std::vector<Coordinate> *vector);

        [[nodiscard]] const Coordinate &getConvection() const;

        static std::vector<PointConvectionDiffusion> *getPtsCd();

        static void setPtsCd(std::vector<PointConvectionDiffusion> *ptsCdt);

        void setConvection(const Coordinate &coefficient);

        void findElement(std::vector<Point> *src) override;

        void calcRepresentationFormula_cross() override;

        void calcRepresentationFormula_neumann() override;

        void calcRepresentationFormula_interface() override;

        void makeDifferentiation_cross() override;

        void makeDifferentiation_boundary() override;

        void makeDifferentiation_interface() override;

        void calcDifferentiation() override;

        void calcContinuity_cross(std::vector<Point> *uvel, std::vector<Point> *vvel) override;

        void calcContinuity_boundary(std::vector<Point> *uvel, std::vector<Point> *vvel) override;

        void calcContinuity_interface(std::vector<Point> *uvel, std::vector<Point> *vvel) override;

        void makePressureTerm_cross() override;

        void makePressureTerm_neumann() override;

        void makePressureTerm_interface() override;

        void updateRb_cross() override;

        void updateRb_neumann() override;

        void updateRb_interface() override;
    };
}

#endif //AGM_POINTCONVECTIONDIFFUSION_H
