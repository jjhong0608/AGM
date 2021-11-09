//
// Created by NIMS-JUNHONG on 2021/01/25.
//

#ifndef AGM_POINTHEAT_H
#define AGM_POINTHEAT_H

#include "Point.h"

namespace AGM {
    class PointHeat : public Point {
    protected:
        static std::vector<Value> *previous_values;
        static double time, delta, CN, BDF;
        static std::function<double (int)> rhsParallelToX, rhsParallelToY;
        static std::vector<PointHeat> *ptsH;

    public:
        static double getTime();

        static void setTime(double d);

        static double getDelta();

        static void setDelta(double d);

        static double getCn();

        static void setCn(double cn);

        static double getBdf();

        static void setBdf(double bdf);

        static std::vector<PointHeat> *getPtsH();

        static void setPtsH(std::vector<PointHeat> *pVector);

        static const std::function<double(int)> &getRhsParallelToX();

        static void setRhsParallelToX(const std::function<double(int)> &function);

        static const std::function<double(int)> &getRhsParallelToY();

        static void setRhsParallelToY(const std::function<double(int)> &function);

        static std::vector<Value> *getPreviousValues();

        static void setPreviousValues(std::vector<Value> *previousValues);

        void findElement(std::vector<Point> *src) override;

        void findElement1(std::vector<Point> *src) override;

        void calcRepresentationFormula_cross() override;

        void calcRepresentationFormula_cross_E();

        void calcRepresentationFormula_cross_RD();

        matrix_row calcRepresentationFormula_cross_E(std::string &string);

        matrix_row calcRepresentationFormula_cross_RD(std::string &string);

        void calcRepresentationFormula_neumann() override;

        void calcRepresentationFormula_neumann_D_E();

        void calcRepresentationFormula_neumann_D_RD();

        void calcRepresentationFormula_neumann_NN_E();

        void calcRepresentationFormula_neumann_NN_RD();

        matrix_row calcRepresentationFormula_neumann_D_RD(std::string &string);

        matrix_row calcRepresentationFormula_neumann_NN_E(std::string &string);

        matrix_row calcRepresentationFormula_neumann_NN_RD(std::string &string);

        void calcRepresentationFormula_interface() override;

        void calcRepresentationFormula_interface_E();

        void calcRepresentationFormula_interface_RD();

        matrix_row calcRepresentationFormula_interface_E(std::string &string);

        matrix_row calcRepresentationFormula_interface_RD(std::string &string);

        void calcRepresentationFormula_boundary() override;

        void makeDifferentiation_cross() override;

        void makeDifferentiation_cross_E();

        void makeDifferentiation_cross_RD();

        void makeDifferentiation_boundary() override;

        void makeDifferentiation_boundary_E();

        void makeDifferentiation_boundary_RD();

        void makeDifferentiation_boundary(const std::string &string) override;

        void makeDifferentiation_boundary_E(const std::string &string);

        void makeDifferentiation_boundary_E_N(const std::string &string);

        void makeDifferentiation_boundary_RD(const std::string &string);

        void makeDifferentiation_boundary_RD_N(const std::string &string);

        void makeDifferentiation_interface() override;

        void makeDifferentiation_interface_E();

        void makeDifferentiation_interface_RD();

        void calcDifferentiation() override;

        void calcDifferentiation(const std::function<double(int)> &f, const std::function<double(int)> &g) override;

        void calcDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                 const std::function<double(int)> &gx, const std::function<double(int)> &gy) override;

        void calcDifferentiationWithPressure(std::vector<Point> *pressure, char i) override;

        void updateRb_cross() override;

        void updateRb_dirichlet() override;

        void updateRb_neumann() override;

        void updateRb_interface() override;

        void makePressureTerm_cross() override;

        void makePressureTerm_neumann() override;

        void makePressureTerm_interface() override;
    };
}


#endif //AGM_POINTHEAT_H
