//
// Created by NIMS-JUNHONG on 2020/12/23.
//

#ifndef AGM_POINT_H
#define AGM_POINT_H

#include "value.h"

namespace AGM {
    class Point;

    class matrix_row;

    class Element {
    private:
        std::array<Point *, 12> element{nullptr,};
    public:
        Element();

        virtual ~Element();

        [[nodiscard]] const std::array<Point *, 12> &getElement() const;

        void setElement(const std::array<Point *, 12> &array);

        Point *&operator[](int i);

        Point *&operator[](EWNS ewns);

        Element &operator=(const Element &rhs);
    };

    class Point {
    protected:
        int idx{};
        Coordinate coordinate{};
        Coordinate normal{};
        double mp{};
        char condition{};
        Element element{};
        Element element1{};
        Value value{};
        std::array<matrix_row, 2> matrixRow{}, rhsMatrixRow{}, phiMatrixRow{};
        std::array<matrix_row, 2> dMatrixRow{}, pMatrixRow{}, pdMatrixRow{};
        std::array<double, 2> rb{}, dv{};
        std::array<AxialLine *, 2> axialLine{};
        static std::vector<Point> *pts;
        static int ptsnum;
    public:
        Point();

        Point(const Point &src);

        explicit Point(int idx);

        explicit Point(const Coordinate &coordinate);

        Point(int idx, const Coordinate &coordinate);

        Point(const Coordinate &coordinate, double mp);

        Point(int idx, const Coordinate &coordinate, double mp);

        virtual ~Point();

        [[nodiscard]] int getIdx() const;

        void setIdx(int i);

        [[nodiscard]] static int getPtsnum();

        static void setPtsnum(int i);

        [[nodiscard]] const Coordinate &getCoordinate() const;

        void setCoordinate(const Coordinate &src);

        [[nodiscard]] const Coordinate &getNormal() const;

        void setNormal(const Coordinate &n);

        [[nodiscard]] const Element &getElement() const;

        void setElement(const Element &src);

        [[nodiscard]] const Element &getElement1() const;

        void setElement1(const Element &src);

        [[nodiscard]] double getMp() const;

        void setMp(double d);

        [[nodiscard]] char getCondition() const;

        void setCondition(char i);

        [[nodiscard]] const Value &getValue() const;

        [[nodiscard]] Value &getValue();

        void setValue(const Value &v);

        [[nodiscard]] const std::array<matrix_row, 2> &getMatrixRow() const;

        void setMatrixRow(const std::array<matrix_row, 2> &row);

        [[nodiscard]] const std::array<matrix_row, 2> &getRhsMatrixRow() const;

        void setRhsMatrixRow(const std::array<matrix_row, 2> &row);

        [[nodiscard]] const std::array<matrix_row, 2> &getPhiMatrixRow() const;

        void setPhiMatrixRow(const std::array<matrix_row, 2> &row);

        [[nodiscard]] const std::array<matrix_row, 2> &getDMatrixRow() const;

        void setDMatrixRow(const std::array<matrix_row, 2> &row);

        [[nodiscard]] const std::array<matrix_row, 2> &getPMatrixRow() const;

        void setPMatrixRow(const std::array<matrix_row, 2> &row);

        [[nodiscard]] const std::array<matrix_row, 2> &getPdMatrixRow() const;

        void setPdMatrixRow(const std::array<matrix_row, 2> &pdMatrixRow);

        [[nodiscard]] const std::array<double, 2> &getRb() const;

        void setRb(const std::array<double, 2> &array);

        [[nodiscard]] const std::array<AxialLine *, 2> &getAxialLine() const;

        AxialLine *&getAxialLine(char i);

        void setAxialLine(const std::array<AxialLine *, 2> &array);

        void setAxialLine(AxialLine *line, char i);

        static std::vector<Point> *getPts();

        static void setPts(std::vector<Point> *pVector);

        void findElement(std::vector<AxialLine> *xline, std::vector<AxialLine> *yline);

        void findElement1(std::vector<AxialLine> *xline, std::vector<AxialLine> *yline);

        virtual void findElement(std::vector<Point> *src);

        virtual void findElement1(std::vector<Point> *src);

        void calcRepresentationFormula();

        virtual void calcRepresentationFormula_cross();

        void calcRepresentationFormula_dirichlet();

        virtual void calcRepresentationFormula_neumann();

        virtual void calcRepresentationFormula_neumann_D();

        virtual void calcRepresentationFormula_neumann_N();

        virtual matrix_row calcRepresentationFormula_neumann_D(std::string &string);

        virtual matrix_row calcRepresentationFormula_neumann_N(std::string &string);

        virtual void calcRepresentationFormula_interface();

        void calcRepresentationFormula_dirichlet_and_Neumann(AGM::Point *pt, const EWNS &ewns, int order);

        void calcRepresentationFormula_phi();

        void calcRepresentationFormula_phi_plane(Point *pt0, Point *pt1, Point *pt2);

        virtual void calcRepresentationFormula_boundary();

        void makeDifferentiation();

        virtual void makeDifferentiation_cross();

        virtual void makeDifferentiation_boundary();

        virtual void makeDifferentiation_boundary(const std::string &string);

        virtual void makeDifferentiation_boundary_D(const std::string &string);

        virtual void makeDifferentiation_boundary_N(const std::string &string);

        virtual void makeDifferentiation_interface();

        virtual void calcDifferentiation();

        virtual void calcDifferentiation(const std::function<double(int)> &f, const std::function<double(int)> &g);

        virtual void calcDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                         const std::function<double(int)> &gx, const std::function<double(int)> &gy);

        virtual void calcDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                         const std::function<double(int)> &gx, const std::function<double(int)> &gy,
                                         const std::function<double(int)> &hx, const std::function<double(int)> &hy);

        virtual void calcDifferentiation_dirichlet_and_Neumann(AGM::Point *pt, const EWNS &ewns, int order,
                                                               const std::string &string);

        virtual void calcDifferentiationWithPressure(std::vector<Point> *pressure, char i);

        void updateRb();

        virtual void updateRb_cross();

        virtual void updateRb_dirichlet();

        virtual void updateRb_neumann();

        virtual void updateRb_interface();

        void updateRb(const std::function<double(int)> &f);

        void updateRb_cross(const std::function<double(int)> &f);

        void updateRb_dirichlet(const std::function<double(int)> &f);

        void updateRb_neumann(const std::function<double(int)> &f);

        void updateRb_interface(const std::function<double(int)> &f);

        void updateRb_t(const std::function<double(int)> &f);

        void updateRb_t_cross(const std::function<double(int)> &f);

        static void updateRb_t_dirichlet(const std::function<double(int)> &f);

        void updateRb_t_neumann(const std::function<double(int)> &f);

        void updateRb_t_interface(const std::function<double(int)> &f);

        void updateRb(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_cross(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_dirichlet(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_neumann(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_interface(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_t(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_t_cross(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_t_dirichlet(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_t_neumann(const std::function<double(int)>& f, const std::function<double(int)>& g);

        void updateRb_t_interface(const std::function<double(int)>& f, const std::function<double(int)>& g);

        double checkRepresentationFormula();

        double checkDifferentiation(const std::function<double(int)> &fx, const std::function<double(int)> &fy,
                                    const std::function<double(int)> &gx, const std::function<double(int)> &gy);

        EWNS findPoint_dirichlet_and_Neumann(Point *&pt);

        void findPoint_non_axial_lines(Point *&pt0, Point *&pt1, Point *&pt2);

        double &operator[](int i);

        const double &operator[](int i) const;

        AGM::Point *&operator[](EWNS ewns);

        double &operator[](const std::string &string);

        const double &operator[](const std::string &string) const;

        double operator-(const Point &rhs);

        Point &operator=(const Point &rhs);

        void calculate_finite_difference(char i);

        void calculate_cross_difference();

        double integrate();

        void translate_rhs(double d);

        void calcContinuity(std::vector<Point> *uvel, std::vector<Point> *vvel);

        virtual void calcContinuity_cross(std::vector<Point> *uvel, std::vector<Point> *vvel);

        virtual void calcContinuity_boundary(std::vector<Point> *uvel, std::vector<Point> *vvel);

        virtual void calcContinuity_interface(std::vector<Point> *uvel, std::vector<Point> *vvel);

        void
        calcContinuity(std::vector<Point *> *uvel, std::vector<Point *> *vvel, std::vector<Coordinate> *convections);

        void calcContinuity_cross(std::vector<Point *> *uvel, std::vector<Point *> *vvel,
                                  std::vector<Coordinate> *convections);

        void calcContinuity_boundary(std::vector<Point *> *uvel, std::vector<Point *> *vvel,
                                     std::vector<Coordinate> *convections);

        void calcContinuity_interface(std::vector<Point *> *uvel, std::vector<Point *> *vvel,
                                      std::vector<Coordinate> *convections);

        void makePressureTerm();

        virtual void makePressureTerm_cross();

        virtual void makePressureTerm_neumann();

        virtual void makePressureTerm_interface();

        void makePressureDifferenceTerm();

        virtual void makePressrueDifferenceTerm_cross();

        virtual void makePressrueDifferenceTerm_boundary();

        virtual void makePressrueDifferenceTerm_interface();

        void calcPressureTerm(std::vector<Point> *pressure, char i);

        void calcPressureTerm_cross(std::vector<Point> *pressure, char i);

        void calcPressureTerm_neumman(std::vector<Point> *pressure, char i);

        void calcPressureTerm_interface(std::vector<Point> *pressure, char i);

        void roundMatrix();

        void printInformation();
    };
}

#endif //AGM_POINT_H
