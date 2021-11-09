//
// Created by NIMS-JUNHONG on 2020/12/24.
//

#ifndef AGM_FUNCTION_H
#define AGM_FUNCTION_H

#include "PointConvectionDiffusionT.h"

namespace AGM {
    class function {
    protected:
        double eps{}, x{}, y{}, time{};
    public:
        function();

        explicit function(const Point &pt);

        explicit function(const PointHeat &pt);

        function(double eps, double x, double y);

        function(double x, double y);

        virtual ~function();

        void setPoint(const Point &pt);

        void setPoint(const PointHeat &pt);

        [[nodiscard]] double getEps() const;

        void setEps(double d);

        [[nodiscard]] double getX() const;

        void setX(double d);

        [[nodiscard]] double getY() const;

        void setY(double d);

        [[nodiscard]] double getTime() const;

        void setTime(double d);

        [[nodiscard]] virtual double u() const;

        [[nodiscard]] virtual double phi() const;

        [[nodiscard]] double f() const;

        [[nodiscard]] double dx() const;

        [[nodiscard]] double dy() const;

        [[nodiscard]] double dxx() const;

        [[nodiscard]] double dyy() const;

        [[nodiscard]] double grad() const;
    };

}

#endif //AGM_FUNCTION_H
