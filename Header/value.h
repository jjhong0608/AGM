//
// Created by JUNHONG-NIMS on 2021/01/09.
//

#ifndef AGM_VALUE_H
#define AGM_VALUE_H

#include "coordinate.h"

namespace AGM {
    class Value {
    private:
        double sol{}, phi{}, rhs{}, bdv{}, dx{}, dy{}, dxx{}, dyy{}, dxy{};
    public:
        Value();

        virtual ~Value();

        double &operator[](const std::string &string);

        const double &operator[](const std::string &string) const;

        Value &operator=(const Value &value);

        [[nodiscard]] double getSol() const;

        void setSol(double d);

        [[nodiscard]] double getPhi() const;

        void setPhi(double d);

        [[nodiscard]] double getRhs() const;

        void setRhs(double d);

        [[nodiscard]] double getBdv() const;

        void setBdv(double d);

        [[nodiscard]] double getDx() const;

        void setDx(double d);

        [[nodiscard]] double getDy() const;

        void setDy(double d);

        [[nodiscard]] double getDxx() const;

        void setDxx(double d);

        [[nodiscard]] double getDyy() const;

        void setDyy(double d);

        [[nodiscard]] double getDxy() const;

        void setDxy(double d);
    };
}



#endif //AGM_VALUE_H
