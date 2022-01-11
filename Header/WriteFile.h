//
// Created by NIMS-JUNHONG on 2021/01/04.
//

#ifndef AGM_WRITEFILE_H
#define AGM_WRITEFILE_H

#include "matrixEIGEN.h"

namespace AGM {
    template<typename T>
    class WriteFile {
    protected:
        T *pt{};
        std::vector<T> *pts{};
    public:
        WriteFile();

        explicit WriteFile(std::vector<T> *pts);

        T *getPt() const;

        std::vector<T> *getPts() const;

        void setPoint(T *pPoint);

        void setPoints(std::vector<T> *points);

        auto calculateErrorAtPoint();

        double calculateError();

        virtual std::array<double, 2> calculateErrorAtPoint(const std::string &string);

        virtual double calculateError(const std::string &string);

        auto calculateIterationErrorAtPoint(AGM::Value *value);

        double calculateIterationError(std::vector<Value> *values);

        virtual void writeResult(const std::string &string);
    };

    template<typename T1, typename T2>
    class WriteFileNS : public WriteFile<T1> {
    private:
        T2 *ptU{}, *ptV{};
        std::vector<T2> *ptsU{}, *ptsV{};
    public:
        WriteFileNS();

        WriteFileNS(std::vector<T1> *pts, std::vector<T2> *ptsU, std::vector<T2> *ptsV);

        void writeResult(const std::string &string) override;

        std::array<double, 2> calculateErrorAtPoint(const std::string &string) override;

        double calculateError(const std::string &string) override;
    };
}


#endif //AGM_WRITEFILE_H
