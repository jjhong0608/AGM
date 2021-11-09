//
// Created by NIMS-JUNHONG on 2021/02/01.
//

#include "matrixOthers.h"
#include "mkl.h"

template<typename pt>
AGM::matrixTwoDimension<pt>::matrixTwoDimension() = default;

template<typename pt>
AGM::matrixTwoDimension<pt>::matrixTwoDimension(std::vector<pt> *pts, std::vector<pt> *pts0) : matrix<pt>(pts),
                                                                                                 pts0(pts0) {}

template<typename pt>
AGM::matrixTwoDimension<pt>::~matrixTwoDimension() = default;

template<typename pt>
std::vector<pt> *AGM::matrixTwoDimension<pt>::getPts0() const {
    return pts0;
}

template<typename pt>
void AGM::matrixTwoDimension<pt>::setPts0(std::vector<pt> *pVector) {
    matrixTwoDimension::pts0 = pVector;
}

template<typename pt>
void AGM::matrixTwoDimension<pt>::calculateMatrix() {
    mkl_set_num_threads(NT);
    mkl_set_dynamic(0);

    int size = matrix<pt>::pts->size();
    int n = size * 2;
    auto *rb = new double[2 * n];
    int rb_idx{}, nrhs{2};
//    for (auto &i : *matrix<pt>::pts) {
//        rb[rb_idx++] = i.getRb()[0];
//    }
//    for (auto &i : *matrix<pt>::pts) {
//        rb[rb_idx++] = i.getRb()[1];
//    }
//    for (auto &i : *pts0) {
//        rb[rb_idx++] = i.getRb()[0];
//    }
//    for (auto &i : *pts0) {
//        rb[rb_idx++] = i.getRb()[1];
//    }
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        rb[i] = matrix<pt>::pts->at(i).getRb()[0];
        rb[i + size] = matrix<pt>::pts->at(i).getRb()[1];
        rb[i + 2 * size] = pts0->at(i).getRb()[0];
        rb[i + 3 * size] = pts0->at(i).getRb()[1];
    }

    void *ppt[64];
    int iparm[64];
    int mtype, maxfct, mnum, phase, error{}, msglvl, idum;
    double x[2 * n];

    for (auto &i : x) i = ZEROVALUE;
    for (auto &i : iparm) i = 0;

    iparm[0] = 1;         /* No solver default */
    iparm[1] = 3;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
    mtype = 11;
    iparm[60] = 1;

    for (auto &i : ppt) i = nullptr;
    phase = 13;
    pardiso(ppt, &maxfct, &mnum, &mtype, &phase, &n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja, &idum, &nrhs,
            iparm, &msglvl, rb, x, &error);
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

//    for (int i = 0; i < size; ++i) {
//        matrix<pt>::pts->at(i)["sol"] = x[i];
//        matrix<pt>::pts->at(i)["phi"] = x[i + size];
//    }
//    for (int i = 2 * size; i < 3 * size; ++i) {
//        pts0->at(i - 2 * size)["sol"] = x[i];
//        pts0->at(i - 2 * size)["phi"] = x[i + size];
//    }
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        matrix<pt>::pts->at(i)["sol"] = x[i];
        matrix<pt>::pts->at(i)["phi"] = x[i + size];
        pts0->at(i)["sol"] = x[i + 2 * size];
        pts0->at(i)["phi"] = x[i + 3 * size];
    }

    phase = -1;
    pardiso(ppt, &maxfct, &mnum, &mtype, &phase, &n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja, &idum, &nrhs,
            iparm, &msglvl, rb, x, &error);

    delete[] rb;
}

template<typename pt>
void AGM::matrixTwoDimension<pt>::calculateMatrixWithPressure(std::vector<Point> *pressure) {
    mkl_set_num_threads(NT);
    mkl_set_dynamic(0);

    int size = matrix<pt>::pts->size();
    int n = size * 2;
    auto *rb = new double[2 * n];
    int rb_idx{}, nrhs{2};

    for (auto &i : *matrix<pt>::pts) {
        rb[rb_idx] = i.getRb()[0];
        for (const auto &item : i.getPMatrixRow()[0]) {
            if (item.idx < Point::getPtsnum()) {
                rb[rb_idx] -= item.value * pressure->at(item.idx)["sol"];
            } else if (item.idx < 2 * Point::getPtsnum()) {
                rb[rb_idx] -= item.value * pressure->at(item.idx - Point::getPtsnum())["sol"];
            }
        }
        ++rb_idx;
    }
    for (auto &i : *matrix<pt>::pts) {
        rb[rb_idx] = i.getRb()[1];
        for (const auto &item : i.getPMatrixRow()[0]) {
            if (item.idx < Point::getPtsnum()) {
                rb[rb_idx] -= item.value * pressure->at(item.idx)["sol"];
            } else if (item.idx < 2 * Point::getPtsnum()) {
                rb[rb_idx] += item.value * pressure->at(item.idx - Point::getPtsnum())["sol"];
            }
        }
//        if (i.getCondition() == 'D' || i.getCondition() == 'N') {
//            for (const auto &item : i.getPdMatrixRow()[0]) {
//                rb[rb_idx] -= item.value * pressure->at(item.idx)["sol"];
//            }
//        }
        ++rb_idx;
    }
    for (auto &i : *pts0) {
        rb[rb_idx] = i.getRb()[0];
        for (const auto &item : i.getPMatrixRow()[1]) {
            if (item.idx < Point::getPtsnum()) {
                rb[rb_idx] -= item.value * pressure->at(item.idx)["sol"];
            } else if (item.idx < 2 * Point::getPtsnum()) {
                rb[rb_idx] -= item.value * pressure->at(item.idx - Point::getPtsnum())["sol"];
            }
        }
        ++rb_idx;
    }
    for (auto &i : *pts0) {
        rb[rb_idx] = i.getRb()[1];
        for (const auto &item : i.getPMatrixRow()[1]) {
            if (item.idx < Point::getPtsnum()) {
                rb[rb_idx] += item.value * pressure->at(item.idx)["sol"];
            } else if (item.idx < 2 * Point::getPtsnum()) {
                rb[rb_idx] -= item.value * pressure->at(item.idx - Point::getPtsnum())["sol"];
            }
        }
//        if (i.getCondition() == 'D' || i.getCondition() == 'N') {
//            for (const auto &item : i.getPdMatrixRow()[1]) {
//                rb[rb_idx] -= item.value * pressure->at(item.idx)["sol"];
//            }
//        }
        ++rb_idx;
    }

    void *ppt[64];
    int iparm[64];
    int mtype, maxfct, mnum, phase, error{}, msglvl, idum;
    double x[2 * n];

    for (auto &i : x) i = ZEROVALUE;
    for (auto &i : iparm) i = 0;

    iparm[0] = 1;         /* No solver default */
    iparm[1] = 3;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
    mtype = 11;
    iparm[60] = 1;

    for (auto &i : ppt) i = nullptr;
    phase = 13;
    pardiso(ppt, &maxfct, &mnum, &mtype, &phase, &n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja, &idum, &nrhs,
            iparm, &msglvl, rb, x, &error);
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

    for (int i = 0; i < size; ++i) {
        matrix<pt>::pts->at(i)["sol"] = x[i];
        matrix<pt>::pts->at(i)["phi"] = x[i + size];
    }
    for (int i = 2 * size; i < 3 * size; ++i) {
        pts0->at(i - 2 * size)["sol"] = x[i];
        pts0->at(i - 2 * size)["phi"] = x[i + size];
    }

    phase = -1;
    pardiso(ppt, &maxfct, &mnum, &mtype, &phase, &n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja, &idum, &nrhs,
            iparm, &msglvl, rb, x, &error);

    delete[] rb;
}

template class AGM::matrixTwoDimension<AGM::Point>;
template class AGM::matrixTwoDimension<AGM::PointHeat>;
template class AGM::matrixTwoDimension<AGM::PointConvectionDiffusion>;
template class AGM::matrixTwoDimension<AGM::PointConvectionDiffusionT>;