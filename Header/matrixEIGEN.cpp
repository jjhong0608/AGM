//
// Created by NIMS-JUNHONG on 2021/04/07.
//

#include "matrixEIGEN.h"
#include "mkl.h"

template<typename pt>
AGM::matrixEIGEN<pt>::matrixEIGEN() = default;

template<typename pt>
AGM::matrixEIGEN<pt>::matrixEIGEN(std::vector<pt> *pts) : matrix<pt>(pts) {}

template<typename pt>
AGM::matrixEIGEN<pt>::matrixEIGEN(std::vector<pt> *pts, int fixedPtrIdx) : matrix<pt>(pts), fixedPtrIdx(fixedPtrIdx) {}

template<typename pt>
AGM::matrixEIGEN<pt>::~matrixEIGEN() = default;

template<typename pt>
int AGM::matrixEIGEN<pt>::getFixedPtrIdx() const {
    return fixedPtrIdx;
}

template<typename pt>
void AGM::matrixEIGEN<pt>::setFixedPtrIdx(int i) {
    matrixEIGEN::fixedPtrIdx = i;
}

template<typename pt>
void AGM::matrixEIGEN<pt>::makeMatrix() {
    auto a = Eigen::SparseMatrix<double, Eigen::RowMajor>(2 * matrix<pt>::pts->size(), 2 * matrix<pt>::pts->size() - 1);

    for (const auto &item : *matrix<pt>::pts) {
        for (const auto &row : item.getMatrixRow()[0]) {
            if (!iszero(row.value)) {
                if (row.idx < fixedPtrIdx) {
                    a.insert(item.getIdx(), row.idx) = row.value;
                } else if (row.idx > fixedPtrIdx) {
                    a.insert(item.getIdx(), row.idx - 1) = row.value;
                }
            }
        }
    }
    for (const auto &item : *matrix<pt>::pts) {
        for (const auto &row : item.getMatrixRow()[1]) {
            if (!iszero(row.value)) {
                if (row.idx < fixedPtrIdx) {
                    a.insert(item.getIdx() + matrix<pt>::pts->size(), row.idx) = row.value;
                } else if (row.idx > fixedPtrIdx) {
                    a.insert(item.getIdx() + matrix<pt>::pts->size(), row.idx - 1) = row.value;
                }
            }
        }
    }
    A = a.transpose() * a;
    A.makeCompressed();

    matrix<pt>::ia = new int[A.outerSize() + 1];
    matrix<pt>::ja = new int[A.nonZeros()];
    matrix<pt>::ent = new double[A.nonZeros()];

    auto oPtr = A.outerIndexPtr();
    for (int j = 0; j < A.outerSize() + 1; ++j) {
        matrix<pt>::ia[j] = *oPtr + 1;
        ++oPtr;
    }

    auto iPtr = A.innerIndexPtr();
    auto vPtr = A.valuePtr();
    for (int j = 0; j < A.nonZeros(); ++j) {
        matrix<pt>::ja[j] = *iPtr + 1;
        matrix<pt>::ent[j] = *vPtr;
        ++iPtr;
        ++vPtr;
    }
    n_entry = A.nonZeros();

    A = a.transpose();
    A.makeCompressed();
    ia_t = new int[A.outerSize() + 1];
    ja_t = new int[A.nonZeros()];
    ent_t = new double[A.nonZeros()];

    oPtr = A.outerIndexPtr();
    for (int j = 0; j < A.outerSize() + 1; ++j) {
        ia_t[j] = *oPtr;
        ++oPtr;
    }

    iPtr = A.innerIndexPtr();
    vPtr = A.valuePtr();
    for (int j = 0; j < A.nonZeros(); ++j) {
        ja_t[j] = *iPtr;
        ent_t[j] = *vPtr;
        ++iPtr;
        ++vPtr;
    }
}

template<typename pt>
void AGM::matrixEIGEN<pt>::calculateMatrix() {
    mkl_set_num_threads(NT);
    mkl_set_dynamic(0);

    int size = int(matrix<pt>::pts->size());
    int n = size * 2 - 1;
    auto *rb = new double[n], *rb0 = new double[n + 1];

//    for (const auto &item : *matrix<pt>::pts) {
//        rb0[item.getIdx()] = item.getRb()[0];
//        rb0[item.getIdx() + size] = item.getRb()[1];
//    }
    #pragma omp parallel
    {
        #pragma omp for
        for (int j = 0; j < size; ++j) {
            rb0[j] = matrix<pt>::pts->at(j).getRb()[0];
            rb0[j + size] = matrix<pt>::pts->at(j).getRb()[1];
        }

        #pragma omp for
        for (int j = 0; j < n; ++j) {
            rb[j] = ZEROVALUE;
            for (int k = ia_t[j]; k < ia_t[j + 1]; ++k) {
                rb[j] += ent_t[k] * rb0[ja_t[k]];
            }
        }
    }

    int nrhs{1};
    void *ppt[64];
    int iparm[64];
    double dparm[64];
    int mtype, maxfct, mnum, phase, error{}, msglvl, idum;
    double x[n];

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

    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        if (i < fixedPtrIdx) {
            matrix<pt>::pts->at(i)["sol"] = x[i];
        } else if (i > fixedPtrIdx) {
            matrix<pt>::pts->at(i)["sol"] = x[i - 1];
        } else {
            matrix<pt>::pts->at(i)["sol"] = ZEROVALUE;
        }
        matrix<pt>::pts->at(i)["phi"] = x[i + size - 1];
    }
//    for (int i = 0; i < size; ++i) {
//        if (i < fixedPtrIdx) {
//            matrix<pt>::pts->at(i)["sol"] = x[i];
//        } else if (i > fixedPtrIdx) {
//            matrix<pt>::pts->at(i)["sol"] = x[i - 1];
//        } else {
//            matrix<pt>::pts->at(i)["sol"] = ZEROVALUE;
//        }
//        matrix<pt>::pts->at(i)["phi"] = x[i + size - 1];
//    }

    phase = -1;
    pardiso(ppt, &maxfct, &mnum, &mtype, &phase, &n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja, &idum, &nrhs,
            iparm, &msglvl, rb, x, &error);

    delete[] rb;
    delete[] rb0;
}

template<typename pt>
void AGM::matrixEIGEN<pt>::exportTranseposedMatrix() {
    std::ofstream ia_output("ia_t");
    std::ofstream ja_output("ja_t");
    std::ofstream ent_output("ent_t");

    for (int i = 0; i < A.outerSize() + 1; ++i) {
        ia_output << i << " " << ia_t[i] << std::endl;
    }
    ent_output.precision(16);
    for (int i = 0; i < A.nonZeros(); ++i) {
        ja_output << i << " " << ja_t[i] << std::endl;
        ent_output << i << " " << std::scientific << ent_t[i] << std::endl;
    }

    ia_output.close();
    ja_output.close();
    ent_output.close();
}

template<typename pt>
void AGM::matrixEIGEN<pt>::exportNormatMatrix() {
    std::ofstream ia_output("ia");
    std::ofstream ja_output("ja");
    std::ofstream ent_output("ent");

    for (int i = 0; i < 2 * matrix<pt>::pts->size(); ++i) {
        ia_output << i << " " << matrix<pt>::ia[i] << std::endl;
    }
    ent_output.precision(16);
    for (int i = 0; i < n_entry; ++i) {
        ja_output << i << " " << matrix<pt>::ja[i] << std::endl;
        ent_output << i << " " << std::scientific << matrix<pt>::ent[i] << std::endl;
    }

    ia_output.close();
    ja_output.close();
    ent_output.close();
}

template class AGM::matrixEIGEN<AGM::Point>;
template class AGM::matrixEIGEN<AGM::PointHeat>;
template class AGM::matrixEIGEN<AGM::PointConvectionDiffusion>;
template class AGM::matrixEIGEN<AGM::PointConvectionDiffusionT>;