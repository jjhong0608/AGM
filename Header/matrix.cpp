//
// Created by NIMS-JUNHONG on 2020/12/31.
//

#include "matrix.h"
#include "mkl.h"

template<typename pt>
AGM::matrix<pt>::matrix() = default;

template<typename pt>
AGM::matrix<pt>::matrix(std::vector<pt> *pts) : pts(pts) {}

template<typename pt>
AGM::matrix<pt>::~matrix() = default;

template<typename pt>
std::vector<pt> *AGM::matrix<pt>::getPts() const {
    return pts;
}

template<typename pt>
void AGM::matrix<pt>::setPts(std::vector<pt> *pVector) {
    matrix::pts = pVector;
}

template<typename pt>
int *AGM::matrix<pt>::getIa() const {
    return ia;
}

template<typename pt>
int *AGM::matrix<pt>::getJa() const {
    return ja;
}

template<typename pt>
double *AGM::matrix<pt>::getEnt() const {
    return ent;
}

template<typename pt>
void AGM::matrix<pt>::makeMatrix() {

    if (!ia) {
        auto row = matrix_row{};
        int n_entry{};
        for (const auto &i : *pts) {
            for (const auto &j : i.getMatrixRow()) {
                row = j;
                n_entry += int(row.size());
            }
        }
        ia = new int[2 * pts->size() + 1];
        ja = new int[n_entry];
        ent = new double[n_entry];

        int ia_idx{}, ja_idx{}, ia_value{1};
        ia[ia_idx++] = ia_value;
        for (const auto &i : *pts) {
            row = i.getMatrixRow()[0];
            ia_value += int(row.size());
            ia[ia_idx++] = ia_value;
            for (const auto &j : row) {
                ja[ja_idx] = j.idx + 1;
                if (std::isnan(j.value)) {
                    std::cout << "i = " << i.getIdx() << '\n';
                    std::cout << "bd = " << i.getCondition() << '\n';
                    exit(123);
                }
                ent[ja_idx++] = j.value;
            }
        }
        for (const auto &i : *pts) {
            row = i.getMatrixRow()[1];
            ia_value += int(row.size());
            ia[ia_idx++] = ia_value;
            for (const auto &j : row) {
                ja[ja_idx] = j.idx + 1;
                if (std::isnan(j.value)) {
                    std::cout << "i = " << i.getIdx() << '\n';
                    std::cout << "bd = " << i.getCondition() << '\n';
                    exit(123);
                }
                ent[ja_idx++] = j.value;
            }
        }
    } else {
        printError("AGM::matrix<pt>::makeMatrix()");
    }
}

template<typename pt>
void AGM::matrix<pt>::calculateMatrix() {
    mkl_set_num_threads(NT);
    mkl_set_dynamic(0);

    int size = int(pts->size());
    int n = size * 2;
    auto *rb = new double[n];
    int rb_idx{}, nrhs{1};
//    for (auto &i : *pts) {
//        rb[rb_idx++] = i.getRb()[0];
//    }
//    for (auto &i : *pts) {
//        rb[rb_idx++] = i.getRb()[1];
//    }

    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        rb[i] = pts->at(i).getRb()[0];
        rb[i + size] = pts->at(i).getRb()[1];
    }

    void *ppt[64];
    int iparm[64];
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
    pardiso(ppt, &maxfct, &mnum, &mtype, &phase, &n, ent, ia, ja, &idum, &nrhs, iparm, &msglvl, rb, x, &error);
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        pts->at(i)["sol"] = x[i];
        pts->at(i)["phi"] = x[i + size];
    }

    phase = -1;
    pardiso(ppt, &maxfct, &mnum, &mtype, &phase, &n, ent, ia, ja, &idum, &nrhs, iparm, &msglvl, rb, x, &error);

    delete[] rb;
}

template<typename pt>
void AGM::matrix<pt>::exportMatrix() {
    auto row = matrix_row{};
    int n_entry{};
    for (const auto &i : *pts) {
        for (const auto &j : i.getMatrixRow()) {
            row = j;
            while (row.back().idx >= 2 * pts->size()) {
                row.pop_back();
            }
            n_entry += int(row.size());
        }
    }

    std::ofstream ia_output("ia");
    std::ofstream ja_output("ja");
    std::ofstream ent_output("ent");
    std::ofstream rb_output("rb");

    rb_output.precision(16);
    for (int i = 0; i < 2 * pts->size() + 1; ++i) {
        ia_output << i << " " << ia[i] << std::endl;
        if (i < pts->size()) {
            rb_output << i << " " << std::scientific << pts->at(i).getRb()[0] << std::endl;
        } else if (i < 2 * pts->size()) {
            rb_output << i << " " << std::scientific << pts->at(i - pts->size()).getRb()[1] << std::endl;
        }
    }
    ent_output.precision(16);
    for (int i = 0; i < n_entry; ++i) {
        ja_output << i << " " << ja[i] << std::endl;
        ent_output << i << " " << std::scientific << ent[i] << std::endl;
    }

    ia_output.close();
    ja_output.close();
    ent_output.close();
    rb_output.close();

    system("cp ia ../docker/AGM_test/");
    system("cp ja ../docker/AGM_test/");
    system("cp ent ../docker/AGM_test/");
    system("cp rb ../docker/AGM_test/");
}

template<typename pt>
void AGM::matrix<pt>::deleteMatrix() {
    delete[] ia;
    delete[] ja;
    delete[] ent;
    ia = nullptr;
    ja = nullptr;
    ent = nullptr;
}

template class AGM::matrix<AGM::Point>;
template class AGM::matrix<AGM::PointHeat>;
template class AGM::matrix<AGM::PointConvectionDiffusion>;
template class AGM::matrix<AGM::PointConvectionDiffusionT>;