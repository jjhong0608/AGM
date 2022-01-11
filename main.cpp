#include "Header/WriteFile.h"
#include "Header/ReadFile.h"
#include "Header/solver.h"
#include "Header/PointConvectionDiffusion.h"
#include <sys/resource.h>

std::ofstream ElapsedTime("ElapsedTime");

int main() {
    struct rlimit rlim{};
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = -1;
    rlim.rlim_max = -1;
    setrlimit(RLIMIT_STACK, &rlim);

    auto start_total = std::chrono::high_resolution_clock::now();
    auto pts = std::vector<AGM::Point>();
    auto xline = std::vector<AGM::AxialLine>();
    auto yline = std::vector<AGM::AxialLine>();
    AGM::ReadFile::loadAxialData("ALG_output", &pts, &xline, &yline);
    auto f = AGM::function();

    std::for_each(pts.begin(), pts.end(), [&](AGM::Point &pt) -> void { pt.findElement(&xline, &yline); });
    std::for_each(pts.begin(), pts.end(), [&](AGM::Point &pt) -> void { pt.findElement1(&xline, &yline); });

    /* ---- */

//    std::ofstream xaxial("/home/jjhong0608/docker/AGM_test/xaxialWavy");
//    std::ofstream yaxial("/home/jjhong0608/docker/AGM_test/yaxialWavy");
//
//    for (auto &item: xline) {
//        xaxial << item[0] << "\t" << item[1] << "\t" << item[2] << "\n";
//    }
//    for (auto &item: yline) {
//        yaxial << item[0] << "\t" << item[1] << "\t" << item[2] << "\n";
//    }
//
//    xaxial.close();
//    yaxial.close();
//
//    exit(1);

    /* ---- */

    auto solver = AGM::solver(&pts, &xline, &yline);

//    auto velocities = solver.NavierStokesSolver();
//    solver.streamSolver();
    solver.ellipticSolver();
//    solver.heatSolver();

//    auto wf = AGM::WriteFileNS<AGM::Point, AGM::PointHeat>(&pts, &velocities.first, &velocities.second);
    auto wf = AGM::WriteFile<AGM::Point>(&pts);
//    auto wf = AGM::WriteFile<AGM::PointHeat>(&velocities.first);

    wf.writeResult("/home/jjhong0608/docker/AGM_test/AGM_Result");

    auto printSingleError = [&]() -> void {
        std::cout.precision(16);
        std::cout << std::setw(32) << std::left << "Relative error of solution = " << std::scientific
                  << wf.calculateError("sol") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of phi = " << std::scientific
                  << wf.calculateError("phi") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of dx = " << std::scientific
                  << wf.calculateError("dx") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of dy = " << std::scientific
                  << wf.calculateError("dy") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of grad = " << std::scientific
                  << wf.calculateError("grad") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of dxx = " << std::scientific
                  << wf.calculateError("dxx") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of dyy = " << std::scientific
                  << wf.calculateError("dyy") << std::endl;
        std::cout << std::setw(32) << std::left << "The number of points = " << pts.size() << std::endl;
        std::cout << std::setw(32) << std::left << "The number of x-axial lines = " << xline.size() << std::endl;
        std::cout << std::setw(32) << std::left << "The number of y-axial lines = " << yline.size() << std::endl;
    };
    auto printMultipleError = [&]() -> void {
        std::cout.precision(16);
        std::cout << std::setw(32) << std::left << "Relative error of u-velocity = " << std::scientific
                  << wf.calculateError("uvel") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of v-velocity = " << std::scientific
                  << wf.calculateError("vvel") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of pressure = " << std::scientific
                  << wf.calculateError("pressure") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of ux = " << std::scientific
                  << wf.calculateError("ux") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of uy = " << std::scientific
                  << wf.calculateError("uy") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of vx = " << std::scientific
                  << wf.calculateError("vx") << std::endl;
        std::cout << std::setw(32) << std::left << "Relative error of vy = " << std::scientific
                  << wf.calculateError("vy") << std::endl;
        std::cout << std::setw(32) << std::left << "The number of points = " << pts.size() << std::endl;
        std::cout << std::setw(32) << std::left << "The number of x-axial lines = " << xline.size() << std::endl;
        std::cout << std::setw(32) << std::left << "The number of y-axial lines = " << yline.size() << std::endl;
    };

//    AGM::PointHeat::setPtsH(&velocities.first);
//    for (auto &item : velocities.first) {
//        item.findElement(&pts);
//    }
//    AGM::PointHeat::setPtsH(&velocities.second);
//    for (auto &item : velocities.second) {
//        item.findElement(&pts);
//    }

    printSingleError();
//    printMultipleError();

    auto end_total = std::chrono::high_resolution_clock::now();

    std::chrono::nanoseconds elapsedNS = end_total - start_total;
    std::chrono::microseconds elapsedMS = std::chrono::duration_cast<std::chrono::microseconds>(elapsedNS);
    std::chrono::seconds elpasedSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsedNS);
    std::time_t end_time = std::chrono::system_clock::to_time_t(end_total);
    ElapsedTime << "finished computation at " << std::ctime(&end_time)
                << "elapsed time: " << elapsedNS.count() << " (ns)" << std::endl
                << "elapsed time: " << elapsedMS.count() << " (ms)" << std::endl
                << "elpased time: " << elpasedSeconds.count() << " (s)" << std::endl;
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsedNS.count() << " (ns)" << std::endl
              << "elapsed time: " << elapsedMS.count() << " (ms)" << std::endl
              << "elpased time: " << elpasedSeconds.count() << " (s)" << std::endl;
    ElapsedTime.close();

    return 0;
}
