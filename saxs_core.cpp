#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cmath>
#include <vector>
#include <iostream>

namespace py = pybind11;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double form_factor_sphere(double q, double r) {
    double qr = q * r;
    if (std::abs(qr) < 1e-8) {
        return 1.0;
    }
    double val = 3.0 * (std::sin(qr) - qr * std::cos(qr)) / (qr * qr * qr); // eq 7
    return val * val;
}

// compute Schulz Sphere Model Intensity
// q: input q array
// r_avg: average radius
// z: width parameter
// scale: scaling factor
// bg: background
py::array_t<double> model_intensity(
    py::array_t<double> q_array,
    double r_avg,
    double z,
    double scale,
    double bg
) {
    // input validation
    if (r_avg <= 0.1 || z <= 0.1) {
        // return array filled with bg
        auto q_buf = q_array.request();
        py::array_t<double> result = py::array_t<double>(q_buf.size);
        auto result_buf = result.request();
        double *ptr_res = (double *)result_buf.ptr;
        for (size_t i = 0; i < q_buf.size; i++) {
            ptr_res[i] = bg;
        }
        return result;
    }

    auto q_buf = q_array.request();
    double *ptr_q = (double *)q_buf.ptr;
    size_t n_q = q_buf.size;

    // prepare result array
    py::array_t<double> result = py::array_t<double>(n_q);
    auto result_buf = result.request();
    double *ptr_res = (double *)result_buf.ptr;

    py::gil_scoped_release release;

    // schulz distribution parameters
    double sigma_x = 1.0 / std::sqrt(z + 1.0);
    double x_min = std::max(0.01, 1.0 - 5.0 * sigma_x);
    double x_max = 1.0 + 8.0 * sigma_x;
    
    int n_grid = 200;
    double dx = (x_max - x_min) / (n_grid - 1);
    
    std::vector<double> x(n_grid);
    std::vector<double> Px_R6(n_grid); // Stores P(x) * R^6
    std::vector<double> R(n_grid);

    // precompute P(x) * R^6
    double term1 = (z + 1.0) * std::log(z + 1.0);
    double term2 = -std::lgamma(z + 1.0);
    
    for (int i = 0; i < n_grid; i++) {
        x[i] = x_min + i * dx;
        double xi = x[i];
        
        // log P(x) = term1 + term2 + z*ln(x) - (z+1)*x
        double log_Px = term1 + term2 + z * std::log(xi) - (z + 1.0) * xi;
        double Px = std::exp(log_Px);
        
        R[i] = xi * r_avg;
        Px_R6[i] = Px * std::pow(R[i], 6);
    }


    for (size_t i = 0; i < n_q; i++) {
        double q_val = ptr_q[i];
        double integral = 0.0;
        
        // Inner Loop: Integration
        for (int j = 0; j < n_grid; j++) {
            double P_qr = form_factor_sphere(q_val, R[j]);
            integral += Px_R6[j] * P_qr;
        }
        integral *= dx;
        
        ptr_res[i] = scale * integral + bg;
    }

    return result;
}

PYBIND11_MODULE(saxs_core, m) {
    m.doc() = "Optimized Calculation in C++";
    m.def("model_intensity", &model_intensity, "Calculate Schulz Sphere Intensity",
          py::arg("q"), py::arg("r_avg"), py::arg("z"), py::arg("scale"), py::arg("bg"));
}
