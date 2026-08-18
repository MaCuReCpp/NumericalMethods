#include "ode.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>

namespace numerical {

Grid1D linear_shooting_bvp(
    const ScalarFunction& p,
    const ScalarFunction& q,
    const ScalarFunction& r,
    double a,
    double b,
    double alpha,
    double beta,
    std::size_t steps)
{
    if (steps == 0) {
        throw std::invalid_argument("steps must be positive");
    }
    const double h = (b - a) / static_cast<double>(steps);

    auto integrate = [&](double y0, double dy0, bool include_r) {
        Vector y(steps + 1, 0.0);
        Vector dy(steps + 1, 0.0);
        y[0] = y0;
        dy[0] = dy0;
        for (std::size_t i = 0; i < steps; ++i) {
            const double t = a + h * static_cast<double>(i);
            auto f1 = [](double velocity) { return velocity; };
            auto f2 = [&](double x, double value, double velocity) {
                return p(x) * velocity + q(x) * value + (include_r ? r(x) : 0.0);
            };
            const double k1y = f1(dy[i]);
            const double k1v = f2(t, y[i], dy[i]);
            const double k2y = f1(dy[i] + 0.5 * h * k1v);
            const double k2v = f2(t + 0.5 * h, y[i] + 0.5 * h * k1y, dy[i] + 0.5 * h * k1v);
            const double k3y = f1(dy[i] + 0.5 * h * k2v);
            const double k3v = f2(t + 0.5 * h, y[i] + 0.5 * h * k2y, dy[i] + 0.5 * h * k2v);
            const double k4y = f1(dy[i] + h * k3v);
            const double k4v = f2(t + h, y[i] + h * k3y, dy[i] + h * k3v);
            y[i + 1] = y[i] + h * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
            dy[i + 1] = dy[i] + h * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0;
        }
        return y;
    };

    const Vector y1 = integrate(alpha, 0.0, true);
    const Vector y2 = integrate(0.0, 1.0, false);
    const double scale = (beta - y1.back()) / y2.back();
    Grid1D grid;
    grid.points.resize(steps + 1);
    grid.values.resize(steps + 1);
    for (std::size_t i = 0; i <= steps; ++i) {
        grid.points[i] = a + h * static_cast<double>(i);
        grid.values[i] = y1[i] + scale * y2[i];
    }
    return grid;
}

std::vector<OdePoint> runge_kutta_fehlberg45(
    const OdeFunction& f,
    double t0,
    double y0,
    double t_end,
    double initial_step,
    double tolerance,
    double min_step,
    double max_step)
{
    if (initial_step <= 0.0 || min_step <= 0.0 || max_step <= 0.0 || tolerance <= 0.0) {
        throw std::invalid_argument("steps and tolerance must be positive");
    }

    std::vector<OdePoint> solution = {{t0, y0}};
    double t = t0;
    double y = y0;
    double h = std::min(initial_step, max_step);

    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }

        const double k1 = h * f(t, y);
        const double k2 = h * f(t + h / 4.0, y + k1 / 4.0);
        const double k3 = h * f(t + 3.0 * h / 8.0, y + 3.0 * k1 / 32.0 + 9.0 * k2 / 32.0);
        const double k4 = h * f(t + 12.0 * h / 13.0, y + 1932.0 * k1 / 2197.0 - 7200.0 * k2 / 2197.0 + 7296.0 * k3 / 2197.0);
        const double k5 = h * f(t + h, y + 439.0 * k1 / 216.0 - 8.0 * k2 + 3680.0 * k3 / 513.0 - 845.0 * k4 / 4104.0);
        const double k6 = h * f(t + h / 2.0, y - 8.0 * k1 / 27.0 + 2.0 * k2 - 3544.0 * k3 / 2565.0 + 1859.0 * k4 / 4104.0 - 11.0 * k5 / 40.0);

        const double y4 = y + 25.0 * k1 / 216.0 + 1408.0 * k3 / 2565.0 + 2197.0 * k4 / 4104.0 - k5 / 5.0;
        const double y5 = y + 16.0 * k1 / 135.0 + 6656.0 * k3 / 12825.0 + 28561.0 * k4 / 56430.0 - 9.0 * k5 / 50.0 + 2.0 * k6 / 55.0;
        const double error = std::abs(y5 - y4);

        if (error <= tolerance || h <= min_step) {
            t += h;
            y = y5;
            solution.push_back({t, y});
        }

        const double factor = (error == 0.0) ? 4.0 : 0.84 * std::pow(tolerance / error, 0.25);
        h = std::clamp(factor * h, min_step, max_step);
    }

    return solution;
}

std::vector<OdePoint> adams_bashforth_moulton4(
    const OdeFunction& f,
    double t0,
    double y0,
    double step,
    std::size_t steps)
{
    if (step <= 0.0) {
        throw std::invalid_argument("step must be positive");
    }
    if (steps <= 3) {
        return runge_kutta4(f, t0, y0, step, steps);
    }

    std::vector<OdePoint> solution = runge_kutta4(f, t0, y0, step, 3);
    for (std::size_t n = 3; n < steps; ++n) {
        const double fn = f(solution[n].t, solution[n].y);
        const double fn1 = f(solution[n - 1].t, solution[n - 1].y);
        const double fn2 = f(solution[n - 2].t, solution[n - 2].y);
        const double fn3 = f(solution[n - 3].t, solution[n - 3].y);
        const double predictor = solution[n].y + step * (55.0 * fn - 59.0 * fn1 + 37.0 * fn2 - 9.0 * fn3) / 24.0;
        const double next_t = solution[n].t + step;
        const double f_predictor = f(next_t, predictor);
        const double corrector = solution[n].y + step * (9.0 * f_predictor + 19.0 * fn - 5.0 * fn1 + fn2) / 24.0;
        solution.push_back({next_t, corrector});
    }
    return solution;
}

} // namespace numerical
