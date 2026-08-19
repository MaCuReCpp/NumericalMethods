#include "errors.hpp"

#include <cmath>
#include <stdexcept>

namespace numerical {

double absolute_error(double exact, double approximate)
{
    return std::abs(exact - approximate);
}

double relative_error(double exact, double approximate)
{
    if (exact == 0.0) {
        throw std::invalid_argument("relative error is undefined when exact value is zero");
    }
    return std::abs((exact - approximate) / exact);
}

double l2_norm(const Vector& values)
{
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

} // namespace numerical
