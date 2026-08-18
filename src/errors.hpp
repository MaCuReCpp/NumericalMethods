#pragma once

#include "types.hpp"

namespace numerical {

double absolute_error(double exact, double approximate);
double relative_error(double exact, double approximate);
double l2_norm(const Vector& values);

} // namespace numerical
