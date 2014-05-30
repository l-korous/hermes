#include "definitions.h"

double CustomExactSolution::value (double x, double y) const 
{
  double slope2 = slope * 3.;
  return Hermes::atan(slope * (Hermes::sqrt(Hermes::sqr(x - 1.25) + Hermes::sqr(y + 0.25)) - M_PI / 3.)) +
    10. * Hermes::atan(slope2 * (Hermes::sqrt(Hermes::sqr(x - 0.25) + Hermes::sqr(y + 0.15)) - M_PI / 3.));
}

void CustomExactSolution::derivatives (double x, double y, double& dx, double& dy) const 
{
  double slope2 = slope * 3.;
  double t = Hermes::sqrt(Hermes::sqr(x - 1.25) + Hermes::sqr(y + 0.25));
  double u = t * (Hermes::sqr(slope) * Hermes::sqr(t - M_PI/3.) + 1);
  dx = slope * (x - 1.25) / u;
  dy = slope * (y + 0.25) / u;
  t = Hermes::sqrt(Hermes::sqr(x - 0.25) + Hermes::sqr(y + 0.25));
  u = t * (Hermes::sqr(slope2) * Hermes::sqr(t - M_PI / 3.) + 1);
  dx += 10. * slope2 * (x - 0.25) / u;
  dy += 10. * slope2 * (y + 0.15) / u;
}

Ord CustomExactSolution::ord(double x, double y) const 
{
  return Ord(20);
}


double CustomFunction::value(double x, double y) const 
{
  double slope2 = slope * 3.;
  double t2 = Hermes::sqr(y + 0.25) + Hermes::sqr(x - 1.25);
  double t = Hermes::sqrt(t2);
  double u = (Hermes::sqr(M_PI - 3.0*t) * Hermes::sqr(slope) + 9.0);

  double to_return = (27.0 / 2.0 * Hermes::sqr(2.0*y + 0.5) * (M_PI - 3.0*t) * Hermes::pow(slope, 3.0) / (Hermes::sqr(u) * t2)
    + 27.0 / 2.0 * Hermes::sqr(2.0*x - 2.5) * (M_PI - 3.0*t) * Hermes::pow(slope, 3.0) / (Hermes::sqr(u) * t2)
    - 9.0 / 4.0 * Hermes::sqr(2.0*y + 0.5) * slope / (u * Hermes::pow(t, 3.0))
    - 9.0 / 4.0 * Hermes::sqr(2.0*x - 2.5) * slope / (u * Hermes::pow(t, 3.0)) + 18.0 * slope / (u * t)
    );

  double t21 = Hermes::sqr(y + 0.15) + Hermes::sqr(x - 0.25);
  double t1 = Hermes::sqrt(t2);
  double u1 = (Hermes::sqr(M_PI - 3.0*t) * Hermes::sqr(slope2) + 9.0);

  to_return += 10. * (27.0 / 2.0 * Hermes::sqr(2.0*y + 0.3) * (M_PI - 3.0*t1) * Hermes::pow(slope2, 3.0) / (Hermes::sqr(u1) * t21)
    + 27.0 / 2.0 * Hermes::sqr(2.0*x - 0.5) * (M_PI - 3.0*t1) * Hermes::pow(slope2, 3.0) / (Hermes::sqr(u1) * t21)
    - 9.0 / 4.0 * Hermes::sqr(2.0*y + 0.3) * slope2 / (u1 * Hermes::pow(t1, 3.0))
    - 9.0 / 4.0 * Hermes::sqr(2.0*x - 0.5) * slope2 / (u1 * Hermes::pow(t1, 3.0)) + 18.0 * slope2 / (u1 * t1)
    );

  return to_return + std::min(1e5, (1. / Hermes::sqrt(Hermes::sqr(Hermes::sqr(x - 0.25)) + Hermes::sqr(Hermes::sqr(y - 0.6)))))
    + std::min(1e5, (1. / Hermes::sqrt(Hermes::sqr(Hermes::sqr(x - 0.23)) + Hermes::sqr(Hermes::sqr(y - 0.5)))))
    + std::min(1e5, (1. / Hermes::sqrt(Hermes::sqr(Hermes::sqr(x - 0.22)) + Hermes::sqr(Hermes::sqr(y - 0.4)))))
    + std::min(1e5, (1. / Hermes::sqrt(Hermes::sqr(Hermes::sqr(x - 0.21)) + Hermes::sqr(Hermes::sqr(y - 0.3)))));
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(20);
}
