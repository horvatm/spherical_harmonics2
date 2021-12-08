#pragma once

/*
  Header providing frequently used mathematical constants.

  Author: Martin Horvat, March 2021
*/

#include <cmath>


template<typename T> struct math_constants{};

//specialization to double precision
template<>
struct math_constants<double>
{
  static constexpr double pi =          3.14159265358979323846;
  static constexpr double twopi =       6.28318530717958647692;
  static constexpr double fourpi =      12.56637061435917295385;
  static constexpr double sqrt2 =       1.41421356237309504880;
  static constexpr double sqrt3 =       1.73205080756887729352;
  static constexpr double sqrt3div2 =   1.22474487139158904909;
  static constexpr double sqrt1div2pi = 0.3989422804014326779399;
  static constexpr double sqrt1div4pi = 0.2820947917738781434740;
  static constexpr double sqrt4pi     = 3.54490770181103205459;
  static constexpr double sqrt2pi     = 2.50662827463100050241;

  static constexpr double one_half = 1.0/2.0;
};

//specialization to long double precision
template<>
struct math_constants<long double>
{
  static constexpr long double pi = 3.1415926535897932384626433832795028842L;
  static constexpr long double twopi = 6.28318530717958647692528676655900576839L;
  static constexpr long double fourpi = 12.56637061435917295385057353311801679935L;
  static constexpr long double sqrt2 = 1.41421356237309504880168872420969807857L;
  static constexpr long double sqrt3 = 1.73205080756887729352744634150587236694L;
  static constexpr long double sqrt3div2 = 1.22474487139158904909864203735294569598L;
  static constexpr long double sqrt1div2pi = 0.39894228040143267793994605993437944938L;
  static constexpr long double sqrt1div4pi = 0.28209479177387814347403972578038615569L;
  static constexpr long double sqrt4pi     = 3.54490770181103205459633496668229110787L;
  static constexpr long double sqrt2pi     = 2.50662827463100050241576528481104577787L;

  static constexpr long double one_half = 1.0/2.0L;
};

//specializations to float precision
template<>
struct math_constants<float>
{
  static constexpr float pi =          3.14159265358979323846f;
  static constexpr float twopi =       6.28318530717958647692f;
  static constexpr float fourpi =      12.566370614359172953850f;
  static constexpr float sqrt2 =       1.41421356237309504880f;
  static constexpr float sqrt3 =       1.73205080756887729352f;
  static constexpr float sqrt3div2 =   1.22474487139158904909f;
  static constexpr float sqrt1div2pi = 0.3989422804014326779399f;
  static constexpr float sqrt1div4pi = 0.2820947917738781434740f;
  static constexpr float sqrt4pi     = 3.54490770181103205459f;
  static constexpr float sqrt2pi     = 2.50662827463100050241f;

  static constexpr float one_half = 1.0/2.0f;
};
