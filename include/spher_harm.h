#pragma once

/*
  Library to generate spherical harmonics and its derivative

  Author: Martin Horvat, December 2020
*/

#include <cmath>
#include <complex>

#include "constants.h"
#include "matrix.h"

/*
  Calculate normalized associated Legendre function (ALF)

    bar P_l^m(theta) = N_l^m P_l^(cos(theta))

    N_l^m = sqrt( (2l + 1) (l - m)!/ ( 4 pi (l+m)!))

    (this is our normalization and other are used in literature)

  and its derivatives wrt theta and use them to determine
  spherical harmonics (SH)

    Y_l^m(theta,phi) = bar P_l^m(theta) exp(I m phi)

  and its derivatives. Associated Legendre polynomials P_l^m are

    P_0^0 (x) = 1   P_1^0 (x) =  x    P_1^1(x) = - (1-x^2)^(1/2) ....

  or more generally

    P_l^m(x) = (-1)^m (1-x^2)^(m/2) d^m/dx^m P_l

    P_l = 1/(2^l l!) d^l/dx^l (x^2-1)^l

  Ref:
  * https://en.wikipedia.org/wiki/Spherical_harmonics
  * https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
  * Holmes, S. & Featherstone, Will. (2002). Journal of Geodesy (2002) 76: 279–299
  * Xing, Z., Li, S., Tian, M. et al. J Geod 94, 2 (2020).
  * T Limpanuparb, J Milthorpe. (2014) Proceedings of The 40th Congress on Science
    and Technology of Thailand; 2014 Dec 2-4, Khon Kaen, Thailand. Bangkok:
    Science Society of Thailand; 2014. P. 233-241
*/

template <class T>
struct TSpherHarm {

  using C = std::complex<T>;

  int n;

  T **A, **B, *a, *b;

  void init(int n) {

    if (this -> n != n) destroy();

    if (n) {

      A = tmatrix<T>(n);
      B = tmatrix<T>(n);

      a = new T [n+1];
      b = new T [n+1];

      for (int l = 2, l_ = 1; l <= n; ++l, ++l_) {

        int ls = l * l, lm1s = l_ * l_;

        for (int m = 0; m < l_; ++m) {

          int ms = m * m;

          A[l][m] = std::sqrt( T(4*ls - 1) / (ls - ms));
          B[l][m] = -std::sqrt( T(lm1s - ms) / (4*lm1s - 1));
        }

        a[l] = std::sqrt(T(2*l + 1));
        b[l] = std::sqrt(1 + T(0.5)/l);
      }
    }

    this-> n = n;
  }

  void destroy(){
    if (this -> n) {

      delete [] a;
      delete [] b;

      a = b = 0;

      free_tmatrix(A);
      free_tmatrix(B);

      this -> n = 0;
    }
  }

  /*
    Precompute coefficients and reserve space

    Input:
      n - maximal degree
  */
  TSpherHarm(int n) : n(0), A(0), B(0), a(0), b(0) { init(n);}

  TSpherHarm(): n(0), A(0), B(0), a(0), b(0) { }

  ~TSpherHarm() { destroy(); }


  /*
    Compute an entire set of ALFs

      bar P_l^m(theta)

    for degree l = 0, .., n, order m = 0,...,l  and store them in P.

    Using column-wise recurrence formulas in P[l][m] table, where l is on vertical,
    m on horizontal axis.

    Input:
      theta in [0,pi]

    Output:
      P - triangle matrix of ALFs up to degree n
  */

  void calc_P(const T & t, const T &u , T **P) {

    T temp = math_constants<T>::sqrt1div4pi;

    P[0][0] = temp;

    if (n > 0) {

      P[1][0] = t * math_constants<T>::sqrt3 * temp;

      temp *= -math_constants<T>::sqrt3div2 * u;

      P[1][1] = temp;

      for (int l = 2, l_ = 1; l <= n; ++l, ++l_) {

        for (int m = 0; m < l_; ++m)
          P[l][m] = A[l][m]* (t *P[l_][m] + B[l][m]* P[l_-1][m]);

        P[l][l_] = t * a[l] * temp;
        P[l][l] = (temp *= -b[l]*u);
      }
    }
  }

  inline void calc_P(const T & theta, T **P) {

    T t = std::cos(theta),
      u = std::sqrt(1 - t*t); // = sin(theta)

    calc_P(t, u, P);
  }


  /*
    Compute an entire set of derivatives of normalized ALFs

      d/dtheta bar P_l^m(theta)

    for degree l = 0, .., n, order m = 0,...,l and store them in dP.

    Input:
      theta in [0,pi]
      P -  tringle matrix of normalized ALFs up to degree n

    Output:
      dP - triangle matrix of derivative of normalized ALFs up to degree n
  */
  void calc_dP(const T & t, const T & u,  T **P, T **dP) {

    T temp = math_constants<T>::sqrt1div4pi;

    dP[0][0] = 0;

    if (n > 0) {

      dP[1][0] = -u * math_constants<T>::sqrt3 * temp;
      dP[1][1] = -t * math_constants<T>::sqrt3div2 * temp;

      for (int l = 2, l_ = 1; l <= n; ++l, ++l_) {

        for (int m = 0; m < l_; ++m)
          dP[l][m] = A[l][m]* (-u*P[l_][m] + t*dP[l_][m] + B[l][m]* dP[l_-1][m]);

        dP[l][l_] = a[l]*(-u*P[l_][l_] + t*dP[l_][l_]);
        dP[l][l] = -b[l]*(t*P[l_][l_] + u*dP[l_][l_]);
      }
    }
  }

  inline void calc_dP(const T & theta, T **P, T **dP) {
    T t = std::cos(theta),
      u = std::sqrt(1 - t*t); // = sin(theta)

    calc_dP(t, u, P, dP);
  }


  /*
    Compute an entire set of spherical harmonics (SH)

      Y_l^m(x)

    for degree l = 0, .., n, order m = 0,...,l and store them in Y.

    Input:
      phi in [0,pi]
      P - triangle matrix of normalized ALFs up to degree n

    Output:
      Y - triangle matrix of complex SHs up to degree n
  */
  void calc_Y(const T &t, const T &u, T **P, C **Y) {

    for (int l = 0; l <= n; ++l) Y[l][0] = P[l][0];

    T c1 = 1, c2 = t,
      s1 = 0, s2 = -u,
      tc = 2* c2, s, c;

    for (int m = 1; m <= n; ++m) {
      s = tc * s1 - s2;
      c = tc * c1 - c2;

      s2 = s1,
      s1 = s;
      c2 = c1;
      c1 = c;

      for (int l = m; l <= n; ++l) Y[l][m] = P[l][m]*C(c,s);
    }
  }

  inline void calc_Y(const T & phi, T **P, C **Y) {
    calc_Y(std::cos(phi), std::sin(phi), P, Y);
  }

  /*
    Compute an entire set of spherical harmonics (SH) and derivatives wrt phi and theta

      Y_l^m(phi,theta)   d/dtheta Y_l^m(phi,theta) d/dphi Y_l^m(phi,theta)

    for degree l = 0, .., n order m = 0,...,l and store them in Y, dYdphi and dYdth.

    Input:
      theta in [0,pi]
      phi in [0,2pi]

      P - triangle matrix of normalied ALFs up to degree n
      dP - triangle matrix of derivatives ALFs up to degree n

    Output:
      Y - triangle matrix of complex SHs up to degree n
      dYdth - triangle matrix of derivatives of SHs wrt theta up to degree n
      dYdphi - triangle matrix of derivatives of SHs wrt phi up to degree n
  */

  void calc_Y(const T &t, const T &u, T **P, T **dP, C **Y, C **dYdth, C **dYdphi) {

    for (int l = 0; l <= n; ++l)  {
      Y[l][0] = P[l][0];
      dYdth[l][0] = dP[l][0];
      dYdphi[l][0] = 0;
    }

    T c1 = 1, c2 = t,
      s1 = 0, s2 = -u,
      tc = 2*c2, s, c, tmp;

    for (int m = 1; m <= n; ++m) {
      s = tc * s1 - s2;
      c = tc * c1 - c2;

      s2 = s1,
      s1 = s;
      c2 = c1;
      c1 = c;

      for (int l = m; l <= n; ++l) {
        Y[l][m] = (tmp = P[l][m])*C(c,s);
        dYdth[l][m] = dP[l][m]*C(c,s);
        dYdphi[l][m] = m*tmp*C(-s, c);
      }
    }
  }

  inline void calc_Y(const T &phi, T **P, T **dP, C **Y, C **dYdth, C **dYdphi) {
    calc_Y(std::cos(phi), std::sin(phi), P, dP, Y, dYdth, dYdphi);
  }
};


