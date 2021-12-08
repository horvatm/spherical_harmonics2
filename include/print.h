#pragma once

/*
  Auxillary routines

  Author: Martin Horvat, December 2021
*/

template <class T> T Linf_norm(const std::complex<T> & x){
  return std::max(std::abs(x.real()), std::abs(x.imag()));
}

template <class T>
void calc_error(std::complex<T> Y, std::complex<T> Y_e, T & de_abs, T & de_rel){

  T dt = Linf_norm(Y_e - Y),
    t = std::max(Linf_norm(Y_e), Linf_norm(Y)),
    eps = std::numeric_limits<T>::epsilon();

  de_rel = (t != 0 && dt*(1 + eps) <= t?  dt/t : 0);
  de_abs = dt;
}

template <class T> void check_values(
  std::string filename_values,
  std::string filename_full_report,
  std::string filename_short_report,
  std::string filename_summary){

  using C = std::complex<T>;

  std::ifstream values_file(filename_values);
  std::ofstream
    full_report(filename_full_report),
    short_report(filename_short_report),
    summary(filename_summary);

  // configure reports
  full_report << std::scientific;
  full_report.precision(std::numeric_limits<T>::digits10);

  short_report << std::scientific;
  short_report.precision(std::numeric_limits<T>::digits10);

  summary << std::scientific;
  summary.precision(std::numeric_limits<T>::digits10);

  full_report << "#filename_values:" << filename_values << '\n';
  short_report <<  "#filename_values:" << filename_values << '\n';
  summary << "#filename_values:" << filename_values << '\n';

  // reading maximal order and number of points that we will test
  int max_order, nr_points;
  values_file >> max_order >> nr_points;

  // initialize spherical harmonics
  TSpherHarm<T> sp(max_order);

  T **P = tmatrix<T>(max_order),
    **dP = tmatrix<T>(max_order);

  C **Y = tmatrix<C>(max_order),
    **dYdth = tmatrix<C>(max_order),
    **dYdphi = tmatrix<C>(max_order);


  // main loop
  T theta, phi, y[2],
    de_abs[3], de_rel[3],
    e_abs[3], e_rel[3],
    e_abs1[3], e_rel1[3];


  C Y_e, dYdth_e, dYdphi_e, dY;

  summary << "#theta phi e_abs_max e_rel_max\n";

  for (int i = 0; i < nr_points; ++i) {

    values_file >> theta >> phi;

    full_report << "#point:" << theta << '\t' << phi << '\n';
    short_report << "#point:" << theta << '\t' << phi << '\n';

    // calcualate values of the spherical harmonics
    auto start = std::chrono::high_resolution_clock::now();

    T t = std::cos(theta), u = std::sqrt(1 - t*t);
    sp.calc_P(t, u, P);
    sp.calc_dP(t, u, P, dP);
    sp.calc_Y(phi, P, dP, Y, dYdth, dYdphi);

    auto end = std::chrono::high_resolution_clock::now();

    // report on calculation and check deviation
    for (int i = 0; i < 3; ++i) e_abs1[i] = e_rel1[i] = 0;

    for (int l = 0; l <= max_order; ++l) {

      for (int i = 0; i < 3; ++i) e_abs[i] = e_rel[i] = 0;

      for (int m = 0; m <= l; ++m) {

        // read correct values
        values_file >> y[0] >> y[1]; Y_e = C(y[0], y[1]);
        values_file >> y[0] >> y[1]; dYdth_e = C(y[0], y[1]);
        values_file >> y[0] >> y[1]; dYdphi_e = C(y[0], y[1]);

        calc_error(Y[l][m], Y_e, de_abs[0], de_rel[0]);
        calc_error(dYdth[l][m], dYdth_e, de_abs[1], de_rel[1]);
        calc_error(dYdphi[l][m], dYdphi_e, de_abs[2], de_rel[2]);

        for (int i = 0; i < 3; ++i){
          e_abs[i] = std::max(e_abs[i], de_abs[i]);
          e_rel[i] = std::max(e_rel[i], de_rel[i]);
        }

        full_report
          << l << '\t' << m << '\t'
          << Y_e << '\t' <<  Y[l][m]  << '\t'
          << dYdth_e << '\t' <<  dYdth[l][m]  << '\t'
          << dYdphi_e << '\t' <<  dYdphi[l][m];
        for (int i = 0; i < 3; ++i) full_report << '\t' << de_abs[i] << '\t' << de_rel[i];
        full_report << '\n';
      }

      short_report << l;
      for (int i = 0; i < 3; ++i) short_report << '\t' << e_abs[i] << '\t' <<  e_rel[i];
      short_report << '\n';

      for (int i = 0; i < 3; ++i){
        e_abs1[i] = std::max(e_abs1[i], e_abs[i]);
        e_rel1[i] = std::max(e_rel1[i], e_rel[i]);
      }
    }

    summary << theta << '\t' << phi;
    for (int i = 0; i < 3; ++i) summary << '\t' << e_abs1[i] << '\t' << e_rel1[i];
    summary
      << '\t'
      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
      << '\n';

    short_report << '\n';
    full_report  << '\n';
  }

  free_tmatrix(P);
  free_tmatrix(dP);
  free_tmatrix(Y);
  free_tmatrix(dYdth);
  free_tmatrix(dYdphi);
}


bool cmp_str(std::string str1, std::string str2){
 return str1.compare(str2)==0;
}
