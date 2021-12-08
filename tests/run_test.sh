#!/bin/bash


echo "Values of SH..."
gzip -d Y_n100.dat.gz

echo "  Running for double precision..."
./test_spher_harm.e Y_n100.dat Y_full_report_d_n100.dat Y_short_report_d_n100.dat Y_summary_d_n100.dat double

echo "  Running for long double precision..."
./test_spher_harm.e Y_n100.dat Y_full_report_ld_n100.dat Y_short_report_ld_n100.dat Y_summary_ld_n100.dat long_double

gzip Y_n100.dat




echo "Values and derivatives of SH..."
gzip -d dY_n100.dat.gz

echo "  Running for double precision..."
./test_dspher_harm.e dY_n100.dat dY_full_report_d_n100.dat dY_short_report_d_n100.dat dY_summary_d_n100.dat double

echo "  Running for long double precision..."
./test_dspher_harm.e dY_n100.dat dY_full_report_ld_n100.dat dY_short_report_ld_n100.dat dY_summary_ld_n100.dat long_double

gzip dY_n100.dat

