#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <complex>
#include <random>
#include <set>
#include <vector>
#include <iomanip>

const double PI = 3.14159265358979323846;

using complexD = std::complex<double>;

inline double nearestPowerOfTwo(double num) { return pow(2, ceil(log2(num))); }

template <class ElemT> class computeFFT {
public:
  static int Resize(std::vector<std::complex<ElemT>> &a) {
    int n = a.size();
    int N = nearestPowerOfTwo(n);
    if (N > n)
      a.resize(N, complexD(0, 0));
    return N;
  }

  static std::vector<std::complex<ElemT>> FFT(std::vector<std::complex<ElemT>> &a, bool reverse = false) {

    int N = a.size();

    if (N == 1) return a;
    
    std::complex<ElemT> w_n = exp(std::complex<ElemT>(0, pow(-1, reverse + 1) * 2 * PI / N));
    std::complex<ElemT> w = 1;

    std::vector<std::complex<ElemT>> a_even(N / 2), a_odd(N / 2);
    for (int i = 0; i < N / 2; i++) {
      a_even[i] = a[2 * i];
      a_odd[i] = a[2 * i + 1];
    }

    std::vector<std::complex<ElemT>> y_even = FFT(a_even, reverse);
    std::vector<std::complex<ElemT>> y_odd = FFT(a_odd, reverse);

    std::vector<std::complex<ElemT>> y(N);
    for (int k = 0; k < N / 2; k++) {
      auto t = w * y_odd[k];
      y[k] = y_even[k] + t;
      y[k + N / 2] = y_even[k] - t;
      w *= w_n;
    }

    return y;
  }
};