#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include <complex>
#include <random>
#include <set>
#include <vector>
#include <iomanip>
#include <numbers>


using complexD = std::complex<double>;

/**
* @brief Вычисляет ближайшую степень двойки, не меньшую заданного числа.
* 
* @param num Входное число.
* @return Ближайшая степень двойки, не меньшая num.
*/
inline double nearestPowerOfTwo(double num) { return pow(2, ceil(log2(num))); }

class computeFFT {
public:
/**
 * @brief Увиличивает размер исходного вектора до степени двойки.
 * 
 * @param a Исходный вектор.
 * @return Ближайшая степень двойки, не меньшая размера исходного вектора.
 */
  static double Resize(std::vector<complexD> &a) {
    int n = a.size();
    double N = nearestPowerOfTwo(n);
    if (N > n)
      a.resize(N, complexD(0, 0));
    return N;
  }

/**
 * @brief Вычисляет быстрое преобразование Фурье.
 * 
 * @param a Исходный вектор значений.
 * @param reverse Если false — выполняется прямое БПФ; если true — обратное БПФ.
 * @return Новый вектор комплексных чисел — результат прямого или обратного БПФ.
 */
  static std::vector<complexD> FFT(std::vector<complexD> &a, bool reverse = false) {

    int N = a.size();

    if (N == 1) return a;
    
    complexD w_n = exp(complexD(0, pow(-1, reverse + 1) * 2 * std::numbers::pi / N));
    complexD w = 1;

    std::vector<complexD> a_even(N / 2), a_odd(N / 2);
    for (int i = 0; i < N / 2; i++) {
      a_even[i] = a[2 * i];
      a_odd[i] = a[2 * i + 1];
    }

    std::vector<complexD> y_even = FFT(a_even, reverse);
    std::vector<complexD> y_odd = FFT(a_odd, reverse);

    std::vector<complexD> y(N);
    for (int k = 0; k < N / 2; k++) {
      auto t = w * y_odd[k];
      y[k] = y_even[k] + t;
      y[k + N / 2] = y_even[k] - t;
      w *= w_n;
    }

    return y;
  }
};