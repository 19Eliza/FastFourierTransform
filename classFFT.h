#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include <complex>
#include <random>
#include <vector>
#include <iomanip>
#include <numbers>

using complexD = std::complex<double>;
const double PI = std::numbers::pi;

/**
* @brief Вычисляет прямое DFT
* 
* @param input Входной вектор.
* @return Новый вектор комплексных чисел — результат прямого ДПФ.
*/
std::vector<complexD> DFT(const std::vector<complexD>& input) {
  size_t N = input.size();
  std::vector<complexD> output(N);

  for (size_t k = 0; k < N; ++k) {
    complexD sum = 0;
    for (size_t n = 0; n < N; ++n) {
      double angle = -2 * PI * k * n / N;
      complexD w(std::cos(angle), std::sin(angle));
      sum += input[n] * w;
    }
    output[k] = sum;
  }

  return output;
}

/**
* @brief Вычисляет обратное DFT
* 
* @param input Входной вектор.
* @return Новый вектор комплексных чисел — результат обратного ДПФ.
*/
std::vector<complexD> IDFT(const std::vector<complexD>& input) {
  size_t N = input.size();
  std::vector<complexD> output(N);

  for (size_t n = 0; n < N; ++n) {
    complexD sum = 0;
    for (size_t k = 0; k < N; ++k) {
      double angle = 2 * PI * k * n / N;
      complexD w(std::cos(angle), std::sin(angle)); 
      sum += input[k] * w;
    }
    output[n] = sum;
  }

  return output;
}

/**
* @brief Вычисляет разложение числа на простые множители.
* 
* @param N Входное число.
* @return Вектор разложения на множители числа N.
*/
std::vector<int> factorize(int N) {
    std::vector<int> factors;
    for (int i = 2; i * i <= N; ++i)
        while (N % i == 0) {
            factors.push_back(i);
            N /= i;
        }
    if (N > 1) factors.push_back(N);
    return factors;
}

class computeFFT {
public:
/**
 * @brief Вычисляет быстрое преобразование Фурье.
 * 
 * @param input Исходный вектор значений.
 * @param reverse Если false — выполняется прямое БПФ; если true — обратное БПФ.
 * @return Новый вектор комплексных чисел — результат прямого или обратного БПФ.
 */
  static std::vector<complexD> FFT(const std::vector<complexD>& input, bool reverse = false) {
    int N = input.size();
    auto factors = factorize(N);
    
    if (factors.size() < 2) {/// вычисление ДПФ, если число простое
        if(!reverse) return DFT(input); 
        return IDFT(input);
    }

    int N1 = factors[0];
    int N2 = N / N1;

    std::vector<std::vector<complexD>> v(N1, std::vector<complexD>(N2));
    for (int i1 = 0; i1 < N1; ++i1)
        for (int i2 = 0; i2 < N2; ++i2)
            v[i1][i2] = input[i1 + N1 * i2];

    for (int i1 = 0; i1 < N1; ++i1)
        v[i1] = FFT(v[i1], reverse);

    double sign = reverse ? 1 : -1;
    for (int i1 = 0; i1 < N1; ++i1)
        for (int i2 = 0; i2 < N2; ++i2) {
            double angle = sign * 2 * PI * i1 * i2 / N;
            complexD w(std::cos(angle), std::sin(angle));
            v[i1][i2] *= w;
        }

    std::vector<std::vector<complexD>> vTranspose(N2, std::vector<complexD>(N1));
    for (int i1 = 0; i1 < N1; ++i1)
        for (int i2 = 0; i2 < N2; ++i2)
            vTranspose[i2][i1] = v[i1][i2];

    for (int i2 = 0; i2 < N2; ++i2)
        vTranspose[i2] = FFT(vTranspose[i2], reverse);

    std::vector<complexD> output(N);
    for (int i1 = 0; i1 < N1; ++i1)
        for (int i2 = 0; i2 < N2; ++i2)
            output[i2 + N2 * i1] = vTranspose[i2][i1];

    return output;
  }
};