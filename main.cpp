#include "classFFT.h"
#include <iostream>
#include <fstream>
#include <set>

const std::string absError = "AbsoluteError.txt";
const std::string relativeError = "RelativeError.txt";

std::vector<complexD> GenerateRandomPoints(int, int);/// Генерирует случайным образом n равномерно распределённых 
                                                        ///пар значений {x,y} в квадрате [-l,l]*[-l,l].

double EuclideanNormVector(const std::vector<complexD> &);/// Евклидова норма вектора.
double AbsoluteError(const std::vector<complexD> &, const std::vector<complexD> &);/// Абсолютная ошибка.
double RelativeError(const std::vector<complexD> &, const std::vector<complexD> &);/// Относительная ошибка.
template<class T1,class T2>
void print(std::ofstream& ,const T1&,const T2 &);/// Форматированный вывод.

int main() {
  
  std::set<int> transformationLengths{2, 3, 4, 5, 8, 9, 16, 25, 27, 32, 64, 81, 125, 128, 243, 256, 
                                  512, 625, 729, 1024, 2048, 2187, 3125, 4096, 6561}; /// множество возможных
                                                                                ///значений длины преобразования.

  std::ofstream fout1{absError};/// Файловый поток для абсолютной погрешности.
  std::ofstream fout2{relativeError};/// Файловый поток для относительной погрешности.

  print(fout1,"Length","Absolute Error");
  print(fout2,"Length","Relative Error");

  for (auto &length : transformationLengths) {

    std::vector<complexD> input = GenerateRandomPoints(length, 10); /// входные данные-комплексные значения (x,y), 
                                                                /// где (x,y) принадлежит квадрату [-10,10]*[-10,10].

    int n = input.size();

    std::vector<complexD> a(input);

    int N = computeFFT::Resize(a);///Ближайшая степень двойки, не меньшая размера исходного вектора a.

    std::vector<complexD> resultDirect = computeFFT::FFT(a); /// прямое преобразование Фурье.

    std::vector<complexD> resultReverse = computeFFT::FFT(resultDirect, true); /// обратное преобразование Фурье.

    for (auto &elem : resultReverse) elem /= N;

    std::vector<complexD> output(resultReverse.begin(),resultReverse.begin() + n); /// выходные данные.

    print(fout1,length,AbsoluteError(input, output)); /// длина преобразования-абсолютная ошибка.

    print(fout2,length,RelativeError(input, output)); /// длина преобразования-относительная ошибка.
  }
  
  return 0;
}


/**
 * @brief Генерирует случайным образом n равномерно распределённых пар значений {x,y} в квадрате [-l,l]*[-l,l] .
 * 
 * @param n Количество комплексных значений.
 * @param l Размер квадрата [-l,l]*[-l,l],в котором генерируются пары значений {x,y}.
 * @return Вектор комплексных значений.
 */
std::vector<complexD> GenerateRandomPoints(int n, int l) {
  std::vector<complexD> complexPoints;
  std::set<std::pair<double, double>> used;
  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<double> dist(-l, l);

  for (int i = 0; i < n; ++i) {/// Создаётся ровно n различных пар значений
    double x, y;
    do {
      x = dist(rng);
      y = dist(rng);
    } while (!used.insert({x, y}).second);/// Условие на отсутствие дубликатов
    complexPoints.emplace_back(x, y);
  }

  return complexPoints;
}

double EuclideanNormVector(const std::vector<complexD> &v) {
  double sum = 0.0;
  for (const auto &x : v)
    sum += norm(x);
  return sqrt(sum);
}

double AbsoluteError(const std::vector<complexD> &v1, const std::vector<complexD> &v2) {
  assert(v1.size() == v2.size());
  int n = v1.size();
  std::vector<complexD> diffVectors(n);
  for (int i = 0; i < n; ++i)
    diffVectors[i] = v1[i] - v2[i];
  return EuclideanNormVector(diffVectors);
}

double RelativeError(const std::vector<complexD> &v1, const std::vector<complexD> &v2) {
  return AbsoluteError(v1, v2) / EuclideanNormVector(v1);
}

template<class T1,class T2>
void print(std::ofstream& of,const T1& val1,const T2 & val2){
    of<< std::setw(6) << val1 << "\t" << std::scientific << std::setprecision(5) << val2<< "\n";
}
