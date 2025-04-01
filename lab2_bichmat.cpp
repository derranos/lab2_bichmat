#include <iostream>
#include <cmath>
#include <vector>
#include <gmp.h>
#include <iomanip>

// Функция для поиска первого отличающегося разряда и вывода разницы
void find_first_diff(double d1, double d2) {
    double diff = std::abs(d1 - d2);

    if (diff == 0.0) {
        std::cout << "Числа идентичны.\n";
        return;
    }

    std::cout << std::fixed << std::setprecision(20);
    std::cout << "Разница между числами: " << diff << "\n";
}

// Наивное вычисление евклидова расстояния
double euclidean_distance_naive(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

// Вычисление евклидова расстояния с суммированием по Кахану
double euclidean_distance_kahan(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    double c = 0.0; // компенсация ошибки
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        double term = diff * diff;
        double y = term - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return std::sqrt(sum);
}

// Вычисление евклидова расстояния с парным суммированием
double pairwise_sum(const std::vector<double>& values, size_t start, size_t end) {
    if (end - start == 1) return values[start];
    if (end - start == 0) return 0.0;

    size_t mid = start + (end - start) / 2;
    return pairwise_sum(values, start, mid) + pairwise_sum(values, mid, end);
}

double euclidean_distance_pairwise(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> squared_diffs(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        squared_diffs[i] = diff * diff;
    }
    double sum = pairwise_sum(squared_diffs, 0, squared_diffs.size());
    return std::sqrt(sum);
}

// Вычисление евклидова расстояния с высокой точностью (GMP)
double euclidean_distance_gmp(const std::vector<double>& a, const std::vector<double>& b) {
    mpf_set_default_prec(256);

    mpf_t sum, diff, term;
    mpf_inits(sum, diff, term, NULL);
    mpf_set_ui(sum, 0);

    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        mpf_set_d(diff, d);
        mpf_mul(term, diff, diff);
        mpf_add(sum, sum, term);
    }

    mpf_sqrt(sum, sum);
    double result = mpf_get_d(sum);

    mpf_clears(sum, diff, term, NULL);
    return result;
}

// Функция проверки симметрии d(A,B) == d(B,A)
void check_symmetry(const std::vector<double>& a, const std::vector<double>& b) {
    double d1 = euclidean_distance_naive(a, b);
    double d2 = euclidean_distance_naive(b, a);
    if (d1 == d2)
        std::cout << "Симметрия пройдена (наивное суммирование).\n";
    else
        std::cout << "Ошибка симметрии (наивное суммирование): разница = " << std::abs(d1 - d2) << "\n";

    d1 = euclidean_distance_kahan(a, b);
    d2 = euclidean_distance_kahan(b, a);
    if (d1 == d2)
        std::cout << "Симметрия пройдена (Кахан).\n";
    else
        std::cout << "Ошибка симметрии (Кахан): разница = " << std::abs(d1 - d2) << "\n";

    d1 = euclidean_distance_pairwise(a, b);
    d2 = euclidean_distance_pairwise(b, a);
    if (d1 == d2)
        std::cout << "Симметрия пройдена (Pairwise).\n";
    else
        std::cout << "Ошибка симметрии (Pairwise): разница = " << std::abs(d1 - d2) << "\n";

    d1 = euclidean_distance_gmp(a, b);
    d2 = euclidean_distance_gmp(b, a);
    if (d1 == d2)
        std::cout << "Симметрия пройдена (GMP).\n";
    else
        std::cout << "Ошибка симметрии (GMP): разница = " << std::abs(d1 - d2) << "\n";
}

void test_multiplication_error() {
    const size_t N = 1000000;
    std::vector<double> a(N), b(N);

    for (size_t i = 0; i < N; ++i) {
        a[i] = 1e5 + i * 1e-3;
        b[i] = 1e5 + (i + 0.1) * 1e-3;
    }

    double d_naive = euclidean_distance_naive(a, b);
    double d_kahan = euclidean_distance_kahan(a, b);
    double d_pairwise = euclidean_distance_pairwise(a, b);
    double d_gmp = euclidean_distance_gmp(a, b);

    std::cout.precision(20);
    std::cout << std::fixed;
    std::cout << "Euclidean distance (наивное суммирование): " << d_naive << "\n";
    std::cout << "Euclidean distance (суммирование по Кахану): " << d_kahan << "\n";
    std::cout << "Euclidean distance (Pairwise): " << d_pairwise << "\n";
    std::cout << "Euclidean distance (GMP, высокая точность): " << d_gmp << "\n";

    std::cout << "\nСравнение (наивное vs Кахан):\n";
    find_first_diff(d_naive, d_kahan);

    std::cout << "\nСравнение (наивное vs Pairwise):\n";
    find_first_diff(d_naive, d_pairwise);

    std::cout << "\nСравнение (наивное vs GMP):\n";
    find_first_diff(d_naive, d_gmp);

    check_symmetry(a, b);
}

// Тест на малых числах
void test_small_multiplication_error() {
    const size_t N = 1000000;
    std::vector<double> a(N), b(N);

    for (size_t i = 0; i < N; ++i) {
        a[i] = 1e-5 + i * 1e-7;
        b[i] = 1e-5 + (i + 0.1) * 1e-7;
    }

    double d_naive = euclidean_distance_naive(a, b);
    double d_kahan = euclidean_distance_kahan(a, b);
    double d_pairwise = euclidean_distance_pairwise(a, b);
    double d_gmp = euclidean_distance_gmp(a, b);

    std::cout.precision(20);
    std::cout << std::fixed;
    std::cout << "\nEuclidean distance (наивное суммирование): " << d_naive << "\n";
    std::cout << "Euclidean distance (суммирование по Кахану): " << d_kahan << "\n";
    std::cout << "Euclidean distance (Pairwise): " << d_pairwise << "\n";
    std::cout << "Euclidean distance (GMP, высокая точность): " << d_gmp << "\n";

    std::cout << "\nСравнение (наивное vs Кахан):\n";
    find_first_diff(d_naive, d_kahan);

    std::cout << "\nСравнение (наивное vs Pairwise):\n";
    find_first_diff(d_naive, d_pairwise);

    std::cout << "\nСравнение (наивное vs GMP):\n";
    find_first_diff(d_naive, d_gmp);

    check_symmetry(a, b);
}

int main() {
    std::cout << "Тест с большими числами:\n";
    test_multiplication_error();

    std::cout << "\nТест с малыми числами:\n";
    test_small_multiplication_error();

    return 0;
}
