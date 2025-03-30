#include <iostream>
#include <cmath>
#include <gmp.h>
#include <bitset>
// Евклидово расстояние в double (IEEE 754)
double euclidean_distance_double(double x1, double y1, double z1, 
                                 double x2, double y2, double z2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Евклидово расстояние с высокой точностью (GMP)
void euclidean_distance_gmp(mpf_t result, double x1, double y1, double z1, 
                                             double x2, double y2, double z2) {
    mpf_t dx, dy, dz, sum, temp;
    mpf_inits(dx, dy, dz, sum, temp, NULL);
    
    // Устанавливаем точность GMP (больше 53 бит для double)
    mpf_set_default_prec(128);

    // Вычисляем разности
    mpf_set_d(dx, x2 - x1);
    mpf_set_d(dy, y2 - y1);
    mpf_set_d(dz, z2 - z1);

    // Квадраты разностей
    mpf_mul(temp, dx, dx);
    mpf_set(sum, temp);

    mpf_mul(temp, dy, dy);
    mpf_add(sum, sum, temp);

    mpf_mul(temp, dz, dz);
    mpf_add(sum, sum, temp);

    // Квадратный корень
    mpf_sqrt(result, sum);

    mpf_clears(dx, dy, dz, sum, temp, NULL);
}
// Функция для сравнения битов двух double
bool compare_bits(double d1, double d2) {
    // Получаем представление чисел в виде битов
    std::bitset<64> bits1(*reinterpret_cast<unsigned long long*>(&d1));
    std::bitset<64> bits2(*reinterpret_cast<unsigned long long*>(&d2));
    return bits1 == bits2;
}

// Сравнение результатов double vs GMP
void compare_distances(double x1, double y1, double z1, 
                       double x2, double y2, double z2) {
    // Вычисление в double
    double d_double = euclidean_distance_double(x1, y1, z1, x2, y2, z2);

    // Вычисление в GMP
    mpf_t d_gmp;
    mpf_init(d_gmp);
    euclidean_distance_gmp(d_gmp, x1, y1, z1, x2, y2, z2);

    // Перевод результата GMP в double
    double d_gmp_double = mpf_get_d(d_gmp);

    // Проверка равенства (сравниваем с точностью до младшего значащего бита)
    if (compare_bits(d_double, d_gmp_double)) {
        std::cout << "Результаты совпадают с точностью до последнего значащего бита!\n";
    } else {
        std::cout << "Ошибка: расхождение в значащих битах!\n";
        std::cout << "Double: " << d_double << "\n";
        std::cout << "GMP: " << d_gmp_double << "\n";
    }

    // Проверка симметрии расстояния d(A, B) == d(B, A)
    double d_swapped = euclidean_distance_double(x2, y2, z2, x1, y1, z1);
    if (d_double == d_swapped) {
        std::cout << "Проверка симметрии пройдена: d(A, B) == d(B, A)\n";
    } else {
        std::cout << "Ошибка: d(A, B) != d(B, A)\n";
    }

    mpf_clear(d_gmp);
}

int main() {
    double x1 = 1.123456789, y1 = 2.987654321, z1 = -3.141592653;
    double x2 = -1.987654321, y2 = 0.123456789, z2 = 2.718281828;

    compare_distances(x1, y1, z1, x2, y2, z2);

    return 0;
}
