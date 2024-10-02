
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


long double factorial(int);
long double count_c(int, int);
void count_bn(long double [], int);
long double sum_of_series(long double, int, long double []);
long double get_pow(long double, int);

void count_bn_wof(long double [], int);
long double sum_of_series_wof(long double, long double, long double []);

int main() {
    long double arg = 0;
    long double number = 0;

    printf("Введите аргумент косеканса и кол_во членов ряда Тейлора для расчёта:\n");
    scanf("%Lf %Lf", &arg, &number);

    if ((fabsl(arg) == 0) || (fabsl(arg) >= M_PI)) {
        printf("Модуль аргумента должен соответствовать условию: 0 < |x| < pi\n");
        return 1;
    }
    
    if (number < 1) {
        printf("Введите НАТУРАЛЬНОЕ кол-во членов ряда\n");
        return 2;
    }

    long double *b_list = malloc(sizeof(long double) * (number * 2 + 3)); 
    long double *b_list_wof = malloc(sizeof(long double) * (number * 2 + 3));

    clock_t start_time_1 = clock();
    printf("Через ряд Тейлора без факториала: %.10Lf\n", sum_of_series_wof(arg, number, b_list_wof));
    clock_t time_result_1 = clock() - start_time_1;
    printf("Затраченное время: %lf\n", (double)time_result_1 / CLOCKS_PER_SEC);

    clock_t start_time_2 = clock();
    printf("Через ряд Тейлора с факторилом: %.10Lf\n", sum_of_series(arg,  number, b_list));
    clock_t time_result_2 = clock() - start_time_2;
    printf("Затраченное время: %.10lf\n", (double)time_result_2 / CLOCKS_PER_SEC);

    printf("Через встроенную библиотеку: %.10Lf\n", (long double) 1 / sinh(arg));
/*  
    printf("Через встроенную библиотеку: %.8lf\n", 1 / sinh(arg));
    printf("Через ряды с факториалом: %.8Lf\n", sum_of_series(arg, number, b_list));
    printf("Через ряды без факториала: %.8Lf\n", sum_of_series_wof(arg, number, b_list_wof));
*/
    free(b_list);
    free(b_list_wof);
    return 0;
}


long double factorial(int n){
    long double result = 1;
    for (int i = 1; i < n + 1; ++i) {
        result *= i;
    }

    return result;
}


long double count_c(int n, int k) {
    long double result = factorial(n) / factorial(k) / factorial(n - k);
    return result; 
}


long double get_pow(long double x, int y) {
    long double result = 1;
    if (y > 0){
        for (int i = 0; i < y; ++i) {
        result *= x;
        }
    }

    else {
        for (int i = 0; i < (-y); ++i){
            result /= x;
        }
    } 

    return result;
}


long double sum_of_series(long double arg, int number, long double *b_list) {
    long double result = 0;

    for (int n = 0; n < number + 1; ++n) {
        count_bn(b_list, 2 * n);
        result += -2 * (get_pow(2, 2 * n - 1) - 1) * b_list[2 * n] * get_pow(arg, 2 * n - 1) / factorial(2 * n);
    }

    return result;
}


void count_bn(long double *b_list, int n) {
    if (n == 0) {
        b_list[0] = 1.0;
        b_list[1] = -0.5;
    }

    else {
        long double bn = 0.0;
        
        for (int k = 0; k < n + 1; ++k) {
            bn += b_list[n - k] * count_c(n + 1, k + 1);
        }

        bn *= (-1.0 / (n + 1));
        b_list[n] = bn;
        b_list[n + 1] = 0.0;
    }
}


void count_bn_wof(long double *b_list_wof, int n) {
    if (n == 0) {
        b_list_wof[0] = 1;
        b_list_wof[1] = -0.5;
    }

    else {
        long double bn = 0;
        long double c = ((long double)n * ((long double)n + 1) / 2);

        for (int k = 1; k < n + 1; k++) {
            bn += b_list_wof[n - k] * c;
            c = c * (long double)(n - k) / (long double)(k + 2);

        }

        b_list_wof[n] = bn * (-1 / ((long double)n + 1));
        b_list_wof[n + 1] = 0;
    }
}


long double sum_of_series_wof(long double arg, long double number, long double *b_list_wof) {
    long double current_pow = -2;
    long double previous_pow = 1;
    
    long double result = 1 / arg;
    long double current_step = 1 / arg;

    count_bn_wof(b_list_wof, 0);
    long double k_pow = 0, k_bn = 0, all_k = 0; 
    for (int n = 1; n < number + 1; n++) {
        count_bn_wof(b_list_wof, 2 * n);
        
        k_pow = current_pow / previous_pow;
        k_bn = b_list_wof[2 * n] / b_list_wof[2 * (n - 1)];
        all_k = arg * arg / (2 * n) / (2 * n - 1) * k_bn * k_pow; 
        
        current_step *= all_k;
        result += current_step;
        
        printf("%d     %.10Lf    %.10Lf     %.10Lf     %.10Lf    %.10Lf\n", n,k_pow , k_bn,  all_k, current_step, result);

        previous_pow = current_pow;
        current_pow = (current_pow - 2) * 4 + 2;
    }

    return result;
}

