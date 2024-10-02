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
    long double epsilon = 0;

    printf("Введите аргумент косеканса и точность для расчёта:\n");
    scanf("%Lf %Lf", &arg, &epsilon);

    if ((fabsl(arg) == 0) || (fabsl(arg) >= M_PI)) {
        printf("Модуль аргумента должен соответствовать условию: 0 < |x| < pi\n");
        return 1;
    }
    
    if (epsilon < 0) {
        printf("Введите точность больше нуля\n");
        return 2;
    }

    long double *b_list_wof = malloc(sizeof(long double) * (300));

    clock_t start_time_1 = clock();
    printf("Через ряд Тейлора без факториала: %.10Lf\n", sum_of_series_wof(arg, epsilon, b_list_wof));
    clock_t time_result_1 = clock() - start_time_1;
    printf("Затраченное время: %lf\n", (double)time_result_1 / CLOCKS_PER_SEC);

    printf("Через встроенную библиотеку:      %.10lf\n", 1 / sinh(arg));

    free(b_list_wof);
    return 0;
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


long double sum_of_series_wof(long double arg, long double epsilon, long double *b_list_wof) {
    long double current_pow = -2;
    long double previous_pow = 1;
    
    long double result = 1 / arg;
    long double current_step = 1 / arg;
    long double previous_step = 1 / arg;

    count_bn_wof(b_list_wof, 0);
    long double k_pow = 0, k_bn = 0, all_k = 0; 
    int n = 1;

    while (fabsl(current_step + previous_step) >= epsilon) {
        count_bn_wof(b_list_wof, 2 * n);
        
        k_pow = current_pow / previous_pow;
        k_bn = b_list_wof[2 * n] / b_list_wof[2 * (n - 1)];
        all_k = arg * arg / (2 * n) / (2 * n - 1) * k_bn * k_pow; 
        
        previous_step = current_step;

        current_step *= all_k;
        result += current_step;
        
        printf("%d     %.10Lf    %.10Lf     %.10Lf     %.10Lf    %.10Lf\n", n,k_pow , k_bn,  all_k, current_step, result);

        previous_pow = current_pow;
        current_pow = (current_pow - 2) * 4 + 2;
        
        n++;
    }
    
    printf("Посчитано членов рядов Тейлора: %d\n", n);

    return result;
}
