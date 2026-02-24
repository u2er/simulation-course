#include <stdlib.h>

void solve_heat_step(double* T, int N, double h, double tau, double rho, double c_heat, double lambda_coef, int steps) {
    double* alpha = (double*)malloc(N * sizeof(double));
    double* beta  = (double*)malloc(N * sizeof(double));

    double A = lambda_coef / (h * h);
    double C = lambda_coef / (h * h);
    double B = 2.0 * lambda_coef / (h * h) + (rho * c_heat) / tau;
    double const_F = (rho * c_heat) / tau;

    for (int s = 0; s < steps; s++) {
        alpha[0] = 0.0;
        beta[0]  = T[0]; 

        for (int i = 1; i < N - 1; i++) {
            double F_i = -const_F * T[i];
            double denominator = B - C * alpha[i - 1];
            
            alpha[i] = A / denominator;
            beta[i]  = (C * beta[i - 1] - F_i) / denominator;
        }

        for (int i = N - 2; i > 0; i--) {
            T[i] = alpha[i] * T[i + 1] + beta[i];
        }
    }

    free(alpha);
    free(beta);
}