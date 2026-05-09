#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>

void print_matrix(double *A, int N, int M);
void set_border(double *A, int N, int M, int L, double Tem_0, double Tem_1);
void sweep(double *A, double *S, int N, int M, int Dim);
void save_matrix(double * A, int N, int M, std::string path);

int main() {

    double *A, *B, *S;
    
    int N = 100;
    int M = 210;
    int L = 30;
    int Dim = 2;

    double Tem_0 = 20;
    double Tem_1 = 600;
    double dt = 0.005;
    double dx = 0.1;

    A = new double[N * M];
    S = new double[2 * Dim + 1];

    // starting point
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j) {
            A[i * M + j] = 20;
        }
    }

    // stencil weights
    double alpha = 0.01;
    S[2* Dim] = 1 - 2 * alpha * dt * (Dim / (dx * dx));
    for(int i = 0; i < 2*Dim; ++i) {
        S[i] = (alpha * dt) / (dx * dx);
    }

    set_border(A, N, M, L, Tem_0, Tem_1);
    // print_matrix(A, N, M);
    // printf("-------------------------------------------------------------------\n");

    for(int i = 0; i < 800000; ++i) {
        sweep(A, S, N, M, Dim);
        set_border(A, N, M, L, Tem_0, Tem_1);    
        // print_matrix(A, N, M);
        // printf("-------------------------------------------------------------------\n");
        if(i%4000 == 0)
            save_matrix(A, N, M, "./RESULT/img_" + std::to_string(i/4000) + ".csv");
    }

    printf("\n------------------------------------\n");
    // print_matrix(A, N, M);

    delete [] A;
    return 0;
}



void print_matrix(double *A, int N, int M) {
    for(int i = 0; i < N; ++i) {        
        printf("[%6.2f", A[i*M]);
        for(int j = 1; j < M; ++j) {
            printf(" %6.2f", A[i*M + j]);
        }
        printf("]\n");
    }
}

void sweep(double *A, double *S, int N, int M, int Dim) {
    
    double a_0 = 0.01;
    double a_1 = 0.0005;
    double dx = 0.1;
    double dt = 0.005;
    double Tem_0 = 20;
    double *A_1 = new double[N * M];
    double *temp;

    // Neuman top
    for(int j = 1; j < M-1; ++j) {
        A_1[0 * M + j] = dt * (a_0 - a_1) * (A[j] - Tem_0) / (dx*dx) +
                         a_0 * dt * (Tem_0 - 4*A[j] + A[M + j] + A[j-1] + A[j+1]) / (dx*dx) +
                         A[0 * M + j];
    }
    // Neuman bottom
    for(int j = 1; j < M-1; ++j) {
        A_1[(N-1) * M + j] = dt * (a_1 - a_0) * (Tem_0 - A[(N-1)*M + j]) / (dx*dx) +
                             a_0 * dt * (Tem_0 - 4*A[(N-1)*M + j] + A[(N-2)*M + j] + A[(N-1)*M + j-1] + A[(N-1)*M + j+1]) / (dx*dx) +
                             A[(N-1)*M + j];
    }
    // Neuman sides
    for(int i = 1; i < N-1; ++i) {
        A_1[i*M + 0] = dt * (a_0 - a_1) * (A[i*M] - Tem_0) / (dx*dx) +
                       a_0 * dt * (Tem_0 - 4*A[i*M] + A[i*M + 1] + A[(i-1)*M] + A[(i+1) * M]) / (dx*dx) +
                       A[i*M + 0];

        A_1[i*M + M-1] = dt * (a_1 - a_0) * (Tem_0 - A[i*M + M-1]) / (dx*dx) +
                         a_0 * dt * (Tem_0 - 4*A[i*M + M-1] + A[i*M + M-2] + A[(i-1)*M + M-1] + A[(i+1)*M + M-1]) / (dx*dx) +
                         A[i*M + M-1];
    }
    // Neuman Corners
    A_1[0] = dt * (a_0 - a_1) * (A[0] - Tem_0) / (dx*dx) +
             dt * (a_0 - a_1) * (A[0] - Tem_0) / (dx*dx) +
             dt *  a_0 * (Tem_0 + Tem_0 + A[M] + A[1] - 4 * A[0]) / (dx*dx) +
             A[0];

    A_1[M-1] = dt * (a_0 - a_1) * (A[M-1] - Tem_0) / (dx*dx) +        // y axis
               dt * (a_1 - a_0) * (Tem_0 - A[M-1]) / (dx*dx) +        // x axis
               dt * a_0 * (Tem_0 + Tem_0 + A[M + M-1] + A[M-2] - 4 * A[M-1]) / (dx*dx) + 
               A[M-1];

    A_1[(N-1) * M] = dt * (a_1 - a_0) * (Tem_0 - A[(N-1) * M]) / (dx*dx) +        // y axis
                     dt * (a_0 - a_1) * (A[(N-1) * M] - Tem_0) / (dx*dx) +        // x axis
                     dt * a_0 * (Tem_0 + Tem_0 + A[(N-2) * M] + A[(N-1)*M + 1] - 4 * A[(N-1) * M]) / (dx*dx) + 
                     A[(N-1) * M];

    A_1[(N-1) * M + M-1] = dt * (a_1 - a_0) * (Tem_0 - A[(N-1)*M + M-1]) / (dx*dx) +        // y axis
                           dt * (a_1 - a_0) * (Tem_0 - A[(N-1)*M + M-1]) / (dx*dx) +        // x axis
                           dt * a_0 * (Tem_0 + Tem_0 + A[(N-2) * M + M-1] + A[(N-1)*M + M-2] - 4 * A[(N-1)*M + M-1]) / (dx*dx) + 
                           A[(N-1) * M];

    // Sweeping inside frame
    double value;
    for(int row = 1; row < N-1; ++row) {
        for(int col = 1; col < M-1; ++col) {
            
            value = A[row * M + col] * S[2*Dim];
            for(int i = 0; i < 2*Dim; ++i) {
                value += A[
                            (row + (i > 1)*(i < 4)*(i%2 * 2 - 1)) * M + 
                            (col + (i > -1)*(i < 2)*(i%2 * 2 - 1))
                        ] * S[i];
            }
            A_1[row * M + col] = value;
        }
    }

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j) {
            A[i * M + j] = A_1[i * M + j];
        }
    }

    delete [] A_1;
}


void set_border(double *A, int N, int M, int L, double Tem_0, double Tem_1) {
    // Top dirichlet
    for(int j = 30; j < L + 30; ++j) {
        A[j] = Tem_1;
    }
    for(int j = M-31; j > M-31-L; --j) {
        A[j] = Tem_1;
    }
}

void save_matrix(double * A, int N, int M, std::string path) {
    std::ofstream file;
    file.open(path);
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M-1; ++j) {
            file << A[i * M + j] << ";";
        }
        file << A[i * M + M-1] << "\n";
    }
    file.close();

}