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
    
    int N = 200;
    int M = 200;
    int L = 150;
    int Dim = 2;

    double Tem_0 = 20;
    double Tem_1 = 600;
    double dt = 0.005;
    double dx = 0.1;

    A = new double[N * M];
    B = new double[N * M];
    S = new double[2 * Dim + 1];

    // starting point
    for(int i = 1; i < N-1; ++i) {
        for(int j = 1; j < M-1; ++j) {
            A[i * M + j] = Tem_0;
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

    for(int i = 0; i < 800000; ++i) {
        sweep(A, S, N, M, Dim);
        set_border(A, N, M, L, Tem_0, Tem_1);
        if(i%4000 == 0)
            save_matrix(A, N, M, "./RESULT/img_" + std::to_string(i/4000) + ".csv");
    }

    printf("\n------------------------------------\n");
    // print_matrix(A, N, M);

    delete [] A;
    delete [] B;
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

void set_border(double *A, int N, int M, int L, double Tem_0, double Tem_1) {
    // Hot plate dirichlet
    for(int j = 1; j < L+1; ++j) {
        A[0 * M + j] = Tem_1;
    }
    A[0] = Tem_0;

    // Air dirichlet
    for(int j = L+1; j < M; ++j) {
        A[0 * M + j] = Tem_0;
    }
    for(int j = 0; j < M; ++j) {
        A[(N-1) * M + j] = Tem_0;
    }
    for(int i = 1; i < N-1; ++i) {
        A[i*M + 0] = Tem_0;
        A[i*M + M-1] = Tem_0;
    }
}

void sweep(double *A, double *S, int N, int M, int Dim) {
    
    
    double *A_1 = new double[N * M];
    double *temp;

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

void save_matrix(double * A, int N, int M, std::string path) {
    std::ofstream file;
    file.open(path);
    for(int i = 1; i < N-1; ++i) {
        for(int j = 1; j < M-2; ++j) {
            file << A[i * M + j] << ";";
        }
        file << A[i * M + M-2] << "\n";
    }
    file.close();

}