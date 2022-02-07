#include <stdio.h>
#include "connect.h"
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <immintrin.h>

#pragma GCC optimize(3)

int main(int argc, char **argv) {
    struct timeval a, b, t_strassen, t_strassen_f, tv, f_tv;
    gettimeofday(&a, NULL);
    char *fileA = NULL;
    char *fileB = NULL;
    char *fileC = NULL;
    if (argc > 1) {
        fileA = argv[1];
        fileB = argv[2];
        fileC = argv[3];
    }
    p_Matrix A = createMatrix(fileA);
    p_Matrix B = createMatrix(fileB);
    FILE *fp = fopen(fileC, "w");

    printf("A's transposition cost: \n");
    p_Matrix tran_A = transpose(A->m_matrix, A->m_row, A->m_column);
    if (tran_A != NULL) {
        fileWriter(fp, tran_A);
        deleteMatrix(tran_A);
    }

    printf("B's transposition cost: \n");
    p_Matrix tran_B = transpose(B->m_matrix, B->m_row, B->m_column);
    if (tran_B != NULL) {
        fileWriter(fp, tran_B);
        deleteMatrix(tran_B);
    }

    p_Matrix sub = subtraction(A->m_matrix, B->m_matrix, A->m_row, A->m_column, B->m_row, B->m_column);
    if (sub != NULL) {
        fileWriter(fp, sub);
        deleteMatrix(sub);
    }

    p_Matrix pul = plus(A->m_matrix, B->m_matrix, A->m_row, A->m_column, B->m_row, B->m_column);
    if (pul != NULL) {
        fileWriter(fp, pul);
        deleteMatrix(pul);
    }

    p_Matrix ladder_A = rank(A->m_matrix, A->m_row, A->m_column);
    A->m_rank = ladder_A->m_rank;
    printf("The A's rank is: %d \n", A->m_rank);
    fileWriter(fp, ladder_A);
    deleteMatrix(ladder_A);

    p_Matrix ladder_B = rank(B->m_matrix, B->m_row, B->m_column);
    B->m_rank = ladder_B->m_rank;
    printf("The B's rank is: %d \n", B->m_rank);
    fileWriter(fp, ladder_B);
    deleteMatrix(ladder_B);

    if(A->m_rank != A->m_row)
        A->hasLU = False;
    if(B->m_rank != B->m_row)
        B->hasLU = False;

    printf("A's determinant is: \n");
    det(A);
    printf("B's determinant is: \n");
    det(B);

    if(A->hasLU == True){
        p_Matrix low_A = (Matrix *) malloc(sizeof(Matrix));
        p_Matrix up_A = (Matrix *) malloc(sizeof(Matrix));
        low_A->m_row = A->m_row;
        up_A->m_row = A->m_row;
        low_A->m_column = A->m_column;
        up_A->m_column = A->m_column;
        low_A->m_matrix = (float *) malloc(sizeof(float) * A->m_row * A->m_column);
        up_A->m_matrix = (float *) malloc(sizeof(float) * A->m_row * A->m_column);
        for(int i = 0 ; i < A->m_row ; ++i){
            for(int j = 0 ; j < A->m_column ; ++j){
                *(low_A->m_matrix + i * A->m_column + j) = 0.f;
                *(up_A->m_matrix + i * A->m_column + j) = 0.f;
            }
        }
        printf("A's LU factorization cost: \n");
        LU(A->m_matrix, low_A->m_matrix, up_A->m_matrix, A->m_row);
        fileWriter(fp, low_A);
        fileWriter(fp, up_A);
        deleteMatrix(low_A);
        deleteMatrix(up_A);
    }
    if(B->hasLU == True){
    p_Matrix low_B = (Matrix *) malloc(sizeof(Matrix));
    p_Matrix up_B = (Matrix *) malloc(sizeof(Matrix));
    low_B->m_row = B->m_row;
    up_B->m_row = B->m_row;
    low_B->m_column = B->m_column;
    up_B->m_column = B->m_column;
    low_B->m_matrix = (float *) malloc(sizeof(float) * B->m_row * B->m_column);
    up_B->m_matrix = (float *) malloc(sizeof(float) * B->m_row * B->m_column);
    for(int i = 0 ; i < B->m_row ; ++i){
            for(int j = 0 ; j < B->m_column ; ++j){
                *(low_B->m_matrix + i * B->m_column + j) = 0.f;
                *(up_B->m_matrix + i * B->m_column + j) = 0.f;
        }
    }
    printf("B's LU factorization cost: \n");
    LU(B->m_matrix, low_B->m_matrix, up_B->m_matrix, B->m_row);
    fileWriter(fp, low_B);
    fileWriter(fp, up_B);
    deleteMatrix(low_B);
    deleteMatrix(up_B);
    }

    printf("A's inversion cost: \n");
    p_Matrix inv_A = inverse(A);
    if (inv_A != NULL) {
        fileWriter(fp, inv_A);
        deleteMatrix(inv_A);
    }

    printf("B's inversion cost: \n");
    p_Matrix inv_B = inverse(B);
    if (inv_B != NULL) {
        fileWriter(fp, inv_B);
        deleteMatrix(inv_B);
    }

    p_Matrix mul = multiply(A->m_matrix, B->m_matrix, A->m_row, A->m_column, B->m_row, B->m_column);
    if (mul != NULL) {
        fileWriter(fp, mul);
        deleteMatrix(mul);
    }

    p_Matrix C = (Matrix *) malloc(sizeof(Matrix));
    C->m_column = B->m_column;
    C->m_row = A->m_row;
    C->m_matrix = (float *) malloc(sizeof(float) * C->m_row * C->m_column);
    for(int i = 0 ; i < A->m_row ; ++i){
           for(int j = 0 ; j < B->m_column ; ++j){
               *(C->m_matrix + i * B->m_column + j) = 0.f;
        }
    }
    if (A->m_column == B->m_row) {
        strassen(A->m_matrix, B->m_matrix, C->m_matrix, A->m_row, A->m_column, B->m_column, A->m_row);
        fileWriter(fp, C);
    } else printf("MULTIPLY FAILURE!");

    // gettimeofday(&tv, NULL);
    // multiply_func(A->m_matrix, B->m_matrix, C->m_matrix, A->m_row, A->m_column, B->m_row, B->m_column);
    // gettimeofday(&f_tv, NULL);
    // printf("%.3f \n", (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000);
    for(int i = 0 ; i < A->m_row ; ++i){
        for(int j = 0 ; j < B->m_column ; ++j){
               *(C->m_matrix + i * B->m_column + j) = 0.f;
        }
    }
    mulSIMD(C, A, B);
    fileWriter(fp, C);
    deleteMatrix(C);
    deleteMatrix(A);
    deleteMatrix(B);
    fclose(fp);
    gettimeofday(&b, NULL);
    printf("The running time is %.3f ms\n",
           (b.tv_sec - a.tv_sec) * (double) 1000 + (b.tv_usec - a.tv_usec) / (double) 1000);
    return 0;
}