#include <stdio.h>
#include "connect.h"
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <immintrin.h>
#pragma GCC optimize(3)
struct timeval tv, f_tv;

p_Matrix createMatrix(char *fileName) {
    p_Matrix matrix = (Matrix *) malloc(sizeof(Matrix));
    matrix->m_row = 0;
    matrix->m_column = 0;
    matrix->m_rank = 0;
    matrix->hasInversion = True;
    matrix->hasLU = True;
    int num = 0;
    char c;
    FILE *file = fopen(fileName, "rb");
    while (fscanf(file, "%c", &c) != EOF) {
        if (c == '\n') {
            matrix->m_row++;
            num++;
        }
        if (c == ' ')
            num++;
    }
    num++;
    matrix->m_column = num / matrix->m_row;
    int i = 0;
    matrix->m_matrix = (float *) malloc(sizeof(float) * num);
    FILE *fid = fopen(fileName, "rb");
    while (!feof(fid)) {
        fscanf(fid, "%f", matrix->m_matrix + (i++));
    }
    fclose(file);
    fclose(fid);
    return matrix;
}

void deleteMatrix(p_Matrix matrix) {
    free(matrix->m_matrix);
    free(matrix);
}

p_Matrix copyMatrix(p_Matrix matrix) {
    p_Matrix newMatrix = (Matrix *) malloc(sizeof(Matrix));
    newMatrix->m_row = matrix->m_row;
    newMatrix->m_column = matrix->m_column;
    newMatrix->m_matrix = (float *) malloc(sizeof(float) * matrix->m_column * matrix->m_row);
    const float *a = matrix->m_matrix;
    float *b = newMatrix->m_matrix;
    for (int i = 0; i < matrix->m_row * matrix->m_column; ++i) {
        *(b++) = *(a++);
    }
    return newMatrix;
}

p_Matrix multiply(const float *A, const float *B, int A_row, int A_column, int B_row, int B_column) {
    if (A_column != B_row) {
        printf("MULTIPLY FAILURE!");
        return NULL;
    }
    p_Matrix C = (Matrix *) malloc(sizeof(Matrix));
    C->m_row = A_row;
    C->m_column = B_column;
    C->m_matrix = (float *) malloc(sizeof(float) * A_row * B_column);
    float *tem = C->m_matrix;
    for (int i = 0; i < A_row; ++i) {
        for (int j = 0; j < B_column; ++j) {
            *(tem++) = 0.0F;
        }
    }
    float *C_mat = C->m_matrix;
    gettimeofday(&tv, NULL);
    float c;
    int i, j, k;
#pragma omp parallel for
    for (i = 0; i < A_row; i++) {
        for (k = 0; k < A_column; k++) {
            c = *(A + i * A_column + k);
            for (j = 0; j < B_column; j += 4) {
                *(C_mat + i * B_column + j) += c * *(B + k * B_column + j);
                *(C_mat + i * B_column + j + 1) += c * *(B + k * B_column + j + 1);
                *(C_mat + i * B_column + j + 2) += c * *(B + k * B_column + j + 2);
                *(C_mat + i * B_column + j + 3) += c * *(B + k * B_column + j + 3);
            }
        }
    }
    gettimeofday(&f_tv, NULL);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of multiplication costs %.3f ms\n", period);
    return C;
}

void multiply_func(const float *A, const float *B, float *C, int A_row, int A_column, int B_row, int B_column) {
    if (A_column != B_row) {
        printf("MULTIPLY FAILURE!");
        return;
    }
//    gettimeofday(&tv, NULL);
    float c;
    int i, j, k;
#pragma omp parallel for
    for (i = 0; i < A_row; i++) {
        for (k = 0; k < A_column; k++) {
            c = *(A + i * A_column + k);
            for (j = 0; j < B_column; j += 4) {
                *(C + i * B_column + j) += c * *(B + k * B_column + j);
                *(C + i * B_column + j + 1) += c * *(B + k * B_column + j + 1);
                *(C + i * B_column + j + 2) += c * *(B + k * B_column + j + 2);
                *(C + i * B_column + j + 3) += c * *(B + k * B_column + j + 3);
            }
        }
    }
//    gettimeofday(&f_tv, NULL);
//    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
//    printf("The calculation of multiplication costs %.3f ms\n", period);
}

p_Matrix plus(const float *A, const float *B, int A_row, int A_column, int B_row, int B_column) {
    if (A_column != B_column || A_row != B_row) {
        printf("PLUS FAILURE!");
        return NULL;
    }
    p_Matrix C = (Matrix *) malloc(sizeof(Matrix));
    C->m_row = A_row;
    C->m_column = A_column;
    C->m_matrix = (float *) malloc(sizeof(float) * C->m_row * C->m_column);
    const float *a = A, *b = B;
    float *c = C->m_matrix;
    int i, j;
    gettimeofday(&tv, NULL);
#pragma omp parallel for
    for (i = 0; i < C->m_row; i++) {
        for (j = 0; j < C->m_column; j += 4) {
            *(c++) = *(a++) + *(b++);
            *(c++) = *(a++) + *(b++);
            *(c++) = *(a++) + *(b++);
            *(c++) = *(a++) + *(b++);
        }
    }
    gettimeofday(&f_tv, NULL);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of addition costs %.3f ms\n", period);
    return C;
}

p_Matrix subtraction(const float *A, const float *B, int A_row, int A_column, int B_row, int B_column) {
    if (A_column != B_column || A_row != B_row) {
        printf("SUBTRACTION FAILURE!");
        return NULL;
    }
    p_Matrix C = (Matrix *) malloc(sizeof(Matrix));
    C->m_row = A_row;
    C->m_column = A_column;
    C->m_matrix = (float *) malloc(sizeof(float) * C->m_row * C->m_column);
    const float *a = A, *b = B;
    float *c = C->m_matrix;
    int i, j;
    gettimeofday(&tv, NULL);
#pragma omp parallel for
    for (i = 0; i < C->m_row; i++) {
        for (j = 0; j < C->m_column; j += 4) {
            *(c++) = *(a++) - *(b++);
            *(c++) = *(a++) - *(b++);
            *(c++) = *(a++) - *(b++);
            *(c++) = *(a++) - *(b++);
        }
    }
    gettimeofday(&f_tv, NULL);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of subtraction costs %.3f ms\n", period);
    return C;
}

p_Matrix transpose(const float *A, int A_row, int A_column) {
    p_Matrix t_matrix = (Matrix *) malloc(sizeof(Matrix));
    t_matrix->m_column = A_row;
    t_matrix->m_row = A_column;
    t_matrix->m_matrix = (float *) malloc(sizeof(float) * A_column * A_row);
    const float *a = A;
    float *t = t_matrix->m_matrix;
    int i, j;
    gettimeofday(&tv, NULL);
#pragma omp parallel for
    for (i = 0; i < t_matrix->m_row; ++i) {
        for (j = 0; j < t_matrix->m_column; j += 4) {
            *(t + i * t_matrix->m_column + j) = *(a + j * A_column + i);
            *(t + i * t_matrix->m_column + j + 1) = *(a + (j + 1) * A_column + i);
            *(t + i * t_matrix->m_column + j + 2) = *(a + (j + 2) * A_column + i);
            *(t + i * t_matrix->m_column + j + 3) = *(a + (j + 3) * A_column + i);
        }
    }
    gettimeofday(&f_tv, NULL);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of transposition costs %.3f ms\n", period);
    return t_matrix;
}

void swap(p_Matrix matrix, int i, int j) {
    float *temp = (float *) malloc(sizeof(float) * matrix->m_column);
    float *t = temp;
    float *m = matrix->m_matrix;
    for (int k = 0; k < matrix->m_column; ++k)
        *(t + k) = *(m + i * matrix->m_column + k);
    for (int k = 0; k < matrix->m_column; ++k)
        *(m + i * matrix->m_column + k) = *(m + j * matrix->m_column + k);
    for (int k = 0; k < matrix->m_column; ++k)
        *(m + j * matrix->m_column + k) = *(t++);
    free(temp);
}

void scaleRow(p_Matrix matrix, int i, float mul) {
    float *m = matrix->m_matrix;
    for (int k = 0; k < matrix->m_column; ++k) {
        *(m + i * matrix->m_column + k) *= mul;
    }
}

void det(p_Matrix matrix) {
    if (matrix->m_column != matrix->m_row) {
        matrix->hasInversion = False;
        printf("DET FAILURE!");
        return;
    }
    p_Matrix mat = copyMatrix(matrix);
    int m_swap = 0;
    float max = 0.0F;
    int max_row_num = 0;
    float temp = 0.0F;
    int det = 1;
    float *m = mat->m_matrix;
    gettimeofday(&tv, NULL);
    for (int i = 0; i < mat->m_row - 1; ++i) {
        max = fabsf(*(m + i * mat->m_column + i));
        max_row_num = i;
        for (int j = i + 1; j < mat->m_row; ++j) {
            if (max < fabsf(*(m + j * mat->m_column + i))) {
                max = fabsf(*(m + j * mat->m_column + i));
                max_row_num = j;
            }
        }
        if (max_row_num != i) {
            swap(mat, i, max_row_num);
            m_swap++;
        }
        for (int j = i + 1; j < mat->m_row; ++j) {
            temp = -*(m + j * mat->m_column + i) / *(m + i * mat->m_column + i);
            for (int k = 0; k < mat->m_column; ++k) {
                *(m + j * mat->m_column + k) += *(m + i * mat->m_column + k) * temp;
            }
        }
    }
    if (m_swap % 2 == 1)
        m_swap = -1;
    else m_swap = 1;
    char res[mat->m_row * 12];
    memset(res, -1, sizeof(char) * mat->m_row * 12);
    Bool directPrint = True;
    int figureSum = 0;
    Bool sign = True;
    int noNegative = 0;
    for (int i = 0; i < mat->m_column; ++i) {
        float mem = mat->m_matrix[i * mat->m_column + i];
        int thisMem = 0;
        long long calMem = llabs((long long) (mem * pow(10, 7)));
        figureSum += 7;
        if (mem == 0)
            matrix->hasInversion = False;
        int isOver = 0;
        int t_mem = (int) mem;
        while (abs(t_mem) > 0) {
            t_mem /= 10;
            isOver++;
        }
        figureSum += isOver;
        if (figureSum <= 9) //不会爆
            det *= (int) calMem;
        else {
            if (directPrint == True) {
                if (!((det < 0 && mem < 0) || (det > 0 && mem > 0)))
                    sign = False;
                sprintf((char *) res, "%lld", (long long) abs(det));
            } else {
                if ((sign == False && mem < 0) || (sign == True && mem > 0))
                    sign = True;
                else sign = False;
            }
            directPrint = False;
            char chs[isOver + 7];
            char *memChars = chs;
            sprintf((char *) memChars, "%lld", calMem);
            for (int j = noNegative; j < figureSum; ++j) {
                if ((int) res[j] > 47)
                    noNegative++;
                else break;
            }
            IntMultiply((char *) res, memChars, noNegative, isOver + 7);
        }
    }
    gettimeofday(&f_tv, NULL);
    if (directPrint == True)
        printf("%d", det * m_swap);
    else {
        for (int j = noNegative; j < figureSum; ++j) {
            if ((int) res[j] > 47)
                noNegative++;
            else break;
        }
        int i = 0;
        while (i < noNegative - 7 * mat->m_row - 1) {
            if (i == 0) {
                if ((m_swap == -1 && sign == True) || (m_swap == 1 && sign == False))
                    printf("%c", '-');
            }
            printf("%c", res[i]);
            i++;
        }
        if ((m_swap == -1 && sign == True) || (m_swap == 1 && sign == False)) {
            if ((int) res[i] >= 53)
                printf("%c", res[noNegative - 7 * mat->m_row - 1] - 1);
            else printf("%c", res[noNegative - 7 * mat->m_row - 1]);
        } else {
            if ((int) res[i] >= 53)
                printf("%c", res[noNegative - 7 * mat->m_row - 1] + 1);
            else printf("%c", res[noNegative - 7 * mat->m_row - 1]);
        }
    }
    printf("\n");
    deleteMatrix(mat);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of determinant costs %.3f ms\n", period);
}

void fileWriter(FILE *fp, p_Matrix matrix) {
    const float *m = matrix->m_matrix;
    for (int i = 0; i < matrix->m_row; ++i) {
        for (int j = 0; j < matrix->m_column; ++j) {
            fprintf(fp, "%.4f ", *(m + i * matrix->m_column + j));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

p_Matrix inverse(p_Matrix matrix) {
    if (matrix->hasInversion == False) {
        printf("INVERSE FAILURE!\n");
        return NULL;
    }
    p_Matrix mat = copyMatrix(matrix);
    p_Matrix inv = (Matrix *) malloc(sizeof(Matrix));
    inv->m_row = mat->m_row;
    inv->m_column = mat->m_column;
    inv->m_matrix = (float *) malloc(sizeof(float) * inv->m_column * inv->m_row);
    float *e = inv->m_matrix;
    float *m = mat->m_matrix;
    for (int i = 0; i < inv->m_row; ++i) {
        for (int j = 0; j < inv->m_column; ++j) {
            if (i == j)
                *(e + i * inv->m_column + j) = 1.0f;
            else *(e + i * inv->m_column + j) = 0.0f;
        }
    }
    float max = 0.0f;
    int max_row_num;
    int m_swap = 0;
    float temp;
    gettimeofday(&tv, NULL);
// #pragma omp parallel for
    for (int i = 0; i < mat->m_row - 1; i++) {
        max = fabsf(*(m + i * mat->m_column + i));
        max_row_num = i;
        for (int j = i + 1; j < mat->m_row; j++) {
            if (max < fabsf(*(m + j * mat->m_column + i))) {
                max = fabsf(*(m + j * mat->m_column + i));
                max_row_num = j;
            }
        }
        if (i != max_row_num) {
            swap(mat, i, max_row_num);
            swap(inv, i, max_row_num);
            m_swap++;
        }
        for (int j = i + 1; j < mat->m_row; j++) {
            temp = -*(m + j * mat->m_column + i) / *(m + i * mat->m_column + i);
            for (int k = 0; k < mat->m_column; k++) {
                *(m + j * mat->m_column + k) += *(m + i * mat->m_column + k) * temp;
                *(e + j * inv->m_column + k) += *(e + i * inv->m_column + k) * temp;
            }
        }
    }
    for (int i = 0; i < mat->m_row; i++) {
        temp = 1 / *(m + i * mat->m_column + i);
        scaleRow(mat, i, temp);
        scaleRow(inv, i, temp);
    }
    int i, j, k;
#pragma omp parallel for
    for (i = mat->m_row - 1; i > 0; i--) {
        for (j = i - 1; j >= 0; j--) {
            temp = -*(m + j * mat->m_column + i) / *(m + i * mat->m_column + i);
            for (k = 0; k < mat->m_column; k++) {
                *(m + j * mat->m_column + k) += temp * *(m + i * mat->m_column + k);
                *(e + j * inv->m_column + k) += temp * *(e + i * mat->m_column + k);
            }
        }
    }
    gettimeofday(&f_tv, NULL);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of inversion costs %.3f ms\n", period);
    deleteMatrix(mat);
    return inv;
}

void IntMultiply(char *a, const char *b, int a_len, int b_len) {
    int c[a_len + b_len];
    memset(c, 0, sizeof(int) * (a_len + b_len));
    int i, j;
#pragma omp parallel for
    for (i = 0; i < a_len; ++i) {
        for (j = 0; j < b_len; ++j) {
            c[i + j + 1] += (a[i] - 48) * (b[j] - 48);
        }
    }
    for (i = a_len + b_len - 1; i > 0; --i) {
        if (c[i] > 9) {
            c[i - 1] = c[i - 1] + c[i] / 10;
            c[i] = c[i] % 10;
        }
    }
    int zeroNum = 0;
    j = 0;
    for (i = 0; i < a_len + b_len; ++i) {
        if (c[i] == 0)
            zeroNum++;
        if (zeroNum != i + 1) {
            a[j] = (char) (c[i] + 48);
            j++;
        }
    }
}

void strassenMul(const float *a, const float *b, float *c, int A_row, int A_column, int B_column) {
    for (int i = 0; i < A_row; ++i) {
        for (int j = 0; j < B_column; ++j) {
            *(c + i * B_column + j) = 0.0F;
        }
    }
    float t;
    int i, j, k;
#pragma omp parallel for
    for (i = 0; i < A_row; i++) {
        for (k = 0; k < A_column; k++) {
            t = *(a + i * A_column + k);
            for (j = 0; j < B_column; j += 4) {
                *(c + i * B_column + j) += t * *(b + k * B_column + j);
                *(c + i * B_column + j + 1) += t * *(b + k * B_column + j + 1);
                *(c + i * B_column + j + 2) += t * *(b + k * B_column + j + 2);
                *(c + i * B_column + j + 3) += t * *(b + k * B_column + j + 3);
            }
        }
    }
}

void strassen(const float *a, const float *b, float *c, int a_row, int a_column, int b_column, int out_decided) {
    struct timeval t_strassen, t_strassen_f;
    if ((a_column % 2 != 0 || b_column % 2 != 0 || a_row % 2 != 0) || a_row < 64) {
        gettimeofday(&t_strassen, NULL);
        strassenMul(a, b, c, a_row, a_column, b_column);
        gettimeofday(&t_strassen_f, NULL);
        if (out_decided == 32)
            printf("The calculation of multiplication using strassen costs %.3f ms (%d)\n",
                   (t_strassen_f.tv_sec - t_strassen.tv_sec) * (double) 1000 +
                   (t_strassen_f.tv_usec - t_strassen.tv_usec) / (double) 1000, out_decided);
        return;
    }
    int A_r = a_row;
    int A_c = a_column;
    int B_c = b_column;
    a_row /= 2;
    b_column /= 2;
    a_column /= 2;

    float *M1 = (float *) malloc(sizeof(float) * a_row * b_column);
    float *M2 = (float *) malloc(sizeof(float) * a_row * b_column);
    float *M3 = (float *) malloc(sizeof(float) * a_row * b_column);
    float *M4 = (float *) malloc(sizeof(float) * a_row * b_column);
    float *M5 = (float *) malloc(sizeof(float) * a_row * b_column);
    float *M6 = (float *) malloc(sizeof(float) * a_row * b_column);
    float *M7 = (float *) malloc(sizeof(float) * a_row * b_column);
    float *A11 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A12 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A22 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A21 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A11022 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A21022 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A11012 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A21111 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *A12122 = (float *) malloc(sizeof(float) * a_row * a_column);
    float *B22 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B11 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B21 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B12 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B11022 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B12122 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B21111 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B11012 = (float *) malloc(sizeof(float) * a_column * b_column);
    float *B21022 = (float *) malloc(sizeof(float) * a_column * b_column);

    // clock_t start, stop;
    // start = clock();
    gettimeofday(&t_strassen, NULL);

    for (int i = 0; i < a_row; ++i) {
        for (int j = 0; j < a_column; ++j) {
            *(A11 + i * a_column + j) = *(a + i * A_c + j);
            *(A12 + i * a_column + j) = *(a + i * A_c + j + a_column);
            *(A22 + i * a_column + j) = *(a + (i + a_row) * A_c + j + a_column);
            *(A21 + i * a_column + j) = *(a + (i + a_row) * A_c + j);
            *(A11022 + i * a_column + j) = *(A11 + i * a_column + j) + *(A22 + i * a_column + j);
            *(A21022 + i * a_column + j) = *(A21 + i * a_column + j) + *(A22 + i * a_column + j);
            *(A11012 + i * a_column + j) = *(A11 + i * a_column + j) + *(A12 + i * a_column + j);
            *(A21111 + i * a_column + j) = *(A21 + i * a_column + j) - *(A11 + i * a_column + j);
            *(A12122 + i * a_column + j) = *(A12 + i * a_column + j) - *(A22 + i * a_column + j);
        }
    }
    for (int i = 0; i < a_column; ++i) {
        for (int j = 0; j < b_column; ++j) {
            *(B11 + i * b_column + j) = *(b + i * B_c + j);
            *(B21 + i * b_column + j) = *(b + (i + a_column) * B_c + j);
            *(B22 + i * b_column + j) = *(b + (i + a_column) * B_c + j + b_column);
            *(B12 + i * b_column + j) = *(b + i * B_c + j + b_column);
            *(B11022 + i * b_column + j) = *(B11 + i * b_column + j) + *(B22 + i * b_column + j);
            *(B12122 + i * b_column + j) = *(B12 + i * b_column + j) - *(B22 + i * b_column + j);
            *(B21111 + i * b_column + j) = *(B21 + i * b_column + j) - *(B11 + i * b_column + j);
            *(B11012 + i * b_column + j) = *(B11 + i * b_column + j) + *(B12 + i * b_column + j);
            *(B21022 + i * b_column + j) = *(B21 + i * b_column + j) + *(B22 + i * b_column + j);
        }
    }
    for (int i = 0; i < a_row; ++i) {
        for (int j = 0; j < b_column; ++j) {
            *(M1 + i * b_column + j) = 0.0f;
            *(M2 + i * b_column + j) = 0.0f;
            *(M3 + i * b_column + j) = 0.0f;
            *(M4 + i * b_column + j) = 0.0f;
            *(M5 + i * b_column + j) = 0.0f;
            *(M6 + i * b_column + j) = 0.0f;
            *(M7 + i * b_column + j) = 0.0f;
        }
    }
    strassen(A11022, B11022, M1, a_row, a_column, b_column, out_decided);
    strassen(A21022, B11, M2, a_row, a_column, b_column, out_decided);
    strassen(A11, B12122, M3, a_row, a_column, b_column, out_decided);
    strassen(A22, B21111, M4, a_row, a_column, b_column, out_decided);
    strassen(A11012, B22, M5, a_row, a_column, b_column, out_decided);
    strassen(A21111, B11012, M6, a_row, a_column, b_column, out_decided);
    strassen(A12122, B21022, M7, a_row, a_column, b_column, out_decided);
    for (int i = 0; i < a_row; i++) {
        for (int j = 0; j < b_column; j++) {
            *(c + i * B_c + j) =
                    *(M1 + i * b_column + j) + *(M4 + i * b_column + j) - *(M5 + i * b_column + j) +
                    *(M7 + i * b_column + j);
            *(c + i * B_c + j + b_column) = *(M3 + i * b_column + j) + *(M5 + i * b_column + j);
            *(c + (i + a_row) * B_c + j + b_column) =
                    *(M1 + i * b_column + j) + *(M3 + i * b_column + j) - *(M2 + i * b_column + j) +
                    *(M6 + i * b_column + j);
            *(c + (i + a_row) * B_c + j) = *(M2 + i * b_column + j) + *(M4 + i * b_column + j);
        }
    }
    gettimeofday(&t_strassen_f, NULL);
    if (out_decided == A_r)
        printf("The calculation of multiplication using strassen costs %.3f ms (%d)\n",
               (t_strassen_f.tv_sec - t_strassen.tv_sec) * (double) 1000 +
               (t_strassen_f.tv_usec - t_strassen.tv_usec) / (double) 1000, out_decided);

    free(M1);
    free(M2);
    free(M3);
    free(M4);
    free(M5);
    free(M6);
    free(M7);
    free(A11);
    free(A12);
    free(A22);
    free(A11012);
    free(A11022);
    free(A12122);
    free(A21022);
    free(A21111);
    free(B11);
    free(B12);
    free(B21);
    free(B22);
    free(B11012);
    free(B11022);
    free(B12122);
    free(B21022);
    free(B21111);
}

void exchangeRow(float *a, int i, int j, int col, int c_mark) {
    int k;
    float t;
#pragma omp parallel for
    for (k = 0; k < col - c_mark; ++k) {
        t = *(a + i * col + c_mark + k);
        *(a + i * col + c_mark + k) = *(a + j * col + c_mark + k);
        *(a + j * col + c_mark + k) = t;
    }
}

void mulRow(float *a, int r, float k, int col, int c_mark) {
    for (int i = 0; i < col - c_mark; i++)
        *(a + r * col + c_mark + i) *= k;
}

void unitization(float *a, int i, int col, int c) {
    float fir = *(a + i * col + c);
    for (int j = c; j < col; ++j) {
        *(a + i * col + j) /= fir;
    }
}

void addRow(float *a, int r1, int r2, float k, int col, int c_mark) {
    for (int i = 0; i < col - c_mark; i++)
        *(a + col * r1 + c_mark + i) += *(a + r2 * col + c_mark + i) * k;
}

p_Matrix rank(const float *mat, int m, int n) {
    int r_mark = 0;
    int c_mark = 0;  //行标记与列标记
    int sign_c;  //某列是否全为0的标志，为1表示全为0
    p_Matrix ladder = (Matrix *) malloc(sizeof(Matrix));
    ladder->m_row = m;
    ladder->m_column = n;
    ladder->m_matrix = (float *) malloc(sizeof(float) * m * n);
    float *c = ladder->m_matrix;
    for (int i = 0; i < m * n; ++i) {
        *(c + i) = *(mat + i);
    }
    gettimeofday(&tv, NULL);
    while (c_mark < n) {
        sign_c = 1;
        for (int i = r_mark; i < m; ++i) {
            if (*(c + i * n + c_mark) != 0) {
                if (i != r_mark) {
                    if (sign_c)
                        exchangeRow(c, r_mark, i, n, c_mark);
                    else {
                        float t = *(c + i * n + c_mark);
                        mulRow(c, i, *(c + r_mark * n + c_mark), n, c_mark);
                        addRow(c, i, r_mark, -t, n, c_mark);
                    }
                } else
                    unitization(c, i, n, c_mark);
                sign_c = 0;
            }
        }
        if (!sign_c) r_mark++;
        c_mark++;
    }
    gettimeofday(&f_tv, NULL);
    printf("The calculation of rank costs %.3f ms \n",
           (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000);
    ladder->m_rank = r_mark;
    return ladder;
}

void LU(const float *A, float *Low, float *Up, int side) {
    float *r1 = (float *) malloc(sizeof(float) * side);
    float *r2 = (float *) malloc(sizeof(float) * side);
    for (int i = 0; i < side; ++i) {
        *(r1 + i) = 0.f;
        *(r2 + i) = 0.f;
    }
    float temp = 0.f;
    for (int i = 0; i < side; i++) {
        for (int j = 0; j < side; j++) {
            if (i == j)
                *(Low + i * side + i) = 1.f;
        }
    }
    gettimeofday(&tv ,NULL);
    for (int k = 0; k < side; k++) {
        for (int j = k; j < side; j++) {
            if (k == 0) {
                temp = 0.f;
            } else {
                for (int i = 0; i < k; i++) {
                    *(r1 + i) = *(Low + k * side + i);
                }
                for (int i = 0; i < k; i++) {
                    *(r2 + i) = *(Up + i * side + j);
                }
                temp = 0.f;
                for (int i = 0; i < k; i++) {
                    temp += *(r1 + i) * *(r2 + i);
                }
            }
            *(Up + k * side + j) = *(A + k * side + j) - temp;
        }

        for (int i = 0; i < side; ++i) {
            *(r1 + i) = 0.f;
            *(r2 + i) = 0.f;
        }
        temp = 0.f;
        for (int t = k + 1; t < side; t++) {
            if (k == 0) {
                temp = 0.f;
            } else {
                for (int i = 0; i < k; i++) {
                    *(r1 + i) = *(Low + t * side + i);
                }
                for (int i = 0; i < k; i++) {
                    *(r2 + i) = *(Up + i * side + k);
                }
                temp = 0.f;
                for (int i = 0; i < k; i++) {
                    temp += *(r1 + i) * *(r2 + i);
                }
            }
            *(Low + t * side + k) = (*(A + t * side + k) - temp) / *(Up + k * side + k);
        }
    }
    gettimeofday(&f_tv , NULL);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of LU factorization costs %.3f ms\n", period);
}

void mulSIMD(Matrix * c, const Matrix  * a, const Matrix * b ) {
    if (a->m_column != b->m_row){
        printf("MULTIPLY FAILURE!");
        return;
    }
    gettimeofday(&tv, NULL);
    __m256 num0, s;
    int i, j, k, n;
#pragma omp parallel for
    for (i = 0; i < a->m_row; ++i) {
        for (k = 0; k < a->m_column ; ++k){
            s = _mm256_broadcast_ss(a->m_matrix + i * a->m_column + k);
            for (j = 0; j < b->m_column; j += 8) {
                n = i * c->m_column;
                num0 = _mm256_loadu_ps(c->m_matrix + n + j);
                num0 = _mm256_add_ps(num0, _mm256_mul_ps(s,_mm256_loadu_ps(b->m_matrix + k * b->m_column + j)));
                _mm256_storeu_ps(c->m_matrix + i * c->m_column + j, num0);
            }
        }
    }
    gettimeofday(&f_tv, NULL);
    double period = (f_tv.tv_sec - tv.tv_sec) * (double) 1000 + (f_tv.tv_usec - tv.tv_usec) / (double) 1000;
    printf("The calculation of multiplication using SIMD costs %.3f ms\n", period);
}