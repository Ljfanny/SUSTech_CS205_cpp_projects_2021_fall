#pragma once
#pragma GCC optimize(3)

typedef enum boolean {
    True, False
} Bool;

typedef struct {
    int m_rank;
    int m_row;
    int m_column;
    float *m_matrix;
    Bool hasInversion;
    Bool hasLU;
} Matrix, *p_Matrix;

// struct timeval  
// {   
//     long tv_sec; /*秒*/   
//     long tv_usec; /*微秒*/   
// };      

p_Matrix createMatrix(char *fileName);

void deleteMatrix(p_Matrix matrix);

p_Matrix copyMatrix(p_Matrix matrix);

p_Matrix multiply(const float *A, const float *B, int A_row, int A_column, int B_row, int B_column);

p_Matrix plus(const float *A, const float *B, int A_row, int A_column, int B_row, int B_column);

p_Matrix subtraction(const float *A, const float *B, int A_row, int A_column, int B_row, int B_column);

p_Matrix transpose(const float *A, int A_row, int A_column);

void swap(p_Matrix matrix, int i, int j);

void det(p_Matrix matrix);

void fileWriter(FILE *fp, p_Matrix matrix);

void scaleRow(p_Matrix matrix, int i, float mul);

p_Matrix inverse(p_Matrix matrix);

void IntMultiply(char *a, const char *b, int a_len, int b_len);

void strassenMul(const float *a, const float *b, float *c, int A_row, int A_column, int B_column);

void strassen(const float *a, const float *b, float *c, int a_row, int a_column, int b_column,int out_decided);

void exchangeRow(float *a, int i, int j, int col, int c_mark);

void mulRow(float *a, int r, float k, int col, int c_mark);

void addRow(float *a, int r1, int r2, float k, int col, int c_mark);

p_Matrix rank(const float *mat, int m, int n);

void LU(const float *A, float *Low, float *Up, int side);

void mulSIMD(Matrix * c, const Matrix  * a, const Matrix * b);

void multiply_func(const float *A, const float *B, float *C, int A_row, int A_column, int B_row, int B_column);