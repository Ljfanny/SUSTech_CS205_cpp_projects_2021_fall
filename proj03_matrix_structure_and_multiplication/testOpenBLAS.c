#include <stdio.h>
#include <string.h>
#include <cblas.h>
#include <time.h>
#include <stdlib.h>

typedef struct {
//    int m_rank;
    int m_row;
    int m_column;
    float *m_matrix;
//    Bool hasInversion;
} Matrix, *p_Matrix;

void fileWriter(FILE *fp, p_Matrix matrix) {
    for (int i = 0; i < matrix->m_row; ++i) {
        for (int j = 0; j < matrix->m_column; ++j) {
            fprintf(fp, "%.4f ", matrix->m_matrix[i * matrix->m_column + j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

//using namespace std;
p_Matrix createMatrix(char *fileName) {
    p_Matrix matrix = (Matrix *) malloc(sizeof(Matrix));
    matrix->m_row = 0;
    matrix->m_column = 0;
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

int main(int argc, char **argv) {
    struct timeval start, stop, up, down;
    gettimeofday(&start,NULL);
    char *fileA = NULL;
    char *fileB = NULL;
    char *fileC = NULL;
    if (argc > 1) {
        fileA = argv[1];
        fileB = argv[2];
        fileC = argv[3];
    }
    FILE *fp = fopen(fileC, "w");
    p_Matrix A = createMatrix(fileA);
    p_Matrix B = createMatrix(fileB);

    p_Matrix C = (Matrix *) malloc(sizeof(Matrix));
    C->m_row = A->m_row;
    C->m_column = B->m_column;
    C->m_matrix = (float *)malloc(sizeof(float) * C->m_column * C->m_row);
    
    const float *a;
    const float *b;
    float *c;
    a = A->m_matrix;
    b = B->m_matrix;
    c = C->m_matrix;
    const int M = A->m_row;
    const int N = A->m_column;
    const int K = B->m_column;
    const float alpha = 1;
    const float beta_positive = 1;
    const float beta_negative = -1.0f;
    const float beta = 0;
    const int lda = M;
    const int ldb = K;
    const int ldc = N;

    float * e = (float *)malloc(sizeof(float) * A->m_row * A->m_column);
    for(int i = 0; i < A->m_row; ++i){
        *(e + i * A->m_column + i) = 1.0f;
    }
    const float * E = e;
    double period;

    // gettimeofday(&up, NULL);
    // cblas_cgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, &alpha,(float *) a, lda,(float *) b, ldb, &beta,(float *) c, ldc);
    // gettimeofday(&down,NULL);
    // double period = (down.tv_sec - up.tv_sec)*(double)1000 + (down.tv_usec - up.tv_usec)/(double)1000;
    // printf("OpenBLAS calculates the multiplication (less matrix mul) costs %.3f ms\n",period);
    // fileWriter(fp, C);

    gettimeofday(&up, NULL);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha,(float *) a, lda,(float *) b, ldb, beta,(float *) c, ldc);
    gettimeofday(&down,NULL);
    period = (down.tv_sec - up.tv_sec)*(double)1000 + (down.tv_usec - up.tv_usec)/(double)1000;
    printf("OpenBLAS calculates the multiplication costs %.3f ms\n",period);
    fileWriter(fp, C);

    gettimeofday(&up, NULL);
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, alpha,(float *) a, lda,(float *) E, ldb, beta,(float *) c, ldc);
    gettimeofday(&down,NULL);
    period = (down.tv_sec - up.tv_sec)*(double)1000 + (down.tv_usec - up.tv_usec)/(double)1000;
    printf("OpenBLAS calculates the transposition costs %.3f ms\n",period);
    fileWriter(fp, C);

    for(int i = 0; i < B->m_row; ++i){
        for(int j = 0; j < B->m_column; ++j)
        *(c + i * B->m_column + j) = *(b + i * B->m_column + j);
    }
    // const float * E = e;

    gettimeofday(&up, NULL);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha,(float *) a, lda,(float *) E, ldb, beta_positive,(float *) c, ldc);
    gettimeofday(&down,NULL);
    period = (down.tv_sec - up.tv_sec)*(double)1000 + (down.tv_usec - up.tv_usec)/(double)1000;
    printf("OpenBLAS calculates the addition costs %.3f ms\n",period);
    fileWriter(fp, C);

    for(int i = 0; i < B->m_row; ++i){
        for(int j = 0; j < B->m_column; ++j)
        *(c + i * B->m_column + j) = *(b + i * B->m_column + j);
    }
     
    // const float alph = -1;
    gettimeofday(&up, NULL);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha,(float *) a, lda,(float *) E, ldb, beta_negative,(float *) c, ldc);
    gettimeofday(&down,NULL);
    period = (down.tv_sec - up.tv_sec)*(double)1000 + (down.tv_usec - up.tv_usec)/(double)1000;
    printf("OpenBLAS calculates the subtraction costs %.3f ms\n",period);
    fileWriter(fp, C);
    
    deleteMatrix(C);
    deleteMatrix(A);
    deleteMatrix(B);
    fclose(fp);
 
    gettimeofday(&stop, NULL);
    period = (stop.tv_sec - start.tv_sec)*(double)1000 + (stop.tv_usec - start.tv_usec)/(double)1000;
    printf("The running time is %.3f ms\n",period);
    return 0 ;
}
