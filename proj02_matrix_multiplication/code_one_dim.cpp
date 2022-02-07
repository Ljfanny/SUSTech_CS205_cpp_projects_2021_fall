#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <sstream>

int row1 = 0;
int row2 = 0;
int column1 = 0;
int column2 = 0;

using namespace std;

int getColumn(string file);

int getRow(string file) ;

void simpleMatrixMulF(float *arr1, float *arr2, float *res);

void simpleMatrixMulD(double *arr1, double *arr2, double *res) ;

void smallMatrixMulF(float *arr1, float *arr2, float *res, int ro1, int col1, int col2) ;

void smallMatrixMulD(double *arr1, double *arr2, double *res, int ro1, int col1, int col2) ;

static void strassenF(float *arr1, float *arr2, float *res, int ro1, int col1, int col2);

static void strassenD(double *arr1, double *arr2, double *res, int ro1, int col1, int col2) ;

static void IKJF(float *arr1, float *arr2, float *res);

static void IKJD(double *arr1, double *arr2, double *res);

static void fileOutPutF(const string &out, float *res) ;

static void fileOutPutD(const string &out, double *res);

#pragma GCC optimize(3)

int main(int argc, char **argv) {
    string inp01, inp02, outfileF_IKJ, outfileD_IKJ, outfileF_N3, outfileD_N3, outfileF_N281, outfileD_N281;
    if (argc > 7) {
        inp01 = argv[1];
        inp02 = argv[2];
        outfileF_IKJ = argv[3];
        outfileD_IKJ = argv[4];
        outfileF_N3 = argv[5];
        outfileD_N3 = argv[6];
        outfileF_N281 = argv[7];
        outfileD_N281 = argv[8];
    } else {
        inp01 = argv[1];
        inp02 = argv[2];
        outfileF_IKJ = argv[3];
        outfileD_IKJ = argv[4];
        outfileF_N3 = argv[5];
        outfileD_N3 = argv[6];
    }
    row1 = getRow(inp01);
    column1 = getColumn(inp01);
    row2 = column1;
    column2 = getColumn(inp02);


    double *aD = new double[row1 * column1];
    float *aF = new float[row1 * column1];
    double *bD = new double[row2 * column2];
    float *bF = new float[row2 * column2];

    // 对第一个矩阵的处理
    ifstream infile1, infile2;
    infile1.open(inp01);
    string temp;
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < column1; ++j) {
            infile1 >> temp;
            aF[i * column1 + j] = stof(temp);
            aD[i * column1 + j] = stod(temp);
        }
    }
    // 对第二个矩阵的处理
    infile2.open(inp02);
    for (int i = 0; i < row2; ++i) {
        for (int j = 0; j < column2; ++j) {
            infile2 >> temp;
            bF[i * column1 + j] = stof(temp);
            bD[i * column1 + j] = stod(temp);
        }
    }
    float *resultF_IKJ = new float[row1 * column2];
    double *resultD_IKJ = new double[row1 * column2];
    float *resultF_N3 = new float[row1 * column2];
    double *resultD_N3 = new double[row1 * column2];
    float *resultF_N281 = new float[row1 * column2];
    double *resultD_N281 = new double[row1 * column2];

    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < column2; ++j) {
            resultF_IKJ[i * column2 + j] = 0.0f;
            resultF_N281[i * column2 + j] = 0.0f;
            resultF_N3[i * column2 + j] = 0.0f;
            resultD_IKJ[i * column2 + j] = 0.0;
            resultD_N3[i * column2 + j] = 0.0;
            resultD_N281[i * column2 + j] = 0.0;
        }
    }

    auto start = chrono::steady_clock::now();
    IKJF((float *) aF, (float *) bF, (float *) resultF_IKJ);
    auto end1 = std::chrono::steady_clock::now();
    chrono::duration<double, micro> elapsed = end1 - start;
    cout << "Use float to calculate and the complexity is IKJ: " << elapsed.count() / 1000 << " ms"
         << endl;

    IKJD((double *) aD, (double *) bD, (double *) resultD_IKJ);
    auto end2 = std::chrono::steady_clock::now();
    elapsed = end2 - end1;
    cout << "Use double to calculate and the complexity is IKJ: " << elapsed.count() / 1000 << " ms"
         << endl;

    simpleMatrixMulF((float *) aF, (float *) bF, (float *) resultF_N3);
    auto end3 = std::chrono::steady_clock::now();
    elapsed = end3 - end2;
    cout << "Use float to calculate and the complexity is O(n^3): " << elapsed.count() / 1000 << " ms"
         << endl;

    simpleMatrixMulD((double *) aD, (double *) bD, (double *) resultD_N3);
    auto end4 = std::chrono::steady_clock::now();
    elapsed = end4 - end3;
    cout << "Use double to calculate and the complexity is O(n^3): " << elapsed.count() / 1000 << " ms"
         << endl;


    if (row1 < 1024) {
        strassenF((float *) aF, (float *) bF, (float *) resultF_N281, row1, column1, column2);
        auto end5 = std::chrono::steady_clock::now();
        elapsed = end5 - end4;
        cout << "Use float to calculate and the complexity is O(n^2.81): " << elapsed.count() / 1000 << " ms"
             << endl;

        strassenD((double *) aD, (double *) bD, (double *) resultD_N281, row1, column1, column2);
        auto end6 = std::chrono::steady_clock::now();
        elapsed = end6 - end5;
        cout << "Use double to calculate and the complexity is O(n^2.81): " << elapsed.count() / 1000 << " ms"
             << endl;
    }

    if (row1 < 1024) {
        fileOutPutF(outfileF_IKJ, resultF_IKJ);
        fileOutPutF(outfileF_N3, resultF_N3);
        fileOutPutF(outfileF_N281, resultF_N281);
        fileOutPutD(outfileD_IKJ, resultD_IKJ);
        fileOutPutD(outfileD_N3, resultD_N3);
        fileOutPutD(outfileD_N281, resultD_N281);
    } else {
        fileOutPutF(outfileF_IKJ, resultF_IKJ);
        fileOutPutF(outfileF_N3, resultF_N3);
        fileOutPutD(outfileD_IKJ, resultD_IKJ);
        fileOutPutD(outfileD_N3, resultD_N3);
    }

    delete[] aD;
    delete[] aF;
    delete[] bD;
    delete[] bF;
    delete[] resultD_IKJ;
    delete[] resultF_IKJ;
    delete[] resultD_N3;
    delete[] resultF_N3;
    delete[] resultD_N281;
    delete[] resultF_N281;

    return 0;
}

int getColumn(string file) {
    ifstream goalFile(file);
    string unit, line;
    getline(goalFile, line);
    stringstream input(line);
    int edgeNum = 0;
    while (input >> unit) {
        edgeNum++;
    }
    return edgeNum;
}

int getRow(string file) {
    ifstream goalFile(file);
    string line;
    int nums = 0;
    while (getline(goalFile, line)) {
        nums++;
    }
    return nums;
}

void simpleMatrixMulF(float *arr1, float *arr2, float *res) {
    float c;
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < column2; ++j) {
            c = 0.0f;
            for (int k = 0; k < column1; ++k) {
                c += (float) arr1[i * column1 + k] * (float) arr2[k * column2 + j];
            }
            res[i * column2 + j] = c;
        }
    }
}

void simpleMatrixMulD(double *arr1, double *arr2, double *res) {
    double c;
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < column2; ++j) {
            c = 0.0;
            for (int k = 0; k < column1; ++k) {
                c += arr1[i * column1 + k] * arr2[k * column2 + j];
            }
            res[i * column2 + j] = c;
        }
    }
}

void smallMatrixMulF(float *arr1, float *arr2, float *res, int ro1, int col1, int col2) {
    float c;
    for (int i = 0; i < ro1; ++i) {
        for (int j = 0; j < col2; ++j) {
            c = 0.0f;
            for (int k = 0; k < col1; ++k) {
                c += arr1[i * col1 + k] * arr2[k * col2 + j];
            }
            res[i * col2 + j] = c;
        }
    }
}

void smallMatrixMulD(double *arr1, double *arr2, double *res, int ro1, int col1, int col2) {
    double c;
    for (int i = 0; i < ro1; ++i) {
        for (int j = 0; j < col2; ++j) {
            c = 0.0;
            for (int k = 0; k < col1; ++k) {
                c += arr1[i * col1 + k] * arr2[k * col2 + j];
            }
            res[i * col2 + j] = c;
        }
    }
}

static void strassenF(float *arr1, float *arr2, float *res, int ro1, int col1, int col2) {
    if (col1 % 2 != 0 || col2 % 2 != 0 || ro1 % 2 != 0)
        return smallMatrixMulF(arr1, arr2, res, ro1, col1, col2);

    float *M1 = new float[ro1 / 2 * col2 / 2];
    float *M2 = new float[ro1 / 2 * col2 / 2];
    float *M3 = new float[ro1 / 2 * col2 / 2];
    float *M4 = new float[ro1 / 2 * col2 / 2];
    float *M5 = new float[ro1 / 2 * col2 / 2];
    float *M6 = new float[ro1 / 2 * col2 / 2];
    float *M7 = new float[ro1 / 2 * col2 / 2];
    float *A11 = new float[ro1 / 2 * col1 / 2];
    float *A12 = new float[ro1 / 2 * col1 / 2];
    float *A22 = new float[ro1 / 2 * col1 / 2];
    float *A21 = new float[ro1 / 2 * col1 / 2];
    float *A11022 = new float[ro1 / 2 * col1 / 2];
    float *A21022 = new float[ro1 / 2 * col1 / 2];
    float *A11012 = new float[ro1 / 2 * col1 / 2];
    float *A21111 = new float[ro1 / 2 * col1 / 2];
    float *A12122 = new float[ro1 / 2 * col1 / 2];
    float *B22 = new float[col1 / 2 * col2 / 2];
    float *B11 = new float[col1 / 2 * col2 / 2];
    float *B21 = new float[col1 / 2 * col2 / 2];
    float *B12 = new float[col1 / 2 * col2 / 2];
    float *B11022 = new float[col1 / 2 * col2 / 2];
    float *B12122 = new float[col1 / 2 * col2 / 2];
    float *B21111 = new float[col1 / 2 * col2 / 2];
    float *B11012 = new float[col1 / 2 * col2 / 2];
    float *B21022 = new float[col1 / 2 * col2 / 2];

    for (int i = 0; i < ro1 / 2; ++i) {
        for (int j = 0; j < col1 / 2; ++j) {
            A11[i * col1 / 2 + j] = arr1[i * col1 + j];
            A12[i * col1 / 2 + j] = arr1[i * col1 + j + col1 / 2];
            A22[i * col1 / 2 + j] = arr1[(i + ro1 / 2) * col1 + j + col1 / 2];
            A21[i * col1 / 2 + j] = arr1[(i + ro1 / 2) * col1 + j];
            A11022[i * col1 / 2 + j] = A11[i * col1 / 2 + j] + A22[i * col1 / 2 + j];
            A21022[i * col1 / 2 + j] = A21[i * col1 / 2 + j] + A22[i * col1 / 2 + j];
            A11012[i * col1 / 2 + j] = A11[i * col1 / 2 + j] + A12[i * col1 / 2 + j];
            A21111[i * col1 / 2 + j] = A21[i * col1 / 2 + j] - A11[i * col1 / 2 + j];
            A12122[i * col1 / 2 + j] = A12[i * col1 / 2 + j] - A22[i * col1 / 2 + j];
        }
    }
    for (int i = 0; i < col1 / 2; ++i) {
        for (int j = 0; j < col2 / 2; ++j) {
            B11[i * col2 / 2 + j] = arr2[i * col2 + j];
            B21[i * col2 / 2 + j] = arr2[(i + col1 / 2) * col2 + j];
            B22[i * col2 / 2 + j] = arr2[(i + col1 / 2) * col2 + j + col2 / 2];
            B12[i * col2 / 2 + j] = arr2[i * col2 + j + col2 / 2];
            B11022[i * col2 / 2 + j] = B11[i * col2 / 2 + j] + B22[i * col2 / 2 + j];
            B12122[i * col2 / 2 + j] = B12[i * col2 / 2 + j] - B22[i * col2 / 2 + j];
            B21111[i * col2 / 2 + j] = B21[i * col2 / 2 + j] - B11[i * col2 / 2 + j];
            B11012[i * col2 / 2 + j] = B11[i * col2 / 2 + j] + B12[i * col2 / 2 + j];
            B21022[i * col2 / 2 + j] = B21[i * col2 / 2 + j] + B22[i * col2 / 2 + j];
        }
    }

    for (int i = 0; i < ro1 / 2; ++i) {
        for (int j = 0; j < col2 / 2; ++j) {
            M1[i * col2 / 2 + j] = 0;
            M2[i * col2 / 2 + j] = 0;
            M3[i * col2 / 2 + j] = 0;
            M4[i * col2 / 2 + j] = 0;
            M5[i * col2 / 2 + j] = 0;
            M6[i * col2 / 2 + j] = 0;
            M7[i * col2 / 2 + j] = 0;
        }
    }
    strassenF((float *) A11022, (float *) B11022, (float *) M1, ro1 / 2, col1 / 2, col2 / 2);
    strassenF((float *) A21022, (float *) B11, (float *) M2, ro1 / 2, col1 / 2, col2 / 2);
    strassenF((float *) A11, (float *) B12122, (float *) M3, ro1 / 2, col1 / 2, col2 / 2);
    strassenF((float *) A22, (float *) B21111, (float *) M4, ro1 / 2, col1 / 2, col2 / 2);
    strassenF((float *) A11012, (float *) B22, (float *) M5, ro1 / 2, col1 / 2, col2 / 2);
    strassenF((float *) A21111, (float *) B11012, (float *) M6, ro1 / 2, col1 / 2, col2 / 2);
    strassenF((float *) A12122, (float *) B21022, (float *) M7, ro1 / 2, col1 / 2, col2 / 2);

    for (int i = 0; i < ro1 / 2; i++) {
        for (int j = 0; j < col2 / 2; j++) {
            res[i * col2 + j] =
                    M1[i * col2 / 2 + j] + M4[i * col2 / 2 + j] - M5[i * col2 / 2 + j] + M7[i * col2 / 2 + j];
            res[i * col2 + j + col2 / 2] = M3[i * col2 / 2 + j] + M5[i * col2 / 2 + j];
            res[(i + ro1 / 2) * col2 + j + col2 / 2] =
                    M1[i * col2 / 2 + j] + M3[i * col2 / 2 + j] - M2[i * col2 / 2 + j] + M6[i * col2 / 2 + j];
            res[(i + ro1 / 2) * col2 + j] = M2[i * col2 / 2 + j] + M4[i * col2 / 2 + j];
        }
    }
}

static void strassenD(double *arr1, double *arr2, double *res, int ro1, int col1, int col2) {
    if (col1 % 2 != 0 || col2 % 2 != 0 || ro1 % 2 != 0)
        return smallMatrixMulD(arr1, arr2, res, ro1, col1, col2);
    double *M1 = new double[ro1 / 2 * col2 / 2];
    double *M2 = new double[ro1 / 2 * col2 / 2];
    double *M3 = new double[ro1 / 2 * col2 / 2];
    double *M4 = new double[ro1 / 2 * col2 / 2];
    double *M5 = new double[ro1 / 2 * col2 / 2];
    double *M6 = new double[ro1 / 2 * col2 / 2];
    double *M7 = new double[ro1 / 2 * col2 / 2];
    double *A11 = new double[ro1 / 2 * col1 / 2];
    double *A12 = new double[ro1 / 2 * col1 / 2];
    double *A22 = new double[ro1 / 2 * col1 / 2];
    double *A21 = new double[ro1 / 2 * col1 / 2];
    double *A11022 = new double[ro1 / 2 * col1 / 2];
    double *A21022 = new double[ro1 / 2 * col1 / 2];
    double *A11012 = new double[ro1 / 2 * col1 / 2];
    double *A21111 = new double[ro1 / 2 * col1 / 2];
    double *A12122 = new double[ro1 / 2 * col1 / 2];
    double *B22 = new double[col1 / 2 * col2 / 2];
    double *B11 = new double[col1 / 2 * col2 / 2];
    double *B21 = new double[col1 / 2 * col2 / 2];
    double *B12 = new double[col1 / 2 * col2 / 2];
    double *B11022 = new double[col1 / 2 * col2 / 2];
    double *B12122 = new double[col1 / 2 * col2 / 2];
    double *B21111 = new double[col1 / 2 * col2 / 2];
    double *B11012 = new double[col1 / 2 * col2 / 2];
    double *B21022 = new double[col1 / 2 * col2 / 2];

    for (int i = 0; i < ro1 / 2; ++i) {
        for (int j = 0; j < col1 / 2; ++j) {
            A11[i * col1 / 2 + j] = arr1[i * col1 + j];
            A12[i * col1 / 2 + j] = arr1[i * col1  + j + col1 / 2];
            A22[i * col1 / 2 + j] = arr1[(i + ro1 / 2) * col1  + j + col1 / 2];
            A21[i * col1 / 2 + j] = arr1[(i + ro1 / 2) * col1  + j];
            A11022[i * col1 / 2 + j] = A11[i * col1 / 2 + j] + A22[i * col1 / 2 + j];
            A21022[i * col1 / 2 + j] = A21[i * col1 / 2 + j] + A22[i * col1 / 2 + j];
            A11012[i * col1 / 2 + j] = A11[i * col1 / 2 + j] + A12[i * col1 / 2 + j];
            A21111[i * col1 / 2 + j] = A21[i * col1 / 2 + j] - A11[i * col1 / 2 + j];
            A12122[i * col1 / 2 + j] = A12[i * col1 / 2 + j] - A22[i * col1 / 2 + j];
        }
    }
    for (int i = 0; i < col1 / 2; ++i) {
        for (int j = 0; j < col2 / 2; ++j) {
            B11[i * col2 / 2 + j] = arr2[i * col2 + j];
            B21[i * col2 / 2 + j] = arr2[(i + col1 / 2) * col2  + j];
            B22[i * col2 / 2 + j] = arr2[(i + col1 / 2) * col2  + j + col2 / 2];
            B12[i * col2 / 2 + j] = arr2[i * col2 + j + col2 / 2];
            B11022[i * col2 / 2 + j] = B11[i * col2 / 2 + j] + B22[i * col2 / 2 + j];
            B12122[i * col2 / 2 + j] = B12[i * col2 / 2 + j] - B22[i * col2 / 2 + j];
            B21111[i * col2 / 2 + j] = B21[i * col2 / 2 + j] - B11[i * col2 / 2 + j];
            B11012[i * col2 / 2 + j] = B11[i * col2 / 2 + j] + B12[i * col2 / 2 + j];
            B21022[i * col2 / 2 + j] = B21[i * col2 / 2 + j] + B22[i * col2 / 2 + j];
        }
    }

    for (int i = 0; i < ro1 / 2; ++i) {
        for (int j = 0; j < col2 / 2; ++j) {
            M1[i * col2 / 2 + j] = 0;
            M2[i * col2 / 2 + j] = 0;
            M3[i * col2 / 2 + j] = 0;
            M4[i * col2 / 2 + j] = 0;
            M5[i * col2 / 2 + j] = 0;
            M6[i * col2 / 2 + j] = 0;
            M7[i * col2 / 2 + j] = 0;
        }
    }

    strassenD((double *) A11022, (double *) B11022, (double *) M1, ro1 / 2, col1 / 2, col2 / 2);
    strassenD((double *) A21022, (double *) B11, (double *) M2, ro1 / 2, col1 / 2, col2 / 2);
    strassenD((double *) A11, (double *) B12122, (double *) M3, ro1 / 2, col1 / 2, col2 / 2);
    strassenD((double *) A22, (double *) B21111, (double *) M4, ro1 / 2, col1 / 2, col2 / 2);
    strassenD((double *) A11012, (double *) B22, (double *) M5, ro1 / 2, col1 / 2, col2 / 2);
    strassenD((double *) A21111, (double *) B11012, (double *) M6, ro1 / 2, col1 / 2, col2 / 2);
    strassenD((double *) A12122, (double *) B21022, (double *) M7, ro1 / 2, col1 / 2, col2 / 2);

    for (int i = 0; i < ro1 / 2; i++) {
        for (int j = 0; j < col2 / 2; j++) {
            res[i * col2 + j] =
                    M1[i * col2 / 2 + j] + M4[i * col2 / 2 + j] - M5[i * col2 / 2 + j] + M7[i * col2 / 2 + j];
            res[i * col2 + j + col2 / 2] = M3[i * col2 / 2 + j] + M5[i * col2 / 2 + j];
            res[(i + ro1 / 2) * col2 + j + col2 / 2] =
                    M1[i * col2 / 2 + j] + M3[i * col2 / 2 + j] - M2[i * col2 / 2 + j] + M6[i * col2 / 2 + j];
            res[(i + ro1 / 2) * col2 + j] = M2[i * col2 / 2 + j] + M4[i * col2 / 2 + j];
        }
    }

}


static void IKJF(float *arr1, float *arr2, float *res) {
    float d = 0.0F;
    {
#pragma omp parallel for
        for (int i = 0; i < row1; i++) {
            for (int k = 0; k < column1; k++) {
                d = arr1[i * column1 + k];
                for (int j = 0; j < column2; j++) {
                    res[i * column2 + j] += d * arr2[k * column2 + j];
                }
            }
        }
    }
}

static void IKJD(double *arr1, double *arr2, double *res) {
    double d = 0.0;
    {
#pragma omp parallel for
        for (int i = 0; i < row1; i++) {
            for (int k = 0; k < column1; k++) {
                d = arr1[i * column1 + k];
                for (int j = 0; j < column2; j++) {
                    res[i * column2 + j] += d * arr2[k * column2 + j];
                }
            }
        }
    }
}

static void fileOutPutF(const string &out, float *res) {
    ofstream outfile;
    outfile.open(out);
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < column2; ++j) {
            outfile << to_string(res[i * column2 + j]) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

static void fileOutPutD(const string &out, double *res) {
    ofstream outfile;
    outfile.open(out);
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < column2; ++j) {
            outfile << to_string(res[i * column2 + j]) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}