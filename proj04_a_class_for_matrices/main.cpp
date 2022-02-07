#include "matclass.hpp"
#include <iostream>
#include <chrono>
#pragma GCC optimize(3)

using namespace std;

//默认使用者知道矩阵的规模！
int main(int argc, char **argv) {
    auto start = chrono::steady_clock::now();
    char *fileA = nullptr;
    char *fileB = nullptr;
    char *fileC = nullptr;
    if (argc > 1) {
        fileA = argv[1];
        fileB = argv[2];
        fileC = argv[3];
    }

    ofstream output;
    output.open(fileC);
    Matrix<int> a = Matrix<int>(6, 6);
    a.readFromFile(fileA, 1);
    Matrix<int> c = a.getSubMat(0, 0, 4, 3);
    cout << "c (the submatrix of a) : " << endl;
    cout << c;
    c.changeSize(3, 4);
    cout << "c (c -> change size to 2 * 2) :" << endl;
    cout << c;

//    Matrix<int> a = Matrix<int>(6, 6);
//    a.readFromFile(fileA, 1);
//    Matrix<int> d = a.getSubMat(3, 3, 3, 3);
//    cout << "d (the submatrix of a) : " << endl;
//    cout << d;
//    测试ROI的改值操作
//    d.setMatPoint(0 ,0 , -111);
//    d.setMatPoint(1 ,1 , -222);
//    d.setMatPoint(2 ,2 , -333);
//    cout << "d (change the number of diagonal element) :" << endl;
//    cout << d;

//    类型强制转换测试
//    Matrix<long> a = Matrix<long>(6, 6);
//    a.readFromFile(fileA, 1);
//    long tem = (long) a;
//    cout << tem;

//    重载测试
//    Matrix<long> a = Matrix<long>(6, 6);
//    a.readFromFile(fileA, 2);
//    Matrix<long> b = Matrix<long>(256, 256);
//    b.readFromFile(fileB, 2);
//    a *= b;
//    a.outputToFile(output);
//    cout << "a : " << endl;
//    cout << a;
//    a -= b;
//    cout << "a -= b, a : " <<endl;
//    cout << a;

//    设置channel测试
//    Matrix<long> a = Matrix<long>(6, 6);
//    a.readFromFile(fileA, 1);
//    a.setCurrentChannel(2);
//    cout << "channel 2 :" << endl;
//    cin >> a;
//    a.setCurrentChannel(1);
//    cout << "channel 1 :" << endl;
//    cin >> a;
//    a.setCurrentChannel(3);
//    cout << "channel 3 :" << endl;
//    cin >> a;
//    cout << "a : " << endl;
//    for (int i = 1; i <= a.getChannel(); ++i) {
//        a.setCurrentChannel(i);
//        cout << a;
//    }
//    cout << "a : " << endl;
//    cout << a;

//    与矩阵相关的方法测试
//    Matrix<long> a = Matrix<long>(6, 6);
//    a.readFromFile(fileA, 1);
//    cout << "a.getTrace() = " << a.getTrace() << endl;
//    a.changeSize(4,9);
//    cout << "a.changeSize(4, 9) : " << endl;
//    cout << a;
//    cout << "a.getRowAve(0) = " << a.getRowAve(0) << endl;
//    cout << "a.getColumnAve(0) = " << a.getColumnAve(0) << endl;
//    a.getRankLadder();
//    cout << "a' inverse : " << endl;
//    cout << a.inverse();
//
//    cout << "a's det = " << a.getDet() << endl;
//    cout << "a'unit : " << endl;
//    cout << a.getUnit();
//    cout << "a'trans : " << endl;
//    cout << a.getTranspose();
//
//    cout << "a'ladder : " << endl;
//    cout << a.getRankLadder();
//    cout <<"a'rank = " << a.getRank() << endl;
//    Matrix<int> low = Matrix<int>(8,8);
//    cout << "a'up : " << endl;
//    cout << a.LU(low);
//    cout << "a'low : " << endl;
//    cout << low;
//    cin >> a;
//    cout << a;

//    内存测试-copy & "="
//    Matrix<long> a = Matrix<long>(6, 6);
//    a.readFromFile(fileA, 1);
//    Matrix<int> b = a;
//    cout << "b : " << endl;
//    cout << b;
//    Matrix<int> c = a.getSubMat(0, 0, 3, 3);
//    cout << "c (the submatrix of a) : " << endl;
//    cout << c;
//    Matrix<int> d = a.getSubMat(3, 3, 3, 3);
//    cout << "d (the submatrix of a) : " << endl;
//    cout << d;
//    Matrix<int> d;
//    d = a;
//    cout << "d : " << endl;
//    cout << d;
//    cout << d.getPointTime()[0];

//    构造函数测试
//    Matrix<int> testData = Matrix<int>(10);
//    cout << testData;
//
//    auto *arr = new double[10]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//    Matrix<double> testOpenCV = Matrix<double>(2, 5, arr);
//    cout << "&arr = " << arr << endl;
//    cout << "&testOpenCV.matrix = " << testOpenCV.getMatrix()<< endl;
//    cout << "pointTime = " << testOpenCV.getPointTime()[0] << endl;
//    cout << testOpenCV;
//    testOpenCV.outputToFile(output);

//    output << "a == b ? " << std::to_string(a == b) << endl;
//    output << "a != b ? " << std::to_string(a != b) << endl;

//    auto testEqual = (double)a;
//    output << "matrix a turn to a double, testEqual = " << std::to_string(testEqual) << endl;
//    long tem = 6;
//    cout << "tem = " << tem << endl;
//    a /= tem;
//    cout << "a /= 6 : " << tem << endl;
//    cout << a;

//    friend测试
//    Matrix<long> ap = 100 + a;
//    Matrix<long> as = 100 - a;
//    Matrix<long> am = 100 * a;
//    long ai = 100 / a;
//    cout << "Matrix<long> ap = 100 + a : " << endl ;
//    cout << ap;
//    cout << "Matrix<long> as = 100 - a : " << endl ;
//    cout << as;
//    cout << "Matrix<long> am = 100 * a : " << endl ;
//    cout << am;
//    cout << "long ai = 100 / a : " << endl ;
//    cout << ai;
//    ap.outputToFile(output);
//    as.outputToFile(output);
//    am.outputToFile(output);
//    output << "100 / a = " << std::to_string(ai) << endl;

//    基本加减乘除，除法定义为两个矩阵平均值相除;
//    Matrix<double> aAb = a + b;
//    Matrix<double> aSb = a - b;
//    Matrix<double> aMb = a * b;
//    double aCb = a / b;
//    aAb.outputToFile(output);
//    aSb.outputToFile(output);
//    aMb.outputToFile(output);
//    output << "a / b = " << std::to_string(aCb) << endl;

    output.close();
    auto end = std::chrono::steady_clock::now();
    chrono::duration<double, micro> elapsed = end  - start;
    cout << "running time : " << elapsed.count() / 1000 << " ms" << endl;
    return 0;
}

