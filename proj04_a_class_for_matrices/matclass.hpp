#ifndef _MAT_H
#define _MAT_H

#include <iostream>
#include <fstream>
#pragma GCC optimize(3)
using namespace std;

template<typename T>
class Matrix {
private:
    int currentChannel;
    int rank;
    int column;
    int row;
    int channel;
    int size;
    int step;
    int *pointTime;
    T *matrix;
public:
    int getRow() {
        return this->row;
    }

    int getCurrentChannel() {
        return this->currentChannel;
    }

    void setCurrentChannel(int a) {
        if (a < 1 || a > this->channel) {
            this->currentChannel = 1;
            return;
        } else
            this->currentChannel = a;
    }

    int getSize() {
        return this->size;
    }

    T *getMatrix() {
        return this->matrix;
    }

    int *getPointTime() {
        return this->pointTime;
    }

    int getColumn() {
        return this->column;
    }

    int getStep() {
        return this->step;
    }

    int getChannel() {
        return this->channel;
    }

    int getRank() {
        return this->rank;
    }

    T getValue(int i, int j) const {
        if (i >= 0 && j >= 0 && i < this->row && j < this->column) {
            return this->matrix[(getCurrentChannel() - 1) * this->size + i * this->step + j];
        } else return 0;
    }

    void setMatPoint(int i, int j, T num) {
        if (i >= 0 && j >= 0 && i < this->row && j < this->column) {
            this->matrix[(getCurrentChannel() - 1) * this->size + i * this->step + j] = num;
        } else return;
    }

    bool release() {
        this->column = 0;
        this->row = 0;
        this->step = 0;
        this->rank = 0;
        this->channel = 0;
        this->currentChannel = 0;
        if (this->matrix != nullptr) {
            delete[] this->matrix;
            this->matrix = nullptr;
        }
        return true;
    }

    ~Matrix() {
        cout << "pointer = " << static_cast<void *>(this) 
        << ", matrix pointer = " << static_cast<void *>(this->matrix) 
        << ", pointTime = " << (*this->pointTime) << endl;
        if (this->pointTime[0] - 1 == 0)
            release();
        else this->pointTime[0]--;
    }

    bool create(int r = 0, int c = 0, int s = 0, int ra = 0, int ch = 1) {
        release();
        this->pointTime = new int{1};
        this->column = c;
        this->row = r;
        this->step = s;
        this->rank = ra;
        this->channel = ch;
        this->currentChannel = ch;
        this->size = this->column * this->row;
        if (this->row != 0 && this->column != 0) {
            this->matrix = new T[r * c * ch];
            for (int i = 0; i < this->row * this->column * this->channel; ++i) {
                this->matrix[i] = 0;
            }
        } else this->matrix = nullptr;
        return true;
    }

    Matrix();

    Matrix(int r, int c, int ch);

    Matrix(int r, int c);

    explicit Matrix(T a);

    Matrix(int r, int c, T *a);

    Matrix(int r, int c, int ch, T *a);

    Matrix(const Matrix &);

    Matrix operator+(const Matrix<T> &) const;

    Matrix operator-(const Matrix<T> &) const;

    friend Matrix<T> operator+(const T &a, const Matrix<T> &m) {
        if (m.matrix == nullptr) {
            cerr << "PLUS FAILURE! MATRIX ERROR!" << endl;
            return Matrix<T>();
        }
        Matrix<T> res = Matrix<T>(m.row, m.column);
        for (int i = 0; i < m.row; i++) {
            for (int j = 0; j < m.column; j++) {
                res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                        m.matrix[(m.currentChannel - 1) * m.size + i * m.step + j] + a;
            }
        }
        return res;
    }

    friend Matrix<T> operator-(const T &a, const Matrix<T> &m) {
        if (m.matrix == nullptr) {
            cerr << "SUBTRACT FAILURE! MATRIX ERROR!" << endl;
            return Matrix<T>();
        }
        Matrix<T> res = Matrix<T>(m.row, m.column);
        for (int i = 0; i < m.row; i++) {
            for (int j = 0; j < m.column; j++) {
                res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                        m.matrix[(m.currentChannel - 1) * m.size + i * m.step + j] - a;
            }
        }
        return res;
    }

    Matrix operator+(T) const;

    Matrix operator-(T) const;

    Matrix operator*(const Matrix<T> &) const;

    friend Matrix<T> operator*(const T &a, const Matrix<T> &m) {
        if (m.matrix == nullptr) {
            cerr << "MUL FAILURE! MATRIX ERROR!" << endl;
            return Matrix<T>();
        }
        Matrix<T> res = Matrix<T>(m.row, m.column);
        for (int i = 0; i < m.row; i++) {
            for (int j = 0; j < m.column; j++) {
                res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                        m.matrix[(m.currentChannel - 1) * m.size + i * m.step + j] * a;
            }
        }
        return res;
    }

    friend T operator/(const T &a, const Matrix<T> &b) {
        if (b.getAverage() == 0 || b.matrix == nullptr) {
            cerr << "DIVISION FAILURE!" << endl;
            return 0;
        }
        return a / b.getAverage();
    }

    Matrix operator*(T) const;

    T getMin() const;

    T getMax() const;

    T getAverage() const;

    Matrix operator/(T) const;

    T operator/(const Matrix<T> &) const;

    Matrix &operator+=(const Matrix<T> &);

    Matrix &operator-=(const Matrix<T> &);

    Matrix &operator*=(const Matrix<T> &);

    Matrix &operator/=(const Matrix<T> &);

    Matrix &operator-=(T);

    Matrix &operator/=(T);

    Matrix &operator*=(T);

    Matrix &operator+=(T);

    friend T &operator*=(T &a, const Matrix<T> &b) {
        a = a * b.getAverage();
        return a;
    }

    friend T &operator+=(T &a, const Matrix<T> &b) {
        a = a + b.getAverage();
        return a;
    }

    friend T &operator-=(T &a, const Matrix<T> &b) {
        a = a - b.getAverage();
        return a;
    }

    friend T &operator/=(T &a, const Matrix<T> &b) {
        if (b.getAverage() == 0) {
            cerr << "DIVISION FAILURE!" << endl;
        } else a = a / b.getAverage();
        return a;
    }

    bool operator==(const Matrix<T> &);

    bool operator!=(const Matrix<T> &);

    explicit operator T() const {
        return this->getAverage();
    }

    void printMatrix() const;

    friend std::ostream &operator<<(std::ostream &os, const Matrix<T> &a) {
        if (a.matrix == nullptr) {
            cerr << "OUTPUT ERROR!" << endl;
            return os;
        }
        std::string str = "row = " + std::to_string(a.row) + ", column = " + std::to_string(a.column);
        os << str << endl;
        for (int i = 0; i < a.row; ++i) {
            for (int j = 0; j < a.column; ++j) {
                string t = to_string(a.matrix[a.size * (a.currentChannel - 1) +
                                              i * a.step + j]) + " ";
                os << t;
            }
            os << endl;
        }
        return os;
    }

    friend std::istream &operator>>(std::istream &is, Matrix<T> &a) {
        for (int i = 0; i < a.row; i++) {
            for (int j = 0; j < a.column; j++) {
                is >> a.matrix[a.size * (a.currentChannel - 1) + i * a.step + j];
            }
        }
        return is;
    }

    Matrix getUnit();

    T getTrace() const;

    Matrix getTranspose() const;

    T getRowAve(int i) const;

    T getColumnAve(int i) const;

    Matrix removeRow(int) const;

    Matrix removeColumn(int) const;

    void changeSize(int i, int j);

    void outputToFile(std::ofstream &);

    void readFromFile(const string &, int);

    Matrix LU(Matrix<T> &) const;

    void swapRow(int, int);

    void scaleRow(int, T);

    T getDet() const;

    Matrix inverse();

    Matrix getRankLadder();

    Matrix getAdjoin() const {
        return this->inverse() * this->getDet();
    }

    Matrix &operator=(const Matrix<T> &);

    Matrix getSubMat(int i, int j, int r, int c);

private:

    void exchangeRow(T *, int, int, int, int);

    void mulRow(T *, int, T, int, int);

    void unitization(T *, int, int, int);

    void addRow(T *, int, int, T, int, int);

};

template<typename T>
Matrix<T>::Matrix() {
    cout << "Matrix()" << endl;
    this->matrix = nullptr;
    create();
}

template<typename T>
Matrix<T>::Matrix(const Matrix &a) {
    cout << "Matrix(const Matrix &)" << endl;
    this->matrix = nullptr;
    create(a.row, a.column, a.step, a.rank, a.channel);
    delete[] this->matrix;
    this->matrix = a.matrix;
    this->pointTime = a.pointTime;
    this->pointTime[0]++;
//    cout << "the address is " << this->pointTime << ", the num of point is " 
//    << this->pointTime[0] << endl;
}

template<typename T>
Matrix<T>::Matrix(int r, int c) {
    cout << "Matrix(int, int)" << endl;
    this->matrix = nullptr;
    create(r, c, c);
//    cout << "the address is " << this->pointTime << ", the num of point is " << this->pointTime[0] << endl;
}

template<typename T>
Matrix<T>::Matrix(int r, int c, T *a) {
    cout << "Matrix(int, int, T *)" << endl;
    this->matrix = nullptr;
    create(r, c, c);
    delete[] this->matrix;
    this->matrix = a;
}

template<typename T>
Matrix<T>::Matrix(int r, int c, int ch, T *a) {
    cout << "Matrix(int, int, int, T* &)" << endl;
    this->matrix = nullptr;
    create(r, c, c, 0, ch);
    delete[] this->matrix;
    this->matrix = a;
}

template<typename T>
Matrix<T>::Matrix(int r, int c, int ch) {
    cout << "Matrix(int, int, int)" << endl;
    this->matrix = nullptr;
    create(r, c, c, 0, ch);
    cout << *(this->pointTime) << endl;
    cout << this->pointTime << endl;
}

template<typename T>
Matrix<T>::Matrix(T a) {
    cout << "Matrix(T)" << endl;
    this->matrix = nullptr;
    create(1, 1, 1);
    this->matrix[0] = a;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &a) const {
    if (this->matrix == nullptr || a.matrix == nullptr) {
        cerr << "PIUS FAILURE! MATRIX ERROR!" << endl;
        return Matrix<T>();
    }
    if (a.column != this->column || this->row != a.row) {
        cerr << "PLUS FAILURE!" << endl;
        return Matrix<T>();
    }
    Matrix<T> res = Matrix<T>(a.row, a.column);
    for (int i = 0; i < a.row; i++) {
        for (int j = 0; j < a.column; j++) {
            res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                    a.matrix[(a.currentChannel - 1) * a.size + i * a.step + j] +
                    this->matrix[(this->currentChannel - 1) *
                                 this->size + i * a.step + j];
        }
    }
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &a) const {
    if (this->matrix == nullptr || a.matrix == nullptr) {
        cerr << "SUBTRACT FAILURE! MATRIX ERROR!" << endl;
        return Matrix<T>();
    }
    if (a.column != this->column || this->row != a.row) {
        cerr << "SUBTRACT FAILURE! ROW OR COLUMN NOT MATCH!" << endl;
        return Matrix<T>();
    }
    Matrix<T> res = Matrix<T>(a.row, a.column);
    for (int i = 0; i < a.row; i++) {
        for (int j = 0; j < a.column; j++) {
            res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                    a.matrix[(a.currentChannel - 1) * a.size + i * a.step + j] -
                    this->matrix[(this->currentChannel - 1) *
                                 this->size + i * a.step + j];
        }
    }
    return res;
}


template<typename T>
Matrix<T> Matrix<T>::operator+(T a) const {
    if (this->matrix == nullptr) {
        cerr << "PLUS FAILURE! MATRIX ERROR!" << endl;
        return Matrix<T>();
    }
    Matrix<T> res = Matrix<T>(this->row, this->column);
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                    a + this->matrix[(this->currentChannel - 1) *
                                     this->size + i * a.step + j];
        }
    }
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(T a) const {
    if (this->matrix == nullptr) {
        cerr << "SUBTRACT FAILURE! MATRIX ERROR!" << endl;
        return Matrix<T>();
    }
    Matrix<T> res = Matrix<T>(this->row, this->column);
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                    this->matrix[(this->currentChannel - 1) *
                                 this->size + i * a.step + j] - a;
        }
    }
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(T a) const {
    if (this->matrix == nullptr) {
        cerr << "MUL FAILURE! MATRIX ERROR!" << endl;
        return Matrix<T>();
    }
    Matrix<T> res = Matrix<T>(this->row, this->column);
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                    this->matrix[(this->currentChannel - 1) *
                                 this->size + i * a.step + j] * a;
        }
    }
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &a) const {
    if (this->matrix == nullptr || a.matrix == nullptr) {
        cerr << "MUL FAILURE! MATRIX ERROR!" << endl;
        return Matrix<T>();
    }
    if (this->column != a.row) {
        cerr << "MUL FAILURE! NOT MATCH!" << endl;
        return Matrix<T>();
    }
    Matrix<T> res = Matrix<T>(this->row, a.column);
    {
        T c;
#pragma omp parallel for
        for (int i = 0; i < this->row; i++) {
            for (int k = 0; k < this->column; k++) {
                c = this->matrix[(this->currentChannel - 1) *
                                 this->size + i * this->step + k];
                for (int j = 0; j < a.column; j++) {
                    res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] +=
                            c * a.matrix[(a.currentChannel - 1) * a.size + k * a.step + j];
                }
            }
        }
    }
    return res;
}

template<typename T>
Matrix<T> &Matrix<T>::operator/=(const Matrix<T> &a) {
    if (this->matrix == nullptr || a.matrix == nullptr) {
        cerr << "DIVIDE FAILURE!" << endl;
        return *this;
    }
    if (a.getAverage() == 0) {
        cerr << "DIVIDE FAILURE!" << endl;
        return *this;
    } else {
        T ave = this->getAverage();
        this->release();
        *this = Matrix<T>(ave / a.getAverage());
        return *this;
    }
}

template<typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &a) {
    if (this->matrix == nullptr || a.matrix == nullptr) {
        cerr << "SUBTRACT FAILURE! MATRIX ERROR!" << endl;
        return *this;
    }
    if (a.column != this->column || this->row != a.row) {
        cerr << "SUBTRACT FAILURE! NOT MATCH!" << endl;
        return *this;
    }
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            this->matrix[(this->currentChannel - 1) *
                         this->size + i * this->step + j] -= a.matrix[(a.currentChannel - 1) * a.size + i * a.step +
                                                                      j];
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &a) {
    if (this->matrix == nullptr || a.matrix == nullptr) {
        cerr << "PLUS FAILURE! MATRIX ERROR!" << endl;
        return *this;
    }
    if (a.column != this->column || this->row != a.row) {
        cerr << "PLUS FAILURE! NOT MATCH!" << endl;
        return *this;
    }
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            this->matrix[(this->currentChannel - 1) *
                         this->size + i * this->step + j] += a.matrix[(a.currentChannel - 1) * a.size + i * a.step +
                                                                      j];
        }
    }
    return *this;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &a) {
    if (this->row != a.row || this->column != a.column)
        return false;
    if (this->matrix == nullptr || a.matrix == nullptr)
        return false;
    for (int i = 0; i < a.row; ++i)
        for (int j = 0; j < a.column; ++j) {
            if (this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j] !=
                a.matrix[(a.currentChannel - 1) * a.size + i * a.step + j])
                return false;
        }
    return true;
}

template<typename T>
bool Matrix<T>::operator!=(const Matrix<T> &a) {
    if (this->row != a.row || this->column != a.column)
        return true;
    if (this->matrix == nullptr || a.matrix == nullptr)
        return false;
    for (int i = 0; i < a.row; ++i)
        for (int j = 0; j < a.column; ++j) {
            if (this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j] !=
                a.matrix[(a.currentChannel - 1) * a.size + i * a.step + j])
                return true;
        }
    return false;
}


template<typename T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &a) {
    if (this->row * a.column > this->size) {
        cerr << "MUL FAILURE! MEMORY LEAK!" << endl;
        return *this;
    }
    if (this->matrix == nullptr || a.matrix == nullptr) {
        cerr << "MUL FAILURE! MATRIX ERROR!" << endl;
        return *this;
    }
    if (this->column != a.row || this->column != a.column) {
        cerr << "MUL FAILURE! NOT MATCH!" << endl;
        return *this;
    }
    Matrix<T> res = Matrix<T>(this->row, a.column);
    {
        T c;
#pragma omp parallel for
        for (int i = 0; i < this->row; i++) {
            for (int k = 0; k < this->column; k++) {
                c = this->matrix[(this->currentChannel - 1) * this->size + i * this->step + k];
                for (int j = 0; j < a.column; j++) {
                    res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] +=
                            c * a.matrix[(a.currentChannel - 1) * a.size + k * a.step + j];
                }
            }
        }
    }
    for (int i = 0; i < res.size; ++i) {
        if (i < this->size)
            this->matrix[(this->currentChannel - 1) * this->size + i] = res.matrix[(res.currentChannel - 1) * res.size +
                                                                                   i];
        else this->matrix[(res.currentChannel - 1) * res.size + i] = 0;
    }
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(T a) const {
    if (a == 0 || this->matrix == nullptr) {
        cerr << "DIVISION FAILURE!" << endl;
        return Matrix<T>();;
    }
    Matrix<T> res = Matrix<T>(this->row, this->column);
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            res.matrix[(res.currentChannel - 1) * res.size + i * res.step + j] =
                    this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j] / a;
        }
    }
    return res;
}

template<typename T>
T Matrix<T>::operator/(const Matrix<T> &a) const {
    if (a.getAverage() == 0 || a.matrix == nullptr || this->matrix == nullptr) {
        cerr << "DIVISION FAILURE!" << endl;
        return 0;
    }
    return this->getAverage() / a.getAverage();
}


template<typename T>
Matrix<T> &Matrix<T>::operator/=(T a) {
    if (a == 0 || this->matrix == nullptr) {
        cerr << "DIVISION FAILURE!" << endl;
        return *this;
    }
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j] /= a;
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator*=(T a) {
    if (this->matrix == nullptr) {
        cerr << "MUL FAILURE!" << endl;
        return *this;
    }
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j] *= a;
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator+=(T a) {
    if (this->matrix == nullptr) {
        cerr << "PLUS FAILURE!" << endl;
        return *this;
    }
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j] += a;
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator-=(T a) {
    if (this->matrix == nullptr) {
        cerr << "SUBTRACT FAILURE!" << endl;
        return *this;
    }
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->column; j++) {
            this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j] -= a;
        }
    }
    return *this;
}


template<typename T>
void Matrix<T>::printMatrix() const {
    if (this->matrix == nullptr)
        return;
    cout << "row = " << this->row << ", column = " << this->column << endl;
    int count = 1;
    while (count <= this->channel) {
        cout << "channel " << count << " :" << endl;
        for (int i = 0; i < this->row; ++i) {
            for (int j = 0; j < this->column; ++j) {
                cout << this->getValue(i, j, count) << " ";
            }
            cout << endl;
        }
        count++;
    }
}


template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &a) {
    if (this == &a)
        return *this;
    else {
        if (this->matrix != nullptr)
            delete[] this->matrix;
        this->matrix = a.matrix;
        this->pointTime = a.pointTime;
        this->pointTime[0]++;
        this->column = a.column;
        this->row = a.row;
        this->step = a.step;
        this->size = a.size;
        this->channel = a.channel;
        this->currentChannel = a.currentChannel;
        this->rank = a.rank;
        cout << "the address is " << this->pointTime 
        << ", the num of point is " << this->pointTime[0] << endl;
        return *this;
    }
}

template<typename T>
Matrix<T> Matrix<T>::getSubMat(int i, int j, int r, int c) {
    if (this->matrix == nullptr) {
        cerr << "ERROR" << endl;
        return *this;
    }
    if (i < 0 || j < 0 || i > this->row - 1 || j > this->column - 1) {
        cerr << "POINT OUT OF MATRIX" << endl;
        return *this;
    }
    if (r < 0 || c < 0 || r > this->row || c > this->column || 
    r + i > this->row || j + c > this->column) {
        cerr << "SIZE OUT OF LIMITATION" << endl;
        return *this;
    }
    T *t = this->matrix + this->size * (this->currentChannel - 1) + i * this->step + j;
    Matrix<T> res = Matrix<T>(r, c, t);
    cout << "pointer matrix = " << t << endl;
    res.pointTime[0] = this->pointTime[0] + 1;
    res.step = this->column;
    return res;
}

template<typename T>
T Matrix<T>::getMin() const {
    if (this->matrix == nullptr) {
        cerr << "MIN CAL FAILURE!" << endl;
        return 0;
    }
    T min = 0;
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->column; ++j) {
            if (min > this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j])
                min = this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
        }
    }
    return min;
}

template<typename T>
T Matrix<T>::getMax() const {
    if (this->matrix == nullptr) {
        cerr << "MAX CAL FAILURE!" << endl;
        return 0;
    }
    T max = 0;
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->column; ++j) {
            if (max < this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j])
                max = this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
        }
    }
    return max;
}

template<typename T>
T Matrix<T>::getAverage() const {
    if (this->matrix == nullptr || this->row == 0 || this->column == 0) {
        cerr << "AVE CAL FAILURE!" << endl;
        return 0;
    }
    T sum = 0;
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->column; ++j) {
            sum += this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
        }
    }
    return sum / (this->column * this->row);
}

template<typename T>
Matrix<T> Matrix<T>::getUnit() {
    if (this->matrix == nullptr || this->column == 0 || this->row == 0) {
        cerr << "GET UNIT ERROR!" << endl;
        return Matrix<T>();
    }
    Matrix<T> unit = Matrix<T>(this->row, this->column);
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->column; ++j) {
            if (i == j)
                unit.matrix[i * this->step + j] = 1;
            else unit.matrix[i * this->step + j] = 0;
        }
    }
    return unit;
}

template<typename T>
Matrix<T> Matrix<T>::getTranspose() const {
    if (this->matrix == nullptr) {
        cerr << "GET TRANS ERROR!" << endl;
        return Matrix<T>();
    }
    Matrix<T> tran = Matrix<T>(this->column, this->row);
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->column; ++j) {
            tran.matrix[j * tran.step + i] = this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
        }
    }
    return tran;
}

template<typename T>
Matrix<T> Matrix<T>::LU(Matrix<T> &low) const {
    if (this->matrix == nullptr || this->row != this->column || low.matrix == nullptr || low.row != low.column) {
        cerr << "LU FACTORIZATION FAILURE" << endl;
        return Matrix<T>();
    }
    if (this->rank != this->row) {
        cerr << "LU FACTORIZATION FAILURE" << endl;
        return Matrix<T>();
    }
    Matrix<T> up = Matrix<T>(this->row, this->column);
    T temp;
    int side = this->row;
    T *r1 = new T[side];
    T *r2 = new T[side];
    for (int i = 0; i < side; i++) {
        for (int j = 0; j < side; j++) {
            if (i == j)
                low.matrix[i * side + i] = 1;
        }
    }
    for (int k = 0; k < side; k++) {
        for (int j = k; j < side; j++) {
            if (k == 0) {
                temp = 0;
            } else {
                for (int i = 0; i < k; i++) {
                    *(r1 + i) = low.matrix[k * side + i];
                }
                for (int i = 0; i < k; i++) {
                    *(r2 + i) = up.matrix[i * side + j];
                }
                temp = 0;
                for (int i = 0; i < k; i++) {
                    temp += *(r1 + i) * *(r2 + i);
                }
            }
            up.matrix[k * side + j] = this->matrix[(this->currentChannel - 1) * this->size + k * side + j] - temp;
        }

        for (int i = 0; i < side; ++i) {
            *(r1 + i) = 0;
            *(r2 + i) = 0;
        }
        temp = 0;
        for (int t = k + 1; t < side; t++) {
            if (k == 0) {
                temp = 0;
            } else {
                for (int i = 0; i < k; i++) {
                    *(r1 + i) = low.matrix[t * side + i];
                }
                for (int i = 0; i < k; i++) {
                    *(r2 + i) = up.matrix[i * side + k];
                }
                temp = 0;
                for (int i = 0; i < k; i++) {
                    temp += *(r1 + i) * *(r2 + i);
                }
            }
            if (up.matrix[k * side + k] != 0)
                low.matrix[t * side + k] =
                        (this->matrix[(this->currentChannel - 1) * this->size + t * side + k] - temp) /
                        up.matrix[k * side + k];
            else {
                cerr << "LU FACTORIZATION FAILURE" << endl;
                return Matrix<T>();
            }
        }
    }
    delete[] r1;
    delete[] r2;
    return up;
}

template<typename T>
void Matrix<T>::exchangeRow(T *a, int i, int j, int col, int c_mark) {
    int k;
    T t;
    for (k = 0; k < col - c_mark; ++k) {
        t = a[i * col + c_mark + k];
        a[i * col + c_mark + k] = a[j * col + c_mark + k];
        a[j * col + c_mark + k] = t;
    }
}

template<typename T>
void Matrix<T>::mulRow(T *a, int r, T k, int col, int c_mark) {
    for (int i = 0; i < col - c_mark; i++)
        a[r * col + c_mark + i] *= k;
}

template<typename T>
void Matrix<T>::unitization(T *a, int i, int col, int c) {
    T fir = a[i * col + c];
    for (int j = c; j < col; ++j) {
        a[i * col + j] /= fir;
    }
}

template<typename T>
void Matrix<T>::addRow(T *a, int r1, int r2, T k, int col, int c_mark) {
    for (int i = 0; i < col - c_mark; i++)
        a[col * r1 + c_mark + i] += a[r2 * col + c_mark + i] * k;
}

template<typename T>
Matrix<T> Matrix<T>::getRankLadder() {
    if (this->matrix == nullptr) {
        cerr << "ERROR" << endl;
        return Matrix<T>();
    }
    int r_mark = 0;
    int c_mark = 0;
    int sign_c;
    Matrix ladder = Matrix<T>(this->row, this->column);
    for (int i = 0; i < this->row * this->column; ++i) {
        ladder.matrix[i] = this->matrix[(this->currentChannel - 1) * this->size + i];
    }
    while (c_mark < ladder.column) {
        sign_c = 1;
        for (int i = r_mark; i < ladder.row; ++i) {
            if (ladder.matrix[i * ladder.column + c_mark] != 0) {
                if (i != r_mark) {
                    if (sign_c)
                        exchangeRow(ladder.matrix, r_mark, i, ladder.column, c_mark);
                    else {
                        T t = ladder.matrix[i * ladder.column + c_mark];
                        mulRow(ladder.matrix, i, ladder.matrix[r_mark * ladder.column + c_mark], ladder.column, c_mark);
                        addRow(ladder.matrix, i, r_mark, -t, ladder.column, c_mark);
                    }
                } else
                    unitization(ladder.matrix, i, ladder.column, c_mark);
                sign_c = 0;
            }
        }
        if (!sign_c) r_mark++;
        c_mark++;
    }
    ladder.rank = r_mark;
    this->rank = r_mark;
    return ladder;
}

template<typename T>
void Matrix<T>::swapRow(int i, int j) {
    if (this->matrix == nullptr)
        return;
    if (i < 0 || i > this->row || j < 0 || j > this->column)
        return;
    T *temp = new T[this->column];
    for (int k = 0; k < this->column; ++k)
        temp[k] = this->matrix[(this->currentChannel - 1) * this->size + i * this->step + k];
    for (int k = 0; k < this->column; ++k)
        this->matrix[(this->currentChannel - 1) * this->size + i * this->step + k] = this->matrix[
                (this->currentChannel - 1) * this->size + j * this->step + k];
    for (int k = 0; k < this->column; ++k)
        this->matrix[(this->currentChannel - 1) * this->size + j * this->step + k] = temp[k];
    delete[] temp;
}

template<typename T>
void Matrix<T>::scaleRow(int i, T mul) {
    if (this->matrix == nullptr)
        return;
    if (i < 0 || i > this->row)
        return;
    for (int k = 0; k < this->column; ++k) {
        this->matrix[(this->currentChannel - 1) * this->size + i * this->step + k] *= mul;
    }
}

template<typename T>
T Matrix<T>::getDet() const {
    if (this->column != this->row || this->matrix == nullptr) {
        cerr << "DET FAILURE!" << endl;
        return 0;
    }
    int m_swap = 0;
    int det = 1;
    int max_row_num = 0;
    T max = 0;
    T temp = 0;
    T *m = new T[this->column * this->row];
    int col = this->column;
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < col; ++j) {
            m[i * col + j] = this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
        }
    }
    for (int i = 0; i < this->row - 1; ++i) {
        max = abs(m[i * col + i]);
        max_row_num = i;
        for (int j = i + 1; j < this->row; ++j) {
            if (max < abs(m[j * col + i])) {
                max = abs(m[j * col + i]);
                max_row_num = j;
            }
        }
        if (max_row_num != i) {
            T *p = new T[col];
            for (int k = 0; k < col; ++k)
                p[k] = m[i * col + k];
            for (int k = 0; k < col; ++k)
                m[i * col + k] = m[max_row_num * col + k];
            for (int k = 0; k < col; ++k)
                m[max_row_num * col + k] = p[k];
            delete[] p;
            m_swap++;
        }
        for (int j = i + 1; j < this->row; ++j) {
            temp = -m[j * col + i] / m[i * col + i];
            for (int k = 0; k < col; ++k) {
                m[j * col + k] += m[i * col + k] * temp;
            }
        }
    }
    if (m_swap % 2 == 1)
        m_swap = -1;
    else m_swap = 1;
    det *= m_swap;
    for (int i = 0; i < col; ++i)
        det *= m[i * col + i];
    return det;
}

template<typename T>
Matrix<T> Matrix<T>::inverse() {
    if (this->matrix == nullptr) {
        cerr << "INVERSE FAILURE!" << endl;
        return Matrix<T>();
    }
    if (this->column != this->row || this->rank < this->row) {
        cerr << "INVERSE FAILURE!" << endl;
        return Matrix<T>();
    }
    T *mat = new T[this->column * this->row];
    Matrix<T> inv = Matrix<T>(this->row, this->column);
    for (int i = 0; i < this->row; ++i) {
        inv.matrix[i * this->step + i] = 1;
    }
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->column; ++j) {
            mat[i * this->column + j] = this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
        }
    }
    int m_row = inv.row;
    int m_column = inv.column;
    T max, temp;
    int max_row_num;
    int m_swap = 0;
    for (int i = 0; i < m_row - 1; i++) {
        max = abs(mat[i * m_column + i]);
        max_row_num = i;
        for (int j = i + 1; j < m_row; j++) {
            if (max < abs(mat[j * m_column + i])) {
                max = abs(mat[j * m_column + i]);
                max_row_num = j;
            }
        }
        if (i != max_row_num) {
            T *p = new T[this->column];
            for (int k = 0; k < this->column; ++k)
                p[k] = mat[i * this->column + k];
            for (int k = 0; k < this->column; ++k)
                mat[i * this->column + k] = mat[max_row_num * this->column + k];
            for (int k = 0; k < this->column; ++k)
                mat[max_row_num * this->column + k] = p[k];
            delete[] p;
            inv.swapRow(i, max_row_num);
            m_swap++;
        }
        for (int j = i + 1; j < m_row; j++) {
            temp = -mat[j * m_column + i] / mat[i * m_column + i];
            for (int k = 0; k < m_column; k++) {
                mat[j * m_column + k] += mat[i * m_column + k] * temp;
                inv.matrix[j * m_column + k] += inv.matrix[i * m_column + k] * temp;
            }
        }
    }
    for (int i = 0; i < m_row; i++) {
        temp = 1 / mat[i * m_column + i];
        for (int k = 0; k < inv.column; ++k) {
            mat[i * inv.column + k] *= temp;
        }
        inv.scaleRow(i, temp);
    }
    int i, j, k;
    for (i = m_row - 1; i > 0; i--) {
        for (j = i - 1; j >= 0; j--) {
            temp = -mat[j * m_column + i] / mat[i * m_column + i];
            for (k = 0; k < m_column; k++) {
                mat[j * m_column + k] += temp * mat[i * m_column + k];
                inv.matrix[j * m_column + k] += temp * inv.matrix[i * m_column + k];
            }
        }
    }
    delete[] mat;
    return inv;
}

template<typename T>
T Matrix<T>::getTrace() const {
    if (this->row != this->column || this->matrix == nullptr) {
        cerr << "FAILURE! TRACE CAN NOT BE CALCULATED!" << endl;
        return 0;
    }
    T trace = 0;
    for (int i = 0; i < this->row; ++i) {
        trace += this->matrix[(this->currentChannel - 1) * this->size + i * this->step + i];
    }
    return trace;
}

template<typename T>
Matrix<T> Matrix<T>::removeRow(int r) const {
    if (r >= this->row || r < 0 || this->matrix == nullptr) {
        cerr << "REMOVE FAILURE!" << endl;
        return *this;
    }
    Matrix aft = Matrix<T>(this->row - 1, this->column);
    int count = 0;
    for (int i = 0; i < this->row; ++i) {
        if (r == i)
            continue;
        for (int j = 0; j < this->column; ++j) {
            aft.matrix[(count) * aft.step + j] = this->matrix[(this->currentChannel - 1) * this->size +
                                                              i * this->step + j];
        }
        count++;
    }
    return aft;
}

template<typename T>
Matrix<T> Matrix<T>::removeColumn(int c) const {
    if (c >= this->row || c < 0 || this->matrix == nullptr) {
        cerr << "REMOVE FAILURE!" << endl;
        return *this;
    }
    Matrix aft = Matrix<T>(this->row, this->column - 1);
    int count = 0;
    for (int i = 0; i < this->row; ++i){
        for (int j = 0; j < this->column; ++j) {
            if (j == c)
                continue;
            aft.matrix[i * aft.step + (count++)] = this->matrix[(this->currentChannel - 1) * this->size +
                                                                i * this->step + j];
        }
        count = 0;
    }

    return aft;
}


template<typename T>
T Matrix<T>::getRowAve(int i) const {
    if (i >= this->row || i < 0 || this->matrix == nullptr)
        return 0;
    T sum = 0;
    for (int j = 0; j < this->column; ++j)
        sum += this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
    return sum / this->column;
}

template<typename T>
T Matrix<T>::getColumnAve(int i) const {
    if (i >= this->column || i < 0 || this->matrix == nullptr)
        return 0;
    T sum = 0;
    for (int k = 0; k < this->row; ++k) {
        for (int j = 0; j < this->column; ++j)
            if (j == i)
                sum += this->matrix[(this->currentChannel - 1) * this->size + k * this->step + j];
    }
    return sum / this->row;
}

template<typename T>
void Matrix<T>::outputToFile(ofstream &outfile) {
    try {
        outfile << "row = " << to_string(this->row) << ", column = " << to_string(this->column) << endl;
        int count = 0;
        for (int i = 0; i < this->row; ++i) {
            for (int j = 0; j < this->column; ++j) {
                outfile << to_string(*(this->matrix + (this->currentChannel - 1) * this->size + i * this->step + j))
                        << " ";
            }
            outfile << endl;
        }
        outfile << endl;
    } catch (int exception) {
        cerr << "ERROR" << endl;
    } catch (double) {
        cerr << "ERROR" << endl;
    }
}

template<typename T>
void Matrix<T>::readFromFile(const string &in, int type) {
    try {
        ifstream infile;
        infile.open(in);
        std::string t;
        int count = 0;
        while (count < this->channel) {
            for (int i = 0; i < this->row; ++i)
                for (int j = 0; j < this->column; ++j) {
                    infile >> t;
                    //1-int 2-long 3-unsigned-long 4-l-l 5-un-l-l 6-float 7-double 8-long-double
                    switch (type) {
                        case 1:
                            this->matrix[i * this->column + j] = stoi(t);
                            break;
                        case 2:
                            this->matrix[i * this->column + j] = stol(t);
                            break;
                        case 3:
                            this->matrix[i * this->column + j] = stoul(t);
                            break;
                        case 4:
                            this->matrix[i * this->column + j] = stoll(t);
                            break;
                        case 5:
                            this->matrix[i * this->column + j] = stoull(t);
                            break;
                        case 6:
                            this->matrix[i * this->column + j] = stof(t);
                            break;
                        case 7:
                            this->matrix[i * this->column + j] = stod(t);
                            break;
                        case 8:
                            this->matrix[i * this->column + j] = stold(t);
                            break;
                        default:
                            this->matrix[i * this->column + j] = stod(t);
                            break;
                    }
                }
            count++;
        }
        infile.close();
    } catch (int exception) {
        cerr << "ERROR" << endl;
    } catch (double) {
        cerr << "ERROR" << endl;
    }
}

template<class T>
void Matrix<T>::changeSize(int i, int j) {
    if (this->matrix == nullptr) {
        cerr << "ERROR!" << endl;
        return;
    }
    if (i < 0 || j < 0 || i * j != this->column * this->row)
        return;
    if (this->step == this->column) {
        this->row = i;
        this->column = j;
        this->step = this->column;
    } else {
        this->row = i;
        this->column = j;
    }
}

#endif