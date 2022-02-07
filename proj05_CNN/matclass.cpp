#include <iostream>
#include <fstream>
#include "matclass.hpp"
#include <opencv2/opencv.hpp>
#pragma GCC optimize(3)
using namespace std;

size_t Matrix::getRow() const {
    return this->row;
}

size_t Matrix::getColumn() const {
    return this->column;
}

size_t Matrix::getCurrentChannel() const {
    return this->currentChannel;
}

void Matrix::setCurrentChannel(size_t a) {
    if (a == 0 || a > this->channel) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : a's value is illegal" << endl;
        this->currentChannel = 1;
        return;
    } else
        this->currentChannel = a;
}

float *Matrix::getMatrix() {
    return this->matrix;
}

size_t *Matrix::getPointTime() {
    return this->pointTime;
}

size_t Matrix::getStep() const {
    return this->step;
}

size_t Matrix::getChannel() const {
    return this->channel;
}

float Matrix::getValue(size_t i, size_t j) const {
    if (i < this->row && j < this->column) {
        return this->matrix[(this->currentChannel - 1) * this->size + i * this->step + j];
    } else {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : row's or column's value is illegal" << endl;
        return 0;
    }
}

void Matrix::setMatPoint(size_t i, size_t j, float num) {
    size_t head = (this->currentChannel - 1) * this->size;
    if (i < this->row && j < this->column) {
        this->matrix[head + i * this->step + j] = num;
    } else {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : row's or column's value is illegal" << endl;
        return;
    }
}

bool Matrix::release() {
    this->column = 0;
    this->row = 0;
    this->step = 0;
    this->size = 0;
    this->channel = 0;
    this->currentChannel = 0;
    if (this->matrix != nullptr) {
        delete[] this->matrix;
        this->matrix = nullptr;
    }
    if (this->pointTime != nullptr) {
        delete this->pointTime;
        this->pointTime = nullptr;
    }
    return true;
}

Matrix::~Matrix() {
    if (this->pointTime[0] - 1 == 0)
        release();
    else this->pointTime[0]--;
}


bool Matrix::create(size_t r = 0, size_t c = 0, size_t s = 0, size_t ch = 1) {
    release();
    this->pointTime = new size_t{1};
    this->column = c;
    this->row = r;
    this->step = s;
    this->channel = ch;
    this->currentChannel = 1;
    this->size = this->column * this->row;
    if (this->row != 0 && this->column != 0) {
        this->matrix = new float[r * c * ch]{};
    } else this->matrix = nullptr;
    return true;
}

std::ostream &operator<<(std::ostream &os, const Matrix &a) {
    if (a.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : a.matrix is null" << endl;
        return os;
    }
    std::string str = "row = " + std::to_string(a.row) + ", column = " + std::to_string(a.column);
    os << str << endl;
    size_t he = a.size * (a.currentChannel - 1);
    for (int i = 0; i < a.row; ++i) {
        size_t tem = i * a.step + he;
        for (int j = 0; j < a.column; ++j) {
            string t = to_string(a.matrix[tem + j]) + " ";
            os << t;
        }
        os << endl;
    }
    return os;
}

std::istream &operator>>(std::istream &is, Matrix &a) {
    if (a.row == 0 || a.column == 0 || a.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : a.matrix is null" << endl;
        return is;
    }
    size_t he = a.size * (a.currentChannel - 1);
    for (int i = 0; i < a.row; i++) {
        size_t tem = he + i * a.step;
        for (int j = 0; j < a.column; j++) {
            is >> a.matrix[tem + j];
        }
    }
    return is;
}

Matrix::Matrix() {
    this->matrix = nullptr;
    create();
}

Matrix::Matrix(const Matrix &a) {
    if (a.row == 0 || a.column == 0 || a.matrix == nullptr)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : a.matrix is null" << endl;
    this->matrix = nullptr;
    this->pointTime = nullptr;
    create(a.row, a.column, a.step, a.channel);
    delete[] this->matrix;
    this->matrix = a.matrix;
    delete this->pointTime;
    this->pointTime = a.pointTime;
    this->pointTime[0]++;
}

Matrix::Matrix(size_t r, size_t c) {
    if (r == 0 || c == 0)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : row == 0 || column == 0" << endl;
    this->matrix = nullptr;
    this->pointTime = nullptr;
    create(r, c, c);
}

Matrix::Matrix(size_t r, size_t c, float *a) {
    if (r == 0 || c == 0)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error :row == 0 || column == 0" << endl;
    if (a == nullptr)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : a is null" << endl;
    this->matrix = nullptr;
    this->pointTime = nullptr;
    create(r, c, c);
    delete[] this->matrix;
    this->matrix = a;
}

Matrix::Matrix(size_t r, size_t c, size_t ch, float *a) {
    if (r == 0 || c == 0 || ch == 0)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : row == 0 || column == 0 || channel == 0" << endl;
    if (a == nullptr)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : a is null" << endl;
    this->matrix = nullptr;
    this->pointTime = nullptr;
    create(r, c, c, ch);
    delete[] this->matrix;
    this->matrix = a;
}

Matrix::Matrix(size_t r, size_t c, size_t ch, cv::Mat &img) {
    if (r == 0 || c == 0 || ch == 0)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : row == 0 || column == 0 || channel == 0" << endl;
    if (img.empty())
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : image is empty" << endl;
    this->matrix = nullptr;
    this->pointTime = nullptr;
    create(r, c, c, ch);
    float *ch1 = this->matrix;
    float *ch2 = this->matrix + this->size;
    float *ch3 = ch2 + this->size;
    size_t cnt = 0;
    uchar * cvImg = img.data;
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->column; ++j) {
            ch1[cnt] = (float)cvImg[2] / 255.0f;
            ch2[cnt] = (float)cvImg[1] / 255.0f;
            ch3[cnt] = (float)cvImg[0] / 255.0f;
            cnt++;
            cvImg += 3;
        }
    }
}

Matrix::Matrix(size_t r, size_t c, size_t ch) {
    if (r == 0 || c == 0 || ch == 0)
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : row == 0 || column == 0 || channel == 0" << endl;
    this->matrix = nullptr;
    this->pointTime = nullptr;
    create(r, c, c, ch);
}

Matrix::Matrix(float a) {
    this->matrix = nullptr;
    this->pointTime = nullptr;
    create(1, 1, 1);
    this->matrix[0] = a;
}

Matrix &Matrix::operator=(const Matrix &a) {
    if (this == &a)
        return *this;
    else {
        delete[] this->matrix;
        this->matrix = a.matrix;
        delete this->pointTime;
        this->pointTime = a.pointTime;
        this->pointTime[0]++;
        this->column = a.column;
        this->row = a.row;
        this->step = a.step;
        this->size = a.size;
        this->channel = a.channel;
        this->currentChannel = a.currentChannel;
        return *this;
    }
}

Matrix Matrix::getSubMat(size_t i, size_t j, size_t r, size_t c, size_t ch) {
    if (this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return *this;
    }
    if (i > this->row - 1 || j > this->column - 1) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : point(i or j) is out of range" << endl;
        return *this;
    }
    if (r > this->row || c > this->column || r + i > this->row || j + c > this->column) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : size is out of range" << endl;
        return *this;
    }
    if (ch == 0 || ch > this->channel) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : channel is out of range" << endl;
        return *this;
    }
    float *t = this->matrix + this->size * (ch - 1) + i * this->step + j;
    Matrix res = Matrix(r, c, t);
    res.pointTime[0] = this->pointTime[0] + 1;
    res.step = this->column;
    return res;
}

void Matrix::outputToFile(ofstream &outfile) {
    if (!outfile.is_open()) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : outfile is close" << endl;
        return;
    }
    if (this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return;
    }
    float *head = (this->currentChannel - 1) * this->size + this->matrix;
    for (int i = 0; i < this->row; ++i) {
        float *tem = head + i * this->step;
        for (int j = 0; j < this->column; ++j) {
            outfile << to_string(tem[j]) << " ";
        }
        outfile << endl;
    }
    outfile << endl;
}

void Matrix::readFromFile(const string &in) {
    if (in.empty()) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : input file is not exist" << endl;
        return;
    }
    if (this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return;
    }

    ifstream infile;
    infile.open(in);
    std::string t;
    size_t count = 0;
    size_t cnt = 0;
    while (count < this->channel) {
        size_t x = count * this->size;
        for (int i = 0; i < this->row; ++i) {
            size_t y = x + i * this->column;
            for (int j = 0; j < this->column; ++j) {
                infile >> t;
                this->matrix[y + j] = stof(t);
            }
        }
        count++;
    }
    infile.close();
}

void Matrix::changeSize(size_t i, size_t j) {
    if (this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return;
    }
    if (i * j != this->column * this->row) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : the whole size is changed" << endl;
        return;
    }
    if (this->step == this->column) {
        this->row = i;
        this->column = j;
        this->step = this->column;
    } else {
        this->row = i;
        this->column = j;
    }
}

//Imcol+MEC, convolution once, focus on kernel size = 3, internal function!
bool Matrix::convolutionImcolMec(Matrix &ker, Matrix &re, size_t stride) {
    size_t n = 3 * stride;
    Matrix t = Matrix(re.row, re.row * n + (9 - n), 1);
    float *ch = t.matrix;
    const float *p = this->matrix + (this->currentChannel - 1) * this->size;
    for (int i = 0; i < t.row; ++i) {
        size_t cnt = 0;
        while (cnt < t.column) {
            for (int k = 0; k < this->row; ++k) {
                size_t cnt_k = 0;
                size_t si = k * this->column;
                while (cnt_k < 3) {
                    ch[cnt + i * t.column] = p[i * stride + cnt_k + si];
                    cnt_k++;
                    cnt++;
                }
            }
        }
    }
    size_t head_k = ker.size * (ker.currentChannel - 1);
    size_t head_re = re.size * (re.currentChannel - 1);
    for (int z = 0; z < re.row; ++z) {
        Matrix tem = t.getSubMat(0, z * stride * ker.row, t.row, ker.size, 1);
        float c;
#pragma omp parallel for
        for (int i = 0; i < t.row; i++) {
            for (int k = 0; k < ker.size; k++) {
                c = tem.matrix[i * tem.step + k];
                for (int j = 0; j < 1; j++) {
                    re.matrix[head_re + z * re.column + i] += c * ker.matrix[head_k + k];
                }
            }
        }
    }
    return true;
}

bool Matrix::pad(Matrix &t, size_t padding) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (t.column == 0 || t.row == 0 || t.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : t.matrix is null" << endl;
        return false;
    }
    if (padding == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : padding is zero" << endl;
        return false;
    }
    for (int k = 0; k < t.channel; ++k) {
        size_t head_th = k * this->size;
        size_t head_t = k * t.size;
        size_t cnt = 0;
        if (padding == 1) {
            size_t edgeRow = this->row + 1;
            size_t edgeCol = this->column + 1;
            for (int i = 0; i < t.row; ++i) {
                if (i == 0 || i == edgeRow)
                    continue;
                else {
                    for (int j = 0; j < t.column; ++j) {
                        if (j == 0 || j == edgeCol)
                            continue;
                        else {
                            t.matrix[head_t + i * t.column + j] = this->matrix[head_th + cnt];
                            cnt++;
                        }
                    }
                }
            }
        } else if (padding == 2) {
            size_t edgeRow = this->row + 2;
            size_t edgeCol = this->column + 2;
            size_t edgeRow1 = this->row + 1;
            size_t edgeCol1 = this->column + 1;
            for (int i = 0; i < t.row; ++i) {
                if (i == 0 || i == 1 || i == edgeRow || i == edgeRow1)
                    continue;
                else {
                    for (int j = 0; j < t.column; ++j) {
                        if (j == 0 || j == 1 || j == edgeCol || j == edgeCol1)
                            continue;
                        else {
                            t.matrix[head_t + i * t.column + j] = this->matrix[head_th + cnt++];
                        }
                    }
                }
            }
        }
    }
    return true;
}

bool Matrix::convolutionImcolMecTotal(Matrix &re, Matrix &a, size_t stride, size_t preSize) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (a.column == 0 || a.row == 0 || a.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : a's matrix is null" << endl;
        return false;
    }
    if (re.column == 0 || re.row == 0 || re.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : re's matrix is null" << endl;
        return false;
    }
    if (stride == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : stride is zero" << endl;
        return false;
    }
    if (this->column < stride || this->row < stride) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this->column < stride || this->row < stride cannot calculate" << endl;
        return false;
    }
    if (preSize == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : perSize(inChannel) is zero" << endl;
        return false;
    }
    for (int i = 1; i <= re.channel; ++i) {
        re.setCurrentChannel(i);
        for (int j = 1; j <= preSize; ++j) {
            size_t curCh = (i - 1) * preSize + j;
            a.setCurrentChannel(curCh);
            this->setCurrentChannel(j);
            this->convolutionImcolMec(a, re, stride);
        }
    }
    return true;
}

bool Matrix::maxPooling(Matrix &re) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (re.column == 0 || re.row == 0 || re.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : re's matrix is null" << endl;
        return false;
    }
    size_t r = this->row - 1;
    size_t c = this->column - 1;
    size_t cnt;
    for (int k = 0; k < re.channel; ++k) {
        float *p = this->matrix + k * this->size;
        float *pRe = re.matrix + k * re.size;
        cnt = 0;
        for (int i = 0; i < r; i += 2) {
            for (int j = 0; j < c; j += 2) {
                size_t r1 = i * this->column;
                size_t r2 = (i + 1) * this->column;
                pRe[cnt] = std::max(std::max(std::max(p[r1 + j], p[r1 + j + 1]), p[r2 + j]), p[r2 + j + 1]);
                cnt++;
            }
        }
    }
    return true;
}

bool Matrix::unFold(Matrix &re) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (re.column == 0 || re.row == 0 || re.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : re's matrix is null" << endl;
        return false;
    }
    delete[] re.matrix;
    re.matrix = this->matrix;
    re.pointTime[0] = this->pointTime[0] + 1;
    return true;
}

bool Matrix::fullConnected(Matrix &re, const float *bia, Matrix &para) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (re.column == 0 || re.row == 0 || re.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : re's matrix is null" << endl;
        return false;
    }
    if (bia == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : bias is null" << endl;
        return false;
    }
    if (para.column == 0 || para.row == 0 || para.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : para's matrix is null" << endl;
        return false;
    }
    {
        float c;
#pragma omp parallel for
        for (int i = 0; i < this->row; i++) {
            for (int k = 0; k < this->column; k++) {
                c = this->matrix[i * this->step + k];
                for (int j = 0; j < para.column; j++) {
                    re.matrix[i + j] += c * para.matrix[k + j];
                }
            }
        }
    }
    re.matrix[0] += bia[0];
    re.matrix[1] += bia[1];
    return true;
}

//only padding = 1/2, kernel size = 3; include all the channel;
bool Matrix::changeIntoMulMat(Matrix &re, size_t padding, size_t kernelSize, size_t stride) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (re.column == 0 || re.row == 0 || re.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : re's matrix is null" << endl;
        return false;
    }
    if (kernelSize == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : kernel's size is zero" << endl;
        return false;
    }
    if (stride == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : stride is zero" << endl;
        return false;
    }
    if (kernelSize > this->row + 2 * padding){
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : kernel's size is illegal" << endl;
        return false;
    }
    size_t last = re.column - 1;
    size_t last1 = re.column - 2;
    size_t last2 = re.column - 3;
    size_t last3 = re.column - 4;
    size_t last4 = re.column - 5;
    size_t last5 = re.column - 6;
    for (int k = 1; k <= this->channel; ++k) {
        re.setCurrentChannel(k);
        this->setCurrentChannel(k);
        float *begin = re.matrix + (re.currentChannel - 1) * re.size;
        float *beginThis = this->matrix + (this->currentChannel - 1) * this->size;
        if (padding == 1) {
            for (size_t i = 0; i < re.row; ++i) {
                float *rePresent = begin + i * re.column;
                size_t i1 = 0;
                size_t i2;
                size_t firC = i * stride;
                size_t lastC = firC + kernelSize - 1;
                if (i == 0)
                    i2 = 0;
                else {
                    i2 = firC - 1;
                }
                for (int j = 0; j < re.column; ++j) {
                    if (j == 0 || j == 1 || j == 2)
                        continue;
                    if (j == last || j == last1 || j == last2)
                        continue;
                    if (i == 0) {
                        if (j % kernelSize == 0)
                            continue;
                        size_t availCol = kernelSize - 1;
                        if (i2 < availCol && i1 < this->row && i2 < this->column) {
                            rePresent[j] = beginThis[i1 * this->column + i2];
                            i2++;
                        }
                        if (i2 == availCol || i2 == this->column) {
                            i1++;
                            i2 = 0;
                        }
                    } else if (i > 0 && i < re.row - 1) {
                        size_t availCol = firC - 1 + kernelSize;
                        if (i2 < availCol && i1 < this->row && i2 < this->column) {
                            rePresent[j] = beginThis[i1 * this->column + i2];
                            i2++;
                        }
                        if (i2 == availCol || i2 == this->column) {
                            i1++;
                            i2 = firC - 1;
                        }
                    } else {
                        if (((this->row + 2) - (kernelSize - stride)) % stride == 0) {
                            if ((j + 1) % kernelSize == 0)
                                continue;
                        }
                        size_t availCol = firC - 1 + kernelSize;
                        if (i2 < availCol && i1 < this->row && i2 < this->column) {
                            rePresent[j] = beginThis[i1 * this->column + i2];
                            i2++;
                        }
                        if (i2 == availCol || i2 == this->column) {
                            i1++;
                            i2 = firC - 1;
                        }
                    }
                }
            }
        } else if (padding == 2) {
            for (size_t i = 0; i < re.row; ++i) {
                float *rePresent = begin + i * re.column;
                size_t i1 = 0;
                size_t i2;
                size_t firC = i * stride;
                size_t lastC = firC + kernelSize - 1;
                if (firC <= 2)
                    i2 = 0;
                else {
                    i2 = firC - 2;
                }
                for (int j = 0; j < re.column; ++j) {
                    if (j == 0 || j == 1 || j == 2 || j == 3 || j == 4 || j == 5)
                        continue;
                    if (j == last || j == last1 || j == last2 || j == last3 || j == last4 || j == last5)
                        continue;
                    if (i == 0) {
                        if (j % kernelSize == 0 || j % kernelSize == 1)
                            continue;
                        size_t availCol = kernelSize - 2;
                        if (i2 < availCol && i1 < this->row && i2 < this->column) {
                            rePresent[j] = beginThis[i1 * this->column + i2];
                            i2++;
                        }
                        if ((i2 == availCol || i2 == this->column)) {
                            i1++;
                            i2 = 0;
                        }
                    } else if (i == 1) {
                        if (stride >= 2) {
                            size_t availCol = firC - 2 + kernelSize;
                            if (i2 < availCol && i1 < this->row && i2 < this->column) {
                                rePresent[j] = beginThis[i1 * this->column + i2];
                                i2++;
                            }
                            if ((i2 == availCol || i2 == this->column)) {
                                i1++;
                                i2 = 0;
                            }
                        } else {
                            if (j % kernelSize == 0)
                                continue;
                            size_t availCol = firC - 2 + kernelSize;
                            if (i2 < availCol && i1 < this->row && i2 < this->column) {
                                rePresent[j] = beginThis[i1 * this->column + i2];
                                i2++;
                            }
                            if ((i2 == availCol || i2 == this->column)) {
                                i1++;
                                i2 = 0;
                            }
                        }
                    } else if (i > 1 && i < re.row - 2) {
                        size_t availCol = firC - 2 + kernelSize;
                        if (i2 < availCol && i1 < this->row && i2 < this->column) {
                            rePresent[j] = beginThis[i1 * this->column + i2];
                            i2++;
                        }
                        if (i2 == availCol || i2 == this->column) {
                            i1++;
                            i2 = firC - 2;
                        }
                    } else if (i == re.row - 2) {
                        if (((this->row + 4) - (kernelSize - stride)) % stride != 0) {
                            if ((j + 1) % kernelSize == 0)
                                continue;
                        }
                        size_t availCol = firC - 2 + kernelSize;
                        if (i2 < availCol && i1 < this->row && i2 < this->column) {
                            rePresent[j] = beginThis[i1 * this->column + i2];
                            i2++;
                        }
                        if (i2 == availCol || i2 == this->column) {
                            i1++;
                            i2 = firC - 2;
                        }
                    } else {
                        if (((this->row + 2) - (kernelSize - stride)) % stride == 0) {
                            if (j % kernelSize == 1 || j % kernelSize == 2)
                                continue;
                        } else {
                            if ((j + 1) % kernelSize == 0)
                                continue;
                        }
                        size_t availCol = firC - 2 + kernelSize;
                        if (i2 < availCol && i1 < this->row && i2 < this->column) {
                            rePresent[j] = beginThis[i1 * this->column + i2];
                            i2++;
                        }
                        if (i2 == availCol || i2 == this->column) {
                            i1++;
                            i2 = firC - 2;
                        }
                    }
                }
            }
        } else if (padding == 0) {
            for (size_t i = 0; i < re.row; ++i) {
                float *rePresent = begin + i * re.column;
                size_t i1 = 0;
                size_t i2;
                size_t firC = i * stride;
                i2 = firC;
                for (int j = 0; j < re.column; ++j) {
                    size_t availCol = firC + kernelSize;
                    if (i2 < availCol && i1 < this->row && i2 < this->column) {
                        rePresent[j] = beginThis[i1 * this->column + i2];
                        i2++;
                    }
                    if (i2 == availCol || i2 == this->column) {
                        i1++;
                        i2 = firC;
                    }
                }
            }
        }
    }
    return true;
}

bool Matrix::conLogicPad(Matrix &con, Matrix &re, size_t stride, size_t preSize) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (re.column == 0 || re.row == 0 || re.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : re's matrix is null" << endl;
        return false;
    }
    if (con.column == 0 || con.row == 0 || con.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : con's matrix is null" << endl;
        return false;
    }
    if (preSize == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : the pre size is zero" << endl;
        return false;
    }
    if (stride == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : stride is zero" << endl;
        return false;
    }
    for (int i = 1; i <= re.channel; ++i) {
        float *head_re = re.matrix + re.size * (i - 1);
        for (int j = 1; j <= preSize; ++j) {
            size_t curCh = (i - 1) * preSize + j;
            this->setCurrentChannel(j);
            float *head_k = con.matrix + con.size * (curCh - 1);
            for (int z = 0; z < re.row; ++z) {
                Matrix tem = this->getSubMat(0, z * stride * con.row, this->row, con.size, j);
                float c;
                size_t p2 = z * re.column;
#pragma omp parallel for
                for (int i1 = 0; i1 < this->row; i1++) {
                    size_t p3 = p2 + i1;
                    size_t p1 = i1 * tem.step;
                    for (int k = 0; k < con.size; k++) {
                        c = tem.matrix[p1 + k];
                        for (int j1 = 0; j1 < 1; j1++) {
                            head_re[p3] += c * head_k[k];
                        }
                    }
                }
            }
        }
    }
    return true;
}
// any padding\stride, any kernel size but only row = column;
bool Matrix::changeIntoMulMatGeneral(Matrix &re, size_t padding, size_t kernelSize, size_t stride) {
    if (this->column == 0 || this->row == 0 || this->matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : this matrix is null" << endl;
        return false;
    }
    if (re.column == 0 || re.row == 0 || re.matrix == nullptr) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : re's matrix is null" << endl;
        return false;
    }
    if (kernelSize == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : the kernel size is zero" << endl;
        return false;
    }
    if (stride == 0) {
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : stride is zero" << endl;
        return false;
    }
    if (kernelSize > this->row + 2 * padding){
        cerr << "file : " << __FILE__ << ", line : " << __LINE__ << ", function : " << __FUNCTION__
             << ", error : kernel's size is illegal" << endl;
        return false;
    }
    for (int k = 0; k < this->channel; ++k) {
        float *begin = re.matrix + k * re.size;
        float *beginThis = this->matrix + k * this->size;
        size_t mod = (this->row + 2 * padding - (kernelSize - padding)) % stride;
        size_t rBound = this->row + padding * 2 - padding;
        for (size_t i = 0; i < re.row; ++i) {
            float *rePresent = begin + i * re.column;
            size_t i1 = 0;
            size_t i2;
            size_t firC = i * stride;
            size_t i3 = padding - firC;
            size_t lastC = firC + kernelSize - 1;
            size_t i4 = kernelSize - (lastC - rBound + 1);
            size_t upBound = padding * kernelSize;
            size_t downBound = re.column - upBound;
            if (firC <= padding)
                i2 = 0;
            else {
                i2 = firC - padding;
            }
            for (int j = 0; j < re.column; ++j) {
                if (j < padding * kernelSize)
                    continue;
                if (j >= downBound)
                    continue;
                if (firC < padding) {
                    if (j % kernelSize < i3)
                        continue;
                }
                if (lastC >= rBound) {
                    if (j % kernelSize >= i4)
                        continue;
                }
                size_t availCol = firC - padding + kernelSize;
                if (i2 < availCol && i1 < this->row && i2 < this->column) {
                    rePresent[j] = beginThis[i1 * this->column + i2];
                    i2++;
                }
                if (i2 == availCol || i2 == this->column) {
                    i1++;
                    if (firC <= padding)
                        i2 = 0;
                    else {
                        i2 = firC - padding;
                    }
                }

            }
        }
    }
    return true;
}