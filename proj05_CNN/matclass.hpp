#pragma once

#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>

using namespace std;
#pragma GCC optimize(3)

class Matrix {
private:
    size_t currentChannel;
    size_t column;
    size_t row;
    size_t channel;
    size_t size;
    size_t step;
    size_t *pointTime;
    float *matrix;
public:
    size_t getRow() const;

    size_t getColumn() const;

    size_t getCurrentChannel() const;

    void setCurrentChannel(size_t a);

    float *getMatrix();

    size_t *getPointTime();

    size_t getStep() const;

    size_t getChannel() const;

    float getValue(size_t i, size_t j) const;

    void setMatPoint(size_t i, size_t j, float num);

    ~Matrix();

    bool release();

    bool create(size_t , size_t , size_t , size_t);

    Matrix();

    Matrix(size_t r, size_t c, size_t ch);

    Matrix(size_t r, size_t c);

    explicit Matrix(float a);

    Matrix(size_t r, size_t c, float *a);

    Matrix(size_t r, size_t c, size_t ch, float *a);

    Matrix(size_t r, size_t c, size_t ch, cv::Mat &);

    Matrix(const Matrix &);

    friend std::ostream &operator<<(std::ostream &os, const Matrix &a);

    friend std::istream &operator>>(std::istream &is, Matrix &a);

    void outputToFile(std::ofstream &);

    void readFromFile(const string &);

    Matrix &operator=(const Matrix &);

    Matrix getSubMat(size_t i, size_t j, size_t r, size_t c,size_t ch);

    bool pad(Matrix &, size_t );

    void changeSize(size_t i, size_t j);

    bool convolutionImcolMecTotal(Matrix & re, Matrix &, size_t,size_t);

    bool maxPooling(Matrix &);

    bool unFold(Matrix &);

    bool fullConnected(Matrix & , const float *, Matrix &);

    bool changeIntoMulMat(Matrix & re, size_t padding, size_t kernelSize, size_t stride);

    bool conLogicPad(Matrix & con, Matrix &re , size_t stride, size_t preSize);

    bool changeIntoMulMatGeneral(Matrix &re, size_t padding, size_t kernelSize, size_t stride);

private:

    bool convolutionImcolMec( Matrix &, Matrix & , size_t);
};