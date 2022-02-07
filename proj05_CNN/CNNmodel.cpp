#include <opencv2/opencv.hpp>
#include "matclass.hpp"
#include <chrono>

using namespace std;

//BRG
//3x128x128 (channel, height, width)

//Explanatory variables: 
// kernel1, kernel2 and kernel3 are the data given in GitHub, they are all the kernel core.
// bias1, bias2, bias3 and blasLast are the data given in GitHub, too. They are the bias.
// the used constructor, which needs the arguments are row, column and channel respectively.

#pragma GCC optimize(3)

int main(int argc, char **argv) {
    auto start = chrono::steady_clock::now();
    char *fileCon1 = nullptr;
    char *fileCon2 = nullptr;
    char *fileCon3 = nullptr;
    char *filePara = nullptr;
    char *fileOut = nullptr;
    if (argc > 1) {
        fileCon1 = argv[1];
        fileCon2 = argv[2];
        fileCon3 = argv[3];
        filePara = argv[4];
        fileOut = argv[5];
    }
    // ofstream output;
    // output.open(fileOut);
   //first convolution (have a try, just the first version)
   cv::Mat img = cv::imread("face.jpg");
   Matrix image = Matrix(img.rows, img.cols, img.channels(), img);
   Matrix kernel1 = Matrix(3, 3, 3 * 16);
   kernel1.readFromFile(fileCon1);
   Matrix t = Matrix(image.getRow() + 1, image.getColumn() + 1, image.getChannel());
   image.pad(t, 1);
   size_t num1 = (t.getRow() - (3 - 2)) / 2;
   Matrix resCon1 = Matrix(num1, num1, 16);
   t.convolutionImcolMecTotal(resCon1, kernel1, 2, 3);
   auto *bias1 = new float[16]{0.5312736f, 2.0255265f, 0.5032426f, 1.0871441f, -0.16811907f, -1.4195297f, 1.3647283f,
                               1.4160137f, 2.0942433f, 0.37322155f, -0.98419213f, -1.6288463f, 0.11098604f, 1.7342286f,
                               0.8017651f, 0.11197768};
   size_t s1 = resCon1.getRow() * resCon1.getColumn();
   size_t te1, te2;
   for (int i = 0; i < resCon1.getChannel(); ++i) {
       te1 = i * s1;
       for (int j = 0; j < s1; ++j) {
           te2 = te1 + j;
           resCon1.getMatrix()[te2] += bias1[i];
           resCon1.getMatrix()[te2] = std::max(resCon1.getMatrix()[te2], 0.f);
       }
   }
   //first maxpooling
   Matrix maxpl1 = Matrix(32, 32, 16);
   resCon1.maxPooling(maxpl1);
   //second convolution
   Matrix kernel2 = Matrix(3, 3, 16 * 32);
   kernel2.readFromFile(fileCon2);
   Matrix resCon2 = Matrix(30, 30, 32);
   maxpl1.convolutionImcolMecTotal(resCon2, kernel2, 1, 16);
   auto *bias2 = new float[32]{-0.8984575f, 0.12183002f, -1.2174866f, 0.15396571f, 0.5613631f, 0.39191365f, 1.5804975f,
                               0.92713225f, -0.28752735f, -1.2465117f, 0.54117906f, 0.45883402f, -0.9297666f,
                               0.9465098f, -2.0196316f, 1.9819049f, 0.03240406f, 0.7440181f, 1.6237522f, -2.151311f,
                               1.003879f, -0.7344588f, -0.8972895f, -0.14268021f, -0.28659603f, 0.6522633f,
                               0.69248044f, -2.569141f, 1.0389522f, 0.6716339f, 0.31258696f, -2.8816516};
   size_t s2 = resCon2.getRow() * resCon2.getColumn();
   for (int i = 0; i < resCon2.getChannel(); ++i) {
        te1 = i * s2;
        for (int j = 0; j < s2; ++j) {
           te2 = te1 + j;
           resCon2.getMatrix()[te2] += bias2[i];
           resCon2.getMatrix()[te2] = std::max(resCon2.getMatrix()[te2], 0.f);
       }
   }
   //second maxpooling
   Matrix maxpl2 = Matrix(15, 15, 32);
   resCon2.maxPooling(maxpl2);
   //third convolution
   Matrix kernel3 = Matrix(3, 3, 32 * 32);
   kernel3.readFromFile(fileCon3);
   Matrix p = Matrix(17, 17, 32);
   maxpl2.pad(p, 1);
   Matrix resCon3 = Matrix(8, 8, 32);
   p.convolutionImcolMecTotal(resCon3, kernel3, 2, 32);
   auto *bias3 = new float[32]{-0.9876442f, -0.5730391f, 0.24176481f, -0.1955811f, -0.63301075f, 0.34056598f,
                               -0.9157471f, 0.65030223f, 0.26016375f, -0.46545726f, 0.33029428f, -0.5477017f,
                               -0.8768753f, 0.63696754f, -0.64720476f, -0.28603157f, -1.1642118f, -0.68296003f,
                               0.10889423f, -0.30630988f, 0.44217426f, -0.015639782f, -0.6834114f, 0.2375867f,
                               -0.8528528f, -0.2255232f, -0.5730451f, 0.41303927f, -0.9590689f, -0.083116576f,
                               -0.35299557f, 0.628371};
   size_t s3 = resCon3.getRow() * resCon3.getColumn();
   for (int i = 0; i < resCon3.getChannel(); ++i) {
       te1 = i * s3;
       for (int j = 0; j < s3; ++j) {
           te2 = te1 + j;
           resCon3.getMatrix()[te2] += bias3[i];
           resCon3.getMatrix()[te2] = std::max(resCon3.getMatrix()[te2], 0.f);
       }
   }
   //full connected
   Matrix para = Matrix(2, 2048, 1);
   para.readFromFile(filePara);
   float *blasLast = new float[2]{-0.0040076836f, 0.00010113005f};
   Matrix unfold = Matrix(2048, 1, 1);
   resCon3.unFold(unfold);
   Matrix result = Matrix(2, 1, 1);
   para.fullConnected(result, blasLast, unfold);

//    //first convolution(logic adding zero, the second version)
//    cv::Mat img = cv::imread("face.jpg");
//    Matrix image = Matrix(img.rows, img.cols, img.channels(), img);
//    Matrix kernel1 = Matrix(3, 3, 3 * 16);
//    kernel1.readFromFile(fileCon1);
//    Matrix t = Matrix(64, 390, image.getChannel());
//    image.changeIntoMulMat(t, 1, 3, 2);
//    size_t num1 = image.getRow() / 2;
//    Matrix resCon1 = Matrix(num1, num1, 16);
//    t.conLogicPad(kernel1, resCon1, 2, 3);
//    auto *bias1 = new float[16]{0.5312736f, 2.0255265f, 0.5032426f, 1.0871441f, -0.16811907f, -1.4195297f, 1.3647283f,
//                                1.4160137f, 2.0942433f, 0.37322155f, -0.98419213f, -1.6288463f, 0.11098604f, 1.7342286f,
//                                0.8017651f, 0.11197768};
//    size_t te1, te2;
//    for (int i = 0; i < resCon1.getChannel(); ++i) {
//        te1 = i * s1;
//        for (int j = 0; j < s1; ++j) {
//            te2 = te1 + j;
//            resCon1.getMatrix()[te2] += bias1[i];
//            resCon1.getMatrix()[te2] = std::max(resCon1.getMatrix()[te2], 0.f);
//        }
//    }
//    //first maxpooling
//    Matrix maxpl1 = Matrix(32, 32, 16);
//    resCon1.maxPooling(maxpl1);
//    //second convolution
//    Matrix kernel2 = Matrix(3, 3, 16 * 32);
//    kernel2.readFromFile(fileCon2);
//    Matrix t1 = Matrix(30, 96, maxpl1.getChannel());
//    maxpl1.changeIntoMulMat(t1, 0, 3, 1);
//    Matrix resCon2 = Matrix(30, 30, 32);
//    t1.conLogicPad(kernel2, resCon2, 1, 16);
//    auto *bias2 = new float[32]{-0.8984575f, 0.12183002f, -1.2174866f, 0.15396571f, 0.5613631f, 0.39191365f, 1.5804975f,
//                                0.92713225f, -0.28752735f, -1.2465117f, 0.54117906f, 0.45883402f, -0.9297666f,
//                                0.9465098f, -2.0196316f, 1.9819049f, 0.03240406f, 0.7440181f, 1.6237522f, -2.151311f,
//                                1.003879f, -0.7344588f, -0.8972895f, -0.14268021f, -0.28659603f, 0.6522633f,
//                                0.69248044f, -2.569141f, 1.0389522f, 0.6716339f, 0.31258696f, -2.8816516};
//    size_t s2 = resCon2.getRow() * resCon2.getColumn();
//    for (int i = 0; i < resCon2.getChannel(); ++i) {
//         te1 = i * s2;
//         for (int j = 0; j < s2; ++j) {
//            te2 = te1 + j;
//            resCon2.getMatrix()[te2] += bias2[i];
//            resCon2.getMatrix()[te2] = std::max(resCon2.getMatrix()[te2], 0.f);
//        }
//    }
//    //second maxpooling
//    Matrix maxpl2 = Matrix(15, 15, 32);
//    resCon2.maxPooling(maxpl2);
//    //third convolution
//    Matrix kernel3 = Matrix(3, 3, 32 * 32);
//    kernel3.readFromFile(fileCon3);
//    Matrix p = Matrix(8 , 51, 32);
//    maxpl2.changeIntoMulMat(p,1,3,2);
//    Matrix resCon3 = Matrix(8, 8, 32);
//    p.conLogicPad(kernel3, resCon3, 2, 32);
//    auto *bias3 = new float[32]{-0.9876442f, -0.5730391f, 0.24176481f, -0.1955811f, -0.63301075f, 0.34056598f,
//                                -0.9157471f, 0.65030223f, 0.26016375f, -0.46545726f, 0.33029428f, -0.5477017f,
//                                -0.8768753f, 0.63696754f, -0.64720476f, -0.28603157f, -1.1642118f, -0.68296003f,
//                                0.10889423f, -0.30630988f, 0.44217426f, -0.015639782f, -0.6834114f, 0.2375867f,
//                                -0.8528528f, -0.2255232f, -0.5730451f, 0.41303927f, -0.9590689f, -0.083116576f,
//                                -0.35299557f, 0.628371};
//    size_t s3 = resCon3.getRow() * resCon3.getColumn();
//    for (int i = 0; i < resCon3.getChannel(); ++i) {
//        te1 = i * s3;
//        for (int j = 0; j < s3; ++j) {
//            te2 = te1 + j;
//            resCon3.getMatrix()[te2] += bias3[i];
//            resCon3.getMatrix()[te2] = std::max(resCon3.getMatrix()[te2], 0.f);
//        }
//    }
//    //full connected
//    Matrix para = Matrix(2, 2048, 1);
//    para.readFromFile(filePara);
//    auto *blasLast = new float[2]{-0.0040076836f, 0.00010113005f};
//    Matrix unfold = Matrix(2048, 1, 1);
//    resCon3.unFold(unfold);
//    Matrix result = Matrix(2, 1, 1);
//    para.fullConnected(result, blasLast, unfold);


    // //first convolution(logic adding zero general, the final version)
    // cv::Mat img = cv::imread("face.jpg");
    // Matrix image = Matrix(img.rows, img.cols, img.channels(), img);
    // Matrix kernel1 = Matrix(3, 3, 3 * 16);
    // kernel1.readFromFile(fileCon1);
    // Matrix t = Matrix(64, 390, image.getChannel());
    // image.changeIntoMulMatGeneral(t, 1, 3, 2);
    // size_t num1 = image.getRow() / 2;
    // Matrix resCon1 = Matrix(num1, num1, 16);
    // t.conLogicPad(kernel1, resCon1, 2, 3);
    // auto *bias1 = new float[16]{0.5312736f, 2.0255265f, 0.5032426f, 1.0871441f, -0.16811907f, -1.4195297f, 1.3647283f,
    //                             1.4160137f, 2.0942433f, 0.37322155f, -0.98419213f, -1.6288463f, 0.11098604f, 1.7342286f,
    //                             0.8017651f, 0.11197768};
    // size_t s1 = resCon1.getRow() * resCon1.getColumn();
//    size_t te1, te2;
//    for (int i = 0; i < resCon1.getChannel(); ++i) {
//        te1 = i * s1;
//        for (int j = 0; j < s1; ++j) {
//            te2 = te1 + j;
//            resCon1.getMatrix()[te2] += bias1[i];
//            resCon1.getMatrix()[te2] = std::max(resCon1.getMatrix()[te2], 0.f);
//        }
//    }
    // //fiest maxpooling
    // Matrix maxpl1 = Matrix(32, 32, 16);
    // resCon1.maxPooling(maxpl1);
    // //second convolution
    // Matrix kernel2 = Matrix(3, 3, 16 * 32);
    // kernel2.readFromFile(fileCon2);
    // Matrix t1 = Matrix(30, 96, maxpl1.getChannel());
    // maxpl1.changeIntoMulMatGeneral(t1, 0, 3, 1);
    // Matrix resCon2 = Matrix(30, 30, 32);
    // t1.conLogicPad(kernel2, resCon2, 1, 16);
    // auto *bias2 = new float[32]{-0.8984575f, 0.12183002f, -1.2174866f, 0.15396571f, 0.5613631f, 0.39191365f, 1.5804975f,
    //                             0.92713225f, -0.28752735f, -1.2465117f, 0.54117906f, 0.45883402f, -0.9297666f,
    //                             0.9465098f, -2.0196316f, 1.9819049f, 0.03240406f, 0.7440181f, 1.6237522f, -2.151311f,
    //                             1.003879f, -0.7344588f, -0.8972895f, -0.14268021f, -0.28659603f, 0.6522633f,
    //                             0.69248044f, -2.569141f, 1.0389522f, 0.6716339f, 0.31258696f, -2.8816516};
    // size_t s2 = resCon2.getRow() * resCon2.getColumn();
 //    for (int i = 0; i < resCon2.getChannel(); ++i) {
//         te1 = i * s2;
//         for (int j = 0; j < s2; ++j) {
//            te2 = te1 + j;
//            resCon2.getMatrix()[te2] += bias2[i];
//            resCon2.getMatrix()[te2] = std::max(resCon2.getMatrix()[te2], 0.f);
//        }
//    }
    // //second maxpooling
    // Matrix maxpl2 = Matrix(15, 15, 32);
    // resCon2.maxPooling(maxpl2);
    // //third convolution
    // Matrix kernel3 = Matrix(3, 3, 32 * 32);
    // kernel3.readFromFile(fileCon3);
    // Matrix p = Matrix(8, 51, 32);
    // maxpl2.changeIntoMulMatGeneral(p, 1, 3, 2);
    // Matrix resCon3 = Matrix(8, 8, 32);
    // p.conLogicPad(kernel3, resCon3, 2, 32);
    // auto *bias3 = new float[32]{-0.9876442f, -0.5730391f, 0.24176481f, -0.1955811f, -0.63301075f, 0.34056598f,
    //                             -0.9157471f, 0.65030223f, 0.26016375f, -0.46545726f, 0.33029428f, -0.5477017f,
    //                             -0.8768753f, 0.63696754f, -0.64720476f, -0.28603157f, -1.1642118f, -0.68296003f,
    //                             0.10889423f, -0.30630988f, 0.44217426f, -0.015639782f, -0.6834114f, 0.2375867f,
    //                             -0.8528528f, -0.2255232f, -0.5730451f, 0.41303927f, -0.9590689f, -0.083116576f,
    //                             -0.35299557f, 0.628371};
    // size_t s3 = resCon3.getRow() * resCon3.getColumn();
//    for (int i = 0; i < resCon3.getChannel(); ++i) {
//        te1 = i * s3;
//        for (int j = 0; j < s3; ++j) {
//            te2 = te1 + j;
//            resCon3.getMatrix()[te2] += bias3[i];
//            resCon3.getMatrix()[te2] = std::max(resCon3.getMatrix()[te2], 0.f);
//        }
//    }
    // //full connected
    // Matrix para = Matrix(2, 2048, 1);
    // para.readFromFile(filePara);
    // auto *blasLast = new float[2]{-0.0040076836f, 0.00010113005f};
    // Matrix unfold = Matrix(2048, 1, 1);
    // resCon3.unFold(unfold);
    // Matrix result = Matrix(2, 1, 1);
    // para.fullConnected(result, blasLast, unfold);

    float dividend = std::exp(result.getMatrix()[0]) + std::exp(result.getMatrix()[1]);
    float res1 = std::exp(result.getMatrix()[0]) / dividend;
    float res2 = std::exp(result.getMatrix()[1]) / dividend;
    cout << "[c0, c1] = " << "[" << res1 << ", " << res2 << "]" << endl;
    cout << "probability of face: " << res2 << endl;
    
    // output.close();
    auto end = std::chrono::steady_clock::now();
    chrono::duration<double, micro> elapsed = end - start;
    cout << "running time : " << elapsed.count() / 1000 << " ms" << endl;
//    cv::waitKey();
    return 0;
}
