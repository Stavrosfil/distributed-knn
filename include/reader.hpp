#ifndef __READER_H__
#define __READER_H__

#include <iostream>
#include <stdlib.h>

namespace rdCorel {

void colorHist(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 68040;
    d = 32;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nCorel ColorHistogram\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/corel/ColorHistogram.asc";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            int count = 0;
            while (std::getline(iss, num, ' ')) {
                if (count++)
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

void colorMom(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 68040;
    d = 9;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nCorel ColorMoments\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/corel/ColorMoments.asc";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            int count = 0;
            while (std::getline(iss, num, ' ')) {
                if (count++)
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

void coocTex(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 68040;
    d = 16;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nCorel CoocTexture\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/corel/CoocTexture.asc";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            int count = 0;
            while (std::getline(iss, num, ' ')) {
                if (count++)
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

} // namespace rdCorel

namespace rdFma {

void features(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 106574;
    d = 518;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nFMA features\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/fma/features.csv";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            int count = 0;
            while (std::getline(iss, num, ',')) {
                if (count++)
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

} // namespace rdFma

// namespace rdMiniboone {

// void mnbPid(int& n, int& d, std::vector<double>& X, int process_rank)
// {
//     n = 130065;
//     d = 50;

//     if (process_rank == 0) {

//         X.resize(n * d);

//         std::cout << "\nMiniBooNE_PID\t" << n << ", " << d << std::endl << std::endl;
//         std::string filePath = "./dataset/miniboone/MiniBooNE_PID.txt";
//         std::ifstream myfile(filePath);
//         std::ifstream input(filePath);
//         std::string s;

//         for (int i = 0; i < 1; i++) {
//             std::getline(input, s);
//             std::istringstream iss(s);
//             std::string num;
//             int j = 0;
//             int count = 0;
//             while (std::getline(iss, num, ' ')) {
//                 if (count++)
//                     X[i * d + j++] = std::stof(num);
//             }
//         }
//         std::cout << X[0] << std::endl;
//         std::cout << "data reading by process 0 completed\n";
//     }
// }

// } // namespace rdMiniboone

// namespace rdTvNewsCom {

// void features(int& n, int& d, std::vector<double>& X, int process_rank)
// {
//     n = 106574;
//     d = 518;

//     if (process_rank == 0) {

//         X.resize(n * d);

//         std::cout << "\nFMA features\t" << n << ", " << d << std::endl << std::endl;
//         std::string filePath = "./dataset/fma/features.csv";
//         std::ifstream myfile(filePath);
//         std::ifstream input(filePath);
//         std::string s;

//         for (int i = 0; i < n; i++) {
//             std::getline(input, s);
//             std::istringstream iss(s);
//             std::string num;
//             int j = 0;
//             int count = 0;
//             while (std::getline(iss, num, ',')) {
//                 if (count++)
//                     X[i * d + j++] = std::stof(num);
//             }
//         }

//         std::cout << "data reading by process 0 completed\n";
//     }
// } // namespace rdTvNewsCom

#endif // __READER_H__