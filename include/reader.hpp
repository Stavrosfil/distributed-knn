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
            int j     = 0;
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
            int j     = 0;
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
            int j     = 0;
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

        for (int i = -4; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j     = 0;
            int count = 0;
            if (i < 0)
                continue;
            while (std::getline(iss, num, ',')) {
                if (count++)
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

} // namespace rdFma

namespace rdMiniboone {

void mnbPid(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 130064;
    d = 50;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nMiniBooNE_PID\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/miniboone/MiniBooNE_PID.txt";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            
            while (std::getline(iss, num, ' ')) {
                X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

} // namespace rdMiniboone

namespace rdTvNewsCom {

void BBC(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 17720;
    d = 17;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nBBC\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/tv_news_com/BBC.txt";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            while (std::getline(iss, num, ' ')) {
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

void CNN(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 22545;
    d = 17;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nCNN\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/tv_news_com/CNN.txt";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            while (std::getline(iss, num, ' ')) {
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

void CNNIBN(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 33117;
    d = 17;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nCNNIBN\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/tv_news_com/CNNIBN.txt";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            while (std::getline(iss, num, ' ')) {
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

void NDTV(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 17051;
    d = 17;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nNDTV\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/tv_news_com/NDTV.txt";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            while (std::getline(iss, num, ' ')) {
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

void TIMESNOW(int& n, int& d, std::vector<double>& X, int process_rank)
{
    n = 39252;
    d = 17;

    if (process_rank == 0) {

        X.resize(n * d);

        std::cout << "\nTIMESNOW\t" << n << ", " << d << std::endl << std::endl;
        std::string filePath = "./dataset/tv_news_com/TIMESNOW.txt";
        std::ifstream myfile(filePath);
        std::ifstream input(filePath);
        std::string s;

        for (int i = 0; i < n; i++) {
            std::getline(input, s);
            std::istringstream iss(s);
            std::string num;
            int j = 0;
            while (std::getline(iss, num, ' ')) {
                    X[i * d + j++] = std::stof(num);
            }
        }

        std::cout << "data reading by process 0 completed\n";
    }
}

}// namespace rdTvNewsCom

#endif // __READER_H__