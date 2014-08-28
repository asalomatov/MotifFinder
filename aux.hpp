#if !defined(_AUX_HPP)
#define _AUX_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <map>
#include <utility>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <queue>

const double EPS = 0.000000001;
std::string TimeStamp(void);
int OpenFile(const char *fname, const int appendYN, std::fstream &myfile);
void CloseFile(std::fstream &myfile);
template <typename T>
int WriteSomethingToFile(std::fstream &myfile, const T &s)
{
    if(!myfile.is_open()) 
        return 0;
    myfile << s;
    myfile.flush();
    return 1;
}
template <typename T>
void PrintVectorSizeToConsole (T& a) {
    std::cout << "\nvector size = " << a.size() <<std::endl;
}
template <typename T>
int PrintVectorToConsole (T& a, const unsigned int i_first, const unsigned int i_last) {
    if(a.empty()){std::cout << "PrintVectorToConsole: supplied vector is empty" <<std::endl;return 0;}
    if(i_first < 0 || i_first > a.size()-1){std::cout << "PrintVectorToConsole: _first is out of range" <<std::endl;return 0;}
    if(i_last < 0 || i_last > a.size()-1 || i_first > i_last){std::cout << "PrintVectorToConsole: i_last is out of range" <<std::endl;return 0;}
    for(unsigned int j = i_first; j <= i_last; j++){
        std::cout <<  *(a.begin() + j) << "\n";
    }
    std::cout <<std::endl;
    return 1;
}
bool Double1IsLessThanDouble2(const double d1, const double d2, const double eps);
bool Double1IsGreaterThanDouble2(const double d1, const double d2, const double eps);
bool Double1IsEqualToDouble2(const double d1, const double d2, const double eps);
double Round(double x, int i);
#endif		// _AUX_HPP
