#include "aux.hpp"

std::string TimeStamp(void)
{
    time_t tt = time(NULL);
    struct tm tm;
    char buf[32];

    tm = *localtime(&tt);
    strftime(buf, 31, "%Y-%m-%d %H:%M:%S", &tm);
//    printf("%s\n", buf);
    return std::string(buf);
}
int OpenFile(const char *fname, const int appendYN, std::fstream &myfile)
{
    if(appendYN == 1)
	    myfile.open(fname, std::fstream::out | std::fstream::app);
    else
	    myfile.open(fname, std::fstream::out);
    if(!myfile.is_open()) 
        return 0;
//    myfile << std::fixed;
    return 1;
}
void CloseFile(std::fstream &myfile)
{
    if(myfile.is_open()) 
        myfile.close();
}
bool Double1IsLessThanDouble2(const double d1, const double d2, const double eps)
{
	if(d1 < d2 - eps) return true;
	else return false;
}
bool Double1IsGreaterThanDouble2(const double d1, const double d2, const double eps)
{
	bool ret = false;
	if(d1 > d2 + eps) 
		ret = true;
	return ret;
}
bool Double1IsEqualToDouble2(const double d1, const double d2, const double eps)
{
	if(d1 > d2 - eps && d1 < d2 + eps) return true;
	else return false;
}
double Round(double x, int i)
{
    return pow(10.0,-i)*std::trunc(pow(10.0,i)*x+.5);
}
