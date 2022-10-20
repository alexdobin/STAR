// system functions
#include <string>
#include <fstream>
#include <sstream>

std::string linuxProcMemory()
{
    std::ifstream t("/proc/self/status");
    std::stringstream buffer;
    buffer << t.rdbuf();


    std::string outString;
    while (buffer.good()) {
        std::string str1;
        std::getline(buffer,str1);
        if ( (str1.rfind("VmPeak",0) == 0) || 
             (str1.rfind("VmSize",0) == 0) ||
             (str1.rfind("VmHWM",0) == 0)  ||
             (str1.rfind("VmRSS",0) == 0) ) {
                 outString += str1+"; ";
             };
    };
    outString += '\n';

    return outString;
};