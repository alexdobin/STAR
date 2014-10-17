#include <string>
#include <time.h>

std::string timeMonthDayTime() {
    time_t rawTime;
    char timeChar[100];
    time(&rawTime);
    strftime(timeChar,80,"%b %d %H:%M:%SS",localtime(&rawTime));
    std::string timeString=timeChar;
    timeString.erase(timeString.end()-1,timeString.end());
    return timeString;
};

std::string timeMonthDayTime(time_t &rawTime) {
    char timeChar[100];
    strftime(timeChar,80,"%b %d %H:%M:%SS",localtime(&rawTime));
    std::string timeString=timeChar;
    timeString.erase(timeString.end()-1,timeString.end());
    return timeString;
};
