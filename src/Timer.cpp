#include <string>

#include "Timer.h"

Timer::Timer() {
    started = false;
}

void Timer::start() {
    if (started) {
        std::cout << "Error: Timer was already started!" << std::endl;
        return;
    }
    
    started = true;
    begin = clock.now();
}

std::string Timer::timeStamp() {
    time_t t = time(0);
    struct tm *now;
    now = localtime(&t);

    int years   = now->tm_year + 1900;
    int months  = now->tm_mon + 1;
    int days    = now->tm_mday;
    int hours   = now->tm_hour;
    int minutes = now->tm_min;

    std::string str_years   = std::to_string(years);
    std::string str_months  = std::to_string(months);
    std::string str_days    = std::to_string(days);
    std::string str_hours   = std::to_string(hours);
    std::string str_minutes = std::to_string(minutes);

    if (months < 10)
       str_months.insert(0, 1, '0');

    if (days < 10)
        str_days.insert(0, 1, '0');

    if (hours < 10)
       str_hours.insert(0, 1, '0');

    if (minutes < 10)
        str_minutes.insert(0, 1, '0');

    return str_years + "-" + str_months + "-" + str_days + "-" + str_hours + "-" + str_minutes;
}
