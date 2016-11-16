#include <string>

#include <chrono>
#include "Timer.h"

#ifdef _WIN32 
Timer::Timer() {
	started = false;
	start();
}

void Timer::start() {
	if (started) {
		std::cout << "Error: Timer was already started!" << std::endl;
		return;
	}

	started = true;
	begin = clock.now();
}

double Timer::elapsed() {
	end = clock.now();

	if (!started) {
		std::cout << "Error: Clock was not started!" << std::endl;
		return 0;
	}

	started = false;
	return (double)(std::chrono::duration_cast<nanoseconds>(end - begin).count() / 1000000000.);
}

std::string Timer::timeStamp() {
	time_t t = time(0);
	struct tm* now;
	now = localtime(&t);

	int years = now->tm_year + 1900;
	int months = now->tm_mon + 1;
	int days = now->tm_mday;
	int hours = now->tm_hour;
	int minutes = now->tm_min;

	std::string str_years = std::to_string(static_cast<long long>(years));
	std::string str_months = std::to_string(static_cast<long long>(months));
	std::string str_days = std::to_string(static_cast<long long>(days));
	std::string str_hours = std::to_string(static_cast<long long>(hours));
	std::string str_minutes = std::to_string(static_cast<long long>(minutes));

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
#endif
