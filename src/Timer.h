#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>
#include <ctime>

typedef std::chrono::nanoseconds nanoseconds;
typedef std::chrono::microseconds microseconds;
typedef std::chrono::milliseconds milliseconds;
typedef std::chrono::seconds seconds;

class Timer {
    public:
        Timer();
        void start();
        template<class T> float stop();
        static std::string timeStamp();
    private:
        std::chrono::steady_clock clock;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        bool started;
};

template<class T>
inline float Timer::stop() {
    end = clock.now();

    if (!started) {
        std::cout << "Error: Clock was not started!" << std::endl;
        return 0;
    }

    started = false;    
    return (float) std::chrono::duration_cast<T>(end - begin).count();
}
#endif //TIMER_H
