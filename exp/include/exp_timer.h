#ifndef _EXP_TIMER_H
#define _EXP_TIMER_H

#include <iostream>
#include <chrono>

class Timer {
public:
    Timer(std::string msg) : msg(std::move(msg)) {
        t1 = std::chrono::steady_clock::now();
    }

    ~Timer() {
        auto t2 = std::chrono::steady_clock::now();

        std::cout << msg << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
    }

    std::string msg;
    std::chrono::steady_clock::time_point t1;
};

#endif // _EXP_TIMER_H
