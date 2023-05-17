#ifndef _ASSERTION_H
#define _ASSERTION_H

#include <iostream>

#ifdef __FILE_NAME__
#undef __FILE_NAME__
#endif

#ifndef __FILE_NAME__
constexpr const char* __file_name(const char* path) {
    const char* file = path;
    while (*path) {
        if (*path++ == '/') {
            file = path;
        }
    }
    return file;
}
#define __FILE_NAME__ __file_name(__FILE__)
#endif // __FILE_NAME__

#ifdef USE_ASSERT
#define Assert(expr) \
    do{ \
        if (!static_cast<bool>(expr)) { \
            std::cerr << "Assertion Failed: " << #expr << '\n' \
                      << "              at: " << __extension__ __PRETTY_FUNCTION__ << '\n' \
                      << "                  [" << __FILE_NAME__ << " : " << __LINE__ << ']' << std::endl; \
            std::exit(1); \
        } \
    } while (0)

#else // !USE_ASSERT
#define Assert(expr) static_const<void>(0)
#endif // USE_ASSERT

#define Assert_endl "\n                  "
#define Assert_msg(expr, msg) \
    do{ \
        if (!static_cast<bool>(expr)) { \
            std::cerr << "Assertion Failed: " << #expr << '\n' \
                      << "         Message: " << msg << '\n' \
                      << "              at: " << __extension__ __PRETTY_FUNCTION__ << Assert_endl \
                      << "[" << __FILE_NAME__ << " : " << __LINE__ << ']' << std::endl; \
            std::exit(1); \
        } \
    } while (0)

#define Warn_endl "\n         "
#define Warn(msg) \
    std::cerr << "Warning: " << msg << '\n' \
              << "     at: " << __extension__ __PRETTY_FUNCTION__ << Warn_endl \
              << "[" << __FILE_NAME__ << " : " << __LINE__ << ']' << std::endl

#define Error_endl "\n       "
#define Error(msg) \
    do { \
        std::cerr << "Error: " << msg << '\n' \
                  << "   at: " << __extension__ __PRETTY_FUNCTION__ << Error_endl \
                  << "[" << __FILE_NAME__ << " : " << __LINE__ << ']' << std::endl; \
        std::exit(1); \
    } while(0)

#endif // _ASSERTION_H
