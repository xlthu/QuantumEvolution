#ifndef _EXP_PRINT_H
#define _EXP_PRINT_H

#include <iostream>
#include <vector>
#include <array>
#include <utility>

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    if (vec.empty()) return out << "[]";

    out << "[" << vec[0];
    for(auto it = vec.begin() + 1; it != vec.end(); ++it) out << ", " << *it;
    out << "]";

    return out;
}

template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& vec) {
    if (vec.empty()) return out << "[]";

    out << "[" << vec[0];
    for(auto it = vec.begin() + 1; it != vec.end(); ++it) out << ", " << *it;
    out << "]";

    return out;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& pair) {
    return out << "(" << pair.first << ", " << pair.second << ")";
}

#endif // _EXP_PRINT_H
