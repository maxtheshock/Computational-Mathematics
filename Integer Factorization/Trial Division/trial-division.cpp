#include <iostream>
#include <chrono>
#include <map>
#include <algorithm>

typedef unsigned __int128 u128;

std::istream& operator>>(std::istream& is, u128& val) {
    std::string s;
    is >> s;
    val = 0;
    for (char c : s) {
        if (isdigit(c)) val = val * 10 + (c - '0');
    }
    return is;
}

std::ostream& operator<<(std::ostream& os, u128 val) {
    if (val == 0) {
        return os << "0";
    }
    std::string s;
    while (val > 0) {
        s += static_cast<char>('0' + (val % 10));
        val /= 10;
    }
    std::reverse(s.begin(), s.end());
    return os << s;
}

//
//

std::map<u128, int> factorize(u128 x) {
    std::map<u128, int> factors;
    u128 d = 2;
    while (d <= x / d) {
        if (x % d == 0) {
            while (x % d == 0) {
                factors[d]++;
                x /= d;
            }
        }
        ++d;
    }
    if (x > 1) {
        factors[x]++;
    }
    return factors;
}

void printMap(const std::map<u128, int>& obj) {
    auto last_it = std::prev(obj.end());
    for (auto it = obj.begin(); it != obj.end(); ++it) {
        (it->second == 1) ? (std::cout << it->first) : (std::cout << it->first << "^" << it->second);
        if (it != last_it) {
            std::cout << " * ";
        }
    }
    std::cout << std::endl;
}

int main() {
    u128 num;
    std::cin >> num;

    auto start = std::chrono::high_resolution_clock::now();
    std::map<u128, int> get = factorize(num);
    auto end = std::chrono::high_resolution_clock::now();
    
    printMap(get);
    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "Time: " << ms << " ms" << std::endl;
    return 0;
}
