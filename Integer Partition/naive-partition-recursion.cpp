// Unoptimized version of integer partition realization
// Recursion is used here, but actually, can be improved by DP

#include <iostream>
#include <cstdint>

int64_t part(int n, int k) {
    if (n == 0 || n == 1) { return 1; }
    if (k == 1) { return 1; }
    if (n < 0) { return 0; }
    return (part(n, k-1) + part(n-k, k));
}

int main() {
    int n;
    std::cin >> n;
    std::cout << part(n, n) << std::endl;
    return 0;
}
