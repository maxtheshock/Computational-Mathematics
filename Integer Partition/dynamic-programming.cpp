// Dynamic programming approach, we reduce 2D-matrix
// to many-times-overlapping vector dp

#include <iostream>
#include <cstdint>
#include <vector>

int64_t partition(int n) {
    std::vector<int64_t> dp(n+1, 0);
    dp[0] = 1; // base case
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j <= n; ++j) {
            dp[j] += dp[j-i];
        }
    }
    return dp[n];
}

int main() {
    int n;
    std::cin >> n;
    std::cout << partition(n) << std::endl;
    return 0;
}