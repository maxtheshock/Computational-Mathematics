        ---------------- Integer partition problem ----------------

This part is dedicated to a mathematical problem of counting the number of integer
partitions of a natural number, i.e. the quantity of ways to represent some natural
number as a sum of smaller non-negative integers (without respecting their order).
For example, say, for n=5:

5 = 1+1+1+1+1 =
  = 2+1+1+1 =
  = 2+2+1 =                 // In this case, we say: p(5) = 7 
  = 3+1+1 =                 // (there are 7 ways to express 5 as a sum of smaller integers)
  = 3+2 =
  = 4+1 =
  = 5;

Here I want to show you 2 possible ways to compute these numbers: (1) using recursion,
which is closer to pure mathematics and theory, (2) dynamic programming - a way more practical


(1) More fundamental approach, better for understanding:

We introduce a function P(n,k): N^2 -> N, defined as the number of partitions of n with
all parts being no more than k. That is, to compute the number of all integer partitions
for a known n, we need to find P(n,n).

Using some combinatorial knowledge we derive the recursive formula:
P(n,k) = P(n,k-1) + P(n-k,k)        (*)

Now, together with some base cases took into account, we can theoretically compute the
respective number of integer partitions for any natural number we would like.

Complexity: O(2^n)


(2) Boost with using only 1D-array:

Instead of computing every intermediate step, we can iteratively increase allowed
upper bound for the parts, thus, becoming able to store then "overlapped" in array cells.

Under the matrix understanding we transform formula (*) to:
dp[j] = dp[j] + dp[j-i] <=> dp[j] += dp[j-i] (**)

Complexity: O(n^2)
Memory: O(n)