Skeleton solver code.

Compile with "g++ -Wall -std=c++11 -O2 skeleton_solver.cpp"

Input is a m x n binary matrix where m, n <= 8. Specify m and n first in that order. Then the matrix. An example is below.

3 3

1 0 1

0 1 0

0 1 1

The output is a all possible 1-Dollo completions to the input.


The algorithm allows you to input multiple matrices in the same file. Simply add the matrix one after the other. An example is below.


3 3

1 0 1

0 1 0

0 1 1

3 3

0 1 0

1 1 0

1 1 1


The output is all possible 1-Dollo phylogenies to the first input followed by phylogenies of the second input.
