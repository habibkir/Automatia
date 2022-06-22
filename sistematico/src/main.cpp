#include <iostream>
#include <vector>
#include "matrice.hpp"

int main() {
    Matrix m ({{0, 2, 3, 4},
	       {2, 3, 6, 0},
	       {2, 7, 5, 8},
	       {1, 2, 3, 4}});

    Matrix n ({{0, 2, 3, 4},
	       {2, 0, 6, 0},
	       {2, 7, 0, 8},
	       {1, 2, 3, 1}});
	
    Matrix mn = m+n;

    for(int i = 0; i<4; ++i) {
	for(int j = 0; j<4; ++j) {
	    std::cout<<mn[i][j]<<' ';
	}
	std::cout<<std::endl;
    }
    return 0;
}
