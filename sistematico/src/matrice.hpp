#include<vector>

class Matrix {
public:
    Matrix(std::vector<std::vector<int>> values): values(values){}
    Matrix(int nRows, int nCols) {
	values = std::vector<std::vector<int>>(nRows, std::vector<int>(nCols));
    }

    int getNRows() {
	return values.size();
    }
    int getNCols() {
	return values[0].size();
    }

    std::vector<int>& operator[](int index) {
	return values[index];
    }

    Matrix operator+(Matrix& right) {
	int nRows = this->getNRows();
	int nCols = this->getNCols();
	if((right.getNRows() != nRows) || (right.getNCols() != nCols)) {
	    std::cout<<"Questa somma non s'ha da fare\n";
	    return std::vector<std::vector<int>>{{0}};
	}

	Matrix res(nRows, nCols);
	for(int i = 0; i<nRows; ++i) {
	    for(int j = 0; j<nCols; ++j) {
		res[i][j] = values[i][j] + right[i][j];
	    }
	}
	return res;
    }

    Matrix operator*(Matrix& right) {
	int nRows = this->getNRows();
	int nCols = this->getNCols();
	/* non mi ricordo come si fa la moltiplicazione tra matrici */
	return right;
    }

private:
    std::vector<std::vector<int>> values;
};
