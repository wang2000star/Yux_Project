#include <vector>
#include <iostream>

using namespace std;
// Function to find the largest square sub-matrix with all 1s
vector<pair<int, int>> largestSquareSubMatrices(const vector<vector<int>>& matrix) {
    if (matrix.empty()) return {};
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    vector<vector<int>> dp(rows, vector<int>(cols, 0));
    int maxSize = 0;

    // Fill the dp matrix and find the maximum size of the sub-matrix
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (i == 0 || j == 0) {
                dp[i][j] = matrix[i][j];
            } else if (matrix[i][j] == 1) {
                dp[i][j] = min(dp[i-1][j], min(dp[i][j-1], dp[i-1][j-1])) + 1;
            }
            maxSize = max(maxSize, dp[i][j]);
        }
    }

    // Collect all top-left coordinates of sub-matrices of the maximum size
    vector<pair<int, int>> result;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (dp[i][j] == maxSize) {
                result.emplace_back(i - maxSize + 1, j - maxSize + 1);
            }
        }
    }

    return result;
}

// Function to print a matrix
void printMatrix(const vector<vector<int>>& matrix, const string& name) {
    cout << name << ":\n";
    for (const auto& row : matrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << "\n";
    }
}

// Function to find the largest sub-matrix with all 1s and return row and column indices
pair<vector<int>, vector<int>> largestSubMatrixIndices(const vector<vector<int>>& matrix) {
    if (matrix.empty()) return {{}, {}};
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    vector<vector<int>> dp(rows, vector<int>(cols, 0));
    int maxSize = 0;
    int maxRow = 0, maxCol = 0;

    // Fill the dp matrix and find the maximum size of the sub-matrix
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (i == 0 || j == 0) {
                dp[i][j] = matrix[i][j];
            } else if (matrix[i][j] == 1) {
                dp[i][j] = min(dp[i-1][j], min(dp[i][j-1], dp[i-1][j-1])) + 1;
            }
            if (dp[i][j] > maxSize) {
                maxSize = dp[i][j];
                maxRow = i;
                maxCol = j;
            }
        }
    }

    // Collect row and column indices of the largest sub-matrix
    vector<int> rowIndices, colIndices;
    for (int i = maxRow - maxSize + 1; i <= maxRow; ++i) {
        rowIndices.push_back(i);
    }
    for (int j = maxCol - maxSize + 1; j <= maxCol; ++j) {
        colIndices.push_back(j);
    }

    return {rowIndices, colIndices};
}

int main() {
    vector<vector<int>> A = {
        {1, 0, 0, 1, 1, 0, 0, -1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 2, 1, 1, -1, 0, 1, 1, 0, 1, 1, 1, 0},
        {0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, -1, 0, 1, 1, 0, 1, 1, 1, 0, -1, 0, 0, -1, -1, 0, 0, 1},
        {1, 1, 0, 1, 0, 1, -1, 0, 0, 1, 1, 1, 1, 0, -1, -1, 1, 2, 1, 2, 1, 1, -2, -1, 0, 1, 1, 1, 1, 0, -1, -1},
        {0, 1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 1, 0, 1, -1, 0, 0, 1, 1, 1, 1, 0, -1, -1, -1, -1, 0, -1, 0, -1, 1, 0},
        {1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, -1, 1, 2, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, -1, 1},
        {1, 1, 1, 0, 1, 0, -1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, -1, 1, -1, 0, -1, -1, 0, -1, -1, 0},
        {-1, 0, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 0, 1, 0, -2, 1, 0, 2, 2, 1, 0, 1, -1, 1, -1, 1, 1, 0, 1, 0},
        {-1, 1, -1, 1, 1, 0, 1, 0, -1, 0, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 0, 1, 0, 1, 0, -1, -1, -1, -1, 1, -1},
        {1, 0, 0, 1, 1, 0, 0, -1, 1, 1, 1, 1, 2, 1, 1, -1, 0, 1, 1, 0, 1, 1, 1, 0, -1, -1, -1, -1, -2, -1, -1, 1},
        {1, 1, 1, 1, 2, 1, 1, -1, 1, 0, 0, 1, 1, 0, 0, -1, 1, 0, 0, 1, 1, 0, 0, -1, 0, -1, -1, 0, -1, -1, -1, 0},
        {1, 1, 0, 1, 0, 1, -1, 0, 1, 2, 1, 2, 1, 1, -2, -1, 0, 1, 1, 1, 1, 0, -1, -1, -1, -2, -1, -2, -1, -1, 2, 1},
        {1, 2, 1, 2, 1, 1, -2, -1, 1, 1, 0, 1, 0, 1, -1, 0, 1, 1, 0, 1, 0, 1, -1, 0, 0, -1, -1, -1, -1, 0, 1, 1},
        {1, 0, 1, 1, 0, 1, 1, 0, 2, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, -1, 1, -2, -1, -2, -1, -1, -1, 0, -1},
        {2, 1, 2, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, -1, -1, -1, 0, -1, 0, 1, -1},
        {-1, 0, 1, 1, 1, 1, -1, 1, -2, 1, 0, 2, 2, 1, 0, 1, -1, 1, -1, 1, 1, 0, 1, 0, 2, -1, 0, -2, -2, -1, 0, -1},
        {-2, 1, 0, 2, 2, 1, 0, 1, -1, 0, 1, 1, 1, 1, -1, 1, -1, 0, 1, 1, 1, 1, -1, 1, 1, -1, 1, -1, -1, 0, -1, 0},
        {1, 1, 1, 1, 2, 1, 1, -1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, -1, -1, 1, 0, -1, -1, -1},
        {0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 2, 1, 1, -1, 1, 0, 0, 1, 1, 0, 0, -1, 0, 1, 1, 0, 1, 1, 1, 0},
        {1, 2, 1, 2, 1, 1, -2, -1, 0, 1, 1, 1, 1, 0, -1, -1, 0, 1, 1, 1, 1, 0, -1, -1, 1, 0, -1, 0, -1, 1, 0, 1},
        {0, 1, 1, 1, 1, 0, -1, -1, 1, 2, 1, 2, 1, 1, -2, -1, 1, 1, 0, 1, 0, 1, -1, 0, 0, 1, 1, 1, 1, 0, -1, -1},
        {2, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, -1, 1, 1, 1, 1, 0, 1, 0, -1, 1, 0, -1, 0, 1, -1, 1, 2, -1},
        {1, 1, 1, 0, 1, 0, -1, 1, 2, 1, 2, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, -1, 1},
        {-2, 1, 0, 2, 2, 1, 0, 1, -1, 1, -1, 1, 1, 0, 1, 0, -1, 1, -1, 1, 1, 0, 1, 0, 0, -1, 2, 0, 0, 1, -2, 1},
        {-1, 1, -1, 1, 1, 0, 1, 0, -2, 1, 0, 2, 2, 1, 0, 1, -1, 0, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 0, 1, 0},
        {1, -1, -1, 1, 0, -1, -1, -1, 1, -1, -1, 1, 0, -1, -1, -1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0},
        {0, -1, -1, 0, -1, -1, -1, 0, 1, 1, 1, 1, 2, 1, 1, -1, 1, 1, 1, 1, 2, 1, 1, -1, 1, -1, -1, 1, 0, -1, -1, -1},
        {1, 0, -1, 0, -1, 1, 0, 1, 1, 0, -1, 0, -1, 1, 0, 1, 0, 1, 1, 1, 1, 0, -1, -1, 0, 1, 1, 1, 1, 0, -1, -1},
        {0, -1, -1, -1, -1, 0, 1, 1, 1, 2, 1, 2, 1, 1, -2, -1, 1, 2, 1, 2, 1, 1, -2, -1, 1, 0, -1, 0, -1, 1, 0, 1},
        {0, -1, 0, 1, -1, 1, 2, -1, 0, -1, 0, 1, -1, 1, 2, -1, 1, 1, 1, 0, 1, 0, -1, 1, 1, 1, 1, 0, 1, 0, -1, 1},
        {-1, -1, -1, 0, -1, 0, 1, -1, 2, 1, 2, 1, 1, 1, 0, 1, 2, 1, 2, 1, 1, 1, 0, 1, 0, -1, 0, 1, -1, 1, 2, -1},
        {0, -1, 2, 0, 0, 1, -2, 1, 0, -1, 2, 0, 0, 1, -2, 1, -1, 1, -1, 1, 1, 0, 1, 0, -1, 1, -1, 1, 1, 0, 1, 0},
        {1, -1, 1, -1, -1, 0, -1, 0, -2, 1, 0, 2, 2, 1, 0, 1, -2, 1, 0, 2, 2, 1, 0, 1, 0, -1, 2, 0, 0, 1, -2, 1}
    };

    vector<vector<int>> A1(A.size(), vector<int>(A[0].size()));
    vector<vector<int>> A2(A.size(), vector<int>(A[0].size()));
    vector<vector<int>> A11(A.size(), vector<int>(A[0].size()));
    vector<vector<int>> A12(A.size(), vector<int>(A[0].size()));
    vector<vector<int>> A21(A.size(), vector<int>(A[0].size()));
    vector<vector<int>> A22(A.size(), vector<int>(A[0].size()));

    // Step 1: Decompose A into A1 and A2
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            if (A[i][j] == 2) {
                A1[i][j] = 1;
                A2[i][j] = 1;
            } else if (A[i][j] == -2) {
                A1[i][j] = -1;
                A2[i][j] = -1;
            } else {
                A1[i][j] = A[i][j];
                A2[i][j] = 0;
            }
        }
    }

    // Step 2: Decompose A1 into A11 and A12
    for (size_t i = 0; i < A1.size(); ++i) {
        for (size_t j = 0; j < A1[i].size(); ++j) {
            if (A1[i][j] == 1) {
                A11[i][j] = 1;
                A12[i][j] = 0;
            } else if (A1[i][j] == -1) {
                A11[i][j] = 0;
                A12[i][j] = 1;
            } else {
                A11[i][j] = 0;
                A12[i][j] = 0;
            }
        }
    }

    // Step 3: Decompose A2 into A21 and A22
    for (size_t i = 0; i < A2.size(); ++i) {
        for (size_t j = 0; j < A2[i].size(); ++j) {
            if (A2[i][j] == 1) {
                A21[i][j] = 1;
                A22[i][j] = 0;
            } else if (A2[i][j] == -1) {
                A21[i][j] = 0;
                A22[i][j] = 1;
            } else {
                A21[i][j] = 0;
                A22[i][j] = 0;
            }
        }
    }

    // Print the matrices for verification
    auto printMatrix = [](const vector<vector<int>>& matrix, const string& name) {
        cout << name << ":\n";
        for (const auto& row : matrix) {
            for (int val : row) {
                cout << val << " ";
            }
            cout << "\n";
        }
    };

    printMatrix(A1, "A1");
    printMatrix(A2, "A2");
    printMatrix(A11, "A11");
    printMatrix(A12, "A12");
    printMatrix(A21, "A21");
    printMatrix(A22, "A22");
    // 检验A = A11-A12+A21-A22
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            if (A[i][j] != A11[i][j] - A12[i][j] + A21[i][j] - A22[i][j]) {
                cout << "Decomposition failed\n";
                return 1;
            }
        }
    }
    cout << "Decomposition successful\n";

    // 对于A11,我想知道最大的全1子矩阵，这个子矩阵就是任取m行任取m列交叉形成，怎么找到它
    {
        auto [rowIndices, colIndices] = largestSubMatrixIndices(A11);
        cout << "Row indices of the largest sub-matrix with all 1s: ";
        for (int row : rowIndices) {
            cout << row << " ";
        }
        cout << "\n";
        cout << "\nColumn indices of the largest sub-matrix with all 1s: ";
        for (int col : colIndices) {
            cout << col << " ";
        }
        cout << "\n";
    }
    // 求A11每行1的个数
    {
        vector<int> rowOnes(A11.size(), 0);
        for (size_t i = 0; i < A11.size(); ++i) {
            for (size_t j = 0; j < A11[i].size(); ++j) {
                if (A11[i][j] == 1) {
                    rowOnes[i]++;
                }
            }
        }
        cout << "Number of 1s in each row of A11:\n";
        for (int val : rowOnes) {
            cout << val-1 << " ";
        }
        cout << "\n";
    }

return 0;
}
