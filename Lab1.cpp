#include <iostream>
#include <vector>
#include <limits>
#include <chrono>
#include <omp.h>
#include <random>

class MatrixProcessor {
private:
    std::vector<std::vector<int>> matrix;
    int rows;
    int cols;

public:
    MatrixProcessor(int r, int c) : rows(r), cols(c) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, 100);

        matrix.resize(rows, std::vector<int>(cols));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = distrib(gen);
            }
        }
    }

    void printMatrix() const {
        for (const auto& row : matrix) {
            for (int val : row) {
                std::cout << val << "\t";
            }
            std::cout << std::endl;
        }
    }

    int sequentialMaxOfMin() const {
        int max_of_min = std::numeric_limits<int>::min();

        for (int j = 0; j < cols; ++j) {
            int min_in_col = std::numeric_limits<int>::max();
            for (int i = 0; i < rows; ++i) {
                if (matrix[i][j] < min_in_col) {
                    min_in_col = matrix[i][j];
                }
            }
            if (min_in_col > max_of_min) {
                max_of_min = min_in_col;
            }
        }

        return max_of_min;
    }

    int parallelMaxOfMin() const {
        int max_of_min = std::numeric_limits<int>::min();

#pragma omp parallel for reduction(max:max_of_min)
        for (int j = 0; j < cols; ++j) {
            int min_in_col = std::numeric_limits<int>::max();
            for (int i = 0; i < rows; ++i) {
                if (matrix[i][j] < min_in_col) {
                    min_in_col = matrix[i][j];
                }
            }
            if (min_in_col > max_of_min) {
                max_of_min = min_in_col;
            }
        }

        return max_of_min;
    }

    template<typename Func>
    double measureTime(Func func) const {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        return elapsed.count();
    }
};

int main() {
    const int rows = 9;
    const int cols = 9;

    MatrixProcessor processor(rows, cols);

    std::cout << "Matrix:" << std::endl;
    processor.printMatrix();
    std::cout << std::endl;

    double seq_time = processor.measureTime([&]() {
        int result = processor.sequentialMaxOfMin();
        std::cout << "Sequential Max of Min: " << result << std::endl;
        });
    std::cout << "Sequential Time: " << seq_time << " seconds" << std::endl;

    double par_time = processor.measureTime([&]() {
        int result = processor.parallelMaxOfMin();
        std::cout << "Parallel Max of Min: " << result << std::endl;
        });
    std::cout << "Parallel Time: " << par_time << " seconds" << std::endl;

    return 0;
}