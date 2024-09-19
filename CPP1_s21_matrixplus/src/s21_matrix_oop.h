#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  // attributes
  int rows_, cols_;  // rows and columns attributes
  double** matrix_;  // pointer to the memory where the matrix will be allocated

 public:
  int GetRows() const;
  int GetCols() const;
  void SetRows(int num);
  void SetCols(int num);
  void EditPointer(double** matrix);
  void EditPointers(double* matrix, int i);
  double** GetPointer() const;
  void PrintMatrix();

  S21Matrix();  // A basic constructor that initialises a matrix of some
                // predefined dimension.
  S21Matrix(int rows, int cols);  //Параметризированный конструктор с
                                  //количеством строк и столбцов.
  S21Matrix(const S21Matrix& other);  //Конструктор копирования.
  ~S21Matrix();                       // destructor
  S21Matrix(S21Matrix&& other);

  void SumMatrix(
      const S21Matrix& other);  // Adds the second matrix to the current one
  bool EqMatrix(const S21Matrix& other)
      const;  // Checks matrices for equality with each other.
  void SubMatrix(
      const S21Matrix& other);  // Subtracts another matrix from the current one
  void MulMatrix(const S21Matrix& other);  // Multiplies the current matrix by
                                           // the second matrix.
  void CreateMatrix(int rows, int cols);   //
  void MulNumber(const double num);
  void CopyData(const S21Matrix& other);
  S21Matrix CopyDataCompiment(const S21Matrix& other, int rows, int cols);
  void AlocateMem(int rows, int cols);
  bool MatrixMulConditions(const S21Matrix& other) const;
  void MutatorReSize(int rows, int cols);

  S21Matrix Transpose();  // Creates a new transposed matrix from the current
                          // one and returns it.
  S21Matrix CalcComplements();  // Calculates the algebraic addition matrix of
                                // the current one and returns it.
  double Determinant();         // Calculates and returns the determinant of the
                                // current matrix.
  S21Matrix InverseMatrix();    // Calculates and returns the inverse matrix.

  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  bool operator==(const S21Matrix& other) const;
  void operator=(const S21Matrix& other);
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  double operator()(int i, int j) const;
  double& operator()(int i, int j);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double num);
  void operator*=(const S21Matrix& other);
  void operator*=(const double num);
  void SwapCols(int place);
};
