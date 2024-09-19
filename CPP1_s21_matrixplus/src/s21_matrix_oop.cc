#include "s21_matrix_oop.h"

// geting rows from matrix
int S21Matrix::GetRows() const { return rows_; };
// geting cols from matrix
int S21Matrix::GetCols() const { return cols_; };
// editing rows from matrix
void S21Matrix::SetRows(int num) {
  if (num < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  if (this->GetRows() != 0) {
    int max_rows = 0;
    this->GetRows() > num ? max_rows = num : max_rows = this->GetRows();
    S21Matrix buf(num, this->GetCols());
    for (int i = 0; i < max_rows; i++) {
      for (int j = 0; j < this->GetCols(); j++) {
        buf(i, j) = (*this)(i, j);
      }
    }
    *this = buf;
  }
  rows_ = num;
};
// editing cols from matrix
void S21Matrix::SetCols(int num) {
  if (num < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  if (this->GetCols() != 0) {
    int max_cols = 0;
    this->GetCols() > num ? max_cols = num : max_cols = this->GetCols();
    S21Matrix buf(this->GetRows(), num);
    for (int i = 0; i < this->GetRows(); i++) {
      for (int j = 0; j < max_cols; j++) {
        buf(i, j) = (*this)(i, j);
      }
    }
    *this = buf;
  }
  cols_ = num;
};
// set **matrix_
void S21Matrix::EditPointer(double** matrix) { matrix_ = matrix; };
void S21Matrix::EditPointers(double* matrix, int i) { matrix_[i] = matrix; };
// return **matrix_
double** S21Matrix::GetPointer() const { return matrix_; };
//создание матриц
void S21Matrix::CreateMatrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0)
    throw std::out_of_range("Incorrect input, index is out of range");
  cols_ = cols;
  rows_ = rows;
  EditPointer(new double*[rows]);
  double** matrix = GetPointer();

  if (matrix == NULL) throw std::underflow_error("Out of memory");

  for (int i = 0; i < rows; i++) {
    EditPointers(new double[cols], i);
  }

  for (int i = 0; i < rows; i++) {
    if (matrix[i] == NULL) throw std::underflow_error("Out of memory");
  }
};
// A basic constructor that initialises a matrix of some predefined dimension.
S21Matrix::S21Matrix() { CreateMatrix(1, 1); };
//Параметризированный конструктор с количеством строк и столбцов.
S21Matrix::S21Matrix(int rows, int cols) { CreateMatrix(rows, cols); };
//конструктор присваивания
S21Matrix::S21Matrix(S21Matrix&& other) {
  this->AlocateMem(other.GetRows(), other.GetCols());
  this->CopyData(other);
  other.cols_ = 0;
  other.rows_ = 0;
  other.~S21Matrix();
  other.EditPointer(nullptr);
};
//Конструктор копирования.
S21Matrix::S21Matrix(const S21Matrix& other) {
  this->CreateMatrix(other.GetRows(), other.GetCols());
  this->CopyData(other);
};
// destructor
S21Matrix::~S21Matrix() {
  if (GetPointer()) {
    for (int i = 0; i < GetRows(); i++) {
      delete[] GetPointer()[i];
    }
    delete[] GetPointer();
  }
};
//принт матриц
// void S21Matrix::PrintMatrix() {
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       printf("%lf ", GetPointer()[i][j]);
//     }
//     std::cout << std::endl;
//   }
// };
//сравнение матриц если они равны возвращает 1
bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  bool result = 1;

  if (this->GetPointer() == NULL || other.GetPointer() == NULL)
    throw std::underflow_error("Out of memory");
  if (this->GetCols() != other.GetCols())
    throw std::out_of_range("cols dont equal");
  if (this->GetRows() != other.GetRows())
    throw std::out_of_range("rows dont equal");

  for (int i = 0; i < GetRows(); i++)
    for (int j = 0; j < GetCols(); j++)
      if (fabs((*this)(i, j) - other(i, j)) > 1e-7) result = 0;

  return result;
};
//сумма матриц одного размера
void S21Matrix::SumMatrix(const S21Matrix& other) {
  this->EqMatrix(other);
  this->EqMatrix(other);
  for (int i = 0; i < GetRows(); i++)
    for (int j = 0; j < GetCols(); j++)
      (*this)(i, j) = (*this)(i, j) + other(i, j);
};
//вычитание матриц
void S21Matrix::SubMatrix(const S21Matrix& other) {
  this->EqMatrix(other);
  for (int i = 0; i < GetRows(); i++)
    for (int j = 0; j < GetCols(); j++)
      (*this)(i, j) = (*this)(i, j) - other(i, j);
};
//умножение матриц
bool S21Matrix::MatrixMulConditions(const S21Matrix& other) const {
  int res = 1;
  if (this->GetRows() < 0 && this->GetCols() < 0) res = 0;
  if ((this->GetCols() != other.GetRows())) res = 0;
  if (!this->GetPointer() || !other.GetPointer()) res = 0;
  return res;
};

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->MatrixMulConditions(other)) {
    S21Matrix res(this->GetRows(), other.GetCols());
    for (int i = 0; i < this->GetRows(); i++)
      for (int j = 0; j < other.GetCols(); j++)
        for (int k = 0; k < this->GetCols(); k++)
          res(i, j) += (*this)(i, k) * other(k, j);
    *this = res;
  } else
    throw std::out_of_range("inappropriate matrix sizes for multiplication");
};
//умножение на число
void S21Matrix::MulNumber(const double num) {
  if (this->GetPointer() == NULL) throw std::underflow_error("Out of memory");
  if (GetCols() <= 0 || GetRows() <= 0)
    throw std::out_of_range("Incorrect input, index is out of range");

  for (int i = 0; i < GetRows(); i++)
    for (int j = 0; j < GetCols(); j++) (*this)(i, j) = (*this)(i, j) * num;
};

//копирование данных из одной матрицы в другую
void S21Matrix::CopyData(const S21Matrix& other) {
  this->EqMatrix(other);
  for (int i = 0; i < this->GetRows(); i++)
    for (int j = 0; j < this->GetCols(); j++) (*this)(i, j) = other(i, j);
};
//алокация памяти
void S21Matrix::AlocateMem(int rows, int cols) {
  this->CreateMatrix(rows, cols);
};
//транспонирование матриц
S21Matrix S21Matrix::Transpose() {
  S21Matrix matrix_t(GetCols(), GetRows());
  for (int i = 0; i < this->GetRows(); i++)
    for (int j = 0; j < this->GetCols(); j++) matrix_t(j, i) = (*this)(i, j);
  *this = matrix_t;
  return *this;
};
//оператор присваивания
void S21Matrix::operator=(const S21Matrix& other) {
  if (this->GetCols() != other.GetCols() ||
      this->GetRows() != other.GetRows()) {
    this->AlocateMem(other.GetRows(), other.GetCols());
  }
  this->CopyData(other);
};
bool S21Matrix::operator==(const S21Matrix& other) const {
  return this->EqMatrix(other);
};
S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
};
S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
};
double& S21Matrix::operator()(int i, int j) {
  if (i <= GetRows() && j <= GetCols())
    return this->GetPointer()[i][j];
  else
    throw std::out_of_range("Incorrect input, index is out of range");
};
double S21Matrix::operator()(int i, int j) const {
  if (i <= GetRows() && j <= GetCols())
    return this->GetPointer()[i][j];
  else
    throw std::out_of_range("Incorrect input, index is out of range");
};
void S21Matrix::operator+=(const S21Matrix& other) { *this = *this + other; };
void S21Matrix::operator-=(const S21Matrix& other) { *this = *this - other; };
S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
};
S21Matrix S21Matrix::operator*(const double num) {
  this->MulNumber(num);
  return *this;
};
void S21Matrix::operator*=(const S21Matrix& other) { *this = *this * other; };
void S21Matrix::operator*=(const double num) { *this = *this * num; };
void S21Matrix::SwapCols(int place) {
  for (int ii = place; ii < this->GetCols(); ii++) {
    if (!(*this)(place, place)) {
      for (int i = 0; i < GetRows(); i++) {
        double rem = (*this)(i, place);
        (*this)(i, place) = (*this)(i, ii);
        (*this)(i, ii) = rem;
      }
    } else
      break;
  }
};
double S21Matrix::Determinant() {
  if (this->GetCols() != this->GetRows()) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  int n = this->GetCols();
  double determinant = 1;
  S21Matrix matrix_det(*this);
  for (int k = 0; k < n; k++) {
    if (!matrix_det(k, k)) determinant *= -1;
    matrix_det.SwapCols(k);
    if (matrix_det(k, k)) {
      for (int i = k + 1; i < n; i++) {
        double K = matrix_det(i, k) / matrix_det(k, k);
        for (int j = 0; j < n; j++)
          matrix_det(i, j) = matrix_det(i, j) - matrix_det(k, j) * K;
      }
      determinant *= matrix_det(k, k);
    } else {
      determinant = 0;
      break;
    }
  }
  return determinant;
};
//вычеркивает из матрицы ряд и строчку
S21Matrix S21Matrix::CopyDataCompiment(const S21Matrix& other, int rows,
                                       int cols) {
  int ii = 0, jj = 0;
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      if (rows != i && cols != j) (*this)(jj, ii++) = other(i, j);
      if (ii == this->GetCols()) ii = 0, jj++;
    }
  }
  return *this;
};
//матрица алгебраических дополнений
S21Matrix S21Matrix::CalcComplements() {
  if (this->GetCols() != this->GetRows()) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  S21Matrix res(this->GetRows(), this->GetCols());
  S21Matrix buf(this->GetRows() - 1, this->GetCols() - 1);
  for (int i = 0; i < res.GetRows(); i++) {
    for (int j = 0; j < res.GetCols(); j++) {
      res(i, j) =
          buf.CopyDataCompiment(*this, i, j).Determinant() * pow(-1, i + j);
    }
  }
  *this = res;
  return res;
};
S21Matrix S21Matrix::InverseMatrix() {
  if (this->Determinant() == 0)
    throw std::out_of_range("Determinant = 0,error");
  double K = (1.0 / this->Determinant());
  this->CalcComplements();
  this->Transpose();
  *this *= K;
  return *this;
};
