#include "matrix.hpp"
#include <expected>
#include <iomanip>

#pragma region VectorImplementaion

// Vector operator+
template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator+(const Vector<T>& other) const {

    if (k != other.k)
        return std::unexpected("Dimensions did not match");

    std::vector<T> add(k);
    for (int i=0; i<k; i++)
        add[i] = data[i] + other.data[i];

    return Vector(add, k);
}

// Vector operator-
template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator-(const Vector<T>& other) const {

    if (k != other.k)
        return std::unexpected("Dimensions did not match");

    std::vector<T> sub(k);
    for (int i=0; i<k; i++)
        sub[i] = data[i] - other.data[i];

    return Vector(sub, k);
}

// Vector operator*
template<Numeric T> Vector<T>
Vector<T>::operator*(Numeric auto& other) const {

    std::vector<T> mul(k);
    for (int i=0; i<k; i++)
        mul[i] = other*data[i];

    return Vector(mul);
}

// Vector operator/
template<Numeric T> Vector<T>
Vector<T>::operator/(Numeric auto& other) const {

    std::vector<T> div(k);
    for (int i=0; i<k; i++)
        div[i] = other*data[i];

    return Vector(div);
}

// Matrix identity
// returns an s x s identity matrix
template<Numeric T> Matrix<T>
Matrix<T>::identity(int s) {

    std::vector<T> result(s*s, T(0));
    for (int i=0; i<s; i++)
        result[i*s+i] = T(1);

    return Matrix(result, s, s);
}

// Matrix transpose
// returns a new matrix equivalent to this transposed
template<Numeric T> Matrix<T>
Matrix<T>::transpose() const {

    std::vector<T> result(m*n);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            result[j*m+i] = data[i*n+j];

    return Matrix(result, (int) n, (int) m);
}

// Matrix operator()
// provides mutable element access using (row, col)
template<Numeric T> T&
Matrix<T>::operator()(int i, int j) {
    return data[i*n+j];
}

// Matrix operator()
// provides read only element access using (row, col)
template<Numeric T> const T&
Matrix<T>::operator()(int i, int j) const {
    return data[i*n+j];
}

// Vector to.col_mat
// converts vector into matrix with sufficient dimensions
template<Numeric T> Matrix<T>
Vector<T>::col_mat() const {
    // return std::unexpected("NOT IMPLEMENTED");
    return Matrix(data, (int) k, 1);
}

// Vector to.row_mat
// converts vector into matrix with sufficient dimentions
template<Numeric T> Matrix<T>
Vector<T>::row_mat() const {
    // return std::unexpected("NOT IMPLEMENTED");
    return Matrix(data, 1, (int) k);
}

// Vector operator~
// returns a new vector equivalent to nomalized vector
template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator~() const {

    T norm = 0;
    for (int i=0; i<k; i++)
        norm += data[i]*data[i];
    norm = std::sqrt(norm);

    if (norm == 0)
        std::unexpected("Normalization not possible as norm=0");

    std::vector<T> result(k);
    for (int i=0; i<k; i++)
        result[i] = data[i]/norm;

    return Vector(result, k);
}

// Vector operator|
// returns a new matrix equivalent to this augment other
template<Numeric T> std::expected<Matrix<T>, std::string>
Vector<T>::operator|(const Vector<T>& other) const {
    // return std::unexpected("NOT IMPLEMENTED");
    auto mat1 = this->col_mat();
    auto mat2 = other.col_mat();
    return mat1 | mat2;
}

// Vector operator|
// returns new matrix equivalent to this augment other
template<Numeric T> std::expected<Matrix<T>, std::string>
Vector<T>::operator|(const Matrix<T>& other) const {
    // return std::unexpected("NOT IMPLEMENTED");
    auto mat1 = this->col_mat();
    auto mat2 = other;
    return mat1 | mat2;
}

// Vector operator|
// check wether two vectors are equal in value
template<Numeric T> bool
Vector<T>::operator==(const Vector<T>& other) const {
    if (k != other.k)
        return false;
    bool equal = true;
    for (int i=0; i<k; i++)
        equal = std::abs(data[i] - other.data[i]) < 1e-6;
    return equal;
}

// Vector operator<<
// simply displays the vector into terminal
template<Numeric T> std::ostream&
operator<<(std::ostream& os, const Vector<T>& vec) {

    auto k    = vec.k;
    auto data = vec.data;

    os << TER_FORE
       << "[";
    for (int i=0; i<k; i++)
        os  << SEC_BACK
            << " "
            << std::setprecision(PRECISION)
            << data[i] << ""
            << (i==k-1 ? " ": "," );
    os << RESET_BACKGROUND
       << "]"
       << RESET_FOREGROUND;

    return os;
}

#pragma endregion

#pragma region MatrixImplementation

// Matrix operator+
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator+(const Matrix<T>& other) const {

    if (m != other.m || n != other.n)
        return std::unexpected("Cannot + on matrices of incompatible sizes");

    std::vector<T> result(m*n);
    for (int i=0; i<m*n; i++)
        result[i] = data[i] + other.data[i];

    return Matrix(result, m, n);
}

// Matrix operator-
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator-(const Matrix<T>& other) const {

    if (m != other.m || n != other.n)
        return std::unexpected("Cannot + on matrices of incompatible sizes");

    std::vector<T> result(m*n);
    for (int i=0; i<m*n; i++)
        result[i] = data[i] - other.data[i];

    return Matrix(result, m, n);
}

// Matrix operator*
template<Numeric T> Matrix<T>
Matrix<T>::operator*(Numeric auto& other) const {

    std::vector<T> result(m*n);
    for (int i=0; i<m*n; i++)
        result[i] = data[i] * other;

    return Matrix(result, m, n);
}

// Matrix operator/
template<Numeric T> Matrix<T>
Matrix<T>::operator/(Numeric auto& other) const {

    std::vector<T> result(m*n);
    for (int i=0; i<m*n; i++)
        result[i] = data[i] / other;

    return Matrix(result, m, n);
}

// Matrix operator+
// Provides support for matrix multiplication if possible
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator*(const Matrix<T>& other) const {
    // return std::unexpected("NOT IMPLEMENTED");
    auto n1 = this->n;
    auto m1 = this->m;
    auto n2 = other.n;
    auto m2 = other.m;
    auto data1 = this->data;
    auto data2 = other.data;

    if (n1 != m2)
        return std::unexpected("Improper dimensions for matrix multiplication");

    std::vector<T> result(m1*n2, 0);
    for (int i=0; i<m1; i++)
        for (int j=0; j<n2; j++)
            for (int k=0; k<n1; k++)
                result[i*n2+j] += data1[i*n1+k]*data2[k*n2+j];

    return Matrix(result, (int) m1, (int) n2);
}

// Matrix operator+
// Provides support for matrix multiplication if possible
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator*(const Vector<T>& other) const {
    auto mat1 = *this;
    auto mat2 = other.col_mat();
    return mat1 * mat2;
}

// Matrix operator[]
// Provides point into matrix using indexing
template<Numeric T> std::expected<std::vector<T>, std::string>
Matrix<T>::operator[](Numeric auto idx) {
    if (idx < 0 || idx > m)
        return std::unexpected("Invalid index into the matrix");
    auto span = std::span(data).subspan(idx*n, n);
    return std::vector<T>(span.begin(), span.end());
}

// Matrix operator|
// Provides support for augmenting matrix with a column vector
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator|(const Vector<T>& other) const {
    // return std::unexpected("NOT IMPLEMENTED");
    auto mat1 = *this;
    auto mat2 = other.col_mat();
    return mat1 | mat2;
}

// Matrix operator|
// Provides support for augmenting matrix with another matrix
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator|(const Matrix<T>& other) const {
    // return std::unexpected("NOT IMPLEMENTED");
    auto k = other.n;
    auto s = n + k;

    if (m != other.m)
        return std::unexpected("cannot augment different sized columns");

    std::vector<T> result((n + other.n)*m);
    for (int i=0; i<m;i++) {
        for (int j=0; j<n; j++)
            result[i*s+j] = data[i*n+j];
        for (int j=0; j<k; j++)
            result[i*s+n+j] = other.data[i*k+j];
    }
    return Matrix(result, (int) m, (int) s);
}

// Matrix operator|
// check wether two matricess are equal in value
template<Numeric T> bool
Matrix<T>::operator==(const Matrix<T>& other) const {
    if (n != other.n || m != other.m)
        return false;
    bool equal = true;
    for (int i=0; i<m*n; i++)
        equal = std::abs(data[i] - other.data[i]) < 1e-6;
    return equal;
}

// Matrix operator<<
// simply displays matrix in ostream, padded so columns line up
template<Numeric T> std::ostream&
operator<<(std::ostream& os, const Matrix<T>& mat) {

    auto m    = mat.m;
    auto n    = mat.n;
    auto data = mat.data;

    std::vector<std::string> cells(m*n);
    size_t width = 0;
    for (int i=0; i<m*n; i++) {
        T value = data[i];

        if (std::abs(value) < 1/PRECISION_FACTOR)
            value = 0;

        std::ostringstream ss;
        ss << std::setprecision(PRECISION) << value;
        cells[i] = ss.str();
        width = std::max(width, cells[i].size());
    }

    os << TER_FORE
       << "["
       << RESET_FOREGROUND
       << std::endl;

    for (int i=0; i<m; i++) {
        os  << "\t"
            << TER_FORE
            << "[";
        for (int j=0; j<n; j++)
            os  << SEC_BACK
                << TER_FORE
                << " "
                << std::right
                << std::setw(width)
                << cells[i*n+j]
                << (j==n-1 ? " ": ",");
        os  << TER_FORE
            << RESET_BACKGROUND
            << "]"
            << (i==m-1 ? " ": ",")
            << RESET_FOREGROUND
            << std::endl;
    }
    os << TER_FORE
       << "]"
       << RESET_FOREGROUND;
    return os;
}

#pragma endregion

#pragma region HelperImplementation (QOL)

// Matrix operator|
// provides support for error and result propogation during augmentation
template<Numeric T> std::expected<Matrix<T>, std::string> 
operator|(const std::expected<Matrix<T>, std::string>& lhs, const Vector<T>& rhs) {
    if (!lhs.has_value())
        return lhs;
    return lhs.value() | rhs;
}

// Matrix operator|
// provides support for error and result propogation during augmentation
template<Numeric T> std::expected<Matrix<T>, std::string> 
operator|(const std::expected<Matrix<T>, std::string>& lhs, const Matrix<T>& rhs) {
    if (!lhs.has_value())
        return lhs;
    return lhs.value() | rhs;
}

// Matrix operator|
// provides support for error and result propogation during augmentation
template<Numeric T> std::expected<Matrix<T>, std::string>
operator|(const Vector<T>& lhs, const std::expected<Matrix<T>, std::string>& rhs) {
    if (!rhs.has_value())
        return rhs;
    return lhs | rhs.value();
}

// Matrix operator|
// provides support for error and result propogation during augmentation
template<Numeric T> std::expected<Matrix<T>, std::string>
operator|(const Matrix<T>& lhs, const std::expected<Matrix<T>, std::string>& rhs) {
    if (!rhs.has_value())
        return rhs;
    return lhs | rhs.value();
}

// Matrix operator*
// provides support fot error and result propogation during multiplication
template<Numeric T> std::expected<Matrix<T>, std::string>
operator*(const std::expected<Matrix<T>, std::string>& lhs, const Vector<T>& rhs) {
    if (!lhs.has_value())
        return lhs;
    return lhs.value() * rhs.col_mat();
}

// Matrix operator*
// provides support fot error and result propogation during multiplication
template<Numeric T> std::expected<Matrix<T>, std::string>
operator*(const std::expected<Matrix<T>, std::string>& lhs, const Matrix<T>& rhs) {
    if (!lhs.has_value())
        return lhs;
    return lhs.value() * rhs;
}

// Matrix operator|
// provides support for error and result propogation during augmentation
template<Numeric T> std::expected<Matrix<T>, std::string>
operator*(const Vector<T>& lhs, const std::expected<Matrix<T>, std::string>& rhs) {
    if (!rhs.has_value())
        return rhs;
    return lhs.col_mat() * rhs.value();
}

// Matrix operator|
// provides support for error and result propogation during augmentation
template<Numeric T> std::expected<Matrix<T>, std::string>
operator*(const Matrix<T>& lhs, const std::expected<Matrix<T>, std::string>& rhs) {
    if (!rhs.has_value())
        return rhs;
    return lhs * rhs.value();
}

// Matrix operator[]
// provides support for chaining and refernce into a expected span
// template<Numeric T> std::expected<T&, std::string>
// operator[](std::expected<std::span<T>, std::string> indexable, Numeric auto idx) {
//     if (!indexable.has_value())
//        return indexable;
//    return &indexable.value()[idx];
// }

#pragma endregion

#pragma region AlgorithmsImplementation

template<Numeric T> std::expected<Vector<T>, std::string>
gaussian(Matrix<T>& coeffs, Vector<T>& value) {

    auto augment = coeffs | value;
    if (!augment.has_value())
        return std::unexpected<std::string>("Cannot apply gaussian to inconsistant coeff and value");
    auto matrix = augment.value();

    std::vector<std::vector<T>> rows;
    for (int i=0; i<matrix.m; i++)
        if (matrix[i].has_value())
            rows.push_back(matrix[i].value());

    int _m = (int) rows.size();
    int _n = (int) rows[0].size();
    if (_m != _n - 1)
        return std::unexpected<std::string>("Gaussian elimination needs a square coefficient matrix");

    for (int i=0; i<_m; i++) {
        int p = i;
        for (int j=i+1; j<_m; j++)
            if (std::abs(rows[j][i]) > std::abs(rows[p][i]))
                p = j;
        std::swap(rows[p], rows[i]);

        T pivot = rows[i][i];
        if (std::abs(pivot) < 1/PRECISION_FACTOR)
            return std::unexpected<std::string>("Singular system, pivot=0");

        for (int j=i; j<_n; j++)
            rows[i][j] /= pivot;

        for (int j=i+1; j<_m; j++) {
            T factor = rows[j][i];
            for (int k=i; k<_n; k++)
                rows[j][k] -= factor * rows[i][k];
        }
    }

    std::vector<T> x(_m);
    for (int i=_m-1; i>=0; i--) {
        T sum = rows[i][_n-1];
        for (int j=i+1; j<_m; j++)
            sum -= rows[i][j] * x[j];
        x[i] = sum;
    }
    return Vector<T>(x, (int) _m);
}

template<Numeric T> std::expected<Matrix<T>, std::string>
cholesky(const Matrix<T>& A) {
    if (A.m != A.n)
        return std::unexpected("Cholesky requires a square matrix");

    int n = (int) A.n;
    auto L = Matrix<T>(std::vector<T>(n*n, T(0)), n, n);

    for (int i=0; i<n; i++) {
        for (int j=0; j<=i; j++) {
            T sum = A(i, j);
            for (int k=0; k<j; k++)
                sum -= L(i, k) * L(j, k);

            if (i == j) {
                if (sum <= 1/PRECISION_FACTOR)
                    return std::unexpected("Matrix is not positive definite");
                L(i, i) = std::sqrt(sum);
            } else {
                L(i, j) = sum / L(j, j);
            }
        }
    }
    return L;
}

template<Numeric T> std::expected<MatPair<T>, std::string>
LU(const Matrix<T>& A) {
    if (A.m != A.n)
        return std::unexpected("LU requires a square matrix");

    int n = (int) A.n;
    auto L = Matrix<T>::identity(n);
    auto U = A;

    for (int k=0; k<n; k++) {
        if (std::abs(U(k, k)) < 1/PRECISION_FACTOR)
            return std::unexpected("Zero pivot encountered, LU without pivoting failed");

        for (int i=k+1; i<n; i++) {
            T factor = U(i, k) / U(k, k);
            L(i, k) = factor;
            for (int j=k; j<n; j++)
                U(i, j) -= factor * U(k, j);
        }
    }
    return MatPair<T>(L, U);
}

template<Numeric T> std::expected<MatPair<T>, std::string>
GS(const Matrix<T>& A) {
    int m = (int) A.m, n = (int) A.n;
    if (m < n)
        return std::unexpected("Gram-Schmidt requires rows >= cols");

    auto Q = Matrix<T>(std::vector<T>(m*n, T(0)), m, n);
    auto R = Matrix<T>(std::vector<T>(n*n, T(0)), n, n);

    for (int j=0; j<n; j++) {
        std::vector<T> v(m);
        for (int r=0; r<m; r++) v[r] = A(r, j);

        for (int i=0; i<j; i++) {
            T dot = 0;
            for (int r=0; r<m; r++) dot += Q(r, i) * v[r];
            R(i, j) = dot;
            for (int r=0; r<m; r++) v[r] -= dot * Q(r, i);
        }

        T norm = 0;
        for (int r=0; r<m; r++) norm += v[r] * v[r];
        norm = std::sqrt(norm);

        if (norm < 1/PRECISION_FACTOR)
            return std::unexpected("Columns are linearly dependent");

        R(j, j) = norm;
        for (int r=0; r<m; r++) Q(r, j) = v[r] / norm;
    }
    return MatPair<T>(Q, R);
}

template<Numeric T> std::expected<MatPair<T>, std::string>
HH(const Matrix<T>& A) {
    int m = (int) A.m, n = (int) A.n;
    auto R = A;
    auto Q = Matrix<T>::identity(m);

    for (int k=0; k<std::min(m-1, n); k++) {
        T norm = 0;
        for (int i=k; i<m; i++) norm += R(i, k) * R(i, k);
        norm = std::sqrt(norm);
        if (norm < 1/PRECISION_FACTOR) continue;

        T alpha = (R(k, k) >= 0) ? -norm : norm;

        std::vector<T> v(m-k);
        for (int i=k; i<m; i++) v[i-k] = R(i, k);
        v[0] -= alpha;

        T vv = 0;
        for (auto x : v) vv += x * x;
        if (vv < 1/PRECISION_FACTOR) continue;

        // R = H R
        for (int j=0; j<n; j++) {
            T s = 0;
            for (int i=k; i<m; i++) s += v[i-k] * R(i, j);
            for (int i=k; i<m; i++) R(i, j) -= 2 * s * v[i-k] / vv;
        }
        // Q = Q H
        for (int i=0; i<m; i++) {
            T s = 0;
            for (int l=k; l<m; l++) s += Q(i, l) * v[l-k];
            for (int l=k; l<m; l++) Q(i, l) -= 2 * s * v[l-k] / vv;
        }
    }

    for (int j=0; j<n; j++)
        for (int i=j+1; i<m; i++)
            R(i, j) = 0;

    return MatPair<T>(Q, R);
}

template<Numeric T> std::expected<MatPair<T>, std::string>
GR(const Matrix<T>& A) {
    int m = (int) A.m, n = (int) A.n;
    auto R = A;
    auto Q = Matrix<T>::identity(m);

    for (int j=0; j<std::min(m-1, n); j++) {
        for (int i=m-1; i>j; i--) {
            T a = R(i-1, j);
            T b = R(i,   j);
            if (std::abs(b) < 1/PRECISION_FACTOR) continue;

            T r = std::hypot(a, b);
            T c = a / r;
            T s = b / r;

            // R = G R
            for (int l=0; l<n; l++) {
                T t1 = R(i-1, l), t2 = R(i, l);
                R(i-1, l) =  c*t1 + s*t2;
                R(i,   l) = -s*t1 + c*t2;
            }
            // Q = Q G^T
            for (int l=0; l<m; l++) {
                T q1 = Q(l, i-1), q2 = Q(l, i);
                Q(l, i-1) =  c*q1 + s*q2;
                Q(l, i)   = -s*q1 + c*q2;
            }
            R(i, j) = 0;
        }
    }
    return MatPair<T>(Q, R);
}

#pragma endregion

#define MENU_LISTING(i, name, description)          \
    std::cout << "\t"                               \
              << SEC_BACK                           \
              << TER_FORE                           \
              << "[" << i << "]"                    \
              << RESET_FOREGROUND                   \
              << RESET_BACKGROUND                   \
              << " : "                              \
              << PRI_FORE                           \
              << name                               \
              << RESET_FOREGROUND                   \
              << " "                                \
              << TER_FORE                           \
              << description                        \
              << RESET_FOREGROUND                   \
              << std::endl

#define PRINT_ERROR(msg)                            \
    std::cout << ERR_BACK                           \
              << " ERROR "                          \
              << RESET_BACKGROUND                   \
              << " "                                \
              << ERR_FORE                           \
              << (msg)                              \
              << RESET_FOREGROUND                   \
              << std::endl

#define PRINT_MATRIX(label, mat)                    \
    std::cout << PRI_FORE                           \
              << (label)                            \
              << RESET_FOREGROUND                   \
              << " = "                              \
              << (mat)                              \
              << std::endl

#define PRINT_PAIR(a, b, res)                       \
    do {                                            \
        auto&& _pair = (res);                       \
        if (!_pair.has_value()) {                   \
            PRINT_ERROR(_pair.error());             \
        } else {                                    \
            PRINT_MATRIX((a), _pair.value().first); \
            PRINT_MATRIX((b), _pair.value().second);\
        }                                           \
    } while (0)

int main(int argc, char** argv) {

    long m = 0, n = 0;
    std::cout << SEC_FORE << "M" << RESET_FOREGROUND << " : " << PRI_FORE;
    std::cin >> m;
    std::cout << SEC_FORE << "N" << RESET_FOREGROUND << " : " << PRI_FORE;
    std::cin >> n;
    std::cout << RESET_FOREGROUND;

    if (!std::cin || m <= 0 || n <= 0) {
        PRINT_ERROR("M and N must be positive integers");
        return 1;
    }

    std::cout << TER_FORE
              << "RECIVED: "
              << PRI_FORE << "M=" << SEC_FORE << m << RESET_FOREGROUND
              << " "
              << PRI_FORE << "N=" << SEC_FORE << n << RESET_FOREGROUND
              << std::endl;
    
    char random_choice = 'n';
    std::cout << SEC_FORE << "RANDOM (Y/N)" << RESET_FOREGROUND << " : " << PRI_FORE;
    std::cin >> random_choice;
    std::cout << RESET_FOREGROUND;

    random_choice = (char) std::tolower((unsigned char) random_choice);

    if (!std::cin || (random_choice != 'y' && random_choice != 'n')) {
        PRINT_ERROR("Invalid input, expected Y or N");
        return 1;
    }

    auto data = std::vector<float>(m*n);

    if (random_choice == 'y') {

        // random integers in [-9, 9], stored as floats so they stay easy to read
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(-9, 9);

        for (int i=0; i<m*n; i++)
            data[i] = (float) dist(rng);

    } else {

        for (int i=0; i<m; i++) {

            std::cout << TER_FORE
                      << "ROW: "
                      << PRI_FORE << "i=" << SEC_FORE << i + 1
                      << RESET_FOREGROUND
                      << std::endl;

            for (int j=0; j<n; j++) {

                std::cout << TER_FORE
                          << "\tCOL:"
                          << PRI_FORE << "j=" << SEC_FORE << j + 1
                          << RESET_FOREGROUND
                          << ": "
                          << SEC_FORE;

                std::cin >> data[i*n+j];
                std::cout << RESET_FOREGROUND;

                if (!std::cin) {
                    std::cout << std::endl;
                    PRINT_ERROR("Invalid matrix entry, expected a number");
                    return 1;
                }
            }
        }
    }
    auto matrix = Matrix(data, m, n);

    std::cout << TER_FORE
              << "RECIVED: "
              << matrix
              << std::endl;

    std::cout << TER_FORE
              << "MENU:"
              << RESET_FOREGROUND
              << std::endl;

    MENU_LISTING(0, "Gaussian", "       Solve Ax=b (asks for b)");
    MENU_LISTING(1, "Cholesky", "       Transform A to L (A=L*L^T)");
    MENU_LISTING(2, "LU", "             Transform A to LU (A=LU)");
    MENU_LISTING(3, "Gram-Schmid", "    Transform A to QR (A=QR)");
    MENU_LISTING(4, "House Holder", "   Transform A to QR (A=QR)");
    MENU_LISTING(5, "Givens Rotation", "Transform A to QR (A=QR)");

    int choice = -1;
    std::cout << SEC_FORE << "CHOICE" << RESET_FOREGROUND << " : " << PRI_FORE;
    std::cin >> choice;
    std::cout << RESET_FOREGROUND;

    if (!std::cin) {
        PRINT_ERROR("Invalid input, expected a number from the menu");
        return 1;
    }

    switch (choice) {

        case 0: {
            auto bdata = std::vector<float>(m);
            std::cout << TER_FORE << "VECTOR b (" << m << " entries)"
                      << RESET_FOREGROUND << std::endl;

            for (int i=0; i<m; i++) {
                std::cout << TER_FORE
                          << "\tb"
                          << PRI_FORE << "[" << SEC_FORE << i + 1 << PRI_FORE << "]"
                          << RESET_FOREGROUND
                          << ": "
                          << SEC_FORE;

                std::cin >> bdata[i];
                std::cout << RESET_FOREGROUND;

                if (!std::cin) {
                    std::cout << std::endl;
                    PRINT_ERROR("Invalid vector entry, expected a number");
                    return 1;
                }
            }
            auto b = Vector<float>(bdata, m);

            std::cout << TER_FORE << "RECIVED: " << b << std::endl;

            auto res = gaussian(matrix, b);
            if (res.has_value()) {
                PRINT_MATRIX("x", res.value());
            } else {
                PRINT_ERROR(res.error());
            }
            break;
        }

        case 1: {
            auto res = cholesky(matrix);
            if (res.has_value()) {
                PRINT_MATRIX("L", res.value());
            } else {
                PRINT_ERROR(res.error());
            }
            break;
        }

        case 2:
            PRINT_PAIR("L", "U", LU(matrix));
            break;

        case 3:
            PRINT_PAIR("Q", "R", GS(matrix));
            break;

        case 4:
            PRINT_PAIR("Q", "R", HH(matrix));
            break;

        case 5:
            PRINT_PAIR("Q", "R", GR(matrix));
            break;

        default:
            PRINT_ERROR("Invalid choice, pick a number between 0 and 5");
            return 1;
    }

    return 0;
}
