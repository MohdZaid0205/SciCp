#include "matrix.hpp"
#include <expected>
#include <iomanip>

#pragma region VectorImplementaion

// Vector operator+
template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator+(const Vector<T>& other) const {

    if (k != other.k)
        return std::unexpected<std::string>("Dimensions did not match");

    std::vector<T> add(k);
    for (int i=0; i<k; i++)
        add[i] = data[i] + other.data[i];

    return Vector(add, k);
}

// Vector operator-
template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator-(const Vector<T>& other) const {

    if (k != other.k)
        return std::unexpected<std::string>("Dimensions did not match");

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

// Vector to.col_mat
// converts vector into matrix with sufficient dimensions
template<Numeric T> Matrix<T>
Vector<T>::col_mat() const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
    return Matrix(data, (int) k, 1);
}

// Vector to.row_mat
// converts vector into matrix with sufficient dimentions
template<Numeric T> Matrix<T>
Vector<T>::row_mat() const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
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
        std::unexpected<std::string>("Normalization not possible as norm=0");

    std::vector<T> result(k);
    for (int i=0; i<k; i++)
        result[i] = data[i]/norm;

    return Vector(result, k);
}

// Vector operator|
// returns a new matrix equivalent to this augment other
template<Numeric T> std::expected<Matrix<T>, std::string>
Vector<T>::operator|(const Vector<T>& other) const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
    auto mat1 = this->col_mat();
    auto mat2 = other.col_mat();
    return mat1 | mat2;
}

// Vector operator|
// returns new matrix equivalent to this augment other
template<Numeric T> std::expected<Matrix<T>, std::string>
Vector<T>::operator|(const Matrix<T>& other) const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
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

    os << "[";
    for (int i=0; i<k; i++)
        os  << " "
            << std::setprecision(PRECISION)
            << data[i] << ""
            << (i==k-1 ? "": ",");
    os << "]";

    return os;
}

#pragma endregion

#pragma region MatrixImplementation

// Matrix operator+
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator+(const Matrix<T>& other) const {

    if (m != other.m || n != other.n)
        return std::unexpected<std::string>("Cannot + on matrices of incompatible sizes");

    std::vector<T> result(m*n);
    for (int i=0; i<m*n; i++)
        result[i] = data[i] + other.data[i];

    return Matrix(result, m, n);
}

// Matrix operator-
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator-(const Matrix<T>& other) const {

    if (m != other.m || n != other.n)
        return std::unexpected<std::string>("Cannot + on matrices of incompatible sizes");

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
Matrix<T>::operator*(const Matrix<T>& other) {
    return std::unexpected<std::string>("NOT IMPLEMENTED");
}

// Matrix operator|
// Provides support for augmenting matrix with a column vector
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator|(const Vector<T>& other) const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
    auto mat1 = *this;
    auto mat2 = other.col_mat();
    return mat1 | mat2;
}

// Matrix operator|
// Provides support for augmenting matrix with another matrix
template<Numeric T> std::expected<Matrix<T>, std::string>
Matrix<T>::operator|(const Matrix<T>& other) const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
    auto k = other.n;
    auto s = n + k;

    if (m != other.m)
        return std::unexpected<std::string>("cannot augment different sized columns");

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
// simply displays matrix in ostream
template<Numeric T> std::ostream&
operator<<(std::ostream& os, const Matrix<T>& mat) {

    auto m    = mat.m;
    auto n    = mat.n;
    auto data = mat.data;

    os << "[" << std::endl;
    for (int i=0; i<m; i++) {
        os  << "\t[";
        for (int j=0; j<n; j++)
            os  << " "
                << std::setprecision(PRECISION)
                << data[i*n+j] << ""
                << (j==n-1 ? "": ",");
        os  << "]"
            << (i==m-1 ? "": ",")
            << std::endl;;
    }
    os << "]";
    return os;
}

#pragma endregion

#pragma region HelperImplementation

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

#pragma endregion


int main(int argc, char** argv) {
    auto vec1 = Vector<double>({1, 0, 0}, 3);
    auto vec2 = Vector<double>({0, 1, 0}, 3);
    auto vec3 = Vector<double>({0, 0, 1}, 3);
    std::cout << "VEC 1" << " " << vec1 << std::endl;
    std::cout << "VEC 2" << " " << vec2 << std::endl;
    std::cout << "VEC 3" << " " << vec3 << std::endl;
    auto aug1 = vec1 | vec2 | vec3;
    std::cout << "MAT 1" << " " << aug1.value() << std::endl;
    return 0;
}
