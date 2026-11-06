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
Matrix<T>::operator*(const Matrix<T>& other) const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
    auto n1 = this->n;
    auto m1 = this->m;
    auto n2 = other.n;
    auto m2 = other.m;
    auto data1 = this->data;
    auto data2 = other.data;

    if (n1 != m2)
        return std::unexpected<std::string>("Improper dimensions for matrix multiplication");

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
        return std::unexpected<std::string>("Invalid index into the matrix");
    auto span = std::span(data).subspan(idx*n, n);
    return std::vector<T>(span.begin(), span.end());
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

template<Numeric T> std::expected<Matrix<T>, std::string>
reduction(Matrix<T> &matrix, bool augmented) {

    std::vector<std::vector<T>> rows;
    for (int i=0; i<matrix.m; i++)
        if (matrix[i].has_value())
            rows.push_back(matrix[i].value());

    auto _n = rows[0].size();
    auto _m = rows.size();
    for (int i=0; i<std::min(_m, _n); i++) {
        int m = i;
        for (int j=i+1; j<_m; j++) {
            if (std::abs(rows[j][i]) > std::abs(rows[i][i]))
                m = j;
        }
        auto tm = rows[m];
        rows[m] = rows[i];
        rows[i] = tm;

        T pivot = rows[i][i];
        if (pivot < 1/PRECISION_FACTOR)
            return std::unexpected<std::string>("Singular system (Infinte Soutions), pivot=0");

        if (augmented) {
            bool zero_coeff = true;
            for (int j=0; j<_n-1; j++)
                if (std::abs(rows[i][j]) > 1/PRECISION_FACTOR)
                    zero_coeff = false;
            if (zero_coeff && std::abs(rows[i][_n-1]) > 1/PRECISION_FACTOR) {
                return std::unexpected<std::string>("Augmented matrix is inconsistant");
            }
        }

        for (int j=0; j<_n; j++) {
           rows[i][j] /= pivot;
        }

        for (int j=i+1; j<_m; j++) {
            T factor = rows[j][i];
            for (int k=0; k<_n; k++)
                rows[j][k] += -factor*rows[i][k];
        }
    }

    std::vector<T> result;
    for (auto& vec:rows)
        for (auto& val:vec) {
            auto rounded = std::round(val*PRECISION_FACTOR)/PRECISION_FACTOR;
            result.push_back((rounded == 0) ? 0 : rounded);
        }


    return Matrix<T>(result, (int) matrix.m, (int) matrix.n);
}

#pragma endregion

int main(int argc, char** argv) {
    auto vec1 = Vector<double>({1, 0, 0}, 3);
    auto vec2 = Vector<double>({0, 1, 0}, 3);
    auto vec3 = Vector<double>({0, 0, 1}, 3);
    auto vec4 = Vector<double>({3, 5, 6}, 3);
    std::cout << "VEC 1" << " " << vec1 << std::endl;
    std::cout << "VEC 2" << " " << vec2 << std::endl;
    std::cout << "VEC 3" << " " << vec3 << std::endl;
    auto aug1 = vec2 | vec1 | vec3;
    auto mat1 = aug1.value();
    auto mat2 = (Matrix<double>({
            1, 2, 3,
            4, 5, 6,
            7, 8, 9,
        }, 3, 3) | vec4).value();
    std::cout << "MAT 1" << " " << mat1 << std::endl;
    std::cout << "MAT 2" << " " << mat2 << std::endl;
    // auto mul1 = mat1 * mat2 * vec1;
    // auto mat3 = mul1.value();
    // std::cout << "MAT 3" << " " << mat3 << std::endl;
    auto gaus = reduction(mat2, true);
    if (gaus.has_value())
        std::cout << "GAUSS" << " " << gaus.value() << std::endl;
    else
        std::cout << "ERROR" << " " << gaus.error() << std::endl;
    return 0;
}
