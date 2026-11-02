#include "matrix.hpp"
#include <iomanip>

// VECTOR function implementations

template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator+(const Vector<T>& other) const {

    if (k != other.k)
        return std::unexpected<std::string>("Dimensions did not match");

    std::vector<T> add(k);
    for (int i=0; i<k; i++)
        add[i] = data[i] + other.data[i];

    return Vector(add, k);
}

template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator-(const Vector<T>& other) const {

    if (k != other.k)
        return std::unexpected<std::string>("Dimensions did not match");

    std::vector<T> sub(k);
    for (int i=0; i<k; i++)
        sub[i] = data[i] - other.data[i];

    return Vector(sub, k);
}

template<Numeric T> Vector<T>
Vector<T>::operator*(Numeric auto& other) const {

    std::vector<T> mul(k);
    for (int i=0; i<k; i++)
        mul[i] = other*data[i];

    return Vector(mul);
}

template<Numeric T> Vector<T>
Vector<T>::operator/(Numeric auto& other) const {

    std::vector<T> div(k);
    for (int i=0; i<k; i++)
        div[i] = other*data[i];

    return Vector(div);
}

template<Numeric T> std::expected<Matrix<T>, std::string>
Vector<T>::col_mat() const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
    return Matrix(data, k, 1);
}

template<Numeric T> std::expected<Matrix<T>, std::string>
Vector<T>::row_mat() const {
    // return std::unexpected<std::string>("NOT IMPLEMENTED");
    return Matrix(data, 1, k);
}

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

    return Vector(result, (int) k);
}

template<Numeric T> std::expected<Vector<T>, std::string>
Vector<T>::operator%(const Vector<T>& other) {
    return std::unexpected<std::string>("NOT IMPLEMENTED");
}

template<Numeric T> std::ostream&
operator<<(std::ostream& os, const Vector<T>& vec) {

    auto k    = vec.k;
    auto data = vec.data;

    os << "[";
    for (int i=0; i<k; i++)
        os  << " "
            << std::setprecision(PRINT_PRECISION)
            << data[i] << ""
            << (i==k-1 ? "": ",");
    os << "]";

    return os;
}

int main(int argc, char** argv) {
    return 0;
}
