#include "matrix.hpp"
#include <expected>

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
