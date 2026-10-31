#include <concepts>
#include <iostream>
#include <vector>
#include <expected>

// Concept for constraining precision of required matices and such with help of a generic
// class that can represent number using matissa and exponent without requiring float/int

template<typename T> concept Floating = std::same_as<T, float> || std::same_as<T, double>;
template<typename T> concept Intiger  = std::same_as<T, int>   || std::same_as<T, long>;
template<typename T> concept Numeric  = Floating<T> || Intiger<T>;

template<Numeric T> class Vector;
template<Numeric T> class Matrix;

template<typename T> concept Arithematic = requires(T a, T b) {
    {a + b} -> std::same_as<T>;     // a + b -> T
    {a - b} -> std::same_as<T>;     // a - b -> T
    {a * b} -> std::same_as<T>;     // a * b -> T
    {a / b} -> std::convertible_to<T>;     // a / b -> can be converted to T
};

template<typename T> concept PartialOrd  = requires(T a, T b) {
    {a == b} -> std::same_as<bool>;
    {a  > b} -> std::same_as<bool>;
    {a  < b} -> std::same_as<bool>;
};


template<Numeric T>
class Vector {
    private:
        std::vector<T> data;
        long long int k;

    public:
        Vector(std::vector<T> _d, Intiger auto _k): data(_d), k(_k) {}

        std::expected<Vector<T>, std::string> operator+(const Vector<T>& other) const;
        std::expected<Vector<T>, std::string> operator-(const Vector<T>& other) const;
        Vector<T> operator*(Numeric auto& other) const;
        Vector<T> operator/(Numeric auto& other) const;

        friend std::expected<Matrix<T>, std::string> operator*(const Vector<T>& other);
        friend std::expected<Matrix<T>, std::string> operator*(const Matrix<T>& other);

        Vector<T> operator~() const;    // normalization
        Vector<T> operator!() const;    // transpose

        std::expected<Vector<T>, std::string> operator%(const Vector<T>& other);    // projection
};


template<Numeric T>
class Matrix {
    private:
        std::vector<T> data;
        long long n;
        long long m;

    public:
        Matrix(std::vector<T> _d, Intiger auto _n, Intiger auto _m): data(_d), n(_n), m(_m) {}

        Matrix<T> operator+(const Matrix<T>& other) const;
        Matrix<T> operator-(const Matrix<T>& other) const;
        Matrix<T> operator*(Numeric auto& other) const;
        Matrix<T> operator/(Numeric auto& other) const;

        std::expected<Matrix<T>, std::string> operator*(const Matrix<T>& other) const;

        static Matrix<T> gaussian(const Matrix<T>& matrix);
        static Matrix<T> cholesky(const Matrix<T>& matrix);
        static Matrix<T> LU(const Matrix<T>& matrix);
        static Matrix<T> GS(const Matrix<T>& matrix);
        static Matrix<T> GR(const Matrix<T>& matrix);
        static Matrix<T> HH(const Matrix<T>& matrix);
};
