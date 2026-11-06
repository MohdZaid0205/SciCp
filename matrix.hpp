#include <concepts>
#include <iostream>
#include <ostream>
#include <vector>
#include <span>
#include <cmath>
#include <expected>

#define PRECISION 6
#define PRECISION_FACTOR 1e6

// Concept for constraining precision of required matices and such with help of a generic
// class that can represent number using matissa and exponent without requiring float/int

template<typename T> concept Floating = std::same_as<T, float> || std::same_as<T, double>;
template<typename T> concept Integer  = std::same_as<T, int>   || std::same_as<T, long>;

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

template<typename T> concept Numeric = Floating<T> || Integer<T> || (Arithematic<T> && PartialOrd<T>);

template<Numeric T> class Vector;
template<Numeric T> class Matrix;


template<Numeric T>
class Vector {
    private:
        std::vector<T> data;

    public:
        long long int k;

    public:
        Vector(std::vector<T> _d, Integer auto _k): data(_d), k(_k) {}

        std::expected<Vector<T>, std::string> operator+(const Vector<T>& other) const;
        std::expected<Vector<T>, std::string> operator-(const Vector<T>& other) const;
        Vector<T> operator*(Numeric auto& other) const;
        Vector<T> operator/(Numeric auto& other) const;
        std::expected<Vector<T>, std::string> operator~() const;  // normalization

        Matrix<T> col_mat() const;
        Matrix<T> row_mat() const;
        std::expected<Matrix<T>, std::string> operator|(const Matrix<T>& other) const;
        std::expected<Matrix<T>, std::string> operator|(const Vector<T>& other) const;
        bool operator==(const Vector<T>& other) const;

        friend std::ostream& operator<< <T>(std::ostream& os, const Vector<T>& vec);

    public:
        friend class Matrix<T>;
};


template<Numeric T>
class Matrix {
    private:
        std::vector<T> data;

    public:
        long long m; // rows
        long long n; // cols

    public:
        Matrix(std::vector<T> _d, Integer auto _m, Integer auto _n): data(_d), m(_m), n(_n) {}

        std::expected<Matrix<T>, std::string> operator+(const Matrix<T>& other) const;
        std::expected<Matrix<T>, std::string> operator-(const Matrix<T>& other) const;
        Matrix<T> operator*(Numeric auto& other) const;
        Matrix<T> operator/(Numeric auto& other) const;

        std::expected<Matrix<T>, std::string> operator*(const Matrix<T>& other) const;
        std::expected<Matrix<T>, std::string> operator*(const Vector<T>& other) const;

        std::expected<std::vector<T>, std::string> operator[](Numeric auto idx);
        std::expected<Matrix<T>, std::string> operator|(const Vector<T>& other) const;
        std::expected<Matrix<T>, std::string> operator|(const Matrix<T>& other) const;
        bool operator==(const Matrix<T>& other) const;

        friend std::ostream& operator<< <T>(std::ostream& os, const Matrix<T>& vec);

    public:
        friend class Vector<T>;
};


// Required Functions for implementation
template<Numeric T> std::expected<Matrix<T>, std::string> reduction(Matrix<T>& matrix);
template<Numeric T> Matrix<T> gaussian(Matrix<T>& matrix);
template<Numeric T> Matrix<T> cholesky(Matrix<T>& matrix);
template<Numeric T> Matrix<T> LU(Matrix<T>& matrix);
template<Numeric T> Matrix<T> GS(Matrix<T>& matrix);
template<Numeric T> Matrix<T> GR(Matrix<T>& matrix);
template<Numeric T> Matrix<T> HH(Matrix<T>& matrix);

