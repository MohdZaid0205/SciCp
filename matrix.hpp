#include <concepts>
#include <iostream>
#include <ostream>
#include <vector>
#include <random>
#include <cctype>
#include <span>
#include <cmath>
#include <string>
#include <utility>
#include <expected>

#define PRECISION 3
#define PRECISION_FACTOR 1e3

#define COLOR_FOREGROUND(r, g, b) "\033[38;2;"#r";"#g";"#b"m"
#define COLOR_BACKGROUND(r, g, b) "\033[48;2;"#r";"#g";"#b"m"
#define RESET_FOREGROUND "\033[39m"
#define RESET_BACKGROUND "\033[49m"

#define PRI_BACK COLOR_BACKGROUND(  0, 163, 133)
#define SEC_BACK COLOR_BACKGROUND(215, 215, 168)
#define TER_BACK COLOR_BACKGROUND(  0,  73,  58)

#define PRI_FORE COLOR_FOREGROUND(  0, 163, 133)
#define SEC_FORE COLOR_FOREGROUND(215, 215, 168)
#define TER_FORE COLOR_FOREGROUND(  0,  73,  58)

#define ERR_FORE COLOR_FOREGROUND(220,  20,  60)
#define ERR_BACK COLOR_BACKGROUND(220,  20,  60)

// Concept for constraining precision of required matices and such with help of a generic
// class that can represent number using matissa and exponent without requiring float/int

template<typename T> concept Floating = std::same_as<T, float> || std::same_as<T, double>;
template<typename T> concept Integer  = std::same_as<T, int>   || std::same_as<T, long>   || std::same_as<T, long long>;

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
        static Matrix<T> identity(int s);
        Matrix<T> transpose() const;

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

        T&       operator()(int i, int j);
        const T& operator()(int i, int j) const;

        friend std::ostream& operator<< <T>(std::ostream& os, const Matrix<T>& vec);

    public:
        friend class Vector<T>;
};

template<typename T> using MatPair = std::pair<Matrix<T>, Matrix<T>>;

// Required Functions for implementation
template<Numeric T> std::expected<Vector<T>, std::string> gaussian(Matrix<T>& coeffs, Vector<T>& value);
template<Numeric T> std::expected<Matrix<T>, std::string> cholesky(const Matrix<T>& A);
template<Numeric T> std::expected<MatPair<T>, std::string>  LU(const Matrix<T>& A);         // {L, U}
template<Numeric T> std::expected<MatPair<T>, std::string>  GS(const Matrix<T>& A);         // {Q, R}
template<Numeric T> std::expected<MatPair<T>, std::string>  HH(const Matrix<T>& A);         // {Q, R}
template<Numeric T> std::expected<MatPair<T>, std::string>  GR(const Matrix<T>& A);         // {Q, R}

