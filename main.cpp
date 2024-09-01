#include <iostream>

#include "polynomial.hpp"

int main() {
    Polynomial<float, 4> a{-1, -2, 3, -4, 15, -6};
    Polynomial<float, 3> b{1, 4, -5, 2, -12, 5};
    Polynomial<float, 0> nullPoly;

    std::cout << "a(x) = " << a << "\n";
    std::cout << "b(x) = " << b << "\n";

    std::cout << "a(x) * b(x) = " << a * b << "\n";
    std::cout << "a(x) / b(x) = " << (a / b) << "\n";
    std::cout << "a(x) % b(x) = " << a % b << "\n";
    std::cout << "a(x)b(x) % b(x) = " << (a*b) % b << "\n";
}