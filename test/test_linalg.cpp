#include "../include/math/linalg.hpp"

int main()
{
    Matrix<int> m1 = Matrix<int>(Row(5), Col(5));
    Matrix<int> m2 = m1;
    m2 = m2 + m1;
    std::cout << m2(1,1) << m2(5,5) << std::endl;
    return 0;
}
