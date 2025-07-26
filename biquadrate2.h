
#ifndef INC_5_TIMES_EQUATION_BIQUADRATE2_H
#define INC_5_TIMES_EQUATION_BIQUADRATE2_H

#include "./typesNeed.h"
#include "biquadrate3.h"

vector<double> solve_quintic_(complex<double> a, complex<double> b, complex<double> c,
                              complex<double> d, complex<double> e, complex<double> f_coeff);
std::vector<vec5_> solveEquation_5_(coefficient6_& src);
bool checkSolution5(coefficient6 src, vec5 res);
bool checkSolution5_(coefficient6_ src, vec5_ res);

double get_real(std::complex<double> x);


vector<double> solve_quintic_(complex<double> a, complex<double> b, complex<double> c,
                              complex<double> d, complex<double> e, complex<double> f_coeff){
    vector<double> res = solve_quintic(get_real(a), get_real(b), get_real(c),
                                       get_real(d), get_real(e), get_real(f_coeff));
    return res;
}

std::vector<vec5_> solveEquation_5_(coefficient6_& src){
    assert(abs(src.f - 1.0) < 1e-12); // 保证系数不为1
    // ag⁴ - bg³ + cg² - dg + e - 1= 0
    vector<double> res = solve_quintic_(src.a, -src.b, src.c, -src.d, src.e ,- 1.0);

//    if(std::isnan(res.x1) && std::isnan(res.x2) && std::isnan(res.x3) && std::isnan(res.x4)){
//        std::cerr << "Error: Can't find solution" << std::endl;
//        std::vector<vec5> empty;
//        return empty;
//    }

    std::complex<double> g;
    std::vector<vec5_> res5;
    for(int i=0; i<res.size(); i++){
//        if( !isnan(get_real(res[i])) ){ // 找到第一个非nan解
            g = res[i];
            std::complex<double> a = src.a*g;
            std::complex<double> b = src.b*g - src.a*g*g;
            std::complex<double> c = src.c*g - src.b*g*g + src.a*g*g*g;
            std::complex<double> d = src.d*g - src.c*g*g + src.b*g*g*g - src.a*g*g*g*g;
            std::cout << "(" << a << "x⁴ + " << b << "x³ + " << c << "x² + " << d << "x + 1) * (x + " << g << ") = 0" << std::endl;
        std::complex<double> A = a/g, B = a+b/g, C = b+c/g, D = c+d/g, E=d+1.0/g, F=1.0;
            std::cout <<  A << "x⁵ + " << B << "x⁴ + " << C << "x³ + " << D << "x² + " << E << "x + " << F << " = 0" << std::endl;
            std::complex<double> beforeF_ = src.beforeF;
//            std::cout << "Before F: " << beforeF_ << std::endl;
            std::cout <<  A*beforeF_ << "x⁵ + " << B*beforeF_ << "x⁴ + " << C*beforeF_ <<
            "x³ + " << D*beforeF_ << "x² + " << E*beforeF_ << "x + " << F*beforeF_ << " = 0\n" << std::endl;
            // (ax⁴ + bx³ + cx² + dx + 1)(x + g) = 0
            vec4_ front4 = solveEquation_4_(coefficient5_{a, b, c, d, 1.0});
        std::complex<double> x5 = -g;
            res5.push_back(vec5_{front4.x1, front4.x2, front4.x3, front4.x4, x5});
//        }
    }

    return res5;
}
bool checkSolution5(coefficient6 src, vec5 res){
    return checkSolution5_(coefficient6_{src.a, src.b, src.c, src.d, src.e, src.f},
                           vec5_{res.x1, res.x2, res.x3, res.x4, res.x5});
}
bool checkSolution5_(coefficient6_ src, vec5_ res){
    std::complex<double> check1 = src.a * res.x1 * res.x1 * res.x1 * res.x1 * res.x1 +
                    src.b * res.x1 * res.x1 * res.x1 * res.x1 +
                    src.c * res.x1 * res.x1 * res.x1 +
                    src.d * res.x1 * res.x1 +
                    src.e * res.x1 +
                    src.f;
    double tolerance = 1e-10;
    bool correct1 = abs(check1) < tolerance;
    if(correct1){
        std::cout << "Solution1: " << res.x1 << " → " << check1 <<
        " is correct" << std::endl;
    }else{
        std::cout << "Solution1: " << res.x1 << " → " << check1 <<
        " is incorrect!!!" << std::endl;
    }
    std::complex<double> check2 = src.a * res.x2 * res.x2 * res.x2 * res.x2 * res.x2 +
                    src.b * res.x2 * res.x2 * res.x2 * res.x2 +
                    src.c * res.x2 * res.x2 * res.x2 +
                    src.d * res.x2 * res.x2 +
                    src.e * res.x2 +
                    src.f;
    bool correct2 = abs(check2) < tolerance;
    if(correct2){
        std::cout << "Solution2: " << res.x2 << " → " << check2 <<
        " is correct" << std::endl;
    }else{
        std::cout << "Solution2: " << res.x2 << " → " << check2 <<
        "is incorrect!!!" << std::endl;
    }
    std::complex<double> check3 = src.a * res.x3 * res.x3 * res.x3 * res.x3 * res.x3 +
                    src.b * res.x3 * res.x3 * res.x3 * res.x3 +
                    src.c * res.x3 * res.x3 * res.x3 +
                    src.d * res.x3 * res.x3 +
                    src.e * res.x3 +
                    src.f;
    bool correct3 = abs(check3) < tolerance;
    if(correct3){
        std::cout << "Solution3: " << res.x3 << " → " << check3 <<
        " is correct" << std::endl;
    }else{
        std::cout << "Solution3: " << res.x3 << " → " << check3 <<
        " is incorrect!!!" << std::endl;
    }
    std::complex<double> check4 = src.a * res.x4 * res.x4 * res.x4 * res.x4 * res.x4 +
                    src.b * res.x4 * res.x4 * res.x4 * res.x4 +
                    src.c * res.x4 * res.x4 * res.x4 +
                    src.d * res.x4 * res.x4 +
                    src.e * res.x4 +
                    src.f;
    bool correct4 = abs(check4) < tolerance;
    if(correct4){
        std::cout << "Solution4: " << res.x4 << " → " << check4 <<
        " is correct" << std::endl;
    }else{
        std::cout << "Solution4: " << res.x4 << " → " << check4 <<
        " is incorrect!!!" << std::endl;
    }
    std::complex<double> check5 = src.a * res.x5 * res.x5 * res.x5 * res.x5 * res.x5 +
                    src.b * res.x5 * res.x5 * res.x5 * res.x5 +
                    src.c * res.x5 * res.x5 * res.x5 +
                    src.d * res.x5 * res.x5 +
                    src.e * res.x5 +
                    src.f;
    bool correct5 = abs(check5) < tolerance;
    if(correct5){
        std::cout << "Solution5: " << res.x5 << " → " << check5 <<
        " is correct" << std::endl;
    }else{
        std::cout << "Solution5: " << res.x5 << " → " << check5 <<
        " is incorrect!!!\n" << std::endl;
    }
    return correct1 && correct2 && correct3 && correct4 && correct5;
}
double get_real(std::complex<double> x){
    return (std::abs(x.imag()) < 1e-12) ? x.real() : NAN;
}
#endif //INC_5_TIMES_EQUATION_BIQUADRATE2_H
