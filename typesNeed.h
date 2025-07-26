//
// Created by iceoc0 on 2025/7/26.
//

#ifndef INC_5_TIMES_EQUATION_TYPESNEED_H
#define INC_5_TIMES_EQUATION_TYPESNEED_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <cassert>
#include <complex>



struct vec4 {
    double x1, x2; // 两个解
    double x3, x4; // 两个解
    double operator[](int index) const {
        switch (index) {
            case 0:
                return x1;
            case 1:
                return x2;
            case 2:
                return x3;
            case 3:
                return x4;
            default:
                throw std::invalid_argument("index out of range");
        }
    }
    // 写入版本 (返回引用)
    double& operator[](int index) {
        switch (index) {
            case 0: return x1;
            case 1: return x2;
            case 2: return x3;
            case 3: return x4;
            default: throw std::invalid_argument("index out of range");
        }
    }

};
struct vec4_ {
    std::complex<double> x1, x2; // 两个解
    std::complex<double> x3, x4; // 两个解
    std::complex<double> operator[](int index) const {
        switch (index) {
            case 0:
                return x1;
            case 1:
                return x2;
            case 2:
                return x3;
            case 3:
                return x4;
            default:
                throw std::invalid_argument("index out of range");
        }
    }
    // 写入版本 (返回引用)
    std::complex<double>& operator[](int index) {
        switch (index) {
            case 0: return x1;
            case 1: return x2;
            case 2: return x3;
            case 3: return x4;
            default: throw std::invalid_argument("index out of range");
        }
    }

};
struct vec5 {
    double x1, x2; // 解
    double x3, x4, x5; // 解
    double operator[](int index) const {
        switch (index) {
            case 0:return x1;
            case 1:return x2;
            case 2:return x3;
            case 3:return x4;
            case 4:return x5;
            default:
                throw std::invalid_argument("index out of range");
        }
    }
    // 写入版本 (返回引用)
    double& operator[](int index) {
        switch (index) {
            case 0: return x1;
            case 1: return x2;
            case 2: return x3;
            case 3: return x4;
            case 4: return x5;
            default: throw std::invalid_argument("index out of range");
        }
    }

};
struct vec5_ {
    std::complex<double> x1, x2; // 解
    std::complex<double> x3, x4, x5; // 解
    std::complex<double> operator[](int index) const {
        switch (index) {
            case 0:return x1;
            case 1:return x2;
            case 2:return x3;
            case 3:return x4;
            case 4:return x5;
            default:
                throw std::invalid_argument("index out of range");
        }
    }
    // 写入版本 (返回引用)
    std::complex<double>& operator[](int index) {
        switch (index) {
            case 0: return x1;
            case 1: return x2;
            case 2: return x3;
            case 3: return x4;
            case 4: return x5;
            default: throw std::invalid_argument("index out of range");
        }
    }

};
struct coefficient5 { double a, b, c, d, e; };
struct coefficient5_ {
    std::complex<double> a, b, c, d, e;
};
struct coefficient6 {
    double a, b, c, d, e, f;
    double operator[](int index) const {
        switch (index) {
            case 0:return a;
            case 1:return b;
            case 2:return c;
            case 3:return d;
            case 4:return e;
            case 5:return f;

            default:
                throw std::invalid_argument("index out of range");
        }
    }
    // 写入版本 (返回引用)
    double& operator[](int index) {
        switch (index) {
            case 0:return a;
            case 1:return b;
            case 2:return c;
            case 3:return d;
            case 4:return e;
            case 5:return f;

            default: throw std::invalid_argument("index out of range");
        }
    }


};
struct coefficient6_ {
    std::complex<double> a, b, c, d, e, f;
    std::complex<double> beforeF;
    std::complex<double> operator[](int index) const {
        switch (index) {
            case 0:return a;
            case 1:return b;
            case 2:return c;
            case 3:return d;
            case 4:return e;
            case 5:return f;
            case 6:return beforeF;
            default:
                throw std::invalid_argument("index out of range");
        }
    }
    // 写入版本 (返回引用)
    std::complex<double>& operator[](int index) {
        switch (index) {
            case 0:return a;
            case 1:return b;
            case 2:return c;
            case 3:return d;
            case 4:return e;
            case 5:return f;
            case 6:return beforeF;
            default: throw std::invalid_argument("index out of range");
        }
    }
    void setBeforeF(){
        beforeF = f;
        std::cout<<"beforeF set to "<<beforeF<<std::endl;
    }
    coefficient6_ transform2Back1(){
        coefficient6_ res;

        res.a = this->a/this->f;
        res.b = this->b/this->f;
        res.c = this->c/this->f;
        res.d = this->d/this->f;
        res.e = this->e/this->f;
        res.f = 1.0;
        res.beforeF = this->beforeF;

        return res;
    }

};




vec4_ solveEquation_4_(coefficient5_ src){
    using namespace std::complex_literals;
    std::complex<double> A = src.a, B = src.b, C = src.c, D = src.d, E = src.e;

    // 安全除法：处理分母接近零的情况
    auto safe_divide =
            [](std::complex<double> num, std::complex<double> den) {
                const double epsilon = 1e-12;
                if (std::abs(den) < epsilon) return num / epsilon;
                return num / den;
            };

    // 计算中间变量（复数版本）
    std::complex<double> delta1 = C*C - 3.0*B*D + 12.0*A*E;
    std::complex<double> delta2 = 2.0*C*C*C - 9.0*B*C*D + 27.0*A*D*D + 27.0*B*B*E - 72.0*A*C*E;
    std::complex<double> sqrt_term = std::sqrt(-4.0*delta1*delta1*delta1 + delta2*delta2);

    // 计算 deltaL 和 deltaR（使用安全除法）
    std::complex<double> denom_pow = std::pow(delta2 + sqrt_term, 1.0/3.0);
    std::complex<double> deltaL = safe_divide(std::pow(2.0, 1.0/3.0) * delta1, 3.0*A * denom_pow);
    std::complex<double> deltaR = safe_divide(denom_pow, 3.0 * std::pow(2.0, 1.0/3.0) * A);
    std::complex<double> delta = deltaL + deltaR;

    // 计算 xL/xR（自动处理复数平方根）
    std::complex<double> term_B2_4A2 = (B*B)/(4.0*A*A);
    std::complex<double> sqrt_arg_L = term_B2_4A2 - (2.0*C)/(3.0*A) + delta;
    std::complex<double> sqrt_L = std::sqrt(sqrt_arg_L);
    std::complex<double> xL_minus = -B/(4.0*A) - 0.5 * sqrt_L;
    std::complex<double> xL_plus  = -B/(4.0*A) + 0.5 * sqrt_L;

    // 计算 xR_minus/xR_plus（处理分母安全）
    std::complex<double> numerator = -(B*B*B)/(A*A*A) + (4.0*B*C)/(A*A) - (8.0*D)/A;
    std::complex<double> denom_xR = 4.0 * sqrt_L;
    std::complex<double> term_inside = safe_divide(numerator, denom_xR);

    std::complex<double> sqrt_arg_R1 = (B*B)/(2.0*A*A) - (4.0*C)/(3.0*A) - delta - term_inside;
    std::complex<double> sqrt_arg_R2 = (B*B)/(2.0*A*A) - (4.0*C)/(3.0*A) - delta + term_inside;
    std::complex<double> xR_minus = 0.5 * std::sqrt(sqrt_arg_R1);
    std::complex<double> xR_plus  = 0.5 * std::sqrt(sqrt_arg_R2);

    // 计算最终解并提取实部
//    const double imag_threshold = 1e-6;
//    auto get_real = [imag_threshold](std::complex<double> x) -> double {
//        return (std::abs(x.imag()) < imag_threshold) ? x.real() : NAN;
//    };

    std::complex<double> x1 = xL_minus - xR_minus;
    std::complex<double> x2 = xL_minus + xR_minus;
    std::complex<double> x3 = xL_plus  - xR_plus;
    std::complex<double> x4 = xL_plus  + xR_plus;

    return vec4_{x1, x2, x3, x4};
}
vec4 solveEquation_4(coefficient5 src){
    using namespace std::complex_literals;
    std::complex<double> A = src.a, B = src.b, C = src.c, D = src.d, E = src.e;
    vec4_ res = solveEquation_4_(coefficient5_{A, B, C, D, E});
    vec4 res2;
    for(int i=0;i<4;i++){
        if(abs(res[i].imag())<1e-15){
            res2[i] = res[i].real();
        }else{
            res2[i] = NAN;
        }
    }
    return res2;
}

#endif //INC_5_TIMES_EQUATION_TYPESNEED_H
