//
// Created by iceoc0 on 2025/7/22.
//
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <cassert>
#include <complex>

using namespace std;

const double INF = numeric_limits<double>::infinity();
const double accuracy = 1e-12;

struct coefficient5 { double a, b, c, d, e; };
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
    coefficient6 transform2Back1(){
        coefficient6 res;
        res.a = a/f;
        res.b = b/f;
        res.c = c/f;
        res.d = d/f;
        res.e = e/f;
        res.f = 1.0;
        return res;
    }

};
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
template<typename Func>
static double bisection(Func&& func_, double a, double b, double accurate = 1e-12, int max_iter = 1000);
static vector<double> solve_quintic(double a, double b, double c, double d, double e, double f_coeff);
// 示例四次方程求解函数（需确保能返回实数解）
static vec4 solveEquation_4(coefficient5 src);

/*
 * rannow: 0.086743871
 * rannow: 4.169851147
 * rannow: 20.808962361
 * rannow: 9.237061995
 * rannow: 11.000000000
 * rannow: 3.662050734
 * 凸点[0] x: -34.2954 y: 824169
 * 凸点[1] x: -3.89681 y: -246.704
 * root: -42.4812557145212111
 * root: -5.2246874213288308
 * root: -0.3598134920840291
 */


int main() {
    class std::random_device rd;  // 获取一个随机数，作为种子
    std::mt19937 gen(rd()); // 使用Mersenne Twister算法生成随机数
    // 定义一个均匀分布的范围
    class std::uniform_int_distribution<> distr(1, 35);

    coefficient6 input3{};
    for (int i = 0; i < 6; i++) {
        double _1 = static_cast<double>(distr(gen)) / 10.0;
        double _2 = static_cast<double>(distr(gen)) / 10.0;
        double randnow = pow(_1, _2) / pow(_2, _1) + pow(_1, _2);
        input3[i] = randnow;
//        printf("rannow: %.9f\n", randnow);
    }
//    input3[0] = 1.0;
//    input3[1] = 0.0;
//    input3[2] = -10.0;
//    input3[3] = 0.0;
//    input3[4] = 20.0;
//    input3[5] = 0.0;
    printf("%.9fx⁵ %c %.9fx⁴ %c %.9fx³ %c %.9fx² %c %.9fx %c %.9f = 0\n",
           input3[0],input3[1] > 0.0 ? '+' : '-', abs(input3[1]),
           input3[2] > 0.0 ? '+' : '-', abs(input3[2]),
           input3[3] > 0.0 ? '+' : '-', abs(input3[3]),
           input3[4] > 0.0 ? '+' : '-', abs(input3[4]),
           input3[5] > 0.0 ? '+' : '-', abs(input3[5]));
    // 1.0,0.0,10.0,0.0,20.0,-4.0
    vector<double> roots = solve_quintic(input3[0], input3[1], input3[2], input3[3], input3[4], input3[5]);
    for (double root : roots) {
        printf("root: x=%.16f, y=0\n", root);
    }
    return 0;
}
// 二分法求根函数
template<typename funcPtr>
double bisection(funcPtr&& func_, double a, double b, double accurate, int max_iter) {
    double fa = func_(a);
    double fb = func_(b);

    if (fa * fb > 0) return NAN;
    if (fa == 0) return a;
    if (fb == 0) return b;

    for (int i = 0; i < max_iter && (b - a) > 1.0E-12; ++i) {
        double c = (a + b) / 2;
        double fc = func_(c);

        if (fc == 0) return c;
        if (fa * fc < 0) b = c;
        else a = c;
    }
    return (a + b) / 2;
}

vector<double> solve_quintic(double a, double b, double c, double d, double e, double f_coeff) {
    bool haveNAN = isnan(a) || isnan(b) || isnan(c) || isnan(d) || isnan(e) || isnan(f_coeff);
    assert(!haveNAN);
    vector<double> roots;

    // 1. 计算导数并求解临界点 get acme
    coefficient5 deriv_coeff = {5*a, 4*b, 3*c, 2*d, e};
    vec4 crits = solveEquation_4(deriv_coeff);

    // 2. 收集实数临界点并排序 order acme
    vector<double> acme;int nanCount = 0;
    for (int i = 0; i < 4; ++i) {
        if (!isnan(crits[i])){  // 简单过滤无效解
            cout << "凸点[" << i << "] x: " << crits[i];
            acme.push_back(crits[i]);
            double y = a*pow(crits[i],5.0) + b*pow(crits[i],4.0) +
                    c*pow(crits[i],3.0) + d*pow(crits[i],2.0) + e*crits[i] + f_coeff;
            cout << " y: " << y << endl;
        } else {
            nanCount++;
        }
    }
    if(nanCount == 4){
        cout << "5 次方程 没有 凸点" << endl;
    }
    sort(acme.begin(), acme.end());


    // 3. 检查临界点是否为根 check acme → roots
    auto quinticGetY = [&](double x) {
        return a*pow(x,5) + b*pow(x,4) + c*pow(x,3) + d*pow(x,2) + e*x + f_coeff;
    };

    for (double topBottom : acme) {
        if (abs(quinticGetY(topBottom)) < accuracy)
            roots.push_back(topBottom);
    }

    // 4. 构建区间（包含正负无穷）get intervals
    vector<double> intervals = {-INF};
    intervals.insert(intervals.end(), acme.begin(), acme.end());
    intervals.push_back(INF);

    // 5. 遍历所有区间寻找根 find roots in intervals
    for (size_t i = 0; i < intervals.size()-1; ++i) {
        double left = intervals[i];
        double right = intervals[i+1];

        // 计算端点符号
        auto get_sign = [&](double x) {
            if (x == -INF) return a > 0 ? -1 : 1;
            if (x == INF) return a > 0 ? 1 : -1;
            return quinticGetY(x) > 0 ? 1 : -1;
        };

        int sign_l = get_sign(left);
        int sign_r = get_sign(right);

        if (sign_l * sign_r >= 0) continue;

        // 确定数值边界
        double l = (left == -INF) ? -1e20 : left;
        double r = (right == INF) ? 1e20 : right;

        // 二分法求根
        double root = bisection(quinticGetY, l, r);
        if (!isnan(root))
            roots.push_back(root);
    }

    // 6. 去重排序 remove duplicates and sort
    sort(roots.begin(), roots.end());
    auto last = unique(roots.begin(), roots.end(),
                       [](double a, double b) { return abs(a - b) < 1e-6; });
    roots.erase(last, roots.end());

    // 最多返回5个实根
    if (roots.size() > 5)
        roots.resize(5);

    return roots;
}
vec4 solveEquation_4(coefficient5 src){
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
    const double imag_threshold = 1e-6;
    auto get_real = [imag_threshold](std::complex<double> x) -> double {
        return (std::abs(x.imag()) < imag_threshold) ? x.real() : NAN;
    };

    double x1 = get_real(xL_minus - xR_minus);
    double x2 = get_real(xL_minus + xR_minus);
    double x3 = get_real(xL_plus  - xR_plus);
    double x4 = get_real(xL_plus  + xR_plus);

    return vec4{x1, x2, x3, x4};
}
