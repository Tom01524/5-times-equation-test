#include <iostream>
#include <cmath> // for pow()
#include <cassert> // for assert()
#include<random> // for random number generation

struct coefficient5 {
    double a, b, c, d, e, f;
};
struct coefficient6 {
    double a, b, c, d, e, f, g;
    double operator[](int index) const {
        switch (index) {
            case 0:return a;
            case 1:return b;
            case 2:return c;
            case 3:return d;
            case 4:return e;
            case 5:return f;
            case 6:return g;
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
            case 6:return g;
            default: throw std::invalid_argument("index out of range");
        }
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
static vec4 solveEquation_4(coefficient5 src);

static std::vector<vec5> solveEquation_5(coefficient6 src);
static bool checkSolution5(coefficient6 src, vec5 res);

int main() {
    class std::random_device rd;  // 获取一个随机数，作为种子
    std::mt19937 gen(rd()); // 使用Mersenne Twister算法生成随机数
    // 定义一个均匀分布的范围
    class std::uniform_int_distribution<> distr(1, 35);

    coefficient6 input3{};
    for (int i = 0; i < 5; i++) {
        double _1 = static_cast<double>(distr(gen)) / 10.0;
        double _2 = static_cast<double>(distr(gen)) / 10.0;
        double randnow = pow(_1, _2) / pow(_2, _1) + pow(_1, _2);
        input3[i] = randnow;
        printf("rannow: %.9f\n", randnow);
    }
    input3[5] = 1.0;
    std::vector<vec5> res2 = solveEquation_5(input3);
    for(int i=0;i<res2.size();i++){
        printf("x1: %.9f, x2: %.9f, x3: %.9f, x4: %.9f, x5: %.9f\n",
               res2[i].x1, res2[i].x2, res2[i].x3, res2[i].x4, res2[i].x5);
        checkSolution5(input3,res2[i]);
    }

//    coefficient5 input1 = {210.00, 1100.00, 2153.00, 1872.00, 594.00};
//    coefficient5 input2 = {input3.a, input3.b, input3.c, input3.d, input3.e};
//    vec4 res = solveEquation_4(input2);
//    printf("x1: %.9f, x2: %.9f, x3: %.9f, x4: %.9f\n\n", res.x1, res.x2, res.x3, res.x4);

//    coefficient6 input3 = {2.1, 11.00, 21.53, 18.72, -29.4, 1.0};



    return 0;
}
vec4 solveEquation_4(coefficient5 src){
    double A = src.a, B = src.b, C = src.c, D = src.d, E = src.e;
    // ax^4 + bx^3 + cx^2 + dx + e = 0
    double delta1 = C*C - 3.0*B*D + 12.0*A*E;
    double delta2 = 2.0*C*C*C - 9.0*B*C*D + 27.0*A*D*D + 27.0*B*B*E - 72.0*A*C*E;

    double deltaL = (pow(2.0, 1/3.0)*delta1) /
            (3.0*A * pow( delta2 + sqrt(-4.0*delta1*delta1*delta1 + delta2*delta2) , 1/3.0));

    double deltaR = (pow( delta2 + sqrt(-4.0*delta1*delta1*delta1 + delta2*delta2) , 1/3.0)) /
            (3.0 * pow(2.0, 1/3.0) * A);
    double delta = deltaL + deltaR;

    double xL_minus = -B/(4.0*A) - 0.5 * sqrt((B*B)/(4.0*A*A) - (2*C)/(3*A) + delta);
    double xL_plus = -B/(4.0*A) + 0.5 * sqrt((B*B)/(4.0*A*A) - (2*C)/(3*A) + delta);

    double xR_minus = 0.5 * sqrt(
            (B*B)/(2*A*A) -
            (4*C)/(3*A) -
            delta -
            (
                    (-(B*B*B)/(A*A*A) + (4.0*B*C)/(A*A) - (8.0*D)/A) /
                    (4 * sqrt((B*B)/(4*A*A) - (2*C)/(3*A) + delta))
            )
            );
    double xR_plus = 0.5 * sqrt(
            (B*B)/(2*A*A) -
            (4*C)/(3*A) -
            delta +
            (
                    (-(B*B*B)/(A*A*A) + (4.0*B*C)/(A*A) - (8.0*D)/A) /
                    (4 * sqrt((B*B)/(4*A*A) - (2*C)/(3*A) + delta))
            )
            );
    double x1 = xL_minus - xR_minus;
    double x2 = xL_minus + xR_minus;
    double x3 = xL_plus - xR_plus;
    double x4 = xL_plus + xR_plus;
    return vec4{x1, x2, x3, x4};
}
std::vector<vec5> solveEquation_5(coefficient6 src){
    assert(abs(src.f - 1.0) < 1e-12); // 保证系数不为1
    // ag⁴ - bg³ + cg² - dg + e - 1= 0
    vec4 res = solveEquation_4(coefficient5{src.a, -src.b, src.c, -src.d, src.e-1.0});
    if(std::isnan(res.x1) && std::isnan(res.x2) && std::isnan(res.x3) && std::isnan(res.x4)){
        std::cerr << "Error: Can't find solution" << std::endl;
        std::vector<vec5> empty;
        return empty;
    }

    double g = std::numeric_limits<double>::quiet_NaN();
    std::vector<vec5> res5;
    for(int i=0; i<4; i++){
        if(!isnan(res[i])){ // 找到第一个非nan解
            g = res[i];
            double a = src.a;
            double b = src.b - src.a*g;
            double c = src.c - src.b*g + src.a*g*g;
            double d = src.d - src.c*g + src.b*g*g - src.a*g*g*g;
            std::cout << "(" << a << "x⁴ + " << b << "x³ + " << c << "x² + " << d << "x + 1) * (x + " << g << ") = 0" << std::endl;
            double A = a, B = a*g+b, C = b*g+c, D = c*g+d, E=d*g+1.0;
            std::cout <<  A << "x⁵ + " << B << "x⁴ + " << C << "x³ + " << D << "x² + " << E << "x + 1 = 0\n" << std::endl;
            // (ax⁴ + bx³ + cx² + dx + 1)(x + g) = 0
            vec4 front4 = solveEquation_4(coefficient5{a, b, c, d, 1.0});
            double x5 = -g;
            res5.push_back(vec5{front4.x1, front4.x2, front4.x3, front4.x4, x5});
        }
    }


    return res5;
}
bool checkSolution5(coefficient6 src, vec5 res){
    double check1 = src.a * res.x1 * res.x1 * res.x1 * res.x1 * res.x1 +
                    src.b * res.x1 * res.x1 * res.x1 * res.x1 +
                    src.c * res.x1 * res.x1 * res.x1 +
                    src.d * res.x1 * res.x1 +
                    src.e * res.x1 +
                    src.f;
    bool correct1 = abs(check1) < 1e-6;
    if(correct1){
        std::cout << "Solution1: " << res.x1 << "→ " << check1 << " is correct" << std::endl;
    }else{
        std::cout << "Solution1: " << res.x1 << "→ " << check1 << " is incorrect!!!" << std::endl;
    }
    double check2 = src.a * res.x2 * res.x2 * res.x2 * res.x2 * res.x2 +
                    src.b * res.x2 * res.x2 * res.x2 * res.x2 +
                    src.c * res.x2 * res.x2 * res.x2 +
                    src.d * res.x2 * res.x2 +
                    src.e * res.x2 +
                    src.f;
    bool correct2 = abs(check2) < 1e-6;
    if(correct2){
        std::cout << "Solution2: " << res.x2 << "→ " << check2 << " is correct" << std::endl;
    }else{
        std::cout << "Solution2: " << res.x2 << "→ " << check2 << "is incorrect!!!" << std::endl;
    }
    double check3 = src.a * res.x3 * res.x3 * res.x3 * res.x3 * res.x3 +
                    src.b * res.x3 * res.x3 * res.x3 * res.x3 +
                    src.c * res.x3 * res.x3 * res.x3 +
                    src.d * res.x3 * res.x3 +
                    src.e * res.x3 +
                    src.f;
    bool correct3 = abs(check3) < 1e-6;
    if(correct3){
        std::cout << "Solution3: " << res.x3 << "→ " << check3 << " is correct" << std::endl;
    }else{
        std::cout << "Solution3: " << res.x3 << "→ " << check3 << " is incorrect!!!" << std::endl;
    }
    double check4 = src.a * res.x4 * res.x4 * res.x4 * res.x4 * res.x4 +
                    src.b * res.x4 * res.x4 * res.x4 * res.x4 +
                    src.c * res.x4 * res.x4 * res.x4 +
                    src.d * res.x4 * res.x4 +
                    src.e * res.x4 +
                    src.f;
    bool correct4 = abs(check4) < 1e-6;
    if(correct4){
        std::cout << "Solution4: " << res.x4 << "→ " << check4 << " is correct" << std::endl;
    }else{
        std::cout << "Solution4: " << res.x4 << "→ " << check4 << " is incorrect!!!" << std::endl;
    }
    double check5 = src.a * res.x5 * res.x5 * res.x5 * res.x5 * res.x5 +
                    src.b * res.x5 * res.x5 * res.x5 * res.x5 +
                    src.c * res.x5 * res.x5 * res.x5 +
                    src.d * res.x5 * res.x5 +
                    src.e * res.x5 +
                    src.f;
    bool correct5 = abs(check5) < 1e-6;
    if(correct5){
        std::cout << "Solution5: " << res.x5 << "→ " << check5 << " is correct" << std::endl;
    }else{
        std::cout << "Solution5: " << res.x5 << "→ " << check5 << " is incorrect!!!\n" << std::endl;
    }
    return correct1 && correct2 && correct3 && correct4 && correct5;
}