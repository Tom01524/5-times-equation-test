//
// Created by iceoc0 on 2025/7/26.
//
#include "biquadrate2.h"


void test2(coefficient6_ input3){


//    coefficient6 input4={1.0,0.0,10.0,0.0,20.0,-4.0};
    coefficient6_ after = input3.transform2Back1();
    std::cout << "after a: " << after.a << ", b: " << after.b << ", c: " << after.c <<
              ", d: " << after.d << ", e: " << after.e << ", f: " << after.f << std::endl;

    std::vector<vec5_> res2 = solveEquation_5_(after);
    std::vector<vec5_>::size_type sizeNow = res2.size();
    std::cout << "res2 size: " << sizeNow << std::endl;
    if(sizeNow == 0){
        std::cerr << "Error: Can't find solution" << std::endl;
    }else{
        for(int i=0;i<res2.size();i++){
            std::cout << "x1: " <<   res2[i].x1 << ", x2: " <<   res2[i].x2 <<
                      ", x3: " <<   res2[i].x3 << ", x4: " <<   res2[i].x4 << ", x5: " <<   res2[i].x5 << std::endl;
            checkSolution5_(after,res2[i]);
        }
    }

}

#include "biquadrate3.h"
void test3(coefficient6 input3){


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
    puts("");
    for (double root : roots) {
        printf("root: x=%.16f, ", root);
        double y_ = checkRoot(input3,root);
        if(y_ < 1e-10){
            printf("y=%.16f √\n",y_);
        }else{
            printf("y=%.16f ×\n",y_);
        }
    }
}

int main(){
    class std::random_device rd;  // 获取一个随机数，作为种子
    std::mt19937 gen(rd()); // 使用Mersenne Twister算法生成随机数
    // 定义一个均匀分布的范围
    class std::uniform_real_distribution<> distr(-4.0, 4.0);

    coefficient6 input3{};
    coefficient6_ input3_{};

    for (int i = 0; i < 6; i++) {
        double randnow = static_cast<double>(distr(gen));

        input3[i] = randnow;
        input3_[i] = randnow;
        printf("randnow: %.9f\n", randnow);
    }
    input3_.setBeforeF();
    test2(input3_);
    std::cout << "=================================" <<
    std::endl;
    test3(input3);
    return 0;
}
