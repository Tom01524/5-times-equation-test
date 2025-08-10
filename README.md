# 5-times-equation-test
1. Ensure the 4 times equation get all correct results. Wirte the formal algorithm and imaginary number function.
2. Calculate the critical point of 5 times equation, by 4 times solve.
3. Check the critical points are root or not(mostly not), then get all monotone interval.
4. Use the iteration of bisection for each interval.
5. Get roots.
6. Check roots in 1e-10.

7. 5 times= 0 → (4times)(1times)=0 need a real root first.
8. x⁵+Bx⁴+Cx³+Dx²+Ex+f=0  (x⁴+bx³+cx²+dx+e)(x+g)=0
   g = - (real root)
   b = B - g;
   c = C - B*g + g*g;
   d = D - C*g + B*g*g - g*g*g;
   e = E - D*g + C*g*g - B*g*g*g + g*g*g*g;
9. 1~4 times equation have standard roots formular. We can get all real and imaginary roots.
