#include "lab8.h"



int main() {

    double alpha = 3.52; 
    double beta = 0.01;
    double gamma = 0.2952;
    double delta = 0.0002;

    double x0 = 2303;
    double y0 = 993; 

    double h = 0.0001; 
    double t_end = 365; // год

    Lotka_Volterra obj(alpha, beta, gamma, delta, x0, y0);

    obj.eulerMethod(h, t_end);

    obj.rungeKutta4(h, t_end);

    obj.trapezoidMethod(h, t_end);

    return 0;
}