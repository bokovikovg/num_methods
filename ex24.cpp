#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>
using namespace std;

double f(double x) { return sin(x) - x*x + 1;}
double df(double x) { return cos(x) - 2*x;}
double d2f(double x) {return -sin(x) - 2;}
double phi(double x) {return (sin(x) + 1)/x;}
double phi1(double x) { return sin(x)-x*x+1+x;}

// solve f(x) = 0
double dihotomie(double (*pf)(double), double a, double b, double eps = 0.001) {
    int k = 0;
    double c = (a + b)/2;
    while (b - a > eps) {
        if (pf(a)*pf(c) < 0) {
            a = a;
            b = c;
        }
        if (pf(c)*pf(b) < 0) {
            a = c;
            b = b;
        }
        c = (a + b)/2;
        k++;
    }
    //cout << "\nnumber of iterations k = " << k++ << endl;
    return c;
}

double fixed_point(double (*pf)(double), double a, double b, double eps = 0.001) {
    int k = 0;
    double px = (a + b)/2, x = pf(px);
    cout << "iteration k = " << k + 1 << " x = " << x << " f(x) = " << f(x) << endl;
    while (abs(x - px) > eps && k < 200) {
        k++;
        px = x;
        x = pf(px);
        //cout << "mod diff = " << abs(pf(x));
        cout << "iteration k = " << k + 1 << " x = " << x << " f(x) = " << f(x) << endl;
    }
    //cout << "\nnumber of iterations k = " << k++ << endl;
    return x;
}

double newton_method(double (*pf)(double), double (*pdf)(double), double a, double b, double eps = 0.001) {
    if (a > b) {
      double temp = a;
      a = b;
      b = temp;
    }
    if (pf(a) * pf(b) > 0) cout << "\nError: no roots in this interval\n";
    int k = 0;
    double px = (a + b)/2, x = px - pf(px)/pdf(px);
    cout << "iteration k = " << k + 1 << " x = " << x << " f(x) = " << pf(x) << endl;
    while (abs(x - px) > eps) {
        px = x;
        x = px - pf(px)/pdf(px);
        k++;
        cout << "iteration k = " << k + 1 << " x = " << x << " f(x) = " << pf(x) << endl;
    }
    //cout << "\nnumber of iterations k = " << k++ << endl;
    return x;
}

double last_method(double (*pf)(double), double a, double b, double eps = 0.001) {
    if (pf(a) * pf(b) > 0) cout << "\nError: no roots in this interval\n";
    int k = 0;
    double dx = 0.1;
    double px = (a + b)/2, x = px - pf(px)*dx/(pf(px + dx) - pf(px));
    while (abs(pf(x)) > eps) {
        px = x;
        x = px - pf(px)*dx/(pf(px + dx) - pf(px));
        k++;
    }
    //cout << "\nnumber of iterations k = " << k++ << endl;
    return x;
}


int main(int argc, const char * argv[]) {
    double (*pf)(double) = f;
    double (*pfi)(double) = phi;
    double (*pfi1)(double) = phi1;
    double (*pdf)(double) = df;
    double res1, res2, eps = 0.001;
    double a = -2.0, b = 0.0, c = 2.0;
    cout << "dihotomie method, eps = " << fixed << setprecision(10) << eps << ":\n";
    res1 = dihotomie(pf, a, b);
    res2 = dihotomie(pf, b, c);
    cout << "x1 = " << fixed << setprecision(10) << res1 << " x2 = " << res2 << endl;
    cout << "fixed point method\n";
    res1 = fixed_point(pfi, a, b);
    res2 = fixed_point(pfi, b, c);
    cout << "x1 = " << fixed << setprecision(10) << res1 << " x2 = " << res2 << endl;
    cout << "newton method\n";
    res1 = newton_method(pf, pdf, a, b);
    res2 = newton_method(pf, pdf, b, c);
    cout << "x1 = " << fixed << setprecision(10) << res1 << " x2 = " << res2 << endl;
    cout << "tangential method\n";
    res1 = last_method(pf, a, b);
    res2 = last_method(pf, b, c);
    cout << "x1 = " << fixed << setprecision(10) << res1 << " x2 = " << res2 << endl;
    return 0;
}
