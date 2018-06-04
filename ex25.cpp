#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <limits>
//#include "ex21.cpp"

#define INF 1000000

using namespace std;
//F(x) = 0
double f1(double x1, double x2) { return 3*x1*x1 - x1 + x2*x2 - 1;}
double f2(double x1, double x2) { return x2 - tan(x1);}
//phi(x) = x
double phi1(double x1, double x2) { return sqrt((x1 - x2*x2 + 1)/3);}
double phi2(double x1, double x2) { return tan(x1);}
// Jacoby
double df11(double x1, double x2) { return  6*x1 - 1;}
double df12(double x1, double x2) { return 2*x2;}
double df21(double x1, double x2) { return -1/(cos(x1)*cos(x1));}
double df22(double x1, double x2) { return 1.0;}

double maxdiff(vector<double> &v1, vector<double> &v2) {
    double max = fabs(v1[0] - v2[0]);
    for (int i = 1; i < v1.size(); i++)
        if (fabs(v1[i] - v2[i]) > max) max = fabs(v1[i] - v2[i]);
    return max;
}

vector<double> gauss(vector< vector<double> > A) {
    int n = A.size();

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

vector<double> iter(vector<double (*)(double, double)> &phi, double eps = 0.001) {
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> dist(0.0, nextafter(1.0, numeric_limits<decltype(1.0)>::max()));
    int n = 2;
    int k = 0;
    vector<double> px(n), x(n);
    for (int i = 0; i < n; i++)   px[i] = dist(generator);
    //px[0] = 0.2; px[1] = 0.3;
    //x[0] = phi1(px[0], px[1]); x[1] = phi2(px[0], px[1]);
    for (int i = 0; i < n; i++) x[i] = phi[i](px[0], px[1]);
    while (maxdiff(x, px) > eps) {
        px = x;
        //x[0] = phi1(px[0], px[1]); x[1] = phi2(px[0], px[1]);
        for (int i = 0; i < n; i++) x[i] = phi[i](px[0], px[1]);
        k++;
    }
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << "\nnumber of iterations: " << k++ << endl;
    return x;
}

vector<double> seidel(vector<double (*)(double, double)> &phi, double eps = 0.001) {
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> dist(0.0, nextafter(1.0, numeric_limits<decltype(1.0)>::max()));
    int n = 2;
    int k = 0;
    vector<double> px(n), x(n);
    //for (int i = 0; i < n; i++) px[i] = dist(generator);
    px[0] = 0.2; px[1] = 0.3;
    x[0] = phi1(px[0], px[1]); x[1] = phi2(x[0], px[1]);
    //for (int i = 0; i < n; i++) x[i] = phi[i](x[0], x[1]);
    while (maxdiff(x, px) > eps) {
        px = x;
        x[0] = phi1(px[0], px[1]); x[1] = phi2(x[0], px[1]);
        //for (int i = 0; i < n; i++) x[i] = phi[i](x[0], x[1]);
        k++;
    }
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << "\nnumber of iterations: " << k++ << endl;
    return x;
}

vector<double> newton(vector<double (*)(double, double)> &F, double eps = 0.001) {
    int n = 2, k = 0;
    vector<double> px(n), x(n);
    px[0] = 0.2; px[1] = 0.3;
    vector<vector<double> > W(2, vector<double>(2));
    W[0][0] = df11(px[0], px[1]); W[0][1] = df12(px[0], px[1]);
    W[1][0] = df21(px[0], px[1]); W[1][1] = df22(px[0], px[1]);
    vector<double> f(2); f[0] = -f1(px[0], px[1]);  f[1] = -f2(px[0], px[1]);
    vector<vector<double> > A(2, vector<double>(3));
    A[0][0] = W[0][0]; A[0][1] = W[0][1]; A[0][2] = f[0];
    A[1][0] = W[1][0]; A[1][1] = W[1][1]; A[1][2] = f[1];
    vector<double> dx = gauss(A);
    for (int i = 0; i < n; i++) x[i] = px[i] + dx[i];
    while (maxdiff(px, x) > eps) {
        px = x;
        W[0][0] = df11(px[0], px[1]); W[0][1] = df12(px[0], px[1]);
        W[1][0] = df21(px[0], px[1]); W[1][1] = df22(px[0], px[1]);
        f[0] = -f1(px[0], px[1]);  f[1] = -f2(px[0], px[1]);
        A[0][0] = W[0][0]; A[0][1] = W[0][1]; A[0][2] = f[0];
        A[1][0] = W[1][0]; A[1][1] = W[1][1]; A[1][2] = f[1];
        dx = gauss(A);
        for (int i = 0; i < n; i++) x[i] = px[i] + dx[i];
        k++;
    }
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << "}\nnumber of iterations: " << k++ << endl;
    return x;
}


int main() {
    double (*pf1)(double, double) = f1, (*pf2)(double, double) = f2;
    double (*pphi1)(double, double) = phi1, (*pphi2)(double, double) = phi2;
    vector<double (*)(double, double)> F;
    vector<double (*)(double, double)> Phi;
    F.push_back(pf1); Phi.push_back(pphi1);
    F.push_back(pf2); Phi.push_back(pphi2);
    //cout << Phi[1](M_PI/4, 0.0) << endl;
    cout << "iter method, ans: v = {";
    vector<double> x1 = iter(Phi);
    cout << "seidel method, ans: v = ";
    vector<double> x2 = seidel(Phi);
    cout << "newton method, ans: v = ";
    vector<double> x3 = newton(F);
}
