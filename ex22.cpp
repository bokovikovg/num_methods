#include <iostream>
#include <vector>

using namespace std;

vector<double> shuttleMatrix(vector<vector<double> > &A, vector<double> &Y) {
    //with the power of shitecode
    int n = (int)Y.size();
    vector<double> x(n);
    vector<double> gamma(n);
    vector<double> alpha(n);
    vector<double> beta(n);
    vector<double> b(n);
    vector<double> a(n);
    vector<double> c(n);
    // hue hue hue
    b[0] = A[0][0];
    c[0] = A[0][1];
    for (int i = 1; i < n - 1; i++) {
        a[i] = A[i][i-1];
        b[i] = A[i][i];
        c[i] = A[i][i+1];
    }
    a[n-1] = A[n-1][n-2];
    b[n-1] = A[n-1][n-1];
    //b[0] = A[0][0];
    gamma[0] = b[0]; //c[0] = A[0][1];
    alpha[0] = -c[0]/gamma[0];
    beta[0] = Y[0]/gamma[0];
    for(int i = 1; i < n - 1; i++) {
        gamma[i] = b[i] + a[i]*alpha[i-1];
        alpha[i] = -c[i]/gamma[i];
        beta[i] = (Y[i] - a[i]*beta[i-1])/gamma[i];
    }
    gamma[n-1] = b[n-1] + a[n-1]*alpha[n-2];
    beta[n-1] = (Y[n-1] - a[n-1]*beta[n-2])/gamma[n-1];
    x[n-1] = beta[n-1];
    for (int i = n-2; i >= 0; i--) x[i] = alpha[i]*x[i+1] + beta[i];
    cout << "прогоночные коэффициенты\n";
    cout << "alpha\tbeta\n";
    for (int i = 0; i < n; i++) cout << alpha[i] << "\t" << beta[i] << endl;
    return x;
}

int main() {
    vector<vector<double> > A = {{12, -5, 0, 0, 0},
                                {-3, 18, -8, 0, 0},
                                {0, -2, -16, -9, 0},
                                {0, 0, -4, 18, -7},
                                {0, 0, 0, 4, -9}};
    vector<double> Y = {148, 45, -155, 11, 3};
    vector<double> x = shuttleMatrix(A, Y);
    cout << "the resulting vector:\n";
    for (int i = 0; i < x.size(); i++) cout << x[i] << " ";
    cout << "\ntesting...\n";
    for (int i = 0; i < x.size(); i++) {
      double sum = 0.0;
      for (int j = 0; j < x.size(); j++) {
        sum += A[i][j] * x[j];
      }
      cout << Y[i] - sum << endl;
    }
    cout << endl;
    return 0;
}
