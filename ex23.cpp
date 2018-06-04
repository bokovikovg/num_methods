#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double findMaxDiff(vector<double> &x, vector<double> &y) {
    size_t n = x.size();
    vector<double> v(n);
    for (int i = 0; i < n; i++) v[i] = x[i] - y[i];
    double max = abs(v[0]);
    for (int i = 1; i < n; i++)
        if (abs(v[i]) > max)
            max = abs(v[i]);
    return max;
}

vector<double> iter_method(vector<vector<double> > &a, vector<double> &y) {
    int n = (int)y.size();
    vector<double> x(n);
    vector<vector<double> > alpha = a;
    vector<double> beta(n);
    for(int i = 0; i < n; i++) {
        alpha[i][i] = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            alpha[i][j] = - a[i][j]/a[i][i];
        }
        beta[i] = y[i]/a[i][i];
    }
    double eps = 0.001;
    int k = 0;
    vector<double> prev_x = beta;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) x[i] += alpha[i][j]*prev_x[j];
        x[i] += beta[i];
    }
    cout << "iteration number: " << k << "\n" << "x = ";
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << endl;
    k++;
    while (findMaxDiff(x, prev_x) > eps) {
        prev_x = x;
        for (int i = 0; i <x.size(); i++) x[i] = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                x[i] += alpha[i][j]*prev_x[j];
            }
            x[i] += beta[i];
        }
        cout << "iteration number: " << k << "\n" << "x = ";
        for (int i = 0; i < n; i++) cout << x[i] << " ";
        cout << endl;
        k++;
        //for(int i = 0; i < n; i++) temp[i] = x[i] - prev_x[i];
    }
    cout << "\nnumber of iterations: " << k << endl;
    return x;
}


vector<double> seidel_method(vector<vector<double> > &a, vector<double> &y) {
    int n = (int)y.size();
    vector<double> x(n);
    vector<vector<double> > alpha = a;
    vector<double> beta(n);
    for(int i = 0; i < n; i++) {
        alpha[i][i] = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            alpha[i][j] = - a[i][j]/a[i][i];
        }
        beta[i] = y[i]/a[i][i];
    }
    vector<vector<double> > L = alpha, U = alpha;
    for(int i = 0; i < n; i++)
        for (int j = 0 ; j < n; j++) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;
        }
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j) U[i][j] = alpha[i][j];
            else L[i][j] = alpha[i][j];
        }
    }
    double eps = 0.001;
    int k = 0;
    vector<double> prev_x = beta;
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            x[i] += L[i][j]*x[j] + U[i][j]*prev_x[j];
        }
        x[i] += beta[i];
    }
    cout << "iteration number: " << k << "\n" << "x = ";
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << endl;
    k++;
    while (findMaxDiff(x, prev_x) > eps) {
        prev_x = x;
        for (int i = 0; i <x.size(); i++) x[i] = 0.0;
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i] += L[i][j]*x[j] + U[i][j]*prev_x[j];
            }
            x[i] += beta[i];
        }
        cout << "iteration number: " << k << "\n" << "x = ";
        for (int i = 0; i < n; i++) cout << x[i] << " ";
        cout << endl;
        k++;
    }
    cout << "\nnumber of iterations: " << k << endl;
    return x;
}

int main(int argc, const char * argv[]) {

    vector<vector<double> > a = {{10, 0, 2, 4},
                                 {2, 16, -3, 8},
                                 {1, 5, 11, -4},
                                 {8, 1, 6, -17}};
    vector<double> y = {110, 128, 102, 81};
    cout << "ITERATIONS METHOD:\n";
    vector<double> x = iter_method(a, y);
    for(int i = 0; i < x.size(); i++) cout << x[i] << " ";
    cout << "\nTesting...\n";
    for (int i = 0; i < x.size(); i++) {
      double sum = 0.0;
      for (int j = 0; j < x.size(); j++) {
        sum += a[i][j] * x[j];
      }
      cout << y[i] - sum << endl;
    }
    cout << '\n' << "SEIDEL METHOD:\n";
    vector<double> x1 = seidel_method(a, y);
    for(int i = 0; i < x1.size(); i++) cout << x[i] << " ";
    cout << "\nTesting...\n";
    for (int i = 0; i < x.size(); i++) {
      double sum = 0.0;
      for (int j = 0; j < x.size(); j++) {
        sum += a[i][j] * x[j];
      }
      cout << y[i] - sum << endl;
    }
    cout << endl;
    return 0;
}
