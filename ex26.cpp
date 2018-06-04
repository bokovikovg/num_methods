#include <iostream>
//#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <vector>


using namespace std;

vector<double> find_max(vector<vector<double>> &A) {
      int n = (int) A[0].size();
      double max = fabs(A[0][1]);
      vector<double> res;
      for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++) {
              if (i == j) continue;
              if (fabs(A[i][j]) > max)  {
                  res.push_back(i);
                  res.push_back(j);
              }
          }
      return res;
}

double t(vector<vector<double>> A) {
    double sum = 0.0;
    int n = (int) A[0].size();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            sum += A[i][j] * A[i][j];
    return sqrt(sum);
}

vector<vector<double>> transpose(vector<vector<double>> A) {
    int n = (int) A[0].size();
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double temp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = temp;
        }
    }

    return A;
}

vector<vector<double>> fillU(vector<vector<double>> &U, int l, int m) {
    int n = (int) U[0].size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((i == l && j == m) || (i == m && j == l)
                || (i == l && j == l) ||  (i == m && j == m)) continue;
            if (i == j) U[i][j] = 1.0;
            else U[i][j] = 0.0;
        }
    }
    return U;
}

vector<vector<double>> mult(vector<vector<double>> &A, vector<vector<double>> &B) {
    const int n = (int) A.size(), m = (int) A[0].size(), l = (int) B[0].size();
    vector<vector<double>> C(n, vector<double>(l, 0));
    for (int j = 0; j < l; j++){
        for (int k = 0; k < m; k++) {
            for (int i = 0; i < n; i++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<vector<double>> jacoby_method (vector<vector<double>> A, double eps = 0.01) {
    int k = 0, n = (int) A[0].size();
    vector<double> max = find_max(A);
    vector<vector<double>> U(n, vector<double>(n));
    int i = max[0], j = max[1];
    double phi = 0.5 * atan(2.0*A[i][j]/(A[i][i] - A[j][j]));
    U[i][j] = - sin(phi); U[j][i] = sin(phi); U[i][i] = U[j][j] = cos(phi);
    U = fillU(U, i, j);
    vector<vector<double>> S = transpose(U);
    vector<vector<double>> T = mult(S, A);
    A = mult(T, U);
    while (t(A) > eps && k < 200) {
        vector<double> max = find_max(A);
        vector<vector<double>> U(n, vector<double>(n));
        int i = max[0], j = max[1];
        double phi = 0.5 * atan(2.0*A[i][j]/(A[i][i] - A[j][j]));
        U[i][j] = - sin(phi); U[j][i] = sin(phi); U[i][i] = U[j][j] = cos(phi);
        U = fillU(U, i, j);
        vector<vector<double>> S = transpose(U);
        vector<vector<double>> T = mult(S, A);
        A = mult(T, U);
        k++;
    }
    return A;
}

void print(vector<vector<double>> &A) {
    int n = (int) A.size(), m = A[0].size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    vector<vector<double>> A = { {3, 2, 6},
                                 {2, -3, -7},
                                 {6, -7, 3} };
    vector<vector<double>> L(3, vector<double>(3, 0));
    L = jacoby_method(A);
    print(L);
    return 0;
}
