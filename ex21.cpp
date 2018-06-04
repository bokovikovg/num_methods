#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void outputMatrix(vector<vector<double> > &a, vector<double> &y) {
    cout << "\nMatrix A\n";
    for(int i = 0; i < y.size(); i++) {
        for (int j = 0; j < y.size(); j++) cout << a[i][j] << " ";
        cout << endl;
    }
    cout << "Vector Y\n";
    for (int i = 0; i < y.size(); i++) cout << y[i] << " ";
}

vector<double> gauss(vector<vector<double> > &a, vector<double> &y) {
    
    int n = y.size();
    vector<double> x(n);
    double max;
    int k = 0, index;
    const double eps = 0.000001;
    while (k < n) {
        //line with max a[i][j]
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > max) {
                max = abs(a[i][k]);
                index = i;
            }
        }
        if (max < eps) {
            cout << "column w/ 0s\n";
            cout << index << "matrix A\n";
            return x;
        }
        //line permutation
        for (int j = 0; j < n; j++) {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[k] = temp;
        for (int i = k; i < n; i++) {
            double temp = a[i][k];
            if (abs(temp) < eps) continue;
            for (int j = 0; j < n; j++) a[i][j] /= temp;
            y[i] /= temp;
            if (i == k) continue;
            for (int j = 0; j < n; j++) a[i][j] -= a[k][j];
            y[i] -= y[k];
        }
        outputMatrix(a, y);
        k++;
    }
    //reverse
    
    for (k = n - 1; k >= 0; k--) {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    outputMatrix(a, y);
    return x;
}


int main() {
    
    vector<vector<double> > a = { {9, -7, -1, 1},
        {2, 7, 3, -6},
        {4, 7, -3, -7},
        {-9, -5, -1, -6 } };
    
    vector<double> y = {55, -66, -43, -24};
    outputMatrix(a, y);
    vector<double> x = gauss(a, y);
    outputMatrix(a, y);
    for (int i = 0; i < x.size(); i++) cout << x[i] << " ";
    cout << endl;
    return 0;
}



