import numpy as np
import math as m
import sys
#A = np.zeros((3, 3))

np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

def abs(a):
    if a < 0.0:
        a *= -1.0
        return a
    else:
        return a

def find_max(A):
    max_ind = [0, 1]
    (n, m) = A.shape
    for i in range(0, n):
        for j in range(i + 1, m):
            k = max_ind[0]
            l = max_ind[1]
            val1 = A[k][l]
            val2 = A[i][j]
            if  abs(val2) > abs(val1):
                max_ind = [i, j]
    return max_ind

#input_mtrx = input("enter matrix: ")
#main script
A = np.array([[5., 1., 2.], [1., 4., 1.], [2., 1., 3.]], np.float32)
print(A)
eps = 0.001
k = 0
(n, m) = A.shape
found_L = False
eig_v = np.identity(n, np.float)
while not found_L:
    [i, j] = find_max(A)
    if abs(A[i][j]) < eps:
        found_L = True
        break
    #calculation of rotation's angle phi
    phi = 0.5 * np.arctan(2.0 * A[i][j] / (A[i][i] - A[j][j]))
    #jacobi matrix H creation for a k step
    H = np.zeros((n, m), np.float32)
    for l in range(0, n):
        H[l][l] = 1.0
    H[i][i] = H[j][j] =  np.cos(phi)
    H[i][j] = -1.0 * np.sin(phi)
    H[j][i] = np.sin(phi)
    eig_v = np.dot(eig_v, H)
    #find new A matrix
    A = np.dot(np.dot(H.T, A), H)
    k = k + 1
for i in range(0, n):
    eig_v[:][i] /= eig_v[i][i]
print("eigenvectors:\n{0}".format(eig_v))
print("eigenvalues\n{0}".format(A))
print("calculated in {0} iterations".format(k))
