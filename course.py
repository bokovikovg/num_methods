import math as m
import sys
import numpy as np
import matplotlib.pyplot as plt

eps = 0.0001

#integrate func f(x) on an interval [a, b] with eps precision
def integrate(f, a, b, eps):
    n = 2
    sum = 0.0
    sum2 = 0.0
    while True:
        h = (b - a) / n
        k1 = 0.0
        k2 = 0.0
        for i in range(1, n, 2):
            k1 += f(a + i * h)
            k2 += f(a + (i + 1) * h)
        sum = h/3 * (f(a) + 4*k1 + 2*k2)
        if m.fabs(sum - sum2) < eps:
            break
        sum2 = sum
        n *= 2
    return sum

#declaring functions
f1 = lambda t: 2.0 * (t + 1)
f2 = lambda t: (8 - t)

def f(t):
    if t >= 0.0 and t <= 2.0:
        return f1(t)/m.log(t + 2)
    if t >= 2.0 and t <= 6.0:
        return f2(t)/m.log(t + 2)
    else:
        return 0.0;

k = integrate(f, 0, 6, eps)

def F(x):
    return k * m.exp(-x) * m.sin(2 * x)

#find a min of f(x) on [a, b]

def find_min(f, a, b, eps):
    phi = 0.5 * (1.0 + m.sqrt(5.0))
    while True:
        x1 = b - (b - a)/phi
        x2 = a + (b - a)/phi
        y1 = f(x1)
        y2 = f(x2)
        if y1 >= y2:
            a = x1
        else:
            b = x2
        if m.fabs(b - a) < eps:
            min_x = (a + b)/2
            break
    return min_x

#main script
a = -2.0
b = 0.0
pt = find_min(F, a, b, eps)
print("minimum of given function is at x = {0}".format(pt))

step = 0.001
x = []
x.append(a)
i = 0
while (x[i] <= b):
    x.append(x[i] + step)
    i = i + 1
y = []
for x0 in x:
    y.append(F(x0))
fig = plt.figure()
plt.plot(x, y)
#plt.text(pt, F(pt), 7, 'minimum of F', fontsize=12, bbox=dict(edgecolor='w', color='w'))
plt.grid(True)
plt.show()
