import numpy as np
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy.interpolate import interpn
# x = np.linspace(0, 10, num=11, endpoint=True)
# y = np.cos(-x**2/9.0)
# f = interp1d(x, y)
# f2 = interp1d(x, y, kind='cubic')
# xnew = np.linspace(0, 10, num=41, endpoint=True)
#
# import matplotlib.pyplot as plt
# plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
# plt.legend(['data', 'linear', 'cubic'], loc='best')
# plt.show()

# Multi-dimensional
# x = np.zeros((2, 2))
# x[0, 1] = 1
# x[1, 1] = 1
# print x
# y = np.array(([ 5.222, 6.916], [6.499, 4.102]))
# print y
# xi = np.array((0.098, 0.08))
#
# print(interpn(x, y, xi))

# 2-D
# dim1 = 2
# dim2 = 4
# x = np.zeros((dim1, dim2))
# # x = np.zeros((2, 3))
# # x[:,1] = 1
# # x[:,2] = 2
# for i in range(dim2):
#     x[:, i] = i
# print 'Achsen: \n', x
# # y = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], np.float32)
# n = 1
# y = np.zeros((dim2, dim2))
# for i in range(dim2):
#     for j in range(dim2):
#         y[(i, j)] = n
#         n += 1
#     # print 'Werte: \n', y
# print 'Werte: \n', y
#
# xi = np.array( (0.5, 1.) )
# print 'Wert am Punkt ', xi
# print(interpn(x, y, xi, method='linear'))

# 3-D
# axis1 = np.array([0., 1., 2.])
# dim1 = axis1.size
# axis2 = np.array([0., 1.])
# dim2 = axis2.size
# axis3 = np.array([0., 2., 4., 6.])
# dim3 = axis3.size
# x = [axis1, axis2, axis3]
# print 'Achsen: \n', x
#
# y = np.zeros((dim1, dim2, dim3))
# n = 1
# for i in range(dim1):
#     for j in range(dim2):
#         for k in range(dim3):
#             y[(i, j, k)] = n
#             n += 1
#
# print 'Werte: \n', y
#
# xi = np.array( [(0.5, 0.5, 0.), (1.5, 0., 0.8)] )
# print 'Wert am Punkt ', xi
# print(interpn(x, y, xi, method='linear'))

# 4-D
axis1 = np.array([0., 1., 2.])
dim1 = axis1.size
axis2 = np.array([0., 1.])
dim2 = axis2.size
axis3 = np.array([0., 2., 4., 6.])
dim3 = axis3.size
axis4 = np.array([0., 1., 2., 3., 4., 5.])
dim4 = axis4.size
x = [axis1, axis2, axis3, axis4]
print 'Achsen: \n', x

y = np.zeros((dim1, dim2, dim3, dim4))
n = 1
for i in range(dim1):
    for j in range(dim2):
        for k in range(dim3):
            for l in range(dim4):
                y[(i, j, k, l)] = n
                n += 1

print 'Werte: \n', y

xi = np.array( [(0.5, 0.5, 0., 3.5), (1.5, 0., 0.8, 0.7)] )
print 'Wert am Punkt ', xi
print(interpn(x, y, xi, method='linear'))