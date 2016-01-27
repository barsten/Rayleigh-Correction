import sys
import math
import numpy as np
import readfile
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

rayADF = readfile.readRayADF()
for i in range(12): print("%i3 12%5.3f\n" % (i, rayADF['RayScattCoeffA'][i, :, 0]))
   # for j in range(12): print(i, j, rayADF['RayScattCoeffA'][i,j,0])

xvalues = rayADF['theta']
yvalues = rayADF['theta']
X, Y = np.meshgrid(xvalues, yvalues)


for i in range(3):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    zvalues = rayADF['RayScattCoeffA'][:, :, i]
    Z = zvalues.reshape(X.shape)
    ax.plot_surface(X, Y, Z)

    ax.set_xlabel('theta sun')
    ax.set_ylabel('theta view')
    ax.set_zlabel('RayScattCoeffA[:, :, i]')

    plt.show()

