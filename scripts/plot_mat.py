#!/usr/bin/env python3
# Plots a matrix sparsity pattern, given a numpy-readable file representing a matrix
import numpy as np
import matplotlib.pyplot as plt
import sys
data = np.loadtxt(sys.argv[1])
data[data==0] = np.nan # white where matrix is zero
plt.imshow(data)
plt.show()
