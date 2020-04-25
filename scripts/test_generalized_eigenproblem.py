#!/usr/bin/env python3
import scipy.linalg as lin
import numpy as np
import sys
lhs = np.loadtxt(sys.argv[1])
rhs = np.loadtxt(sys.argv[2])
eigs = lin.eig(lhs, rhs)
eigs1 = lin.eig(lhs)
print('general problem')
print(eigs[0])
print('lhs eigs only')
print(eigs1[0])
