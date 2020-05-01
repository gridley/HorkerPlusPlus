#!/usr/bin/env python3
# Generates the shape functions for a quadrilateral on the reference element
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as lin
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os
import sys

# Max degree of the polynomial
d = int(sys.argv[1])

specialization = '''
template<>
class ReferenceElement<%i>
{
public:
  static constexpr unsigned Npts = %i;
  static constexpr unsigned NEdgeNodes = %i;
  static constexpr unsigned NInternalNodes = %i;
  static constexpr unsigned Degree = %i;
  static constexpr std::array<double, Npts> node_xpoints = {%s};
  static constexpr std::array<double, Npts> node_ypoints = {%s};
  static inline std::array<double, %i> evaluate_basis(double x, double y) {
    return %s
  }
  static inline std::array<double, %i> evaluate_deriv_x(double x, double y) {
    return %s
  }
  static inline std::array<double, %i> evaluate_deriv_y(double x, double y) {
    return %s
  }
};
'''

# Generate x and y point values
pt_dir = d+1 # points along one direction
npts = pt_dir**2

basis_matrix = np.zeros((npts, npts))

# Now determine where the points lie:
points = np.zeros((npts,2))

xvals = np.linspace(0, 1, d+1)

# always have corner points..
points[0,0]=0
points[0,1]=0
points[1,0]=1
points[1,1]=0
points[2,0]=1
points[2,1]=1
points[3,0]=0
points[3,1]=1

# Edge nodes come next
n_edge = pt_dir-2
points[4:4+n_edge,0] = xvals[1:-1]
points[4:4+n_edge,1] = 0
points[4+n_edge:4+2*n_edge,0] = 1
points[4+n_edge:4+2*n_edge,1] = xvals[1:-1]
points[4+2*n_edge:4+3*n_edge,0] = xvals[1:-1]
points[4+2*n_edge:4+3*n_edge,1] = 1
points[4+3*n_edge:4+4*n_edge,0] = 0
points[4+3*n_edge:4+4*n_edge,1] = xvals[1:-1]

# Finally, all internal nodes:
n = 4+4*n_edge
for i in range(n_edge):
    points[n:n+n_edge, 0] = xvals[1:-1]
    points[n:n+n_edge, 1] = xvals[1+i]
    n += n_edge

# Now generate the matrix, which, when inverted, will give the shape functions
def evaluate_basis(x, y):
    result = np.zeros(npts)
    k = 0
    for i in range(pt_dir):
        for j in range(pt_dir):
            result[k] = x**j * y**i
            k += 1
    return result

# Finally, create the basis matrix
for i in range(npts):
    row = evaluate_basis(points[i,0], points[i,1])
    basis_matrix[i,:] = evaluate_basis(points[i,0], points[i,1])
result = lin.inv(basis_matrix)


# make a fine mesh to plot the newly found basis on
n_fine = 200
fine_x = np.linspace(0, 1, n_fine)
fine_points = np.zeros((n_fine**2,3))
xx, yy = np.meshgrid(fine_x, fine_x)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

dirname = 'order_%i'%d
if os.path.exists(dirname):
    os.remove(dirname)
os.mkdir(dirname)

# Plot each basis function
for fcn in range(npts):
    # plt.axis('off')
    k=0
    print(fcn)
    for j in range(n_fine):
        for i in range(n_fine):
            fine_points[k,0] = fine_x[i]
            fine_points[k,1] = fine_x[j]
            fine_points[k,2] = np.dot(evaluate_basis(fine_x[i], fine_x[j]), result[:,fcn])
            k+=1

    # plt.pcolormesh(fine_x, fine_x, fine_points[:,2].reshape(n_fine,n_fine))
    arr = fine_points[:,2].reshape(n_fine,n_fine)
    ax.plot_surface(xx, yy, arr, cmap=cm.viridis)
    ax.plot(points[:,0], points[:,1], 'bs')
    plt.show()
    fig.savefig(os.path.join(dirname, 'basis_%i.png'%fcn))
    ax.clear()

# # Save the basis function coefficients
# np.savetxt(os.path.join(dirname, 'basis_coeffs'), result.T)
# np.savetxt(os.path.join(dirname, 'node_locations'), points)


######
# quit()
#######

# Generate C++ code for reference element basis evaluation

def monom(xi, yi):
    if xi==0 and yi==0:
        result='1'
    elif xi < 0 or yi < 0:
        result='1'
    else:
        yterm = 'y*' + 'y*'*(yi-1) if yi != 0 else ''
        xterm = 'x*' + 'x*' * (xi-1) if xi != 0 else ''
        result = xterm + yterm
    if result[-1] == '*':
        return result[:-1]
    else:
        return result

basis_eval = ['{']
der_x_eval = ['{']
der_y_eval = ['{']
for b in range(npts):
    k = 0
    thisterm = []
    thisx = []
    thisy = []
    for i in range(pt_dir):
        for j in range(pt_dir):
            sign = '' if k==0 else ('-' if result[k,b]<0 else '+')
            if 1-1e-9 < abs(result[k,b]) < 1+1e-9:
                thisterm.append(sign+monom(j, i))
                # Handle derivatives
                if j>0 and i==0:
                    thisx.append(sign+str(1.0))
                elif i>0 and j==0:
                    thisy.append(sign+str(1.0))
                else:
                    thisx.append(sign+str(j)+'*'+monom(j-1,i))
                    thisy.append(sign+str(i)+'*'+monom(j,i-1))
            elif abs(result[k,b]) > 1e-9:
                thisterm.append(sign + str(abs(result[k,b]))+'*' + monom(j, i))
                thisx.append(sign+str(j*abs(result[k,b]))+'*' + monom(j-1, i))
                thisy.append(sign+str(i*abs(result[k,b]))+'*' + monom(j, i-1))
            k += 1
    com = '};' if b+1 == npts else ','
    basis_eval.append(' '*4+''.join(thisterm) + com)
    der_x_eval.append(' '*4+''.join(thisx) + com)
    der_y_eval.append(' '*4+''.join(thisy) + com)
basis = '\n'.join(basis_eval)
derx = '\n'.join(der_x_eval)
dery = '\n'.join(der_y_eval)

xpts = ','.join([str(s) for s in points[:,0]])
ypts = ','.join([str(s) for s in points[:,1]])
this_spec = specialization%(d, npts, d-1, (d-1)**2, d, xpts, ypts, npts, basis, npts, derx, npts, dery)
print(this_spec)
