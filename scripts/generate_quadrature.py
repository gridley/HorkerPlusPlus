#!/usr/bin/env python3
import quadpy as qp
import matplotlib.pyplot as plt
import sys
import numpy as np

frmt = '''
struct %s : Quadrature<%i>
{
  static constexpr std::array<Position, Npts> points = {%s};
  static constexpr std::array<double, Npts> weights = {%s};
};
'''

# This is a particularly nice quadrature: order 15, and not an unreasonable amount of points.
# rule = qp.quadrilateral.rabinowitz_richter_6()
rule = qp.quadrilateral.sommariva_16()
these_points = rule.points
these_weights = rule.weights
these_points += np.array([1.0, 1.0]) # shift so origin is at bottom left  corner
these_points /= 2.0 # scale so top right corner is at (1,1)
these_weights /= 4.0 # area should be one after scaling

cname = rule.name.replace('-','')
cname = cname.replace(' ','')
points = ['constpos(%f, %f, 0.0),'%(x, y) for x,y in zip(these_points[:,0], these_points[:,1])]
weights = [str(s) for s in these_weights]

pointstring = '\n    '.join(points)
weightstr = ',\n    '.join(weights)

result = frmt%(cname, len(points), pointstring, weightstr)
print(result)
