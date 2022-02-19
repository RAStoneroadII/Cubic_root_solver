import numpy as np
import cmath
a = float(input("Enter the x-cubed coeffecient:"))
b = float(input("Enter the x-squared coeffecient:"))
c = float(input("Enter the x coeffecient:"))
d = float(input("Enter the constant's coeffecient:"))
coefs = np.array([d,c,b,a])
def cubRtsFn(coefs):
	#Find e and f for y**3 + ey + f = 0
	eterm = (1/a) * (c - (b **2 ) / ( 3 * a)) 
	fterm = (1/a) * (d + ((2 * b ** 3) / (27 * a **2)) - ((b * c) / (3 * a)))
	sterm = (-eterm) / 3 #Derived from Vieta's sub s.t. s = -e/3
	#Changes the first equation to z**6 + fz - (e**3)/27 = 0
	#let w = z**3 s.t. the equation is now w**2 + fw - (e**3)/27 
	ebasedconstant = sterm ** 3
	determinant = cmath.sqrt((fterm ** 2) - (4 * ebasedconstant)) 
	w = (-fterm + determinant) / 2 #Either value of w can find the roots so only one needs to be solved for
	#Solve for the roots using w = z**3: (w**(1/3)) * (cos(alpha) + j*sin(alpha) for alpha = (theta + 2pik)/3 for k = 0,1,2
	theta = cmath.phase(w) #This is to find theta of the equation
	a1 = theta / 3
	a2 = (theta + 2 * np.math.pi) / 3
	a3 = (theta + 4 * np.math.pi) / 3
	z1 = (w ** (1/3)) * complex(cmath.cos(a1), cmath.sin(a1))
	z2 = (w ** (1/3)) * complex(cmath.cos(a2), cmath.sin(a2))
	z3 = (w ** (1/3)) * complex(cmath.cos(a3), cmath.sin(a3))
	#Backsubstitution to solve for y  = z + (s/z)
	y1 = z1 + (sterm / z1)
	y2 = z2 + (sterm / z2)
	y3 = z3 + (sterm / z3)
	#Backsubstitute again to solve for x = (y - b) / (3 * a)
	x1 = y1 - b / (3 * a)
	x2 = y2 - b / (3 * a)
	x3 = y3 - b / (3 * a)
	xvalues = [x1, x2, x3]
	return xvalues
print(cubRtsFn(coefs))