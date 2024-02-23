from sympy import *
init_printing()
x, y, z, delta_x, delta_y, delta_z, sigma_x, sigma_y, sigma_z = symbols('x y z delta_x delta_y delta_z sigma_x sigma_y sigma_z')
Points = Matrix([[-2,2,2,-2,-2,2,2,-2,-2,2], [-1,-1,1,1,-1,-1,1,1,-0.5,-0.5], [3,3,3,3,-3,-3,-3,-3,-1,-1], [1,1,1,1,1,1,1,1,1,1]])
T =  Matrix([[1,0,0,delta_x], [0,1,0,delta_y], [0,0,1,delta_z], [0,0,0,1]])
S =  Matrix([[sigma_x,0,0,0], [0,sigma_y,0,0], [0,0,sigma_z,0], [0,0,0,1]])
Rx = Matrix([[1,0,0,0], [0,cos(x),-sin(x),0], [0,sin(x),cos(x),0], [0,0,0,1]])
Ry = Matrix([[cos(y),0,sin(y),0], [0,1,0,0], [-sin(y),0,cos(y),0], [0,0,0,1]])
Rz = Matrix([[cos(z),-sin(z),0,0], [sin(z),cos(z),0,0], [0,0,1,0], [0,0,0,1]])
Ryx = Ry * Rx
Rzyx = Rz * Ryx
RzyxEval = Rzyx.subs([(x, 0.1), (y, 0.2), (z, 0.1)]).evalf(3)
M = T * Rzyx * S
MEval = M.subs([(x, 0.1), (y, 0.2), (z, 0.1), (delta_x, 0), (delta_y, 1), (delta_z, 1), (sigma_x, 0.7), (sigma_y, 0.7), (sigma_z, 0.7)]).evalf(3)
TransformedPoints = (M * Points).subs([(x, 0.1), (y, 0.2), (z, 0.1), (delta_x, 0), (delta_y, 1), (delta_z, 1), (sigma_x, 0.7), (sigma_y, 0.7), (sigma_z, 0.7)]).evalf(3)