from sympy import *
from sympy.core.rules import Transform

def roundDP(self, precision):
    """
    Rounds each element of the matrix to the specified number of decimal places.

    Parameters:
        self (MutableDenseMatrix): The matrix to round.
        precision (int): Number of decimal places to round to.

    Returns:
        MutableDenseMatrix: A new matrix with rounded elements.
    """
    return self.applyfunc(lambda e: round(e, precision))

def compute_centroid(points):
    """
    Compute the centroid of a set of 3D points.

    Args:
        points (Matrix): A matrix where each row represents a point in 3D space.

    Returns:
        Matrix: A column matrix representing the centroid coordinates (x, y, z).
    """
    # Calculate the total sum of x, y, and z coordinates
    sum_x = sum(points[0, :])
    sum_y = sum(points[1, :])
    sum_z = sum(points[2, :])

    # Compute the number of points
    num_points = len(Points[0,:])

    # Compute the centroid
    centroid_x = Add(sum_x, 0) / num_points
    centroid_y = Add(sum_y, 0) / num_points
    centroid_z = Add(sum_z, 0) / num_points

    centroid = Matrix([centroid_x, centroid_y, centroid_z])

    return centroid

def custom_evaluate_no_rounding(self):
    return self.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal), (xc, xcVal), (yc, ycVal), (dc, dcVal)])

def custom_evaluate(self):
    return roundDP(custom_evaluate_no_rounding(self),precision)

init_printing()
x, y, z, delta_x, delta_y, delta_z, sigma_x, sigma_y, sigma_z, xc, yc, dc = symbols('theta_x theta_y theta_z delta_x delta_y delta_z sigma_x sigma_y sigma_z theta_c phi_c d_c')
xVal = 0.1
yVal = 0.2
zVal = 0.1
delta_xVal = 0
delta_yVal = 1
delta_zVal = 1
sigma_xVal = 0.7
sigma_yVal = 0.7
sigma_zVal = 0.7
xcVal = 0.15
ycVal = 0.07
dcVal = 12
precision = 2

Points = Matrix([[-2,2,2,-2,-2,2,2,-2,-2,2], [-1,-1,1,1,-1,-1,1,1,-0.5,-0.5], [3,3,3,3,-3,-3,-3,-3,-1,-1], [1,1,1,1,1,1,1,1,1,1]])

# Question 1. a) Translation Matrix
T =  Matrix([[1,0,0,delta_x], [0,1,0,delta_y], [0,0,1,delta_z], [0,0,0,1]])
TEval = custom_evaluate(T)
print("Q1. a) T:")
print_latex(T)
print("-----\nQ1. a) T (Evaluated):")
print_latex(TEval)
print("--------------------------------------------------")

# Question 1. a) Scaling Matrix
S = Matrix([[sigma_x,0,0,0], [0,sigma_y,0,0], [0,0,sigma_z,0], [0,0,0,1]])
SEval = custom_evaluate(S)
print("Q1. a) S:")
print_latex(S)
print("-----\nQ1. a) S (Evaluated):")
print_latex(SEval)
print("--------------------------------------------------")

# Question 1. a) Rotation Matrix
Rx = Matrix([[1,0,0,0], [0,cos(x),-sin(x),0], [0,sin(x),cos(x),0], [0,0,0,1]])
Ry = Matrix([[cos(y),0,sin(y),0], [0,1,0,0], [-sin(y),0,cos(y),0], [0,0,0,1]])
Rz = Matrix([[cos(z),-sin(z),0,0], [sin(z),cos(z),0,0], [0,0,1,0], [0,0,0,1]])
Ryx = Ry * Rx
Rzyx = Rz * Ryx
RzyxEval = custom_evaluate(Rzyx)

print("Q1. a) Rx:")
print_latex(Rx)
print("-----\nQ1. a) Ry:")
print_latex(Ry)
print("-----\nQ1. a) Rz:")
print_latex(Rz)
print("-----\nQ1. a) Ryx:")
print_latex(Ryx)
print("-----\nQ1. a) Rzyx:")
print_latex(Rzyx)
print("-----\nQ1. a) Rzyx (Evaluated):")
print_latex(RzyxEval)
print("--------------------------------------------------")

# Question 1. c) Model Matrix
ST = S * T
Mmodel = Rzyx * ST
MmodelEval = custom_evaluate(Mmodel)

print("Q1. c) ST:")
print_latex(ST)
print("-----\nQ1. c) Mmodel:")
print_latex(Mmodel)
print("-----\nQ1. c) Mmodel (Evaluated):")
print_latex(MmodelEval)
print("--------------------------------------------------")

# Question 1. d) Model Transformed Points (P')
MPoints = (Mmodel * Points)
MPointsEval = custom_evaluate(MPoints)

print("Q1. d) Model Transformed Points (P'):")
print_latex(MPointsEval)
print("--------------------------------------------------")

# Question 2. a) Cartesian Points P1' to P5'
print("Q2. a) Cartesian Point P1'")
print_latex(MPointsEval[0:3,0])
print("-----\nQ2. a) Cartesian Point P2'")
print_latex(MPointsEval[0:3,1])
print("-----\nQ2. a) Cartesian Point P3'")
print_latex(MPointsEval[0:3,2])
print("-----\nQ2. a) Cartesian Point P4'")
print_latex(MPointsEval[0:3,3])
print("-----\nQ2. a) Cartesian Point P5'")
print_latex(MPointsEval[0:3,4])
print("--------------------------------------------------")

# Question 2. b) Vectors v1 (P2'-P1'), v2 (P4'-P1') and v3 (P5'-P1')
v1 = MPoints[0:3,1] - MPoints[0:3,0]
v1Eval = custom_evaluate(v1)
v2 = MPoints[0:3,3] - MPoints[0:3,0]
v2Eval = custom_evaluate(v2)
v3 = MPoints[0:3,4] - MPoints[0:3,0]
v3Eval = custom_evaluate(v3)

print("Q2. b) v1 (P2'-P1'):")
print_latex(v1Eval)
print("-----\nQ2. b) v2 (P4'-P1'):")
print_latex(v2Eval)
print("-----\nQ2. b) v3 (P5'-P1'):")
print_latex(v3Eval)
print("--------------------------------------------------")

# Question 2. c) Dot Products v1Dotv2, v1Dotv3 and v2Dotv3
v1Dotv2 = v1.dot(v2)
v1Dotv2Eval = round(custom_evaluate_no_rounding(v1Dotv2),precision)
v1Dotv3 = v1.dot(v3)
v1Dotv3Eval = round(custom_evaluate_no_rounding(v1Dotv3),precision)
v2Dotv3 = v2.dot(v3)
v2Dotv3Eval = round(custom_evaluate_no_rounding(v2Dotv3),precision)

print("Q2. c) v1 dot v2")
print_latex(v1Dotv2Eval)
print("-----\nQ2. c) v1 dot v3")
print_latex(v1Dotv3Eval)
print("-----\nQ2. c) v2 dot v3")
print_latex(v2Dotv3Eval)
print("--------------------------------------------------")

# Question 2. d) Showing P1', P2', P3' and P4' are coplanar
v1Crossv2 = v1.cross(v2) # This is the normal.
v1Crossv2Eval = custom_evaluate(v1Crossv2)
rx, ry, rz = symbols('x y z')
r = Matrix([rx, ry, rz])
planeEquation = (r.T * v1Crossv2)[0]
# Round all numbers in sympy object - smichr - https://stackoverflow.com/questions/56280007/round-all-numbers-in-sympy-object - Accessed 02.03.2024
planeEquationEval = custom_evaluate_no_rounding(planeEquation).xreplace(Transform(lambda x: x.round(precision), lambda x: isinstance(x, Float)))
planeEquationSubP1 = round(custom_evaluate_no_rounding(planeEquation.subs([(rx, MPoints[0,0]), (ry, MPoints[1,0]), (rz, MPoints[2,0])])), precision)
planeEquationSubP3 = round(custom_evaluate_no_rounding(planeEquation.subs([(rx, MPoints[0,2]), (ry, MPoints[1,2]), (rz, MPoints[2,2])])), precision)

print("Q2. d) v1 cross v2 (Normal):")
print_latex(v1Crossv2Eval)
print("-----\nQ2. d) Plane Equation:")
print_latex(planeEquationEval)
print("-----\nQ2. d) Plane Equation (Substituting P1'):")
print_latex(planeEquationSubP1)
print("-----\nQ2. d) Plane Equation (Substituting P3'):")
print_latex(planeEquationSubP3)
print("--------------------------------------------------")

# Question 3. a) Rotation Matrix
Rxc = Matrix([[1,0,0,0], [0,cos(xc),-sin(xc),0], [0,sin(xc),cos(xc),0], [0,0,0,1]])
Ryc = Matrix([[cos(yc),0,sin(yc),0], [0,1,0,0], [-sin(yc),0,cos(yc),0], [0,0,0,1]])
Ryxc = Ryc * Rxc
RyxcEval = custom_evaluate(Ryxc)

print("Q3. a) Rxc:")
print_latex(Rxc)
print("-----\nQ3. a) Ryc:")
print_latex(Ryc)
print("-----\nQ3. a) Ryxc:")
print_latex(Ryxc)
print("-----\nQ3. a) Ryxc (Evaluated):")
print_latex(RyxcEval)
print("--------------------------------------------------")

# Question 3. a) Translation Matrix
Tc = Matrix([[1,0,0,0], [0,1,0,0], [0,0,1,dc], [0,0,0,1]])
TcEval = custom_evaluate(Tc)

print("Q3. a) Tc:")
print_latex(Tc)
print("-----\nQ3. a) Tc (Evaluated):")
print_latex(TcEval)
print("--------------------------------------------------")

# Question 3. a) Camera Matrix
Mcamera = Ryxc * Tc
McameraEval = custom_evaluate(Mcamera)

print("Q3. a) Camera Matrix:")
print_latex(Mcamera)
print("-----\nQ3. a) Camera Matrix (Evaluated):")
print_latex(McameraEval)
print("--------------------------------------------------")

# Question 3. b) View Matrix
Mview = Mcamera.inv()
MviewEval = custom_evaluate(Mview)

print("Q3. b) View Matrix (Evaluated):")
print_latex(MviewEval)
print("--------------------------------------------------")

# Question 3. c) Model-View Matrix
Mmodelview = Mview * Mmodel
MmodelviewEval = custom_evaluate(Mmodelview)

print("Q3. c) Model-View Matrix (Evaluated):")
print_latex(MmodelviewEval)
print("--------------------------------------------------")

# Question 3. d) Model-View Transformed Points (P'')
MVPoints = (Mmodelview * Points)
MVPointsEval = custom_evaluate(MVPoints)

print("Q3. d) Model-View Transformed Points (P''):")
print_latex(MVPointsEval)
print("--------------------------------------------------")

# Question 4.1 In-place Scaler
MmodelviewInverse = Mmodel.inv() * Mcamera
DoubleXAtOrigin =  Matrix([[2,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])
InPlaceScaler = Mmodelview * DoubleXAtOrigin * MmodelviewInverse
InPlaceScalerEval = custom_evaluate(InPlaceScaler)

print("Q4.1 In-place Scaler:")
print_latex(InPlaceScalerEval)
print("--------------------------------------------------")

# Question 4.2 Scaled Model-View Transformed Points (P''')
SMVPoints = InPlaceScaler * MVPoints
SMVPointsEval = custom_evaluate(SMVPoints)

print("Q4.2 Scaled Model-View Transformed Points (P'''):")
print_latex(SMVPointsEval)
print("--------------------------------------------------")
print("Correct scaling verification:")
print(f"P1: {Points[:,0]}. P4: {Points[:,3]}")
print(f"P1P4Before: {round((Points[:,3] - Points[:,0]).norm(),precision)}. P1P4Model: {round((MPointsEval[:,3] - MPointsEval[:,0]).norm(),precision)}. P1P4ModelView: {round((MVPointsEval[:,3] - MVPointsEval[:,0]).norm(),precision)}. P1P4Scaled: {round((SMVPointsEval[:,3] - SMVPointsEval[:,0]).norm(),precision)}.")
print(f"P1P2Before: {round((Points[:,1] - Points[:,0]).norm(),precision)}. P1P2Model: {round((MPointsEval[:,1] - MPointsEval[:,0]).norm(),precision)}. P1P2ModelView: {round((MVPointsEval[:,1] - MVPointsEval[:,0]).norm(),precision)}. P1P2Scaled: {round((SMVPointsEval[:,1] - SMVPointsEval[:,0]).norm(),precision)}.")
print("--------------------------------------------------")