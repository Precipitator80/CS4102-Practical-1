from sympy import *

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

init_printing()
x, y, z, delta_x, delta_y, delta_z, sigma_x, sigma_y, sigma_z, xc, yc, dc = symbols('x y z delta_x delta_y delta_z sigma_x sigma_y sigma_z x y d_c')
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

T =  Matrix([[1,0,0,delta_x], [0,1,0,delta_y], [0,0,1,delta_z], [0,0,0,1]])

S =  Matrix([[sigma_x,0,0,0], [0,sigma_y,0,0], [0,0,sigma_z,0], [0,0,0,1]])

Rx = Matrix([[1,0,0,0], [0,cos(x),-sin(x),0], [0,sin(x),cos(x),0], [0,0,0,1]])
Ry = Matrix([[cos(y),0,sin(y),0], [0,1,0,0], [-sin(y),0,cos(y),0], [0,0,0,1]])
Rz = Matrix([[cos(z),-sin(z),0,0], [sin(z),cos(z),0,0], [0,0,1,0], [0,0,0,1]])
Ryx = Ry * Rx
Rzyx = Rz * Ryx
RzyxEval = roundDP(Rzyx.subs([(x, xVal), (y, yVal), (z, zVal)]), precision)

M = T * Rzyx * S
MEval = roundDP(M.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal)]),precision)

TransformedPoints = (M * Points)
TransformedPointsEval = roundDP(TransformedPoints.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal)]),precision)
# TransformedPointsRoundingError = roundDP((T.subs([(delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal)]),precision) * roundDP(Rz.subs([(x, xVal), (y, yVal), (z, zVal)]),precision) * roundDP(Ry.subs([(x, xVal), (y, yVal), (z, zVal)]),precision) * roundDP(Rx.subs([(x, xVal), (y, yVal), (z, zVal)]),precision) * roundDP(S.subs([(sigma_x, sigma_xVal), (sigma_y,sigma_yVal), (sigma_z, sigma_zVal)]),precision)) * Points

v1 = TransformedPoints.col(1) - TransformedPoints.col(0)
v1Eval = roundDP(v1.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal)]),precision)
v2 = TransformedPoints.col(3) - TransformedPoints.col(0)
v2Eval = roundDP(v2.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal)]),precision)
v3 = TransformedPoints.col(4) - TransformedPoints.col(0)
v3Eval = roundDP(v3.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal)]),precision)

v1Dotv2 = v1.dot(v2)
v1Dotv2Eval = v1Dotv2.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal)]).round(precision)
v2Dotv3 = v2.dot(v3)
v2Dotv3Eval = v2Dotv3.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal)]).round(precision)

Rxc = Matrix([[1,0,0,0], [0,cos(xc),-sin(xc),0], [0,sin(xc),cos(xc),0], [0,0,0,1]])
Ryc = Matrix([[cos(yc),0,sin(yc),0], [0,1,0,0], [-sin(yc),0,cos(yc),0], [0,0,0,1]])
Ryxc = Ryc * Rxc
RyxcEval = roundDP(Ryxc.subs([(xc, xcVal), (yc, ycVal)]), precision)
Tc = Matrix([[1,0,0,0], [0,1,0,0], [0,0,1,dc], [0,0,0,1]])

Mcamera = Ryxc * Tc
McameraEval = roundDP(Mcamera.subs([(xc, xcVal), (yc, ycVal), (dc, dcVal)]), precision)

Mview = Mcamera.inv()
MviewEval = roundDP(Mview.subs([(xc, xcVal), (yc, ycVal), (dc, dcVal)]), precision)

Mvm = Mview * M
MvmEval = roundDP(Mvm.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal), (xc, xcVal), (yc, ycVal), (dc, dcVal)]),precision)
MvmEvalAlt = roundDP(MviewEval * MEval, precision)
# Double check the symbols work properly!

TransformedPointsMVM = (Mvm * Points)
TransformedPointsMVMEval = roundDP(TransformedPointsMVM.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal), (xc, xcVal), (yc, ycVal), (dc, dcVal)]),precision)
#TransformedPointsMVMEvalAlt = roundDP((Mview * TransformedPoints).subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal), (xc, xcVal), (yc, ycVal), (dc, dcVal)]),precision)
#TransformedPointsMVMEvalAlt2 = roundDP(MviewEval * TransformedPointsEval,precision)

OriginalCentroid = compute_centroid(Points)
TransformedMVMCentroid = compute_centroid(TransformedPointsMVMEval)
CentroidDifference = TransformedMVMCentroid - OriginalCentroid

ReturnMVMToOriginTransform = Matrix([[1,0,0,-CentroidDifference[0]], [0,1,0,-CentroidDifference[1]], [0,0,1,-CentroidDifference[2]], [0,0,0,1]])
ReturnMVMToIdentityRotation = (Ryxc * Rzyx).T
#DoubleXAtOrigin =  Matrix([[2,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])
#ReturnIdentityToMVMRotation = ReturnMVMToIdentityRotation.T
#ReturnOriginToMVMTransform = ReturnMVMToOriginTransform.inv()

#ScaledTransformedPointsMVM = ReturnOriginToMVMTransform * ReturnIdentityToMVMRotation * DoubleXAtOrigin * ReturnMVMToIdentityRotation * ReturnMVMToOriginTransform * TransformedPointsMVM
#ScaledTransformedPointsMVMEval = roundDP(ScaledTransformedPointsMVM.subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal), (xc, xcVal), (yc, ycVal), (dc, dcVal)]),precision)

#compute_centroid(ReturnMVMToOriginTransform * TransformedPointsEval)

TransformedPointsMVMAtOrigin = roundDP((ReturnMVMToIdentityRotation * ReturnMVMToOriginTransform * TransformedPointsMVM).subs([(x, xVal), (y, yVal), (z, zVal), (delta_x, delta_xVal), (delta_y, delta_yVal), (delta_z, delta_zVal), (sigma_x, sigma_xVal), (sigma_y, sigma_yVal), (sigma_z, sigma_zVal), (xc, xcVal), (yc, ycVal), (dc, dcVal)]),precision)
compute_centroid(TransformedPointsMVMAtOrigin)

# Remember print_latex()!