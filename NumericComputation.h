#include "stdafx.h"
#include <math.h>
#include <iostream>
using namespace std;

#ifndef NUMERICCOMPUTATION
#define NUMERICCOMPUTATION


#define COMPUTATION_SMALL_NUMBER 1.0e-5
#define COMPUTATION_DOUBLE_MAX 1.0e+299


enum  Sign
{
	NEGATIVE = -1, ZERO = 0, POSITIVE = 1,

	// Orientation constants:
	RIGHT_TURN = -1, LEFT_TURN = 1,

	CLOCKWISE = -1, COUNTERCLOCKWISE = 1,

	COLLINEAR = 0, COPLANAR = 0, DEGENERATE = 0,

	// Oriented_side constants:
	ON_NEGATIVE_SIDE = -1, ON_ORIENTED_BOUNDARY = 0, ON_POSITIVE_SIDE = 1,

	// Comparison_result constants:
	SMALLER = -1, EQUAL = 0, LARGER = 1
};

typedef Sign Comparison_result;

template<typename T>
class NumericComputation
{
public:

	
	static Comparison_result compare_xyz(const T*p, const T*q) ;
	static Comparison_result compare(const T x, const T y);

	static Comparison_result cmp_dist_to_pointC3(const T &px, const T &py, const T &pz,
		const T &qx, const T &qy, const T &qz,
		const T &rx, const T &ry, const T &rz)
	{
		return  compare(squared_distanceC3(px, py, pz, qx, qy, qz),
			squared_distanceC3(px, py, pz, rx, ry, rz));
	}
	
	static T squared_distanceC3(const T &px, const T &py, const T &pz,
		const T &qx, const T &qy, const T &qz)
	{
		return  (px - qx)*(px-qx) +  (py-qy)*(py - qy) +
			 (pz-qz)*(pz - qz);
	}
	// Description:
	// Compute determinant of 2x2 matrix. Two columns of matrix are input.
	static const T Determinant2x2(const T c1[2], const T c2[2]) {
		return (c1[0] * c2[1] - c2[0] * c1[1]);
	};
	static const T Determinant2x2(const T a, const T b, const T c, const T d) {
		return (a * d - b * c);
	};
	// Description:
	// Given a 3D point x[3], determine the barycentric coordinates of the point.
	// Barycentric coordinates are a natural coordinate system for simplices that
	// express a position as a linear combination of the vertices. For a
	// tetrahedron, there are four barycentric coordinates (because there are
	// four vertices), and the sum of the coordinates must equal 1. If a
	// point x is inside a simplex, then all four coordinates will be strictly
	// positive.  If three coordinates are zero (so the fourth =1), then the
	// point x is on a vertex. If two coordinates are zero, the point x is on an
	// edge (and so on). In this method, you must specify the vertex coordinates
	// x1->x4. Returns 0 if tetrahedron is degenerate.
	//计算点的空间质心，用于查找包含四面体，判断可见性
	static int BarycentricCoords(double x[3], double  x1[3], double x2[3],
		double x3[3], double x4[3], double bcoords[4]);
	// Description:
	// Solve linear equations Ax = b using Crout's method. Input is square
	// matrix A and load vector x. Solution x is written over load vector. The
	// dimension of the matrix is specified in size. If error is found, method
	// returns a 0.
	static int SolveLinearSystem(const T **A, const T *x, int size);
	static void LUSolveLinearSystem(const T **A, int *index,
		const T *x, int size);
	//----------------------------------------------------------------------------
	// Factor linear equations Ax = b using LU decompostion A = LU where L is
	// lower triangular matrix and U is upper triangular matrix. Input is
	// square matrix A, integer array of pivot indices index[0->n-1], and size
	// of square matrix n. Output factorization LU is in matrix A. If error is
	// found, method returns 0.
	static int LUFactorLinearSystem(const T **A, int *index, int size);

	//------------------------------------------------------------------

	// Description:
	// Dot product of two 3-vectors (const T-precision version).
	static const T Dot(const T x[3], const T y[3]) {
		return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
	};
	//计算两个三维向量的向量积（叉乘）
	static T* Cross(const T a[3], const T b[3])
	{
		T* c=new T[3];
		c[0]=a[1]*b[2]-a[2]*b[1];
		c[1]=a[2]*b[0]-a[0]*b[2];
		c[2]=a[0]*b[1]-a[1]*b[0];
		return c;
	}
	//计算两点间距离的平方
	static const T SquareDistance(const T x[3], const T y[3])
	{
		return ((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]));
	};
	static const T Determinant3x3(const T A[3][3]);
	static Sign SignOfDeterminant2x2(const T x1, const T y1, const T x2, const T y2)
	{
		if (long double(x1*y2 - y1*x2) > COMPUTATION_SMALL_NUMBER)
			return POSITIVE;
		else if (long double( y1*x2-x1*y2) > COMPUTATION_SMALL_NUMBER)
			return NEGATIVE;
		return ZERO;
	};
	static const T Determinant3x3(const T& a00, const T& a01, const T& a02,
		const T& a10, const T& a11, const T& a12,
		const T& a20, const T& a21, const T& a22);
	static const T determinant4x4(
		const T& a00, const T& a01, const T& a02, const T& a03,
		const T& a10, const T& a11, const T& a12, const T& a13,
		const T& a20, const T& a21, const T& a22, const T& a23,
		const T& a30, const T& a31, const T& a32, const T& a33);
	static Sign SignOfDeterminant4x4(
		const T& a00, const T& a01, const T& a02, const T& a03,
		const T& a10, const T& a11, const T& a12, const T& a13,
		const T& a20, const T& a21, const T& a22, const T& a23,
		const T& a30, const T& a31, const T& a32, const T& a33)
	{
		long double det = determinant4x4(
			 a00,  a01,  a02,  a03,
			 a10,  a11,  a12,  a13,
			 a20,  a21,  a22,  a23,
			 a30,  a31,  a32,  a33);
		if (det > COMPUTATION_SMALL_NUMBER)
			return POSITIVE;
		else if (-det > COMPUTATION_SMALL_NUMBER)
			return NEGATIVE;
		return ZERO;
	};
		static Sign SignOfDeterminant3x3(
		const T& a00, const T& a01, const T& a02,
		const T& a10, const T& a11, const T& a12,
		const T& a20, const T& a21, const T& a22)
	{
		long double det = Determinant3x3(
			 a00,  a01,  a02,
			 a10,  a11,  a12, 
			 a20,  a21,  a22);
		if (det > COMPUTATION_SMALL_NUMBER)
			return POSITIVE;
		else if (-det > COMPUTATION_SMALL_NUMBER)
			return NEGATIVE;
		return ZERO;
	};
	
};


template<typename T>
const T NumericComputation<T>::
determinant4x4(
const T& a00, const T& a01, const T& a02, const T& a03,
const T& a10, const T& a11, const T& a12, const T& a13,
const T& a20, const T& a21, const T& a22, const T& a23,
const T& a30, const T& a31, const T& a32, const T& a33)
{
	// First compute the det2x2
	const T m01 = a10*a01 - a00*a11;
	const T m02 = a20*a01 - a00*a21;
	const T m03 = a30*a01 - a00*a31;
	const T m12 = a20*a11 - a10*a21;
	const T m13 = a30*a11 - a10*a31;
	const T m23 = a30*a21 - a20*a31;
	// Now compute the minors of rank 3
	const T m012 = m12*a02 - m02*a12 + m01*a22;
	const T m013 = m13*a02 - m03*a12 + m01*a32;
	const T m023 = m23*a02 - m03*a22 + m02*a32;
	const T m123 = m23*a12 - m13*a22 + m12*a32;
	// Now compute the minors of rank 4
	const T m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
	return m0123;
}



template<typename T>
int NumericComputation<T>::
SolveLinearSystem(const T **A, const T *x, int size)
{
	// if we solving something simple, just solve it
	//
	if (size == 2)
	{
		const T det, y[2];

		det = NumericComputation::Determinant2x2(A[0][0], A[0][1], A[1][0], A[1][1]);

		if (det == 0.0)
		{
			// Unable to solve linear system
			return 0;
		}

		y[0] = (A[1][1] * x[0] - A[0][1] * x[1]) / det;
		y[1] = (-A[1][0] * x[0] + A[0][0] * x[1]) / det;

		x[0] = y[0];
		x[1] = y[1];
		return 1;
	}
	else if (size == 1)
	{
		if (A[0][0] == 0.0)
		{
			// Unable to solve linear system
			return 0;
		}

		x[0] /= A[0][0];
		return 1;
	}

	//
	// System of equations is not trivial, use Crout's method
	//

	// Check on allocation of working vectors
	//
	int *index, scratch[10];
	index = (size < 10 ? scratch : new int[size]);

	//
	// Factor and solve matrix
	//
	if (NumericComputation::LUFactorLinearSystem(A, index, size) == 0)
	{
		return 0;
	}
	NumericComputation::LUSolveLinearSystem(A, index, x, size);

	if (size >= 10) delete[] index;
	return 1;
}
template<typename T>
void NumericComputation<T>::
LUSolveLinearSystem(const T **A, int *index,
	const T *x, int size)
{
	int i, j, ii, idx;
	const T sum;
	//
	// Proceed with forward and backsubstitution for L and U
	// matrices.  First, forward substitution.
	//
	for (ii = -1, i = 0; i < size; i++)
	{
		idx = index[i];
		sum = x[idx];
		x[idx] = x[i];

		if (ii >= 0)
		{
			for (j = ii; j <= (i - 1); j++)
			{
				sum -= A[i][j] * x[j];
			}
		}
		else if (sum)
		{
			ii = i;
		}

		x[i] = sum;
	}
	//
	// Now, back substitution
	//
	for (i = size - 1; i >= 0; i--)
	{
		sum = x[i];
		for (j = i + 1; j < size; j++)
		{
			sum -= A[i][j] * x[j];
		}
		x[i] = sum / A[i][i];
	}
}


//----------------------------------------------------------------------------
// Factor linear equations Ax = b using LU decompostion A = LU where L is
// lower triangular matrix and U is upper triangular matrix. Input is
// square matrix A, integer array of pivot indices index[0->n-1], and size
// of square matrix n. Output factorization LU is in matrix A. If error is
// found, method returns 0.
template<typename T>
int NumericComputation<T>::
LUFactorLinearSystem(const T **A, int *index, int size)
{
	const T scratch[10];
	const T *scale = (size<10 ? scratch : new const T[size]);

	int i, j, k;
	int maxI = 0;
	const T largest, temp1, temp2, sum;

	//
	// Loop over rows to get implicit scaling information
	//
	for (i = 0; i < size; i++)
	{
		for (largest = 0.0, j = 0; j < size; j++)
		{
			if ((temp2 = fabs(A[i][j])) > largest)
			{
				largest = temp2;
			}
		}

		if (largest == 0.0)
		{
			//vtkGenericWarningMacro(<<"Unable to factor linear system");
			cout << "Unable to factor linear system" << endl;
			return 0;
		}
		scale[i] = 1.0 / largest;
	}
	//
	// Loop over all columns using Crout's method
	//
	for (j = 0; j < size; j++)
	{
		for (i = 0; i < j; i++)
		{
			sum = A[i][j];
			for (k = 0; k < i; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;
		}
		//
		// Begin search for largest pivot element
		//
		for (largest = 0.0, i = j; i < size; i++)
		{
			sum = A[i][j];
			for (k = 0; k < j; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;

			if ((temp1 = scale[i] * fabs(sum)) >= largest)
			{
				largest = temp1;
				maxI = i;
			}
		}
		//
		// Check for row interchange
		//
		if (j != maxI)
		{
			for (k = 0; k < size; k++)
			{
				temp1 = A[maxI][k];
				A[maxI][k] = A[j][k];
				A[j][k] = temp1;
			}
			scale[maxI] = scale[j];
		}
		//
		// Divide by pivot element and perform elimination
		//
		index[j] = maxI;

		if (fabs(A[j][j]) <= COMPUTATION_SMALL_NUMBER)
		{
			//vtkGenericWarningMacro(<<"Unable to factor linear system");
			cout << "Unable to factor linear system" << endl;
			return 0;
		}

		if (j != (size - 1))
		{
			temp1 = 1.0 / A[j][j];
			for (i = j + 1; i < size; i++)
			{
				A[i][j] *= temp1;
			}
		}
	}

	if (size >= 10) delete[] scale;

	return 1;
}
template<typename T>
const T NumericComputation<T>::
Determinant3x3(const T A[3][3])
{
	return A[0][0] * A[1][1] * A[2][2] + A[1][0] * A[2][1] * A[0][2] +
		A[2][0] * A[0][1] * A[1][2] - A[0][0] * A[2][1] * A[1][2] -
		A[1][0] * A[0][1] * A[2][2] - A[2][0] * A[1][1] * A[0][2];
}

template<typename T>
const T NumericComputation<T>::
Determinant3x3(const T& a00, const T& a01, const T& a02,
	const T& a10, const T& a11, const T& a12,
	const T& a20, const T& a21, const T& a22)
{
	// First compute the det2x2
	const T m01 = a00*a11 - a10*a01;
	const T m02 = a00*a21 - a20*a01;
	const T m12 = a10*a21 - a20*a11;
	// Now compute the minors of rank 3
	const T m012 = m01*a22 - m02*a12 + m12*a02;
	return m012;
}

template<typename T>
Comparison_result NumericComputation<T>::
compare_xyz(const T*p, const T*q) 
{
	//return geom_traits().compare_xyz_3_object()(p, q);
	T px, py, pz, qx, qy, qz;
	px = p[0]; py = p[1]; pz = p[2];
	qx = q[0]; qy = q[1]; qz = q[2];

	Comparison_result c = compare(px, qx);
	if (c != EQUAL) return c;
	c = compare(py, qy);
	if (c != EQUAL) return c;
	return  compare(pz, qz);

}

template<typename T>
Comparison_result NumericComputation<T>::
compare(const T x, const T y)
{
	return (x < y) ? SMALLER : (y < x) ? LARGER : EQUAL;
}



template<typename T>
int NumericComputation<T>::BarycentricCoords(double x[3], double  x1[3], double x2[3],
	double x3[3], double x4[3], double bcoords[4])
{
	double *A[4], p[4], a1[4], a2[4], a3[4], a4[4];
	int i;

	// Homogenize the variables; load into arrays.
	//
	a1[0] = x1[0]; a1[1] = x2[0]; a1[2] = x3[0]; a1[3] = x4[0];
	a2[0] = x1[1]; a2[1] = x2[1]; a2[2] = x3[1]; a2[3] = x4[1];
	a3[0] = x1[2]; a3[1] = x2[2]; a3[2] = x3[2]; a3[3] = x4[2];
	a4[0] = 1.0;   a4[1] = 1.0;   a4[2] = 1.0;   a4[3] = 1.0;
	p[0] = x[0]; p[1] = x[1]; p[2] = x[2]; p[3] = 1.0;

	//   Now solve system of equations for barycentric coordinates
	//
	A[0] = a1;
	A[1] = a2;
	A[2] = a3;
	A[3] = a4;

	if (SolveLinearSystem(A, p, 4))
	{
		for (i = 0; i<4; i++)
		{
			bcoords[i] = p[i];
		}
		return 1;
	}
	else
	{
		return 0;
	}
}


















#endif // !NUMERICCOMPUTATION
