#include"math.h"
#include"DataArray.h"
#include"NumericComputation.h"
#ifndef GEOMETRICTRAITS
#define GEOMETRICTRAITS

typedef Sign Orientation;
typedef Sign Oriented_side;

enum  Bounded_side
{
	ON_UNBOUNDED_SIDE = -1,
	ON_BOUNDARY,
	ON_BOUNDED_SIDE
};
enum  Angle
{
	OBTUSE = -1,
	RIGHT,
	ACUTE
};
enum COLLINEAR_POSITION { BEFORE, SOURCE, MIDDLE, TARGET, AFTER };

template<typename T>
class GeometricTraits 
{
private:

public:
	template <class T>
	inline T opposite( T& t)
	{
		return -t;
	}


	inline Bounded_side opposite(Bounded_side bs)
	{
		return static_cast<Bounded_side>(-static_cast<int>(bs));
	}

	inline Angle opposite(Angle a)
	{
		return static_cast<Angle>(-static_cast<int>(a));
	}


	template < typename T, typename U >
	inline T enum_cast( U& u)
	{
		return static_cast<T>(u);
	}

public:
	static Orientation orientationC2( const T&ux,  const T&uy,  const T&vx,  const T&vy)
	{
		return NumericComputation<T>::SignOfDeterminant2x2(ux, uy, vx, vy);
	}
	static Orientation orientationC2( const T&px,  const T&py,
		 const T&qx,  const T&qy,
		 const T&rx,  const T&ry)
	{
		return NumericComputation<T>::SignOfDeterminant2x2(qx - px, qy - py, rx - px, ry - py);
	}
	static Orientation coplanar_orientationC3( const T&px,  const T&py,  const T&pz,
		 const T&qx,  const T&qy,  const T&qz,
		 const T&rx,  const T&ry,  const T&rz)
	{
		//typedef typename Same_uncertainty_nt<Orientation, T>::type  Ori;
		Orientation oxy_pqr = orientationC2(px, py, qx, qy, rx, ry);
		if (oxy_pqr != COLLINEAR)
			return oxy_pqr;

		Orientation oyz_pqr = orientationC2(py, pz, qy, qz, ry, rz);
		if (oyz_pqr != COLLINEAR)
			return oyz_pqr;

		return orientationC2(px, pz, qx, qz, rx, rz);
	}
	static Orientation coplanar_orientation(const T* p1, const T* p2, const T* p3)
	{
		return coplanar_orientationC3(p1[0],p1[1],p1[2],
							p2[0],p2[1],p2[2],
							p3[0],p3[1],p3[2]);
	}

	static bool collinear(const T* p, const T* q, const T* r) 
	{
		return coplanar_orientation(p, q, r) == COLLINEAR;
	}

	static bool equal(const T* p, const T* q) 
	{
		return NumericComputation<T>::compare_xyz(p, q) == EQUAL;
	}

	
	static Oriented_side Side_of_oriented_sphere_3( const T* p,  const T* q,  const T* r,
		 const T* s,  const T* t) ;
	static Bounded_side coplanar_side_of_bounded_circleC3(const T* p, const T* q, const T* r,
		   const T* t);
	static Bounded_side side_of_bounded_sphereC3(const T* p, const T* q, const T* r,
		const T* t);
	static COLLINEAR_POSITION collinear_position( const T*s,  const T*p,  const T*t) ;

	static Orientation 
		orientation(const T* p, const T* q, const T* r, const T* s);
	static Orientation
		inexact_orientation(const T* p, const T* q, const T* r, const T* s);
	//q->r->s满足右手定则,求点到面的距离
	static T distance_point_to_facet(const T* p, const T* q, const T* r, const T* s)
	{
		T rq[3]={q[0]-r[0],q[1]-r[1],q[2]-r[2]};
		T rs[3]={s[0]-r[0],s[1]-r[1],s[2]-r[2]};
		const T* normal=NumericComputation<T>::Cross(rs,rq);
		T mo=NumericComputation<T>::Dot(normal,normal);
		mo=sqrt(mo);
		T rp[3]={p[0]-r[0],p[1]-r[1],p[2]-r[2]};
		T dot=NumericComputation<T>::Dot(rp,normal);
		delete []normal;
		return dot/mo;
	}
	//面的方向即其法向的方向，面q0,r0,s0与面q1,r1,s1的法向均满足右手定则
	//返回值为法向量夹角的补角的余弦值，所以注意法向量方向
	static T cos_dihedral(const T* q0, const T* r0, const T* s0, const T* q1, const T* r1, const T* s1)
	{
		T rq0[3] = { q0[0] - r0[0], q0[1] - r0[1], q0[2] - r0[2] };
		T rs0[3] = { s0[0] - r0[0], s0[1] - r0[1], s0[2] - r0[2] };
		const T* normal0 = NumericComputation<T>::Cross(rs0, rq0);
		T mo0 = NumericComputation<T>::Dot(normal0, normal0);
		mo0 = sqrt(mo0);

		T rq1[3]={q1[0]-r1[0],q1[1]-r1[1],q1[2]-r1[2]};
		T rs1[3]={s1[0]-r1[0],s1[1]-r1[1],s1[2]-r1[2]};
		const T* normal1=NumericComputation<T>::Cross(rs1,rq1);
		T mo1=NumericComputation<T>::Dot(normal1,normal1);
		mo1=sqrt(mo1);
		if (mo0<=0.000000001||mo1<0.000000001)
		{
			delete []normal0;
			delete []normal1;
			return 1;
		}
		T dot=NumericComputation<T>::Dot(normal0, normal1);
		delete []normal0;
		delete []normal1;
		return -dot / (mo0*mo1);
	}

};

//p,q,r组成面，s点，判断点在面的positive,negative;
//p,q,r,s符合右手定则，POSITIVE（1）；共面，ZERO(0)；否则,NEGATIVE(-1)
//求混合积，向量<pq>叉乘向量<pr>，点乘向量<ps>,p为0，q为1,r为2，右手定则
template<typename T>
Orientation GeometricTraits<T>::
inexact_orientation(const T* p, const T* q, const T* r, const T* s)
{
	double px = p[0], py = p[1], pz = p[2],
		qx = q[0], qy = q[1], qz = q[2],
		rx = r[0], ry = r[1], rz = r[2],
		sx = s[0], sy = s[1], sz = s[2];

	double pqx = qx - px;
	double pqy = qy - py;
	double pqz = qz - pz;
	double prx = rx - px;
	double pry = ry - py;
	double prz = rz - pz;
	double psx = sx - px;
	double psy = sy - py;
	double psz = sz - pz;

	double det = NumericComputation<T>::Determinant3x3(pqx, pqy, pqz,
		prx, pry, prz,
		psx, psy, psz);
	if (det > 0) return POSITIVE;
	if (det < 0) return NEGATIVE;
	return ZERO;
}

//p,q,r组成面，s点，判断点在面的positive,negative
//p,q,r,s符合右手定则，POSITIVE（1,LARGER）；共面，ZERO(0,EQUAL)；否则,NEGATIVE(-1,SMALLER)；即，POSITIVE表示s在有向面pqr法线一侧
//求混合积，向量<pq>叉乘向量<pr>，点乘向量<ps>,p为0，q为1,r为2，右手定则
//该函数在inexact_orientation基础上考虑了一些扰动，更适用浮点数计算
template<typename T>
Orientation GeometricTraits<T>::
orientation(const T* p, const T* q, const T* r, const T* s)
{
	double px=p[0], py=p[1], pz=p[2],
		qx=q[0], qy=q[1], qz=q[2], 
		rx=r[0], ry=r[1], rz=r[2],
		sx=s[0], sy=s[1], sz=s[2];
	double pqx = qx - px;
	double pqy = qy - py;
	double pqz = qz - pz;
	double prx = rx - px;
	double pry = ry - py;
	double prz = rz - pz;
	double psx = sx - px;
	double psy = sy - py;
	double psz = sz - pz;

	double maxx = abs(pqx);
	double maxy = abs(pqy);
	double maxz = abs(pqz);

	double aprx = abs(prx);
	double apsx = abs(psx);

	double apry = abs(pry);
	double apsy = abs(psy);

	double aprz = abs(prz);
	double apsz = abs(psz);

	if (maxx < aprx) maxx = aprx;
	if (maxx < apsx) maxx = apsx;
	if (maxy < apry) maxy = apry;
	if (maxy < apsy) maxy = apsy;
	if (maxz < aprz) maxz = aprz;
	if (maxz < apsz) maxz = apsz;

	double det = NumericComputation<T>::Determinant3x3(pqx, pqy, pqz,
		prx, pry, prz,
		psx, psy, psz);

	double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;  //eps ,epsilon meaning small enough to be insignificant

	if (maxx > maxz)
		std::swap(maxx, maxz);
	if (maxy > maxz)
		std::swap(maxy, maxz);
	else if (maxy < maxx)
		std::swap(maxx, maxy);  //maxx,maxy,maxz 从小到大排列

	// Protect against underflow in the computation of eps.
	if (maxx < 1e-97) /* cbrt(min_double/eps) */ {
		if (maxx == 0)
			return ZERO;
	}
	// Protect against overflow in the computation of det.
	else if (maxz < 1e102) /* cbrt(max_double [hadamard]/4) */ {

		if (det > eps)  return POSITIVE;
		if (det < -eps) return NEGATIVE;
	}
	//等同于求Determinant3x3符合，不过|det|<COMPUTATION_SMALL_NUMBER,认为det为ZERO
	return NumericComputation<T>::SignOfDeterminant3x3(qx-px,rx-px,sx-px,
															qy-py,ry-py,sy-py,
															qz-pz,rz-pz,sz-pz);	
}

//p,q,r,s需要符合右手定则
//t在四面体pqrs的外接球上，返回ON_ORIENTED_BOUNDARY（0）
//t在四面体pqrs的外接球外，返回ON_NEGATIVE_SIDE（-1）
//t在四面体pqrs的外接球内，返回ON_POSITIVE_SIDE（1）
template<typename T>
Oriented_side GeometricTraits<T>::
Side_of_oriented_sphere_3( const T* p,  const T* q,  const T* r,
 const T* s,  const T* t) 
{
	double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz;
	px = p[0]; py = p[1]; pz = p[2];
	qx = q[0]; qy = q[1]; qz = q[2];
	rx = r[0]; ry = r[1]; rz = r[2];
	sx = s[0]; sy = s[1]; sz = s[2];
	tx = t[0]; ty = t[1]; tz = t[2];
	double ptx = px - tx;
	double pty = py - ty;
	double ptz = pz - tz;
	//double pt2 =  square(ptx) +  square(pty)+  square(ptz);
	double pt2 = ptx*ptx + pty*pty + ptz*ptz;
		
	double qtx = qx - tx;
	double qty = qy - ty;
	double qtz = qz - tz;
	//double qt2 =  square(qtx) +  square(qty)+  square(qtz);
	double qt2 = qtx*qtx + qty*qty + qtz*qtz;
		
	double rtx = rx - tx;
	double rty = ry - ty;
	double rtz = rz - tz;
	//double rt2 =  square(rtx) +  square(rty)+  square(rtz);
	double rt2 = rtx*rtx + rty*rty + rtz*rtz;
		
	double stx = sx - tx;
	double sty = sy - ty;
	double stz = sz - tz;
	//double st2 =  square(stx) +  square(sty)+  square(stz);
	double st2 = stx*stx + sty*sty + stz*stz;

	// Compute the semi-static bound.
	double maxx = abs(ptx);
	double maxy = abs(pty);
	double maxz = abs(ptz);

	double aqtx = abs(qtx);
	double artx = abs(rtx);
	double astx = abs(stx);

	double aqty = abs(qty);
	double arty = abs(rty);
	double asty = abs(sty);

	double aqtz = abs(qtz);
	double artz = abs(rtz);
	double astz = abs(stz);


	if (maxx < aqtx) maxx = aqtx;
	if (maxx < artx) maxx = artx;
	if (maxx < astx) maxx = astx;

	if (maxy < aqty) maxy = aqty;
	if (maxy < arty) maxy = arty;
	if (maxy < asty) maxy = asty;

	if (maxz < aqtz) maxz = aqtz;
	if (maxz < artz) maxz = artz;
	if (maxz < astz) maxz = astz;

	double eps = 1.2466136531027298e-13 * maxx * maxy * maxz;

	if (maxx > maxz)
		std::swap(maxx, maxz);
	if (maxy > maxz)
		std::swap(maxy, maxz);
	else if (maxy < maxx)
		std::swap(maxx, maxy); //排序 maxx<maxy<maxz

	double det = NumericComputation<T>::determinant4x4(ptx, pty, ptz, pt2,
		rtx, rty, rtz, rt2,
		qtx, qty, qtz, qt2,
		stx, sty, stz, st2);

	// Protect against underflow in the computation of eps.
	if (maxx < 1e-58) /* sqrt^5(min_double/eps) */ {
		if (maxx == 0)
			return ON_ORIENTED_BOUNDARY;
	}
	// Protect against overflow in the computation of det.
	else if (maxz < 1e61) /* sqrt^5(max_double/4 [hadamard]) */ {
		eps *= (maxz * maxz);
		if (det > eps)  return ON_POSITIVE_SIDE;
		if (det < -eps) return ON_NEGATIVE_SIDE;
	}
	return NumericComputation<T>::SignOfDeterminant4x4(ptx, pty, ptz, pt2,
		rtx, rty, rtz, rt2,
		qtx, qty, qtz, qt2,
		stx, sty, stz, st2);
}

//要求pqrt共面,pqr逆时针排序
template<typename T>
Bounded_side GeometricTraits<T>::
coplanar_side_of_bounded_circleC3( const T* p,  const T* q,  const T* r,
								   const T* t)
{
	// The approach is to compute side_of_bounded_sphere(p,q,r,t+v,t),
	// with v = pq ^ pr.
	// Note : since the circle defines the orientation of the plane, it can not
	// be considered oriented.
	T px, py, pz, qx, qy, qz, rx, ry, rz, tx, ty, tz;
	px = p[0]; py = p[1]; pz = p[2];
	qx = q[0]; qy = q[1]; qz = q[2];
	rx = r[0]; ry = r[1]; rz = r[2];
	tx = t[0]; ty = t[1]; tz = t[2];

	T ptx = px - tx;
	T pty = py - ty;
	T ptz = pz - tz;
	double pt2 = ptx*ptx + pty*pty + ptz*ptz;
	T qtx = qx - tx;
	T qty = qy - ty;
	T qtz = qz - tz;
	double qt2 = qtx*qtx + qty*qty + qtz*qtz;
	T rtx = rx - tx;
	T rty = ry - ty;
	T rtz = rz - tz;
	double rt2 = rtx*rtx + rty*rty + rtz*rtz;
	T pqx = qx - px;
	T pqy = qy - py;
	T pqz = qz - pz;
	T prx = rx - px;
	T pry = ry - py;
	T prz = rz - pz;
	T vx = pqy*prz - pqz*pry;
	T vy = pqz*prx - pqx*prz;
	T vz = pqx*pry - pqy*prx;  //向量v = pq ^ pr,^为叉乘
	double v2 = vx*vx + vy*vy + vz*vz;

	return static_cast<Bounded_side>(NumericComputation<T>::SignOfDeterminant4x4(ptx, pty, ptz, pt2,
		rtx, rty, rtz, rt2,
		qtx, qty, qtz, qt2,
		vx, vy, vz, v2));
}


//点t在pqs最小外接球内,外还是上;
//t is outside of bounded sphere, return ON_UNBOUNDED_SIDE(-1);t is inside of bounded sphere, return ON_BOUNDED_SIDE(1);ON_BOUNDARY(0)
//m表示三角形pqs外心,向量<sm>=m-s;向量<sr>=<sp>^<sq>,其中^表示叉乘
//<sm>*2|sr|^2=|sp|^2 <sq>*<sr> - |sq|^2 <sp>*<sr>
template<typename T>
Bounded_side GeometricTraits<T>::
side_of_bounded_sphereC3(const T* p, const T* q, const T* s,
	const T* t)
{
	T px, py, pz, qx, qy, qz, sx, sy, sz, tx, ty, tz;
	px = p[0]; py = p[1]; pz = p[2];
	qx = q[0]; qy = q[1]; qz = q[2];
	sx = s[0]; sy = s[1]; sz = s[2];
	tx = t[0]; ty = t[1]; tz = t[2];

	T psx = px - sx;
	T psy = py - sy;
	T psz = pz - sz;
	T ps2 = psx*psx + psy*psy + psz*psz;
	T qsx = qx - sx;
	T qsy = qy - sy;
	T qsz = qz - sz;
	T qs2 = qsx*qsx + qsy*qsy + qsz*qsz;
	T rsx = psy*qsz - psz*qsy;
	T rsy = psz*qsx - psx*qsz;
	T rsz = psx*qsy - psy*qsx; //向量sr=向量sp叉乘向量sq
	T tsx = tx - sx;
	T tsy = ty - sy;
	T tsz = tz - sz;

	T num_x = ps2 * NumericComputation<T>::Determinant2x2(qsy, qsz, rsy, rsz)
		- qs2 * NumericComputation<T>::Determinant2x2(psy, psz, rsy, rsz);
	T num_y = ps2 * NumericComputation<T>::Determinant2x2(qsx, qsz, rsx, rsz)
		- qs2 * NumericComputation<T>::Determinant2x2(psx, psz, rsx, rsz);
	T num_z = ps2 * NumericComputation<T>::Determinant2x2(qsx, qsy, rsx, rsy)
		- qs2 * NumericComputation<T>::Determinant2x2(psx, psy, rsx, rsy); //(num_x,-num_y,num_z)是<sm>*2|sr|^2坐标,即<sm>*den2

	T den2 = 2 * NumericComputation<T>::Determinant3x3(psx, psy, psz,
		qsx, qsy, qsz,
		rsx, rsy, rsz);   //向量sp,向量sq,向量sr混合积的2倍;den2即2*sr^2

	// The following could be simplified a bit.
	return static_cast<Bounded_side>(
		NumericComputation<T>::cmp_dist_to_pointC3(num_x, -num_y, num_z,
		psx*den2, psy*den2, psz*den2,
		tsx*den2, tsy*den2, tsz*den2)); //实际上比较|pm|,|tm|
}


template<typename T>
COLLINEAR_POSITION GeometricTraits<T>::
collinear_position( const T*s,  const T*p,  const T*t) 
// (s,t) defines a line, p is on that line.
// Depending on the position of p wrt s and t, returns :
// --------------- s ---------------- t --------------
// BEFORE       SOURCE    MIDDLE    TARGET       AFTER
{

	Comparison_result ps = NumericComputation<T>::compare_xyz(p, s);
	if (ps == EQUAL)
		return SOURCE;
	Comparison_result st = NumericComputation<T>::compare_xyz(s, t);
	if (ps == st)
		return BEFORE;
	Comparison_result pt = NumericComputation<T>::compare_xyz(p, t);
	if (pt == EQUAL)
		return TARGET;
	if (pt == st)
		return MIDDLE;
	return AFTER;
}



#endif // !GEOMETRICTRAITS
