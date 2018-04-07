#ifndef UTILITY
#define UTILITY


template <int i, typename T>
struct Tuple_get;

template <typename T>
struct Tuple_get<0, T>
{
	typedef typename T::first_type result_type;
	static result_type       & get(T       & t) { return t.first; }
	static result_type const & get(T const & t) { return t.first; }
};

template <typename T>
struct Tuple_get<1, T>
{
	typedef typename T::second_type result_type;
	static result_type       & get(T       & t) { return t.second; }
	static result_type const & get(T const & t) { return t.second; }
};

template <typename T>
struct Tuple_get<2, T>
{
	typedef typename T::third_type result_type;
	static result_type       & get(T       & t) { return t.third; }
	static result_type const & get(T const & t) { return t.third; }
};

template <typename T>
struct Tuple_get<3, T>
{
	typedef typename T::fourth_type result_type;
	static result_type       & get(T       & t) { return t.fourth; }
	static result_type const & get(T const & t) { return t.fourth; }
};



//+---------------------------------------------------------------------+
//| Triple class                                                        |
//+---------------------------------------------------------------------+

template <class T1, class T2, class T3>
class Triple
{
	typedef Triple<T1, T2, T3> Self;

public:

	typedef T1 first_type;
	typedef T2 second_type;
	typedef T3 third_type;

	T1 first;
	T2 second;
	T3 third;

	Triple() {}

	Triple(const T1& a, const T2& b, const T3& c)
		: first(a), second(b), third(c)
	{}

	template <class U, class V, class W>
	Triple(const U& a, const V& b, const W& c)
		: first(a), second(b), third(c)
	{}

	template <class U, class V, class W>
	Triple& operator=(const Triple<U, V, W> &t) {
		first = t.first;
		second = t.second;
		third = t.third;
		return *this;
	}

	template < int i >
	typename Tuple_get<i, Self>::result_type const &
		get() const
	{
		return Tuple_get<i, Self>::get(*this);
	}

	template < int i >
	typename Tuple_get<i, Self>::result_type &
		get()
	{
		return Tuple_get<i, Self>::get(*this);
	}

};

template <class T1, class T2, class T3>
inline
Triple<T1, T2, T3> make_triple(const T1& x, const T2& y, const T3& z)
{
	return Triple<T1, T2, T3>(x, y, z);
}

template <class T1, class T2, class T3>
inline
Triple<T1, T2, T3> make_tuple(const T1& x, const T2& y, const T3& z)
{
	return Triple<T1, T2, T3>(x, y, z);
}

template <class T1, class T2, class T3>
inline bool operator==(const Triple<T1, T2, T3>& x,
	const Triple<T1, T2, T3>& y)
{
	return ((x.first == y.first) &&
		(x.second == y.second) &&
		(x.third == y.third));
}

template <class T1, class T2, class T3>
inline bool operator!=(const Triple<T1, T2, T3>& x,
	const Triple<T1, T2, T3>& y)
{
	return !(x == y);
}

template <class T1, class T2, class T3>
inline bool operator<(const Triple<T1, T2, T3>& x,
const Triple<T1, T2, T3>& y)
{
	return (x.first < y.first ||
		(!(y.first < x.first) &&
		(x.second < y.second ||
		(!(y.second < x.second) && x.third < y.third))));
}

#endif // !UTILITY
