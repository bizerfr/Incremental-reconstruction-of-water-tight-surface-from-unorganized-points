

#ifndef TYPELIST_H
#define TYPELIST_H

#ifndef IDTYPE
typedef int Idtype;
#endif

template<typename T>
class TypeList
{
protected:
	T* Types;
	int Size;
	Idtype MaxId;

public:
	TypeList(){ MaxId = -1; Size = 0; Types = NULL; };
	void init(){ delete[] Types; MaxId = -1; Size = 0; Types = NULL; };
	void init(int size){ delete[] Types; MaxId = -1; Size = size; Types = new T[size]; };
	void clear_to_0(){ memset(Types, 0, sizeof(T)*Size); };        
		

	int _allocate(const Idtype sz);

	T* get_data(){return this->Types;};

	Idtype get_max_id(){ return this->MaxId; };

	Idtype get_element(const Idtype i){ return this->Types[i]; };

	void set_element(const Idtype i, const T _type);

	//与set_element()还是有些区别的
	void insert_element(const Idtype i,const T _type);

	Idtype insert_next_element(const T _type);

	T* get_pointer(const Idtype i){ return this->Types + i; };

	void reset(){ this->MaxId = -1; };

	T* _resize(const Idtype sz);

	void _squeeze(){ this->_resize(this->MaxId+1); };

};



template<typename T>
inline void TypeList<T>::
insert_element(const Idtype i, const T _type)
{
	if (i >= this->Size)
	{
		this->_resize(i+1);
	}
	this->Types[i] = _type;
	if (i > this->MaxId)
	{
		this->MaxId = i;
	}
}

template<typename T>
inline void TypeList<T>::
set_element(const Idtype i, const T _type)
{
	if (i > MaxId)
	{
		cout << "error0: TypeList, id exceed its maxid" << endl;
		return;
	}
	this->Types[i] = _type;
	if (i > this->MaxId)
	{
		this->MaxId = i;
	}
}
template<typename T>
inline Idtype TypeList<T>::
insert_next_element(const T _type)
{
	if (this->MaxId >= this->Size - 1)
	{
		this->_resize(this->MaxId+2);
	}
	this->Types[++this->MaxId] = _type;
	return this->MaxId;
}



template<typename T>
int TypeList<T>::
_allocate(const Idtype sz)
{
	if (sz > this->Size)
	{
		this->init();
		this->Size = (sz>0?sz:1);
		if ((this->Types = new T[this->Size]) == NULL)
			return 0;
	}
	this->MaxId = -1;
	return 1;
}

template<typename T>
T* TypeList<T>::
_resize(const Idtype sz)
{
	T* newTypes;
	Idtype newSize;
	if (sz > this->Size)
	{
		newSize = this->Size + sz;
	}
	else if (sz == this->Size)
	{
		return this->Types;
	}
	else
	{
		newSize = sz;
	}
	if (newSize <= 0)
	{
		this->init();
		return 0;
	}
	if ((newTypes = new T[newSize]) == NULL)
		return 0;
	if (this->Types)
	{
		memcpy(newTypes,this->Types,static_cast<size_t>(sz<this->Size?sz:this->Size)*sizeof(T));
		delete[] this->Types;
	}
	this->Size = newSize;
	this->Types = newTypes;
	return this->Types;
}



#endif // !TYPELIST_H
