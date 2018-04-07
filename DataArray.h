
#ifndef DATAARRAY_H
#define DATAARRAY_H


#ifndef IDTYPE
typedef int Idtype;
#endif



template<typename T>
class DataArray 
{
protected:
	T* Data;
	Idtype TupleSize;//numbers of tuples including the redundant parts
	Idtype Tuple;//the size of each tuple
	int MaxId;
	//int AppendSize;
public:

	
	//initialize Data，TupleSize,Tuple
	DataArray(size_t tuple, size_t tupleSize)
	{
		//AppendSize = 32;
		MaxId = -1;
		Tuple = tuple;
		TupleSize = tupleSize;
		Data = new T[Tuple*TupleSize];
	}

	DataArray()
	{
		Tuple = 0;
		TupleSize = 0;
		MaxId = -1;
		Data = NULL;
	}

	~DataArray(){ delete[] this->Data; };
	// Description:
	// Allocate a capacity for sz ids in the list and
	// set the number of stored ids in the list to 0.

	void init()
	{
		delete[] this-> Data;
		Tuple = 0;
		TupleSize = 0;
		MaxId = -1;
		Data = NULL;
	}
	void init(size_t tuple, size_t tupleSize)
	{
		MaxId = -1;
		Tuple = tuple;
		TupleSize = tupleSize;
		Data = new T[Tuple*TupleSize];
	}
	
	void clear_to_0(){ memset(Data, 0, sizeof(T)*TupleSize*Tuple); };        
	int size_of_tuple(){ return Tuple; };
	
	int _allocate(const Idtype sz);

	T* get_data(){return this->Data;};

	// Description:
	// Return the number of id's in the list.
	int get_max_tuple_id(){ return this->MaxId; };

	const T* get_tuple(const Idtype i);
	T* get_tuple_nr(Idtype i);

	T get_tuple_element(const Idtype i,const int j);

	void insert_tuple(const Idtype i,const T* tuple);

	Idtype insert_next_tuple(const T* tuple);
	//插入第i个tuple的第j个元素为ele
	void insert_tuple_element(const Idtype i, const int j, const T ele);

	void set_tuple(const Idtype i,const T* tuple);

	T* _resize(const Idtype sz);
	
	void _squeeze(){ this->_resize(this->MaxId+1); };
};

template<typename T>
inline void DataArray<T>::
set_tuple(const Idtype i, const T* tuple)
{
	if (i > MaxId)
	{
		cout << "error0: TypeList, id exceed its maxid" << endl;
		return;
	}
	else{
		for (int j = 0; j < this->Tuple; j++)
		{
			this->Data[i*Tuple + j] = tuple[j];
		}
		if (i > this->MaxId)
			this->MaxId = i;
	}
	
}

template<typename T>
inline void DataArray<T>::
insert_tuple(const Idtype i, const T* tuple)
{
	if (i >= this->TupleSize)
	{
		this->_resize(i + 1);
	}
	for (int j = 0; j < this->Tuple; j++)
	{
		this->Data[i*Tuple + j] = tuple[j];
	}
	if (i > this->MaxId)
		this->MaxId = i;
}


//i(0,1...)是第几个Tuple,j(0,1,...)是每组Tuple的相对编号,ele插入元素
template<typename T>
inline void DataArray<T>::
insert_tuple_element(const Idtype i, const int j, const T ele)
{
	if (i >= this->TupleSize)
	{
		this->_resize(i + 1);
	}
	this->Data[i*Tuple + j] = ele;
	if (i > this->MaxId)
		this->MaxId = i;
}

template<typename T>
inline Idtype DataArray<T>::
insert_next_tuple(const T* tuple)
{
	if (this->MaxId >= this->TupleSize-1)
	{
		this->_resize(this->MaxId+2);
	}
	for (int j = 0; j < Tuple; j++)
	{
		this->Data[(this->MaxId + 1)*Tuple + j] = tuple[j];
	}
	return this->MaxId++ ;
}

template<typename T>
int DataArray<T>::_allocate(const Idtype sz)
{
	if (sz > this->TupleSize)
	{
		this->init();
		this->TupleSize = (sz > 0 ? sz : 1);
		if ((this->Data = new T[this->TupleSize*Tuple]) == NULL)
			return 0;
	}
	this->MaxId = -1;
	return 1;
}

template<typename T>
T* DataArray<T>::
_resize(const Idtype sz)
{
	T* newData;
	Idtype newSize;

 	if (sz > this->TupleSize)
	{
		newSize = this->TupleSize + sz;
	}
	else if (sz == this->TupleSize)
	{
		return this->Data;
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
	if ((newData = new T[newSize*Tuple]) == NULL)
	{
		return 0;
	}

	if (this->Data)
	{
		memcpy(newData,this->Data,
			static_cast<size_t>(sz<this->TupleSize?sz:this->TupleSize)*Tuple*sizeof(T));
		delete[] this->Data;

	}
	this->TupleSize = newSize;
	this->Data = newData;
	return this->Data;
}

template<typename T>
const T* DataArray<T>::
get_tuple(const Idtype i)
{
	/*T* tuple=new T[Tuple];
	for (int j = 0; j < Tuple; j++)
	{
		tuple[j] = Data[i*Tuple+j];
	}
	return tuple;*/
	if (i<=MaxId)
	{
		return (Data+i*Tuple);
	}
	return NULL;
}


//i范围0，1，...；0<=j<=Tuple-1
//i(0,1...)是第几个Tuple,j(0,1,...)是每组Tuple的相对编号,
template<typename T>
T DataArray<T>::
get_tuple_element(const Idtype i,const int j)
{
	if (i<=MaxId)
	{
		return Data[i*Tuple+j];
	}
	return -1;
}






#endif