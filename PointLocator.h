#include"DataArray.h"
#include"TypeList.h"
#include"NumericComputation.h"
#ifndef POINTLOCATOR_H
#define POINTLOCATOR_H

typedef Idtype Indextype;
typedef Idtype Vertex_Idtype;
typedef Idtype Cell_Idtype;
//point 的坐标值应该放在PointLocator里面的。
template<typename T,typename T_INFO>
class PointLocator 
{
	typedef TypeList<bool> FlagList;
	typedef TypeList<Idtype> IdList;
	typedef DataArray<T> DataList;
	
	

	

private:
	size_t Level;//Level>=1
	DataList Points;
	DataArray<T_INFO> PointAttributes;
	T Bounds[2][3];    //点所在区域的边界  下界---->上界
	FlagList ParentOctant;
	
	std::vector<Vertex_Idtype>* TerminalOctant;
public:
	void get_point_data(DataList& PointSet)
	{PointSet=Points;}
	void init_point_locator(int level,T *bounds[3]);
	const const T* get_point(Vertex_Idtype i);
	void set_point(Vertex_Idtype IdV,const T* p);
	Vertex_Idtype insert_vertex(const T*);
	Vertex_Idtype nearest_inexact(const T* p);
	Vertex_Idtype nearest_exact(const T* p);
	Indextype index(const T* p, int level);

};

template<typename T,typename T_INFO>
void PointLocator<T,typename T_INFO>::
init_point_locator(int level,T *bounds[3])
{
	this->Level = level;

	for (int i = 0;i < 3;  i++)
	{
		this->Bounds[0][i] = bounds[0][i];
		this->Bounds[1][i] = bounds[1][i];
	}
	Points.init(3,0);
	//int NUM = (8 ^ this->Level - 1) / 7;
	//pow计算8^this->Level次幂
	int termin_num = pow(8,this->Level);
	int NUM = (termin_num - 1) / 7;
	this->ParentOctant.init(NUM);
	this->ParentOctant.clear_to_0();
	this->TerminalOctant = new std::vector<Vertex_Idtype>[termin_num];
}

template<typename T,typename T_INFO>
const const T* PointLocator<T,T_INFO>::
get_point(Vertex_Idtype i)
{
	return Points.get_tuple(i);
}

template<typename T,typename T_INFO>
Vertex_Idtype PointLocator<T,T_INFO>::
insert_vertex(const T* p)
{
	Indextype index_temp;
	Indextype index_pre;
	Idtype vertex_id=Points.insert_next_tuple(p);
	if (ParentOctant.get_element(0) == false)
		ParentOctant.insert_element(0,true);
	for (int i = 1; i < Level; i++)
	{
		index_temp = index(p,i);
		index_pre = (pow(8,i) - 1) / 7;
		index_mark = index_temp + index_pre - 1;
		if (ParentOctant.get_element(0) == false)
			ParentOctant.insert_element(0, true);

	}
	/*Indextype index_temp = index(p, i);
	index_pre = (8 ^ i - 1) / 7;
	index_mark = index_temp + index_pre - 1;*/
	index_mark = index(p,Level);
	TerminalOctant[index_mark].push_back(vertex_id);
	return vertex_id;

}

template<typename T,typename T_INFO>
void PointLocator<T,T_INFO>::
set_point(Vertex_Idtype IdV, const T* p)
{

	Points.insert_tuple(IdV,p);

}

template<typename T,typename T_INFO>
Indextype PointLocator<T,T_INFO>::
index(const T* p, int level)
{
	Indextype index_temp[3];
	
	for (int i = 0; i < 3; i++)
		index_temp[i] = floor((p[i]-this->Bounds[0][i])/(this->Bounds[1][i]-this->Bounds[0][i])*(1<<Level));
	return index_temp[0] + index_temp[1] * (1 << Level) + index_temp[2] * (1 << (Level * 2));

}

template<typename T,typename T_INFO>
Vertex_Idtype PointLocator<T,T_INFO>::
nearest_inexact(const T* p)
{
	Vertex_Idtype Nei_Point;
	int OctreeDepth = Level;
	int depth = OctreeDepth + 1;
	int MaxIndex = (1 << (OctreeDepth)) - 1;
	int IndexTemp[3];
	int NeiPointIndex[3][2];
	for (int i = 0; i<3; i++)
	{
		IndexTemp[i] = floor(((p[i] - Bounds[1][i]) / (Bounds[2][i] - Bounds[1][i]))*(1 << depth));
		if ((IndexTemp[i]<0) || (IndexTemp[i]>(MaxIndex << 2)))
		{
			return 0;
		}
		int temp = floor(IndexTemp[i] / 2);
		if (IndexTemp[i] % 2)
		{
			NeiPointIndex[i][0] = temp;
			if (temp<MaxIndex)
			{
				NeiPointIndex[i][1] = temp + 1;
			}
			else
				NeiPointIndex[i][1] = temp;
		}
		else
		{
			if (temp>0)
			{
				NeiPointIndex[i][0] = temp - 1;
			}
			else
				NeiPointIndex[i][0] = temp;
			NeiPointIndex[i][1] = temp;
		}
	}

	//------------------------计算插入点附近的方格索引-----------------------
	vector<Idtype> NeiborPointId;
	//取附近八个方格
	//取Z=TerminalOctant[3][0]的平面上的方格
	for (int i = NeiPointIndex[0][0]; i <= NeiPointIndex[0][1]; i++)
	{
		//int  IndexMarkX=NeiPointIndex[0][i];
		for (int j = NeiPointIndex[1][0]; j <= NeiPointIndex[1][1]; j++)
		{
			if (NeiPointIndex[2][0]<NeiPointIndex[2][1])
			{
				int IndexMark0 = i + (1 << OctreeDepth)*j + (1 << (2 * OctreeDepth))*NeiPointIndex[2][0];
				int IndexMark1 = i + (1 << OctreeDepth)*j + (1 << (2 * OctreeDepth))*NeiPointIndex[2][1];
				if (!TerminalOctant[IndexMark0].empty())
				{
					auto PointIdIt = TerminalOctant[IndexMark0].begin();
					while (PointIdIt != TerminalOctant[IndexMark0].end())
					{
						NeiborPointId.push_back(*PointIdIt);
						PointIdIt++;
					}
				}
				if (!TerminalOctant[IndexMark1].empty())
				{
					auto PointIdIt = TerminalOctant[IndexMark1].begin();
					while (PointIdIt != TerminalOctant[IndexMark1].end())
					{
						NeiborPointId.push_back(*PointIdIt);
						PointIdIt++;
					}
				}
			}
			else
			{
				int IndexMark = i + (1 << OctreeDepth)*j + (1 << (2 * OctreeDepth))*NeiPointIndex[2][1];
				if (!TerminalOctant[IndexMark].empty())
				{
					auto PointIdIt = TerminalOctant[IndexMark].begin();
					while (PointIdIt != TerminalOctant[IndexMark].end())
					{
						NeiborPointId.push_back(*PointIdIt);
						PointIdIt++;
					}
				}
			}
		}
	}
	//若邻域的八个方格不包含点,则向外扩展一排方格
	while (NeiborPointId.empty())
	{
		if ((NeiPointIndex[0][0] == 0) && (NeiPointIndex[1][0] == 0) && (NeiPointIndex[2][0] = 0) &&
			(NeiPointIndex[0][1] == MaxIndex) && (NeiPointIndex[1][1] == MaxIndex) && (NeiPointIndex[2][1] == MaxIndex))
		{
			return 0; cout << "find no point" << endl;
		}
		int NeiPointIndexLast[3][2];
		for (int i = 0; i<3; i++)
		{
			NeiPointIndexLast[i][0] = NeiPointIndex[i][0];
			NeiPointIndexLast[i][1] = NeiPointIndex[i][1];
			if (NeiPointIndex[i][0]>0)
				NeiPointIndex[i][0]--;
			if (NeiPointIndex[i][1]<MaxIndex)
				NeiPointIndex[i][1]++;
		}
		int ZMarked[2] = { 0, 0 };
		int YMarked[2] = { 0, 0 };
		//扩展面X=
		if ((NeiPointIndexLast[2][0]>0) || (NeiPointIndexLast[2][1]<MaxIndex))
		{
			for (int i = NeiPointIndex[0][0]; i <= NeiPointIndex[0][1]; i++)
			{
				//int  IndexMarkX=NeiPointIndex[0][i];
				for (int j = NeiPointIndex[1][0]; j <= NeiPointIndex[1][1]; j++)
				{
					if ((NeiPointIndexLast[2][0]>0) && (NeiPointIndexLast[2][1]<MaxIndex))
					{
						/*NeiPointIndex[2][0]--;
						NeiPointIndex[2][1]++;*/
						ZMarked[0] = 1;
						ZMarked[1] = -1;
						int IndexMark0 = i + (1 << OctreeDepth)*j + (1 << (2 * OctreeDepth))*NeiPointIndex[2][0];
						int IndexMark1 = i + (1 << OctreeDepth)*j + (1 << (2 * OctreeDepth))*NeiPointIndex[2][1];
						if (!TerminalOctant[IndexMark0].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark0].begin();
							while (PointIdIt != TerminalOctant[IndexMark0].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
						if (!TerminalOctant[IndexMark1].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark1].begin();
							while (PointIdIt != TerminalOctant[IndexMark1].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
					}
					else
					{
						int IndexMark;
						if (NeiPointIndexLast[2][0]>0)
						{
							/*--NeiPointIndex[2][0];*/
							ZMarked[0] = 1;
							IndexMark = i + (1 << OctreeDepth)*j + (1 << (2 * OctreeDepth))*NeiPointIndex[2][0];
						}
						else
						{
							/* ++NeiPointIndex[2][1];*/
							ZMarked[1] = -1;
							IndexMark = i + (1 << OctreeDepth)*j + (1 << (2 * OctreeDepth))*NeiPointIndex[2][1];
						}
						if (!TerminalOctant[IndexMark].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark].begin();
							while (PointIdIt != TerminalOctant[IndexMark].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
					}
				}
			}

		}
		//扩展面Y=？
		if ((NeiPointIndexLast[1][0]>0) || (NeiPointIndexLast[1][1]<MaxIndex))
		{
			for (int i = NeiPointIndex[0][0]; i <= NeiPointIndex[0][1]; i++)
			{
				//int  IndexMarkX=NeiPointIndex[0][i];
				for (int j = NeiPointIndex[2][0] + ZMarked[0]; j <= NeiPointIndex[2][1] + ZMarked[1]; j++)
				{
					if ((NeiPointIndexLast[1][0]>0) && (NeiPointIndexLast[1][1]<MaxIndex))
					{
						YMarked[0] = 1;
						YMarked[1] = -1;
						int IndexMark0 = i + (1 << OctreeDepth)*NeiPointIndex[1][0] + (1 << (2 * OctreeDepth))*j;
						int IndexMark1 = i + (1 << OctreeDepth)*NeiPointIndex[1][1] + (1 << (2 * OctreeDepth))*j;
						if (!TerminalOctant[IndexMark0].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark0].begin();
							while (PointIdIt != TerminalOctant[IndexMark0].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
						if (!TerminalOctant[IndexMark1].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark1].begin();
							while (PointIdIt != TerminalOctant[IndexMark1].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
					}
					else
					{
						int IndexMark;
						if (NeiPointIndexLast[1][0]>0)
						{
							/*--NeiPointIndex[2][0];*/
							YMarked[0] = 1;
							IndexMark = i + (1 << OctreeDepth)*NeiPointIndex[1][0] + (1 << (2 * OctreeDepth))*j;
						}
						else
						{
							/* ++NeiPointIndex[2][1];*/
							YMarked[1] = -1;
							IndexMark = i + (1 << OctreeDepth)*NeiPointIndex[1][1] + (1 << (2 * OctreeDepth))*j;
						}
						if (!TerminalOctant[IndexMark].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark].begin();
							while (PointIdIt != TerminalOctant[IndexMark].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
					}
				}
			}

		}
		//扩展面Z=
		if ((NeiPointIndexLast[0][0]>0) || NeiPointIndexLast[1][0]<MaxIndex)
		{
			for (int i = NeiPointIndex[1][0] + YMarked[0]; i <= NeiPointIndex[1][1] + YMarked[1]; i++)
			{
				//int  IndexMarkX=NeiPointIndex[0][i];
				for (int j = NeiPointIndex[2][0] + ZMarked[0]; j <= NeiPointIndex[2][1] + ZMarked[1]; j++)
				{
					if ((NeiPointIndexLast[0][0]>0) && (NeiPointIndexLast[0][1]<MaxIndex))
					{
						int IndexMark0 = NeiPointIndex[0][0] + (1 << OctreeDepth)*i + (1 << (2 * OctreeDepth))*j;
						int IndexMark1 = NeiPointIndex[0][1] + (1 << OctreeDepth)*i + (1 << (2 * OctreeDepth))*j;
						if (!TerminalOctant[IndexMark0].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark0].begin();
							while (PointIdIt != TerminalOctant[IndexMark0].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
						if (!TerminalOctant[IndexMark1].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark1].begin();
							while (PointIdIt != TerminalOctant[IndexMark1].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
					}
					else
					{
						int IndexMark;
						if (NeiPointIndexLast[0][0]>0)
						{

							IndexMark = NeiPointIndex[0][0] + (1 << OctreeDepth)*i + (1 << (2 * OctreeDepth))*j;
						}
						else
						{

							IndexMark = NeiPointIndex[0][1] + (1 << OctreeDepth)*i + (1 << (2 * OctreeDepth))*j;
						}
						if (!TerminalOctant[IndexMark].empty())
						{
							auto PointIdIt = TerminalOctant[IndexMark].begin();
							while (PointIdIt != TerminalOctant[IndexMark].end())
							{
								NeiborPointId.push_back(*PointIdIt);
								PointIdIt++;
							}
						}
					}
				}
			}
		}

	}

	if (NeiborPointId.empty())
	{
		cout << "error" << endl; 
		return 0;
	}
	else
	{
		auto it = NeiborPointId.begin();
		Nei_Point = *it;
		it++;
		while (it != NeiborPointId.end())
		{

			double currentMinDistance = NumericComputation<T>::SquareDistance(get_point(Nei_Point), p);
			double currentDistance = NumericComputation<T>::SquareDistance(get_point(*it), p);
			if ((currentDistance<COMPUTATION_SMALL_NUMBER) || (currentMinDistance<COMPUTATION_SMALL_NUMBER))
				return -1;
			else if (currentMinDistance>currentDistance)
				Nei_Point = *it;
			it++;
		}
	}
	return Nei_Point;
}
























#endif // !POINTLOCATOR_H
