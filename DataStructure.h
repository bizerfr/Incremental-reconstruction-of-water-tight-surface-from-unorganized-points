//变量使用了类模板和tuple设计，使用前应初始化

#include "DataArray.h"
#include "PointLocator.h"
#include "Utility.h"
#include "SurfaceDelaunayUtils.h"
#include <stack>
#include <set>
#include <unordered_set>
#include <memory>
#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

typedef Idtype Facet_Idtype;
typedef std::pair<Cell_Idtype, int>              Facet;
typedef Triple<Vertex_Idtype, Vertex_Idtype, Vertex_Idtype> Vertex_triple;
typedef std::pair<Facet_Idtype,Indextype>        Edge;
typedef std::pair<Facet, int>                    EdgeInFacet;
typedef std::pair<Vertex_Idtype, Vertex_Idtype>  Vertex_pair;

size_t hash_function_pair_int(pair<int, int> p)
{
	return hash<int>()(p.first - p.second);
}
typedef  unordered_set < pair<int, int>, decltype(hash_function_pair_int)* > set_pair_int;

enum Locate_type {
	VERTEX = 0,
	EDGE, //1
	FACET, //2
	CELL, //3
	OUTSIDE_CONVEX_HULL, //4
	OUTSIDE_AFFINE_HULL //5
};

template<typename T,typename T_INFO>
class DataStructure
{

	//HHHHHHHHHHHHHHHHHHHHHHHHH------------------delaunay-------------------HHHHHHHHHHHHHHHHHHHHHHHH
	
	typedef DataArray<Idtype> IdList;
	typedef IdList Vertex_Id_List;
	typedef IdList Cell_Id_List;
	typedef DataArray<T> DataList;
	typedef DataArray<T_INFO> InfoList;
	typedef TypeList<bool> FlagList;


//提取数据
public:
	int K_NN;//取几个最近邻点
	int get_vertex_num()
	{
		return VertexUse.get_max_id();
	}
	void get_vertex(Indextype IdV,double* PCoord,bool& Use)
	{
		const double* pp=point(IdV);
		PCoord[0]=pp[0];
		PCoord[1]=pp[1];
		PCoord[2]=pp[2];
	
		Use=VertexUse.get_element(IdV);
	}
	int get_facet_num()
	{
		return FacetUse.get_max_id();
	}
	void get_facet_vertex(Indextype IdF,const Indextype* IdVerF,bool& Use)
	{
		Vertex_triple vtriF=make_vertex_triple(get_facet_cell_link(IdF));
		*IdVerF++=vtriF.first;
		*IdVerF++=vtriF.second;
		*IdVerF=vtriF.third;
		Use=FacetUse.get_element(IdF);
	}



private:
	int Dimension;
	
	//DataList Points;
	PointLocator<T,T_INFO> PointDataSet;
	InfoList CellAttributes;//---------------（未定义） 几何单元的属性
	Vertex_Id_List Cells;//Cell的初始化均为-1
	//IdList CellLocation;//所有的Cell都是4个点，所以不需要定位
	Cell_Id_List Neighbor;
	TypeList<Cell_Idtype> CellIncident;//与vertex关联的cell
	FlagList VertexUse;//true:正被使用；flase：已经被删除         注意初始化
	FlagList CellUse;//true:正被使用；false：已经被删除                  注意初始化
	TypeList<unsigned int> MarkCellConflictState;  //1:conflict;2:boundary;0:clear
	Vertex_Idtype Infinite;
	
	
	TypeList<unsigned int> VisitedForVertexExtractor;
	

	std::stack<Vertex_Idtype> Hole_Vertices;
	std::stack<Cell_Idtype> Hole_Cells; //删除的cell

	
	
public:
	//internally used for create_star_3(faster than a tuple)
	struct iAdjacency_info{
		int v1;
		Cell_Idtype v2;
		int v3;
		Cell_Idtype v4;
		int v5;
		int v6;
		iAdjacency_info(){}
		iAdjacency_info(int a1, Cell_Idtype a2, int a3, Cell_Idtype a4, int a5, int a6) :
			v1(a1), v2(a2), v3(a3), v4(a4), v5(a5), v6(a6) {}
		void update_variables(int& a1, Cell_Idtype& a2, int& a3, Cell_Idtype& a4, int& a5, int& a6)
		{
			a1 = v1;
			a2 = v2;
			a3 = v3;
			a4 = v4;
			a5 = v5;
			a6 = v6;
		}
	};

	
	

public:
	DataStructure() :Dimension(-2)
	{
		K_NN=5;
		Cells.init(4,0);
		Neighbor.init(4,0);
	
		FacetNeighbor.init(3, 0);
		FacetSurface.init(2, 0);
		CellSurface.init(4, 0);
		NegCos=0;		
	}

	void init_point_set(int level, double* bounds[3]){ this->PointDataSet.init_point_locator(level,bounds); };
	Vertex_Idtype insert_infinite_vertex(){Infinite=create_vertex();return Infinite;};
	Vertex_Idtype nearest_point_inexact(const T* p){ return PointDataSet.nearest_inexact(p); };
	void insert_first_finite_cell();
	Vertex_Idtype insert_first_finite_cell(Vertex_Idtype &v0, Vertex_Idtype &v1, Vertex_Idtype &v2, Vertex_Idtype &v3,
		Vertex_Idtype Infinite);
	Vertex_Idtype insert_increase_dimension(Idtype star = 0);

	Vertex_Idtype insert_in_edge(Cell_Idtype c,int i,int j);


	template<typename CellIt>
	Vertex_Idtype insert_in_hole( CellIt cell_begin, CellIt cell_end,
		Cell_Idtype begin, int i);

	template<typename CellIt>
	Vertex_Idtype insert_in_hole(CellIt cell_begin, CellIt cell_end,
		Cell_Idtype begin, int i,set_pair_int boundEdge,list<Edge> ConflictBoundEdges,
		std::map<Vertex_Idtype,size_t>& VertexDeletedFacet,
		std::set<Vertex_Idtype>& VertexCreatedFacet,vector<Facet>& NewBoundarySurfacets,
		vector<Facet>& NewConflictSurfacets,bool IsIsolate,bool isInside);

	Cell_Idtype create_star_3(Vertex_Idtype v, Cell_Idtype c,
		int li, int prev_ind2 = -1);

	Cell_Idtype create_star_3(Vertex_Idtype v, Cell_Idtype c,
		int li, int prev_ind2 ,set_pair_int boundEdge,list<Edge> ConflictBoundEdges,
		std::map<Vertex_Idtype,size_t>& VertexDeletedFacet,
		std::set<Vertex_Idtype>& VertexCreatedFacet,vector<Facet>& NewBoundarySurfacets,
		vector<Facet>& NewConflictSurfacets,bool IsoIsolate,bool isInside);

	Cell_Idtype recursive_create_star_3(Vertex_Idtype v, Cell_Idtype c, int li,
		int prev_ind2, int depth);

	Cell_Idtype surface_recursive_create_star_3(Vertex_Idtype v, Cell_Idtype c, int li,
		 int depth,set_pair_int boundEdge,
		 std::map<Vertex_Idtype,size_t>& VertexDeletedFacet,
		 std::set<Vertex_Idtype>& VertexCreatedFacet,
		 std::vector<Facet>& NewCreateBoundFacet,
		 bool label,bool SideChange=false,int prev_ind2=-1);


	void update_surface_connection(Vertex_Idtype V,list<Edge> ConflictBound,vector<Facet> NewBoundFacet,
		vector<Facet>& NewBonudarySurfacets);

	Cell_Idtype non_recursive_create_star_3(Vertex_Idtype v, Cell_Idtype c, int li, int prev_ind2);

	Cell_Idtype create_star_2(Vertex_Idtype v, Cell_Idtype c, int li);

	Facet mirror_facet(Facet f) ;
	EdgeInFacet mirror_edge(EdgeInFacet E);

	template<typename IncidentCellIterator, typename IncidentFacetIterator>
	void incident_cells_3(Vertex_Idtype v, Cell_Idtype d,std::pair<IncidentCellIterator,
														IncidentFacetIterator> it) ;
			

	template<typename OutputIteratorFacet, typename OutputIteratorVertex>
	void incident_cells(Vertex_Idtype v,
		OutputIteratorFacet outputFacets, OutputIteratorVertex outputVertex) ;

	template<typename OutputIteratorCell>
	void incident_cells(Vertex_Idtype v, OutputIteratorCell outputCell) ;

	template<typename OutputIteratorCell,typename OutputIteratorVertex,typename OutputIteratorFacet>
	void incident_cells_mark_side(Vertex_Idtype v, OutputIteratorCell outputCell,
		OutputIteratorVertex outputVertex,OutputIteratorFacet outputFacet,bool& IsIso,bool& PInside) ;

	template<typename OutputIteratorVertex, typename OutputIteratorFacet>
	void adjacent_vertices_facets(Vertex_Idtype v,
		OutputIteratorFacet outputFacets, OutputIteratorVertex vertices) ;
	
	//创建vertex,即往后追加
	Vertex_Idtype create_vertex()
	{
		//Vertex_Idtype IdV;
		if (!Hole_Vertices.empty())
		{
			Vertex_Idtype IdV = Hole_Vertices.top();
			Hole_Vertices.pop();
			clear_vertex_deleted(IdV);
			clear_vertex_visited(IdV);
			CellIncident.insert_element(IdV,-1);
			return IdV;
		}
		else
		{
			//clear_vertex_visited(IdV);
			VisitedForVertexExtractor.insert_next_element(false);
			CellIncident.insert_next_element(-1);
			return VertexUse.insert_next_element(true);
		}
	};
	//创建cell,即往后追加或替代被标记删除的
	Cell_Idtype create_cell()
	{
		Cell_Idtype* t=new Cell_Idtype[Cells.size_of_tuple()];
		for (int i = 0; i < Cells.size_of_tuple(); i++)
		{
			t[i] = -1;
		}
		if (!Hole_Cells.empty())
		{
			Cell_Idtype IdC = Hole_Cells.top();
			Hole_Cells.pop();
			clear_cell_deleted(IdC);
			clear(IdC);
			Cells.set_tuple(IdC,t);
			Neighbor.set_tuple(IdC, t);

			if (Surface)
			{
				CellSurface.insert_tuple(IdC, t);
				VisitedForCellExtractor.insert_element(IdC, false);
			}

			delete []t;
			return IdC;
		}
		else
		{
			MarkCellConflictState.insert_next_element(0);
			Neighbor.insert_next_tuple(t);
			Cells.insert_next_tuple(t);

			if (Surface)
			{
				CellSurface.insert_next_tuple(t);
				VisitedForCellExtractor.insert_next_element(false);
			}

			delete []t;
			return CellUse.insert_next_element(true);

		}
	}
	Cell_Idtype create_cell(Vertex_Idtype v0, Vertex_Idtype v1, Vertex_Idtype v2, Vertex_Idtype v3)
	{
		Cell_Idtype i = create_cell();
		set_vertex(i,0,v0);
		set_vertex(i, 1, v1);
		set_vertex(i, 2, v2);
		set_vertex(i, 3, v3);
		return i;
	}
	//取dimension
	int dimension(){ return Dimension; };
	//设定dimension
	void set_dimension(int i){ this->Dimension = i; };
	//判断是否为无限点
	bool is_infinite(Vertex_Idtype IdV){ return IdV == Infinite; };
	
	//取与点IdV关联的Cell
	Cell_Idtype cell(Vertex_Idtype IdV)
	{
		return CellIncident.get_element(IdV);
	};
	//将IdV的Cell设定为IdC
	void set_cell(Vertex_Idtype IdV, Cell_Idtype IdC)
	{
		CellIncident.insert_element(IdV, IdC);//用insert更稳定
	};
	//设定相邻关系，IdC1的第I1个领域是IdC2；IdC2的第I2个邻域是IdC1
	void set_adjacency(Cell_Idtype IdC1, int I1, Cell_Idtype IdC2, int I2)
	{
		//IdC1的第I1个邻域是IdC2
		set_neighbor(IdC1,I1,IdC2);
		set_neighbor(IdC2,I2,IdC1);
	};
	//取IdC的第I个相邻Cell
	Cell_Idtype neighbor(Cell_Idtype IdC, int I)
	{
		return Neighbor.get_tuple(IdC)[I];
	};
	//设定邻域Cell,将IdC1的第I个邻域设定为IdC2
	void set_neighbor(Cell_Idtype IdC1, int I, Cell_Idtype IdC2)
	{
		Neighbor.insert_tuple_element(IdC1,I,IdC2);
	};
	//将Cell（IdC）的第I个点设定为IdV
	void set_vertex(Cell_Idtype IdC, int I, Vertex_Idtype IdV)
	{
		Cells.insert_tuple_element(IdC,I,IdV);
	};
	void set_vertices(Cell_Idtype IdC, Vertex_Idtype IdV0, Vertex_Idtype IdV1, Vertex_Idtype IdV2, Vertex_Idtype IdV3)
	{
		Cells.insert_tuple_element(IdC, 0, IdV0);
		Cells.insert_tuple_element(IdC, 1, IdV1);
		Cells.insert_tuple_element(IdC, 2, IdV2);
		Cells.insert_tuple_element(IdC, 3, IdV3);

	};
	//取IdC的第I个点
	Vertex_Idtype vertex(Cell_Idtype IdC, int I)
	{
		return Cells.get_tuple(IdC)[I];
	};

	Cell_Idtype create_face();
	Cell_Idtype create_face(Vertex_Idtype v0, Vertex_Idtype v1,
		Vertex_Idtype v2);
	Cell_Idtype create_face(Cell_Idtype f0, int i0,
		Cell_Idtype f1, int i1,
		Cell_Idtype f2, int i2);
	Cell_Idtype create_face(Cell_Idtype f0, int i0,
		Cell_Idtype f1, int i1);

	void set_point(Vertex_Idtype IdV, const T* p)
	{
		PointDataSet.set_point(IdV,p);
	};
	//取点坐标
	const const T* point(Vertex_Idtype IdV)
	{
		return PointDataSet.get_point(IdV);
	};
	//求IdV是IdC的第几个点
	int vertex_index(Cell_Idtype IdC, Vertex_Idtype IdV)
	{
		const Vertex_Idtype* vertices = Cells.get_tuple(IdC);
		for (int i = 0; i < Cells.size_of_tuple(); i++)
		{
			if (vertices[i] == IdV)
				return i;
		}
		return -1;
	};
	//取Cell中第一个没有被删除的点
	Cell_Idtype cells_begin()
	{
		for (int i = 0; i<=Cells.get_max_tuple_id(); i++)
		{
			if (!cell_is_deleted(i))
				return i;
		}
		return -1;
	};
	//判断点是否被删，被删则返回true
	bool vertex_is_deleted(Vertex_Idtype IdV)
	{
		if (VertexUse.get_element(IdV) == false)
			return true;
		else
			return false;
	};
	//判断Cell是否被删,被删则返回true
	bool cell_is_deleted(Cell_Idtype IdC)
	{
		if (CellUse.get_element(IdC) == false)
			return true;
		else
			return false;
	};
	//判断IdC中是否有点IdV,有则返回true
	bool has_vertex(Cell_Idtype IdC, Vertex_Idtype IdV)
	{
		const Vertex_Idtype* vertices = Cells.get_tuple(IdC);
		for (int i = 0; i < Cells.size_of_tuple(); i++)
		{
			if (vertices[i] == IdV)
				return true;
		}
		return false;
	};
	//重载，若存在点则将index赋给i
	//四面体Idc是否含有点Idv,并将该点相对索引（0，1，2，3）赋值给I
	bool has_vertex(Cell_Idtype IdC, Vertex_Idtype IdV, int& I)
	{
		const Vertex_Idtype* vertices = Cells.get_tuple(IdC);
		for (int i = 0; i < Cells.size_of_tuple(); i++)
		{
			if (vertices[i] == IdV)
			{
				I = i;
				return true;
			}
				
		}
		return false;
	};
	//求IdC2是IdC1的第几个临近Cell
	Indextype neighbor_index(Cell_Idtype IdC1, Cell_Idtype IdC2)
	{
		const Cell_Idtype* nei = Neighbor.get_tuple(IdC1);
		for (int i = 0; i < Neighbor.size_of_tuple(); i++)
		{
			if (nei[i] == IdC2)
				return i;
		}
	};

	void reorient();
	//计算共多少个vertex，包括已经被删除的,其实返回MaxId
	size_t num_of_vertices(){ return VertexUse.get_max_id()+1; };
	//计算共有多少cell,包括已经被删除的,其实返回MaxId
	size_t num_of_cells(){ return CellUse.get_max_id() + 1; };
	//改变第IdC的方向
	void change_orientation(Cell_Idtype IdC) ;
	//将IdC标记成conflict(1)
	void mark_in_conflict(Cell_Idtype IdC){ MarkCellConflictState.insert_element(IdC,1); };
	//判断IdC是否conflict(1)
	bool is_in_conflict(Cell_Idtype IdC){ return MarkCellConflictState.get_element(IdC)==1; };
	//将IdC标记成boundary(2)
	void mark_on_boundary(Cell_Idtype IdC){ MarkCellConflictState.insert_element(IdC, 2); };
	//判断IdC是否为boundary(2)
	bool is_on_boundary(Cell_Idtype IdC){ return MarkCellConflictState.get_element(IdC)==2; };
	//将conflict_state置clear(0)
	void clear(Cell_Idtype IdC){ MarkCellConflictState.insert_element(IdC, 0); };
	//检查IdC的MarkCellConflictState标记是否被clear(0)
	bool is_clear(Cell_Idtype IdC){ return MarkCellConflictState.get_element(IdC) == 0; };
	//delete vertex
	void delete_vertex(Vertex_Idtype IdV)
	{
		Hole_Vertices.push(IdV);
		mark_vertex_deleted(IdV);
	}

	void delete_cell(Cell_Idtype IdC)
	{
		Hole_Cells.push(IdC);
		mark_cell_deleted(IdC);
	}
	
	//删除cells，实际上只是标记
	template<class CellIt>
	void delete_cells(CellIt cell_begin, CellIt cell_end)
	{
		CellIt cit;
		for (cit = cell_begin; cit != cell_end; cit++)
		{
			//内存标记HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
			Hole_Cells.push(*cit);
			mark_cell_deleted(*cit);
			//end----------------------------------
		}
	};
	
	//将cell标记为删除
	void mark_cell_deleted(Cell_Idtype IdC){ CellUse.insert_element(IdC,false); };
	//将vertex标记为删除
	void mark_vertex_deleted(Vertex_Idtype IdV){ VertexUse.insert_element(IdV, false); };
	//将cell标记为未删除
	void clear_cell_deleted(Cell_Idtype IdC){ CellUse.insert_element(IdC, true); };
	//将vertex标记为未删除
	void clear_vertex_deleted(Vertex_Idtype IdV){ VertexUse.insert_element(IdV, true); };
	
	//检查Vertex是否被访问过
	unsigned int visited_for_vertex_extractor(Vertex_Idtype IdV){ return VisitedForVertexExtractor.get_element(IdV); };
	//将VisitedForVertexExtractor置为true，默认为false
	void mark_vertex_visited(Vertex_Idtype IdV)
	{
		
		unsigned int ii=VisitedForVertexExtractor.get_element(IdV);
		ii++;
		VisitedForVertexExtractor.insert_element(IdV,ii);
	};
	//清掉 VisitedForVertexExtractor
	void clear_vertex_visited(Vertex_Idtype IdV){ VisitedForVertexExtractor.insert_element(IdV,0); };

private:
	bool Surface;
	//bool 类型标注每个cell是否为Inside，tuple=4
	FlagList CellInsideLabel;
	//标注每个cell的四个面是否为表面，若不是则为-1，若是则为surface中facet的标识号,tuple=4
	DataArray<Facet_Idtype> CellSurface; 
	//按标识号顺序存储facet，以（Cell_Idtype,index）的形式,tuple=2
	DataArray<Cell_Idtype> FacetSurface;    
	//按标识号顺序存储facet的邻接facet，tuple=3
	DataArray<Facet_Idtype> FacetNeighbor;
	//标注facet是否被删除
	FlagList FacetUse;
	//存储被删除的facet
	std::stack<Facet_Idtype> Hole_Facets;
	//不mark cell的情况下抽取cell，做标记
	FlagList VisitedForCellExtractor;
	//记录孤立点
	std::list<Vertex_Idtype> IsolateVertices;
	//记录facet是否被visited及被visited的次数
	TypeList<unsigned int> VisitedForFacetExtractor;

	//用于测试的变量
public:
	size_t NegCos;

public:
	void get_surface_data_structure(DataArray<T>& Point,DataArray<Idtype>& FacetVertex)
	{
		for(int i=0;i<=FacetSurface.get_max_tuple_id();i++)
		{
			if(is_facet_to_delete(i))
				continue;
			else
			{
				Facet f=get_facet_cell_link(i);
				Vertex_triple vertexFacet=make_vertex_triple(f);
				Vertex_Idtype v[3]={vertexFacet.first,vertexFacet.second,vertexFacet.third};
				FacetVertex.insert_next_tuple(v);
			}

		}

		PointDataSet.get_point_data(Point);


	}
	void set_surface_extractor(bool surface){ Surface = surface; }

	bool extract_surface(){ return this->Surface; };

	void clear_facet_deleted(Facet_Idtype IdF)
	{
		FacetUse.set_element(IdF,true);
	}
	Facet_Idtype create_surface_facet()
	{
		Facet_Idtype* t = new Facet_Idtype[FacetSurface.size_of_tuple()];
		Facet_Idtype t2[3]={-1,-1,-1};
		for (int i = 0; i < FacetSurface.size_of_tuple(); i++)
		{
			t[i] = -1;
		}

		if (!Hole_Facets.empty())
		{
			Facet_Idtype IdF = Hole_Facets.top();
			Hole_Facets.pop();
			clear_facet_deleted(IdF);			
			FacetSurface.set_tuple(IdF, t);
			clear_facet_visited(IdF);
			FacetNeighbor.set_tuple(IdF, t2);
			delete []t;
			return IdF;
		}
		else
		{
			FacetNeighbor.insert_next_tuple(t2);
			FacetSurface.insert_next_tuple(t);
			VisitedForFacetExtractor.insert_next_element(0);
			delete []t;
			return FacetUse.insert_next_element(true);

		}
		
	}

	Facet_Idtype create_surface_facet(Cell_Idtype IdC, Indextype ii)
	{
		Facet_Idtype idF = create_surface_facet();
		set_surface_facet(idF,IdC,ii);
		return idF;
	}

	Facet_Idtype create_surface_facet(Facet F)
	{
		Facet_Idtype idF=create_surface_facet();
		set_surface_facet(idF,F.first,F.second);
		return idF;
	}

	Vertex_Idtype vertex_facet(Facet F,Indextype I)
	{
		Vertex_triple vtriF=make_vertex_triple(F);
		if(I==0)
			return vtriF.first;
		if(I==1)
			return vtriF.second;
		if(I==2)
			return vtriF.third ;
		return -1;
	}

	//把E从（（cell_id,0/1/2/3）,0/1/2）改成(面id,0/1/2)
	Edge turn_edgeinfacet_to_edge(EdgeInFacet E)
	{
		return Edge(get_facet_index(E.first),E.second);
	}

	//面F为（cell,id）形式，把F的面编号idF从面（F.first,F.second）改成面（IdC,ii）对应
	void change_surface_facet_link(Facet F,Cell_Idtype IdC,Indextype ii)
	{
		//---------for test------//
		//Vertex_triple vtriiF=make_vertex_triple(Facet(IdC,ii));
		//========test end=======//
		Facet_Idtype idF=CellSurface.get_tuple_element(F.first,F.second);
		CellSurface.insert_tuple_element(F.first,F.second,-1);
		set_surface_facet(idF,IdC,ii);//把面编号为IdF的面设置（IdC,ii）（处理FacetSurface），并把面（IdC,ii）编号设置为IdF（处理CellSurface）
	}
	//把面编号为IdF的面设置（IdC,ii），并把面（IdC,ii）编号设置为IdF
	void set_surface_facet(Facet_Idtype IdF, Cell_Idtype IdC, Indextype ii)
	{
		//---------for test------//
		//求面（IdC，ii）中按逆时针实践顶点编号（i,j,k）
		//Vertex_triple vtriiF=make_vertex_triple(Facet(IdC,ii));
		//========test end=======//
		FacetSurface.insert_tuple_element(IdF,0,IdC);
		FacetSurface.insert_tuple_element(IdF, 1, ii);
		CellSurface.insert_tuple_element(IdC,ii,IdF);
		
	}

	Facet get_facet_cell_link(Facet_Idtype IdF)
	{
		Facet f;
		f.first=FacetSurface.get_tuple_element(IdF,0);
		f.second = FacetSurface.get_tuple_element(IdF,1);
		return f;
	}
	Facet_Idtype get_facet_index(Facet F)
	{
		return CellSurface.get_tuple_element(F.first,F.second);
	}
	//求F1是F0的第几个相邻的surface facet
	Indextype neighbor_index_facet(Facet F0,Facet F1)
	{
		Facet_Idtype idF0=CellSurface.get_tuple_element(F0.first,F0.second);
		Facet_Idtype idF1=CellSurface.get_tuple_element(F1.first,F1.second);
		for(int i=0;i<3;i++)
		{
			if(FacetNeighbor.get_tuple_element(idF0,i)==idF1)
				return i;
		}
	}
	//求IdF1是IdF0的第几个相邻的surface facet
	Indextype neighbor_index_facet(Facet_Idtype IdF0,Facet_Idtype IdF1)
	{
		for(int i=0;i<3;i++)
		{
			if(FacetNeighbor.get_tuple_element(IdF0,i)==IdF1)
				return i;
		}
	}
	//设置面IdF0的第ii个相邻的面是IdF1
	void set_facet_neighbor(Facet_Idtype IdF0, Indextype ii, Facet_Idtype IdF1)
	{
		FacetNeighbor.insert_tuple_element(IdF0,ii,IdF1);
	}

	//IdF0的i0个相邻面是IdF1,IdF1的第i1个相邻面是IdF0
	void set_facet_adjacency(Facet_Idtype IdF0, Indextype i0, Facet_Idtype IdF1, Indextype i1)
	{
		set_facet_neighbor(IdF0,i0,IdF1);
		set_facet_neighbor(IdF1,i1,IdF0);
	}
	void set_facet_adjacency(Facet F0,Indextype i0,Facet F1,Indextype i1)
	{
		Facet_Idtype idF0=get_facet_index(F0);
		Facet_Idtype idF1=get_facet_index(F1);
		set_facet_adjacency(idF0,i0,idF1,i1);
	}
	//求点IdV是facet F的第几个点
	Indextype vertex_index_facet(Facet F,Vertex_Idtype IdV)
	{
		Vertex_triple vtriF=make_vertex_triple(F);
		if(vtriF.first==IdV)
			return 0;
		if(vtriF.second==IdV)
			return 1;
		if(vtriF.third==IdV)
			return 2;
		return -1;
	}
	void label_cell_side(Cell_Idtype IdC, bool IsInside)
	{
 		CellInsideLabel.insert_element(IdC,IsInside);
	}
	
	//将IdF标记为非表面，并不处理表面的相邻关系，若以后处理奇异边则应该加入表面相邻关系的处理
	void delete_surface_facet(Facet_Idtype IdF)
	{
		Hole_Facets.push(IdF);
		mark_facet_to_deleted(IdF);
		Facet fTemp = get_facet_cell_link(IdF);
	
		CellSurface.insert_tuple_element(fTemp.first,fTemp.second,-1);
		clear_facet_visited(IdF);
	}
	//将IdF标记为非表面，并不处理表面的相邻关系，若以后处理奇异边则应该加入表面相邻关系的处理
	void delete_surface_facet(Facet F)
	{
		Facet_Idtype idF = CellSurface.get_tuple_element(F.first,F.second);
		Hole_Facets.push(idF); //Hole_Facets存储被删除的facet的id
		mark_facet_to_deleted(idF);
		clear_facet_visited(idF);
		CellSurface.insert_tuple_element(F.first, F.second, -1);
		
		
	}
	//标记facet在conflict region内
	void mark_facet_to_deleted(Facet_Idtype IdF)
	{
		FacetUse.insert_element(IdF,false);
		
	}

	void mark_facet_to_deleted(Facet F)
	{
		Facet_Idtype idF = CellSurface.get_tuple_element(F.first, F.second);
		FacetUse.insert_element(idF, false);
	
	}

	bool is_facet_to_delete(Facet_Idtype IdF)
	{
		return !FacetUse.get_element(IdF);
	}

	bool is_facet_to_delete(Facet F)
	{
		Facet_Idtype idF = CellSurface.get_tuple_element(F.first, F.second);
		return !FacetUse.get_element(idF);
	}
	bool is_facet_in_conflict(Facet_Idtype IdF)
	{
		Facet F=get_facet_cell_link(IdF);
		return is_facet_in_conflict(F);
	}
	bool is_facet_in_conflict(Facet F)
	{
		Cell_Idtype c=F.first;
		Cell_Idtype cOpp=neighbor(c,F.second);
		bool b0=is_in_conflict(c);
		bool b1=is_in_conflict(cOpp);
		return (b0&&b1);
	}
	bool is_facet_on_conflict_boundary(Facet F)
	{
		Cell_Idtype c=F.first;
		Cell_Idtype cOpp=neighbor(c,F.second);
		bool b0=is_in_conflict(c);
		bool b1=is_in_conflict(cOpp);
		bool b2=is_on_boundary(c);
		bool b3=is_on_boundary(cOpp);
		return (b0&&b3)||(b2&&b1);
	}
	Vertex_triple make_vertex_triple(Facet& f);

	Vertex_pair make_vertex_pair(EdgeInFacet E)
	{
		Facet f=E.first;
		Vertex_triple vtrif=make_vertex_triple(f);
		if(E.second==0)
			return make_pair(vtrif.second,vtrif.third);
		if(E.second==1)
			return make_pair(vtrif.third,vtrif.first);
		if(E.second==2)
			return make_pair(vtrif.first,vtrif.second);
	}
	
	std::pair<Indextype,Indextype> make_edge_index(int I)
	{
		if(I==0)
			return make_pair(1,2);
		if(I==1)
			return make_pair(2,0);
		if(I==2)
			return make_pair(0,1);
	}

	
	bool is_label_inside(Cell_Idtype IdC)
	{
		return CellInsideLabel.get_element(IdC);
	}
	
	bool is_surface(Facet F)
	{
		if (CellSurface.get_tuple_element(F.first,F.second) == -1)
			return false;
		else
			return true;
	}
	//判断点是否在表面内
	bool is_vertex_inside_surface(Locate_type lt, Cell_Idtype c, int li, int lj)
	{
		bool pInside=true;
	
		pInside=is_label_inside(c);
		if (lt==FACET&&is_surface(Facet(c,li)))
		{
			pInside=false;
		}
		if (lt==EDGE)
		{
			//检查此边是否为表面上的边
			bool isBegin=true;
			Cell_Idtype initCellEdge=c;
			Cell_Idtype tmpCellEdge=c;
			Vertex_Idtype vj1=vertex(c,li);
			Vertex_Idtype vj2=vertex(c,lj);
			//turn around the edge vj1 vj2
			while (initCellEdge!=tmpCellEdge||isBegin)
			{
				isBegin=false;
				if (pInside!=is_label_inside(tmpCellEdge))
				{
					pInside=false;
					break;
				}
				int zz = Triangulation_utils_3::next_around_edge(vertex_index(tmpCellEdge, vj1), vertex_index(tmpCellEdge, vj2));
				tmpCellEdge=neighbor(tmpCellEdge,zz);
			}
		}
		return pInside;
	}
	//将孤立点V加入孤立点链表集合中,且不重复
	void insert_to_isolate_vertices(Vertex_Idtype V)
	{
		bool is_include=false;
		auto iiv=IsolateVertices.begin();
		while(iiv!=IsolateVertices.end())
		{
			if(*iiv==V)
			{
				is_include=true;
				break;
			}
			iiv++;
		}
		if(!is_include)
			IsolateVertices.push_back(V);
	}
	//将一个孤立点从孤立点链表中删除
	bool remove_from_isolate_vertices(Vertex_Idtype V)
	{
		
		auto iiv=IsolateVertices.begin();
		while(iiv!=IsolateVertices.end())
		{
			if(*iiv==V)
			{
				IsolateVertices.erase(iiv);
				return true;
			}
			iiv++;
		}
		return false;
	}
	//判断vertex是否为isolate vertex，是则给出距离此点最近的正对的surface facet
	bool is_vertex_isolate(Vertex_Idtype V,Facet& SurfaceFacet,bool& Inside);
	//将孤立点插入表面中
	bool insert_vertex_into_surface(Vertex_Idtype V,Facet Surfacet,bool Inside);
	//将孤立点V插入表面中（新方法）,PInside标注孤立点是在表面内还是外
	bool insert_vertex_into_surface(Vertex_Idtype V);
	//求surface facet 某一边相邻的surface facet
	bool neighbor_surface_facet(EdgeInFacet SurfaceEdge,Facet& NeighborFacet,bool& IsToDelete);
	//求surface facet的第i个相邻的surface facet
	Facet neighbor_surfacet(Facet F,Indextype I)
	{
		Facet_Idtype idF0=CellSurface.get_tuple_element(F.first,F.second);
		Facet_Idtype idF1=FacetNeighbor.get_tuple_element(idF0,I);
		//表面idF1换成（cell,0/1/2/3）形式
		return get_facet_cell_link(idF1);
	}
	Facet_Idtype neighbor_surfacet(Facet_Idtype IdF,Indextype I)
	{
		return FacetNeighbor.get_tuple_element(IdF,I);;
	}
	//绕表面三角面片上相应的边旋转四面体，旋转的四面体均为表面外的四面体
	void neighbor_surfacet_around_outside(Facet F0,Indextype I0,Facet& F1,Indextype& I1);
	//检查Cell是否被访问过
	bool visited_for_cell_extractor(Cell_Idtype IdC){ return VisitedForCellExtractor.get_element(IdC); };
	//将VisitedForCellExtractor置为true，默认为false
	void mark_cell_visited(Cell_Idtype IdC){ VisitedForCellExtractor.insert_element(IdC, true); };
	//清掉VisitedForCellExtractor
	void clear_cell_visited(Cell_Idtype IdC){ VisitedForCellExtractor.insert_element(IdC, false); };
	//检查facet是否被访问过
	int visited_for_facet_extractor(Facet_Idtype IdF){return VisitedForFacetExtractor.get_element(IdF);};
	int visited_for_facet_extractor(Facet F)
	{
		Facet_Idtype idF=CellSurface.get_tuple_element(F.first,F.second);
		return VisitedForFacetExtractor.get_element(idF);
	};
	//将facet设置为已经访问，true
	void mark_facet_visited(Facet_Idtype IdF){
		Indextype times=VisitedForFacetExtractor.get_element(IdF)+1;
		VisitedForFacetExtractor.insert_element(IdF,times);};
		
	void mark_facet_visited(Facet F)
	{
		Facet_Idtype idF=CellSurface.get_tuple_element(F.first,F.second);
		Indextype times=VisitedForFacetExtractor.get_element(idF)+1;
		VisitedForFacetExtractor.insert_element(idF,times);
	};
	//将facet设置为未访问，false
	void clear_facet_visited(Facet_Idtype IdF){VisitedForFacetExtractor.insert_element(IdF,0);};
	void clear_facet_visited(Facet F)
	{
		Facet_Idtype idF=CellSurface.get_tuple_element(F.first,F.second);
		VisitedForFacetExtractor.insert_element(idF,0);
	}
	//求某一点在surface上的umbrella，不做标记
	template<class OutputIteratorFacet>
	void incident_surface_facet(Vertex_Idtype IdV, OutputIteratorFacet FacetDelete, OutputIteratorFacet FacetUndelete);
	//求某一点在surface上的umbrella，不做标记
	template<class IteratorFacet>
	void incident_surface_facet(Vertex_Idtype IdV, IteratorFacet Surfacets);
	//求某一点的incident vertex和在surface上的umbrella，不做标记
	template<class IteratorVertex,class IteratorFacet>
	void incident_vertex_and_surface_facet(Vertex_Idtype IdV,IteratorVertex IncidentVertices,IteratorFacet IncidentFacets);
	//不对cell做conflict mark的前提下求incident_cell 3d
	template<typename IncidentCellIterator, class OutputIteratorFacet>
	void incident_cells_3_withoutmark(Vertex_Idtype IdV, Cell_Idtype d, IncidentCellIterator IncidentCell,
		OutputIteratorFacet, OutputIteratorFacet FacetUndelete);
	//不对cell做conflict mark的前提下求incident_cell 3d, 不标记是否被删除
	template<typename IteratorCell, class IteratorFacet>
	void incident_cells_surfacet_3_withoutmark(Vertex_Idtype IdV, Cell_Idtype d,IteratorCell IncidentCells,
		IteratorFacet Surfacets);
	//求距离某一点最近的surface facet
	bool surface_facet_nearest_to_point(const T* p, Facet& NearestFacet, std::vector<Facet>& FacetDelete, std::vector<Facet>& FacetUndelete);
	//求incident surface facet中的某个fact使与point相连的另外的两条边最长的一个最短
	template<class IteratorIncidentFacet>
	Facet incident_facet_with_min_longest_edge(const T* p, IteratorIncidentFacet FacetBegin, IteratorIncidentFacet FacetEnd, Vertex_Idtype VertexIncident);

	//求point到vertex的距离
	double distant_to_vertex(const T* p,Vertex_Idtype IdV);

	//求交界面，使用Gabriel的规则决定替换面的延展（主要针对符合pseudo-concave的surface），凸延展
	void find_conflict_surfacet(const T* P,Facet InitFacet,list<Edge>& BoundEdges,set_pair_int& BoundEdge,Facet& Begin,
		std::set<Vertex_Idtype>& VertexOnBoundary,map<Vertex_Idtype,size_t>& VertexOnConflictFacet,
		bool PInside,vector<Facet_Idtype>& SurfacetsInConflict,bool& quit);

	//求点到面三个点的最大距离
	double longest_edge_to_facet(const T* p,Facet f);
	//求点到三角面片三个顶点的最小距离
	double shortest_edge_to_facet(const T* p,Facet f);
	//求两个矢量面二面角的余弦值
	double cos_dihedral(Facet F0,Facet F1)
	{
		Vertex_triple vf0=make_vertex_triple(F0);
		Vertex_triple vf1=make_vertex_triple(F1);
		return GeometricTraits<T>::cos_dihedral(point(vf0.first),point(vf0.second),point(vf0.third),
			point(vf1.first),point(vf1.second),point(vf1.third));
	}
	//求某一点p到面Facet形成的四面体的另外三个面与Facet组成的最大二面角的余弦值（对应的余弦值最小）
	T cos_min_dihedral_point_to_facet(const T* P,Facet F)
	{
		Vertex_triple vf=make_vertex_triple(F);
		T cosD0=GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.first),P,point(vf.second));
		T cosD1=GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.second),P,point(vf.third));
		T cosD2=GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.third),P,point(vf.first));
		T maxCos=cosD0<cosD1?(cosD0<cosD2?cosD0:cosD2):(cosD1<cosD2?cosD1:cosD2);
		return maxCos;
	}
	//求某一点P到面Facet形成的四面体的另外三个内面与Facet组成的三个二面角的余弦值
	void cos_dihedral_point_to_facet(const T* P,Facet F,T* CosDihe)
	{
		Vertex_triple vf=make_vertex_triple(F);
		
		CosDihe[0]=GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.second),P,point(vf.third));
		CosDihe[1]=GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.third),P,point(vf.first));
		CosDihe[2]=GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.first),P,point(vf.second));
		
	}
	//求某一点P到面Facet关于某个边EdgeInFacet(F,i)形成二面角的余弦值
	double cos_dihedral_point_to_Edge(const T* P,Facet F,Indextype I)
	{
		Vertex_triple vf=make_vertex_triple(F);
		if(I==0)
			return GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.second),P,point(vf.third));
		if(I==1)
			return GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.third),P,point(vf.first));
		if(I==2)
			return GeometricTraits<T>::cos_dihedral(point(vf.first),point(vf.second),point(vf.third),point(vf.first),P,point(vf.second));
	}
	//处理孤立点
	void isolate_vertices(std::map<Vertex_Idtype,size_t> VertexWithSurfacetDeleted,std::set<Vertex_Idtype> VertexWithSurfacetCreated);
	//检查孤立点
	bool has_isolate_vertex(std::map<Vertex_Idtype,size_t> VertexWithSurfacetDeleted,std::set<Vertex_Idtype> VertexWithSurfacetCreated);

	//判断面F的三个顶点,至少有一个在集合Nearest中
	bool is_facet_in_neighbor(Facet F,list<Vertex_Idtype> Nearest)
	{
		Vertex_triple vtriF=make_vertex_triple(F);
		return(is_in_list(vtriF.first,Nearest.begin(),Nearest.end())||is_in_list(vtriF.second,Nearest.begin(),Nearest.end())||
			is_in_list(vtriF.third,Nearest.begin(),Nearest.end()));
	}

	//插入一个点，更新最近邻列表
	void refresh_K_nearest(const T* p,Vertex_Idtype v,list<Vertex_Idtype>& nearest);
	//比较一个点到两个点的最近点
	Vertex_Idtype nearest_vertex(const T* p,Vertex_Idtype v,Vertex_Idtype w);
	//将孤立点连入表面
	bool link_isolate_to_surface(Vertex_Idtype V,Facet SeedFacet,bool PInside);
	//寻找孤立点相关的交界面
	void find_isolate_conflict_surfacet(Vertex_Idtype V, Facet SeedFacet,bool PInside,std::set<Vertex_Idtype>& VertexOnBoundary,
		set<Vertex_Idtype>& VertexOnConflictFacet,vector<Facet>& Conflicts,
		list<Edge>& ConflictBoundEdges,vector<Facet>& NewCreateSurfacets,bool& quit);

	//用迭代的方式将孤立点连入surface
	void recursive_find_isolate_conflict_surfacet(Vertex_Idtype V, Facet SeedFacet,bool PInside,std::set<Vertex_Idtype>& VertexOnBoundary,
		set<Vertex_Idtype>& VertexOnConflictFacet,vector<Facet>& Conflicts,
		list<Edge>& ConflictBoundEdges,vector<Facet>& NewCreateSurfacets,int prev_ind3=-1);
	//求孤立点类交界面的种子三角面片
	bool seed_facet_of_isolate_vertex(Vertex_Idtype V,Facet InitFacet,bool PInside,Facet& SeedFacet);
	//求最近邻的方法，从candidate中选择
	void nearest_surfacet(const T* P,vector<Facet> SurfacetsInternalConflict,
		vector<Facet> SurfacetsIncidentKNN,vector<Facet> SurfacetsIncidentSparse,
		bool PInside,Facet& NearestFacet,vector<Facet_Idtype>& SurfacetsInConflict,
		std::map<Vertex_Idtype,size_t>& verticesWithDeletedSurfacet);
	//求某一点的邻域内的surface facet，按照点到面所形成二面角的最大角最小优先排列
	void nearest_surface_facet(const T* p,list<Vertex_Idtype> NearestVertices,
		vector<Facet> InitFacet,bool PInside,Facet& NearestFacet,bool& IsConcave);
	//求最近邻的新方法，加入几何约束和拓扑约束
	void nearest_surface_facet(const T* p,list<Vertex_Idtype> NearestVertices,
		vector<Facet> InitFacet,Triple<Facet,double,double> InitFacetConflict,bool PInside,Facet& NearestFacet);

	template<typename IterList,typename ValueTypeList>
	static bool is_in_list(ValueTypeList v,IterList IBegin,IterList Iend)
	{
		for(;IBegin!=Iend;IBegin++)
		{
			if(*IBegin==-1)
				break;
			if(*IBegin==v)
				return true;
		}
		return false;
	}
	template<typename IterList,typename ValueTypeList>
	static int rank_in_list(ValueTypeList v,IterList IBegin,IterList Iend)
	{
		int i=0;
		for(;IBegin!=Iend;IBegin++)
		{
			if(*IBegin==v)
				return i;
			i++;
		}
		return -1;
	}
	bool is_edge_in_surface(Cell_Idtype C,Indextype I0,Indextype I1);
	//sculpture操作+表面数据结构更新,C为要sculpture的cell，NumOfSurface为原有的surface facet的个数，SurfaceIndex为这些surface facet对应的标记号。
	void sculpture(Cell_Idtype C,int NumOfSuface,int* SurfaceIndex,vector<Facet>& NewCreateSurfacet);
	//inflate操作+表面数据结构更新
	void inflate(Cell_Idtype C,int NumOfSurface,int* SurfaceIndex,vector<Facet>& NewCreateSurfacet);
	//迭代sculpture，优化表面,flip、连入孤立点仅以E(S)的增量为限制。MarkIso表示是否可以产生孤立点，MarkIso=true时不能产生孤立点如果可以则最多产生一个孤立点
	//且不能将孤立点VertexIso连入表面
	template<typename IteratorFacet>
	bool iterative_sculpture_with_mark(Facet F,Vertex_Idtype VertexIso,IteratorFacet ItSurF,Vertex_Idtype IdV,Facet& FCur,bool MarkIso);
	//迭代inflate，优化表面,flip、连入孤立点仅以E(S)的增量为限制，最多产生一个孤立点
	template<typename IteratorFacet>
	bool iterative_inflate_with_mark(Facet F,Vertex_Idtype VertexIso,IteratorFacet ItSurF,Vertex_Idtype IdV,Facet& FCur,bool MarkIso);
	//迭代更新K-NN的incident surfacet，用inflate、sculpture
	template<typename IteratorFacet>
	void update_surfacets_KNN(bool PInside,Vertex_Idtype VertexIso,list<Vertex_Idtype> NearestVertices,vector<Facet> InitFacets,IteratorFacet ItSurF,IteratorFacet ItSurSparse);


private:
	//存储插入点时候,改变的信息,用于异常处理后复原

	//存储递归过程新生成的cell
	unique_ptr<stack<Cell_Idtype> > unique_pCnews_stack;

	//存储递归过程中新生成和删除的表面
	struct MarkedSurfacet
	{
		Facet f;  //递归过程中新生成和删除的表面
		bool mark; //true为递归中新生成表面,false为递归中删除的表面

		//存储递归中被删除面的邻接面,留恢复时候使用;对递归过程中新生成的表面不需要,
		Facet_Idtype original_nei0; 
		Facet_Idtype original_nei1;
		Facet_Idtype original_nei2;

		MarkedSurfacet(Facet f_, bool mark_,
			Facet_Idtype original_nei0_,Facet_Idtype original_nei1_,Facet_Idtype original_nei2_): 
		    f(f_),mark(mark_),
			original_nei0(original_nei0_),original_nei1(original_nei1_),original_nei2(original_nei2_){}

		MarkedSurfacet(Cell_Idtype c,Indextype i, bool mark_,
			Facet_Idtype original_nei0_,Facet_Idtype original_nei1_,Facet_Idtype original_nei2_): f(c,i),mark(mark_),
		    original_nei0(original_nei0_),original_nei1(original_nei1_),original_nei2(original_nei2_){}
	};
	unique_ptr<stack<MarkedSurfacet> > unique_pCreated_DeletedSurfacets_stack;

	//存储递归过程中改变的(c,cnew,li),三个一组
	struct ChangedSurfacet
	{
		Cell_Idtype c;
		Cell_Idtype cnew;
		Indextype li;
		ChangedSurfacet(Cell_Idtype c_,Cell_Idtype cnew_,Indextype li_):
			c(c_),cnew(cnew_),li(li_){} 
	};
	unique_ptr<stack<ChangedSurfacet> > unique_pChangedSurfacets_stack;

	//存储递归过程中改变的(c,c_li,原来的neighbor_index(c_li, c)),三个一组
	struct ChangedCellAjcacency 
	{
		Cell_Idtype c;
		Cell_Idtype c_li;
		Indextype original_index;
		ChangedCellAjcacency(Cell_Idtype c_,Cell_Idtype c_li_,Indextype original_index_):
			c(c_),c_li(c_li_),original_index(original_index_){}
	};
	unique_ptr<stack<ChangedCellAjcacency> > unique_pChangedCellAjcacency_stack;

	//存储递归过程中改变的(顶点,与之关联cell),两个一组
	struct ChangedCellIncident
	{
		Vertex_Idtype changed_cell_vertex;
		Cell_Idtype original_cell;
		ChangedCellIncident(Vertex_Idtype changed_cell_vertex_,Cell_Idtype original_cell_):
			changed_cell_vertex(changed_cell_vertex_),original_cell(original_cell_){}
	};
	unique_ptr<stack<ChangedCellIncident> > unique_pChangedCellIncident_stack;

	void ini_uniqueptr_of_Changeds();
	void return_to_before_insert();

};



template<typename T, typename T_INFO>
void DataStructure<T,T_INFO>::
ini_uniqueptr_of_Changeds()
{
	unique_pCnews_stack=
		unique_ptr<stack<Cell_Idtype> >(new stack<Cell_Idtype>);
	unique_pCreated_DeletedSurfacets_stack=
		unique_ptr<stack<MarkedSurfacet> >(new stack<MarkedSurfacet>);
	unique_pChangedSurfacets_stack=
		unique_ptr<stack<ChangedSurfacet> >(new stack<ChangedSurfacet>);
	unique_pChangedCellAjcacency_stack=
		unique_ptr<stack<ChangedCellAjcacency> >(new stack<ChangedCellAjcacency>);
	unique_pChangedCellIncident_stack=
		unique_ptr<stack<ChangedCellIncident> >(new stack<ChangedCellIncident>);

}


template<typename T, typename T_INFO>
void DataStructure<T,T_INFO>::
return_to_before_insert()
{
	//递归过程中新生成的面删除,被删除的面重新恢复,并且保证面id跟以前一样
	while (!unique_pCreated_DeletedSurfacets_stack->empty())
	{
		auto tmp=unique_pCreated_DeletedSurfacets_stack->top();
		unique_pCreated_DeletedSurfacets_stack->pop();
		if (tmp.mark)
		{
			delete_surface_facet(tmp.f);
		}
		else
		{
			create_surface_facet(tmp.f);
			Facet_Idtype fid=CellSurface.get_tuple_element((tmp.f).first,(tmp.f).second);

			Facet_Idtype original_nei0=tmp.original_nei0;
			Facet_Idtype original_nei1=tmp.original_nei1;
			Facet_Idtype original_nei2=tmp.original_nei2;
			set_facet_neighbor(fid,0,original_nei0);
			set_facet_neighbor(fid,1,original_nei1);
			set_facet_neighbor(fid,2,original_nei2);
		}
	}
	//递归过程中因change_surface_facet_link(Facet(c,li),cnew,li)改变的状态恢复
	while(!unique_pChangedSurfacets_stack->empty())
	{
		auto tmp=unique_pChangedSurfacets_stack->top();
		unique_pChangedSurfacets_stack->pop();
		Cell_Idtype c=tmp.c;
		Cell_Idtype cnew=tmp.cnew;
		Indextype li=tmp.li;
		change_surface_facet_link(Facet(cnew,li),c,li);
	}
	//递归过程中因set_adjacency(cnew, li, c_li, neighbor_index(c_li, c))改变的状态恢复
	while (!unique_pChangedCellAjcacency_stack->empty())
	{
		auto tmp=unique_pChangedCellAjcacency_stack->top();
		unique_pChangedCellAjcacency_stack->pop();
		Cell_Idtype c=tmp.c;
		Cell_Idtype c_li=tmp.c_li;
		Indextype original_index=tmp.original_index;
		set_neighbor(c_li,original_index,c);
	}
	//递归过程中因set_cell(vertex(cnew, ii), cnew);改变的状态恢复
	while (!unique_pChangedCellIncident_stack->empty())
	{
		auto tmp=unique_pChangedCellIncident_stack->top();
		unique_pChangedCellIncident_stack->pop();
		Vertex_Idtype changed_cell_vertex=tmp.changed_cell_vertex;
		Cell_Idtype original_cell=tmp.original_cell;
		set_cell(changed_cell_vertex,original_cell);
	}
	while (!unique_pCnews_stack->empty())
	{
		Cell_Idtype cnew=unique_pCnews_stack->top();
		unique_pCnews_stack->pop();
		delete_cell(cnew);
	}
}


template<typename T, typename T_INFO>
Vertex_Idtype DataStructure<T,T_INFO>::
insert_increase_dimension(Vertex_Idtype star)
{
	Vertex_Idtype idVertex = create_vertex();
	int dim = dimension();
	set_dimension(dim+1);
	switch (dim)
{
	
	case -2:
	    // insertion of the first vertex
	    // ( geometrically : infinite vertex )
	{
		Infinite=idVertex;
	
		Cell_Idtype c = create_face(idVertex, -1, -1);
		set_cell(idVertex,c);
		break;
	}

	case -1:
		// insertion of the second vertex
		// ( geometrically : first finite vertex )
	{
		Cell_Idtype d = create_face(idVertex, -1, -1);
		set_cell(idVertex,d);
		set_adjacency(d, 0, cell(star), 0);
		break;
	}

	case 0:
		// insertion of the third vertex
		// ( geometrically : second finite vertex )
	{
		Cell_Idtype c = cell(star);
		Cell_Idtype d = neighbor(c,0);

		set_vertex(c,1, vertex(d,0));
		set_vertex(d,1, idVertex);
		set_neighbor(d,1, c);
		Cell_Idtype e = create_face(idVertex, star, -1);
		set_adjacency(e, 0, c, 1);
		set_adjacency(e, 1, d, 0);

		set_cell(idVertex,d);
		break;
	}

	case 1:
		// general case : 4th vertex ( geometrically : 3rd finite vertex )
		// degenerate cases geometrically : 1st non collinear vertex
	{
		Cell_Idtype c = cell(star);
		int i = vertex_index(c,star); // i== 0 or 1
		
		int j = (i == 0) ? 1 : 0;
		Cell_Idtype d = neighbor(c,j);

		set_vertex(c,2, idVertex);

		Cell_Idtype e = neighbor(c,i);
		Cell_Idtype cnew = c;
		Cell_Idtype enew = -1;

		while (e != d){
			enew = create_cell();
			set_vertex(enew,i, vertex(e,j));
			set_vertex(enew,j, vertex(e,i));
			set_vertex(enew,2, star);

			set_adjacency(enew, i, cnew, j);
			// false at the first iteration of the loop where it should
			// be neighbor 2
			// it is corrected after the loop
			set_adjacency(enew, 2, e, 2);
			// neighbor j will be set during next iteration of the loop

			set_vertex(e,2, idVertex);

			e = neighbor(e,i);
			cnew = enew;
		}

		set_vertex(d,2, idVertex);
		set_adjacency(enew, j, d, 2);

		// corrections for star->cell() :
		c = cell(star);
		//set_neighbor(c, 2, c->neighbor(i)->neighbor(2)); 
		set_neighbor(c,2,neighbor(neighbor(c, i), 2));
		set_neighbor(c,j, d);

		set_cell(idVertex,d);

		break;
	}

	case 2:
		// general case : 5th vertex ( geometrically : 4th finite vertex )
		// degenerate cases : geometrically 1st non coplanar vertex
	{
		// used to store the new cells, in order to be able to traverse only
		// them to find the missing neighbors.
		std::vector<Cell_Idtype > new_cells;
		new_cells.reserve(16);

		// allowed since the dimension has already been set to 3
		//Cell_iterator it = cells_begin();
		Cell_Idtype iCell = cells_begin();
		
		// ok since there is at least one ``cell''
		set_cell(idVertex,iCell); 
		for (; iCell <= num_of_cells(); ++iCell) {
			// Here we must be careful since we create_cells in a loop controlled
			// by an iterator.  So we first take care of the cells newly created
			// by the following test :
			if (cell_is_deleted(iCell))
				continue;
			if (neighbor(iCell,0) == -1)
				continue;
			set_neighbor(iCell,3, -1);
			set_vertex(iCell,3, idVertex);
			if (!has_vertex(iCell,star)) {
				Cell_Idtype cnew = create_cell(vertex(iCell,0), vertex(iCell,2),
					vertex(iCell,1), star);
				// The Intel compiler has a problem with passing "it" directly to
				// function "set_adjacency": the adjacency is not changed.
				Cell_Idtype ch_it = iCell;
				set_adjacency(cnew, 3, ch_it, 3);
				set_neighbor(cnew,0, -1);
				new_cells.push_back(cnew);
			}
		}

		// traversal of the new cells only, to add missing neighbors
		for (typename std::vector<Cell_Idtype>::iterator ncit = new_cells.begin();
			ncit != new_cells.end(); ++ncit) {
			Cell_Idtype n = neighbor(*ncit,3); // opposite to star
			for (int i = 0; i<3; i++) {
				int j;
				if (i == 0) j = 0;
				else j = 3 - i; // vertex 1 and vertex 2 are always switched when
				// creating a new cell (see above)
				//Cell_Idtype c = n->neighbor(i)->neighbor(3);
				Cell_Idtype c = neighbor(neighbor(n,i),3);
				if (c != -1) {
					// i.e. star is not a vertex of n->neighbor(i)
					set_neighbor(*ncit,j, c);
					// opposite relation will be set when ncit arrives on c
					// this avoids to look for the correct index
					// and to test whether *ncit already has neighbor i
				}
				else {
					// star is a vertex of n->neighbor(i)
					set_adjacency(*ncit, j, neighbor(n,i), 3);//neighbor opposite to idVertex
				}
			}
		}
	}
}// end switch

return idVertex;
}

template<typename T, typename T_INFO>
Vertex_Idtype DataStructure<T, T_INFO>::
insert_in_edge(Cell_Idtype c, int i, int j)
{
	Vertex_Idtype v = create_vertex();
	Cell_Idtype cnew = create_face(v, vertex(c, 1), -1);
	set_cell(vertex(c, 1), cnew);
	set_vertex(c, 1, v);
	set_adjacency(cnew, 0, neighbor(c, 0), 1);
	set_adjacency(cnew, 1, c, 0);
	set_cell(v, cnew);
	return v;

}


template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
reorient()
{
	for (rsize_t i = 0; i < num_of_cells(); i++)
	{
		if (!cell_is_deleted(i))
			change_orientation(i);
	}
}

template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
change_orientation(Cell_Idtype c) 
{
	Vertex_Idtype tmp_v = vertex(c,0);
	set_vertex(c,0, vertex(c,1));
	set_vertex(c,1, tmp_v);
	Cell_Idtype tmp_c = neighbor(c,0);
	set_neighbor(c,0, neighbor(c,1));
	set_neighbor(c,1, tmp_c);
}

template<typename T, typename T_INFO>
template<typename CellIt>
Vertex_Idtype DataStructure<T, T_INFO>::
insert_in_hole(CellIt cell_begin, CellIt cell_end,
Cell_Idtype begin, int i)
{

	Vertex_Idtype newv = create_vertex();

	Cell_Idtype cnew;
	if (dimension() == 3)
		cnew = create_star_3(newv, begin, i);
	else
		cnew = create_star_2(newv, begin, i);

	set_cell(newv,cnew);
	delete_cells(cell_begin, cell_end);
	return newv;
}

//要求影响域的cells相连，其边界也相连,p在影响域内
//填充Delaunay影响域，维护表面的完整性，以替换面为边界对hole进行内外标注
//cell_begin为迭代器中初始删除的cell;cell_end为最后删除的cell;
//（begin,i）是影响域边界上面，即begin属于影响域，但neighbor(begin,i)不属于影响域
//boundEdge为被替换面边界线，以端点(vertex_id,vertex_id)表示；ConflictBoundEdges也是边界线，以Edge(面id,0/1/2)表示
//VertexDeletedFacet为(被替换面上的顶点id,与之关联的表面三角面片被删除的个数);
//VertexCreatedFacet为边界线上的点;
//NewBoundarySurfacets(cell,0/1/2/3)表示与新形成表面共替换面边界线的表面;NewConflictSurfacets为新形成的非替换面的表面三角面片;
//IsoIsolate表示p的SeedFacet是否在影响域内或边界上；isInside为p关于表面的位置关系
template<typename T, typename T_INFO>
template<typename CellIt>
Vertex_Idtype DataStructure<T, T_INFO>::
insert_in_hole(CellIt cell_begin, CellIt cell_end,
		Cell_Idtype begin, int i, set_pair_int boundEdge,list<Edge> ConflictBoundEdges,
		std::map<Vertex_Idtype,size_t>& VertexDeletedFacet,
		std::set<Vertex_Idtype>& VertexCreatedFacet,vector<Facet>& NewBoundarySurfacets,
		vector<Facet>& NewConflictSurfacets,bool IsoIsolate,bool isInside)
{
	Vertex_Idtype newv = create_vertex();
	ini_uniqueptr_of_Changeds();
	try
	{
		Cell_Idtype cnew;
		if (dimension() == 3)
			cnew = create_star_3(newv, begin, i,-1,boundEdge,ConflictBoundEdges,VertexDeletedFacet,
				VertexCreatedFacet,NewBoundarySurfacets,NewConflictSurfacets,IsoIsolate,isInside);	
		else
			cnew = create_star_2(newv, begin, i);
		set_cell(newv, cnew);
		delete_cells(cell_begin, cell_end);
		update_surface_connection(newv,ConflictBoundEdges,NewConflictSurfacets,NewBoundarySurfacets);
		return newv;
	}
	catch(exception& )
	{
		//恢复插入该点之前状态
		delete_vertex(newv);
		return_to_before_insert();
		return -1;
	}

}


template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
create_star_3(Vertex_Idtype v, Cell_Idtype c,
	int li, int prev_ind2 = -1)
{
	return recursive_create_star_3(v, c, li, prev_ind2, 0);
}

//要求影响域的cells相连，其边界也相连,新顶点v在影响域内
//填充Delaunay影响域，维护表面的完整性，以替换面为边界对hole进行内外标注
//v新插入顶点编号；（c,li）是影响域边界上面，即c属于影响域，但neighbor(c,li)不属于影响域
//prev_ind2即previous index
//boundEdge为被替换面边界线，以端点(vertex_id,vertex_id)表示；ConflictBoundEdges也是边界线，以Edge(面id,0/1/2)表示
//VertexDeletedFacet为(被替换面上的顶点id,与之关联的表面三角面片被删除的个数);
//VertexCreatedFacet为边界线上的点；NewBoundarySurfacets为新形成的替换面上的三角面片;
//NewConflictSurfacets为新形成的非替换面的表面三角面片；IsoIso表示p的SeedFacet是否在影响域内或边界上；isInside为p关于表面的位置关系
template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
create_star_3(Vertex_Idtype v, Cell_Idtype c,
int li, int prev_ind2,set_pair_int boundEdge,list<Edge> ConflictBoundEdges,
std::map<Vertex_Idtype,size_t>& VertexDeletedFacet,
std::set<Vertex_Idtype>& VertexCreatedFacet,vector<Facet>& NewBonudarySurfacets,
vector<Facet>& NewConflictSurfacets,bool IsIso,bool isInside)
{
	Cell_Idtype cnew;

	cnew= surface_recursive_create_star_3(v, c, li, 0,boundEdge,VertexDeletedFacet,
		VertexCreatedFacet,NewConflictSurfacets,true,false, prev_ind2);	

	return cnew;
}

template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
recursive_create_star_3(Vertex_Idtype v, Cell_Idtype c, int li,
						int prev_ind2, int depth)
{
	if (depth == 100) return non_recursive_create_star_3(v, c, li, prev_ind2);


	Cell_Idtype cnew = create_cell(vertex(c,0),
		vertex(c,1),
		vertex(c,2),
		vertex(c,3));
	set_vertex(cnew,li, v);
	Cell_Idtype c_li = neighbor(c,li);
	set_adjacency(cnew, li, c_li, neighbor_index(c_li,c));

	// Look for the other neighbors of cnew.
	for (int ii = 0; ii<4; ++ii) {
		//if (ii == prev_ind2 || neighbor(cnew,ii) != -1)
		if (ii==prev_ind2||ii==li)
			continue;
		//cnew->vertex(ii)->set_cell(cnew);
		set_cell(vertex(cnew,ii),cnew);
		// Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
		Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
		Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
		Cell_Idtype cur = c;
		int zz = ii;
		Cell_Idtype n = neighbor(cur,zz);
		// turn around the oriented edge vj1 vj2
		while (is_in_conflict(n)) {
			cur = n;
			zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
			n = neighbor(cur,zz);
		}
		// Now n is outside region, cur is inside.
		clear(n); // Reset the flag for boundary cells.

		int jj1 = vertex_index(n,vj1);
		int jj2 = vertex_index(n,vj2);
		Vertex_Idtype vvv = vertex(n, Triangulation_utils_3::next_around_edge(jj1, jj2));
		Cell_Idtype nnn = neighbor(n, Triangulation_utils_3::next_around_edge(jj2, jj1));
		int zzz = vertex_index(nnn,vvv);
		if (nnn == cur) {
			// Neighbor relation is reciprocal, ie
			// the cell we are looking for is not yet created.
			nnn = recursive_create_star_3(v, nnn, zz, zzz, depth + 1);
		}

		set_adjacency(nnn, zzz, cnew, ii);
	}

	return cnew;
}


//要求影响域的cells相连，其边界也相连,新顶点v在影响域内
//填充Delaunay影响域，维护表面的完整性，以替换面为边界对hole进行内外标注
//v新插入顶点编号；（c,li）是影响域边界上面，即c属于影响域，但neighbor(c,li)不属于影响域
//depth递归深度；boundEdge为被替换面边界线，以端点(vertex_id,vertex_id)表示；
//VertexDeletedFacet为(被替换面上的顶点id,与之关联的表面三角面片被删除的个数);
//VertexCreatedFacet为边界线上的点；NewCreateBoundFacet为新形成的非替换面上的三角面片;
//label表示新产生cell是否inside,SideChange:当有向边vj1->vj2（vj1、vj2为该函数内变量）为被替换面边界线为true,prev_ind2影响内四面体已经处理过的相对顶点
template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
surface_recursive_create_star_3(Vertex_Idtype v, Cell_Idtype c, int li,
		 int depth,set_pair_int boundEdge,
		 std::map<Vertex_Idtype,size_t>& VertexDeletedFacet,
		 std::set<Vertex_Idtype>& VertexCreatedFacet,std::vector<Facet>& NewCreateBoundFacet,
         bool label,bool SideChange,int prev_ind2)
{


	//---------for test-------//
	/*cout<<"recursive depth: "<<depth<<" starting,"<<endl;*/
	//=========test end=======//

	Cell_Idtype cnew = create_cell(vertex(c, 0),
		vertex(c, 1),
		vertex(c, 2),
		vertex(c, 3));
	//将Cell（cnew）的第li个点设定为v
	set_vertex(cnew, li, v);
    unique_pCnews_stack->push((cnew));

	Cell_Idtype c_li = neighbor(c, li);
	
	//--------------------提取交界面的表面---------------
    if (SideChange)
	{
		if (label)
		{
			create_surface_facet(cnew,prev_ind2);
			unique_pCreated_DeletedSurfacets_stack->push(MarkedSurfacet(cnew,prev_ind2,true,-1,-1,-1));
			//---------for test-------//
			//Vertex_triple vtriFTT=make_vertex_triple(Facet(cnew,prev_ind2));
			//cout<<"!!create facet opposite to pre_ind2: ("<<vtriFTT.first<<","<<vtriFTT.second<<","<<vtriFTT.third<<")"
			//	<<"/("<<cnew<<","<<prev_ind2<<")"<<endl;
			//=========test end=======//
		
		}
		
	}


	label_cell_side(cnew,label);

	//----------for test------------//
	//Vertex_triple vtrT=make_vertex_triple(Facet(c,li));
	//cout<<"create cell "<<cnew<<"("<<vertex(cnew, 0)<<","<<vertex(cnew, 1)<<","<<vertex(cnew, 2)<<","<<vertex(cnew, 3)<<")"
	//	<<"[label true:"<<label<<"] opposite to facet ("<<c<<","<<li<<")/("<<vtrT.first<<","<<vtrT.second<<","<<vtrT.third<<")"
	//	<<"   c:"<<c<<"["<<is_label_inside(c)<<"]("<<vertex(c, 0)<<","<<vertex(c, 1)<<","<<vertex(c, 2)<<","<<vertex(c, 3)<<") "
	//	<<"c_li:"<<c_li<<"["<<is_label_inside(c_li)<<"]("<<vertex(c_li, 0)<<","<<vertex(c_li, 1)<<","<<vertex(c_li, 2)<<","<<vertex(c_li, 3)<<")"<<endl;
	//===========test end===========//

	//--------------判断boundary facet的surface facet状态的变化---------
	if (label == is_label_inside(c_li)) //cnew和c_li都在表面内或表面外
	{
		if (is_surface(Facet(c, li)))
		{
			Facet_Idtype fid=CellSurface.get_tuple_element(c,li);
			Facet_Idtype original_nei0=neighbor_surfacet(fid,0);
			Facet_Idtype original_nei1=neighbor_surfacet(fid,1);
			Facet_Idtype original_nei2=neighbor_surfacet(fid,2);
			unique_pCreated_DeletedSurfacets_stack->push(MarkedSurfacet(c,li,false,
				original_nei0,original_nei1,original_nei2));
			delete_surface_facet(Facet(c,li));		
			//-------------标记被删除表面上的点--------------------//
			Vertex_triple vf=make_vertex_triple(Facet(c,li));
			//cout<<"!!delete facet: ("<<vf.first<<","<<vf.second<<","<<vf.third<<")"
			//	<<"/("<<c<<","<<li<<")"<<endl;
			VertexDeletedFacet[vf.first]++;
			VertexDeletedFacet[vf.second]++;
			VertexDeletedFacet[vf.third]++;
			//===================标记end===========================//
		}
		else if (is_surface(Facet(c_li, neighbor_index(c_li, c))))    //neighbor_index(c_li, c)求c是c_li的第几个临近Cell;找mirror_facet
		{
			int neiI=neighbor_index(c_li, c);
			Facet_Idtype fid=CellSurface.get_tuple_element(c_li,neiI);
			Facet_Idtype original_nei0=neighbor_surfacet(fid,0);
			Facet_Idtype original_nei1=neighbor_surfacet(fid,1);
			Facet_Idtype original_nei2=neighbor_surfacet(fid,2);
			unique_pCreated_DeletedSurfacets_stack->push(MarkedSurfacet(c_li,neiI,false,
				original_nei0,original_nei1,original_nei2));
			delete_surface_facet(Facet(c_li, neiI));
			//-------------标记被删除表面上的点--------------------//
			Vertex_triple vf=make_vertex_triple(Facet(c_li, neiI));
			//cout<<"!!delete facet: ("<<vf.first<<","<<vf.second<<","<<vf.third<<")"
			//	<<"/("<<c_li<<","<<neiI<<")"<<endl;
			VertexDeletedFacet[vf.first]++;
			VertexDeletedFacet[vf.second]++;
			VertexDeletedFacet[vf.third]++;
			//===================标记end===========================//
		}
	}
	else
	{
		if (label) //cnew表面内,c_li表面外
		{
			if (is_surface(Facet(c, li)))
			{
				//delete_surface_facet(Facet(c,li));
				//面（c,li）的编号改成面（cnew,li）对应
				change_surface_facet_link(Facet(c,li),cnew,li);
				unique_pChangedSurfacets_stack->push(ChangedSurfacet(c,cnew,li));
				//-------------标记被删除表面上的点--------------------//
				Vertex_triple vfc=make_vertex_triple(Facet(c,li));
				//cout<<"!!change surface facet link: ("<<vfc.first<<","<<vfc.second<<","<<vfc.third<<")"
				//	<<"/("<<c<<","<<li<<")"<<endl;
				VertexCreatedFacet.insert(vfc.first);
				VertexCreatedFacet.insert(vfc.second);
				VertexCreatedFacet.insert(vfc.third);
				//===================标记end===========================//
			}
			else
			{
				throw exception();
				//create_surface_facet(cnew,li);
				//NewCreateBoundFacet.push_back(Facet(cnew,li));
				////-----------------标记新生成表面上的点-------------------//
				//Vertex_triple vfc=make_vertex_triple(Facet(cnew,li));
				////cout<<"!!create facet without the new point: ("<<vfc.first<<","<<vfc.second<<","<<vfc.third<<")"
				////	<<"/("<<cnew<<","<<li<<")"<<endl;
				//VertexCreatedFacet.insert(vfc.first);
				//VertexCreatedFacet.insert(vfc.second);
				//VertexCreatedFacet.insert(vfc.third);
				////=======================标记end=========================//
			}
		}
		else  //cnew表面外,c_li表面内
		{
			if (!is_surface(Facet(c_li, neighbor_index(c_li, c))))
			{
				throw exception();
				//int neiI=neighbor_index(c_li,c);
				//create_surface_facet(c_li, neiI);
				//NewCreateBoundFacet.push_back(Facet(c_li,neiI));
				////-----------------标记新生成表面上的点-------------------//
				//Vertex_triple vfc=make_vertex_triple(Facet(c_li,neiI));
				////cout<<"!!create facet without the new point: ("<<vfc.first<<","<<vfc.second<<","<<vfc.third<<")"
				////	<<"/("<<c_li<<","<<neiI<<")"<<endl;
				//VertexCreatedFacet.insert(vfc.first);
				//VertexCreatedFacet.insert(vfc.second);
				//VertexCreatedFacet.insert(vfc.third);
				////=======================标记end=========================//
			}
				
		}
	}
	//================================判断end===============================
	////neighbor_index(c_li, c)求c是c_li的第几个临近Cell
	//设定相邻关系，cnew的第li个领域是c_li；c_li的第neighbor_index(c_li, c)个邻域是cnew
	unique_pChangedCellAjcacency_stack->push(ChangedCellAjcacency(c,c_li,neighbor_index(c_li, c)));
	set_adjacency(cnew, li, c_li, neighbor_index(c_li, c));

	// Look for the other neighbors of cnew.
	for (int ii = 0; ii<4; ++ii) {
		//if (ii == prev_ind2 || neighbor(cnew, ii) != -1)
		if (ii==prev_ind2||ii==li)
			continue;
		//cnew->vertex(ii)->set_cell(cnew);
		//将vertex(cnew, ii)的Cell设定为cnew
		unique_pChangedCellIncident_stack->push(ChangedCellIncident(vertex(cnew,ii),cell(vertex(cnew,ii))));
		set_cell(vertex(cnew, ii), cnew);
		// Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
		Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
		Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
		Cell_Idtype cur = c;
		int zz = ii;
		Cell_Idtype n = neighbor(cur, zz);
		// turn around the oriented edge vj1 vj2
		//沿着有向边vj1->vj2依次搜索与vj1->vj2关联的四面体,搜索方向与vj1->vj2方向（右手拇指方向）符合右手螺旋定则
		while (is_in_conflict(n)) {
			cur = n;
			zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
			n = neighbor(cur, zz);
		}
		// Now n is outside region, cur is inside.
		//但cur在影响域，且n=neighbor(cur,zz)不属于影响域，zz在迭代过程中用于替换形参li
		clear(n); // Reset the flag for boundary cells.

		int jj1 = vertex_index(n, vj1);
		int jj2 = vertex_index(n, vj2);
		Vertex_Idtype vvv = vertex(n, Triangulation_utils_3::next_around_edge(jj1, jj2));
		Cell_Idtype nnn = neighbor(n, Triangulation_utils_3::next_around_edge(jj2, jj1));
		//zzz表示当nnn（属于影响域）中已处理过的相对顶点，在迭代过程中用于替换形参pre_ind2
		int zzz = vertex_index(nnn, vvv);
		//-----------------------------判断label的符号是否改变--------------------------
		bool label2=label;
		bool SideChange2=false;
		if (boundEdge.find(make_pair(vj1, vj2)) != boundEdge.end())
		{
			label2 = !label;
			SideChange2 = true;
			//--------------提取交界面中的表面----------------//
			if (label == true)
			{
				create_surface_facet(cnew, ii);
				unique_pCreated_DeletedSurfacets_stack->push(MarkedSurfacet(cnew,ii,true,-1,-1,-1));
				//---------for test-------//
				Vertex_triple vtriFTT=make_vertex_triple(Facet(cnew,ii));
				//cout<<"!!create facet when finding the boundEdge: ("<<vtriFTT.first<<","<<vtriFTT.second<<","<<vtriFTT.third<<")"
				//	<<"/("<<cnew<<","<<ii<<")"<<endl;
				//=========test end=======//
		
			}
		}


		//=====================================判断end==================================
		//当nnn==cur时，表示nnn仍然需要进一步迭代产生新四面体；
		if (nnn == cur) {
			// Neighbor relation is reciprocal, ie
			// the cell we are looking for is not yet created.
			//nnn在影响域，neighbor(nnn,zz)不属于影响域,zz在迭代过程中用于替换形参li；zzz表示nnn（属于影响域）中已处理过的相对顶点，在迭代过程中用于替换形参pre_ind2
			nnn = surface_recursive_create_star_3(v, nnn, zz, depth + 1,boundEdge,
				VertexDeletedFacet,VertexCreatedFacet,NewCreateBoundFacet,label2,SideChange2,zzz);
		}
		//当nnn!=cur时，或者nnn已更新为新生成的cell，只需要更新与cnew的邻接关系
		set_adjacency(nnn, zzz, cnew, ii);
	}

	//---------for test-------//
	/*cout<<"recursive depth: "<<depth<<" ending!"<<endl;*/
	//=========test end=======//

	return cnew;
}

template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
non_recursive_create_star_3(Vertex_Idtype v, Cell_Idtype c, int li, int prev_ind2)
{
	

	Cell_Idtype cnew = create_cell(vertex(c,0),
		vertex(c,1),
		vertex(c,2),
		vertex(c,3));
	set_vertex(cnew,li, v);
	Cell_Idtype c_li = neighbor(c,li);
	set_adjacency(cnew, li, c_li, neighbor_index(c_li,c));

	//---------for test-------//
	//cout<<"[non recursive]create cell "<<cnew<<"("<<vertex(cnew, 0)<<","<<vertex(cnew, 1)
	//	<<","<<vertex(cnew, 2)<<","<<vertex(cnew, 3)<<")"<<endl;
	//=========test end=======//


	std::stack<iAdjacency_info> adjacency_info_stack;

	int ii = 0;
	do
	{
		// Look for the other neighbors of cnew.
		if (!(ii == prev_ind2 || neighbor(cnew,ii) != -1)){  //此处-1由Cell_Idtype()替换
			set_cell(vertex(cnew,ii),cnew);

			// Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
			Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
			Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
			Cell_Idtype cur = c;
			int zz = ii;
			Cell_Idtype n = neighbor(cur,zz);
			// turn around the oriented edge vj1 vj2
			//沿着有向边vj1->vj2依次搜索与vj1->vj2关联的四面体,搜索方向与vj1->vj2方向（右手拇指方向）符合右手螺旋定则
			while (is_in_conflict(n)) {
				cur = n;
				zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
				n = neighbor(cur,zz);
			}
			// Now n is outside region, cur is inside.
			clear(n); // Reset the flag for boundary cells.

			int jj1 = vertex_index(n,vj1);
			int jj2 = vertex_index(n,vj2);
			Vertex_Idtype vvv = vertex(n, Triangulation_utils_3::next_around_edge(jj1, jj2));
			Cell_Idtype nnn = neighbor(n, Triangulation_utils_3::next_around_edge(jj2, jj1));
			int zzz = vertex_index(nnn,vvv);
			//当nnn==cur时，表示nnn仍然需要进一步产生新四面体；
			if (nnn == cur) {
				// Neighbor relation is reciprocal, ie
				// the cell we are looking for is not yet created.
				//re-run the loop
				adjacency_info_stack.push(iAdjacency_info(zzz, cnew, ii, c, li, prev_ind2));
				c = nnn;  li = zz; prev_ind2 = zzz; ii = 0;
				//copy-pasted from beginning to avoid if ii==0
				cnew = create_cell(vertex(c,0), vertex(c,1), vertex(c,2), vertex(c,3));
				set_vertex(cnew,li, v);
				c_li = neighbor(c,li);
				set_adjacency(cnew, li, c_li, neighbor_index(c_li,c));
				continue;
			}
			set_adjacency(nnn, zzz, cnew, ii);//当nnn!=cur时，或者nnn已更新为新生成的cell，只需要更新与cnew的邻接关系
		}
		while (++ii == 4){
			if (adjacency_info_stack.empty()) return cnew;
			Cell_Idtype nnn = cnew;
			int zzz;
			adjacency_info_stack.top().update_variables(zzz, cnew, ii, c, li, prev_ind2);
			adjacency_info_stack.pop();
			set_adjacency(nnn, zzz, cnew, ii);
		}
	} while (true);
}

template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
create_star_2(Vertex_Idtype v, Cell_Idtype c, int li)
{

	Cell_Idtype cnew;

	// i1 i2 such that v,i1,i2 positive
	int i1 = Triangulation_cw_ccw_2::ccw(li);
	// traversal of the boundary of region in Triangulation_cw_ccw_2::ccw order to create all
	// the new facets
	Cell_Idtype bound = c;
	Vertex_Idtype v1 = vertex(c,i1);
	//int ind = c->neighbor(li)->index(c); // to be able to find the
	int ind = neighbor_index(neighbor(c,li),c);
	// first cell that will be created
	Cell_Idtype cur;
	//Cell_Idtype pnew = Cell_Idtype();
	Cell_Idtype pnew = -1;
	do {
		cur = bound;
		// turn around v2 until we reach the boundary of region
		while (is_in_conflict(neighbor(cur,Triangulation_cw_ccw_2::cw(i1)))) {
			// neighbor in conflict
			cur = neighbor(cur,Triangulation_cw_ccw_2::cw(i1));
			i1 = vertex_index(cur,v1);
		}
		//cur->neighbor(Triangulation_cw_ccw_2::cw(i1))->tds_data().clear();
		clear(neighbor(cur,Triangulation_cw_ccw_2::cw(i1)));
		// here cur has an edge on the boundary of region
		cnew = create_face(v, v1, vertex(cur,Triangulation_cw_ccw_2::ccw(i1)));
		set_adjacency(cnew, 0, neighbor(cur,Triangulation_cw_ccw_2::cw(i1)),
			neighbor_index(neighbor(cur,Triangulation_cw_ccw_2::cw(i1)),cur));
		set_neighbor(cnew,1, -1);
		set_neighbor(cnew,2, pnew);
		// pnew is null at the first iteration
		set_cell(v1,cnew);
		//pnew->set_neighbor( Triangulation_cw_ccw_2::cw(pnew->index(v1)), cnew );
		if (pnew != -1) { set_neighbor(pnew,1, cnew); }

		bound = cur;
		i1 = Triangulation_cw_ccw_2::ccw(i1);
		v1 = vertex(bound,i1);
		pnew = cnew;
		//} while ( ( bound != c ) || ( li != Triangulation_cw_ccw_2::cw(i1) ) );
	} while (v1 != vertex(c,Triangulation_cw_ccw_2::ccw(li)));
	// missing neighbors between the first and the last created cells
	cur = neighbor(neighbor(c,li),ind); // first created cell
	set_adjacency(cnew, 1, cur, 2);
	return cnew;
}



//求面（cell，Id）中按逆时针实践顶点编号（i,j,k）,并且该面三点与该cell中剩余一点符合右手螺旋定则
template<typename T, typename T_INFO>
Vertex_triple DataStructure<T, T_INFO>::
make_vertex_triple(Facet& f)
{
	Cell_Idtype ch = f.first;
	int i = f.second;

	//vertex_triple_index(i, j),i为四面体相对顶点编号，j(0/1/2)是i点相对面的编号的点，返回该点在该四面体的相对编号
	//vertex(Idc,I)取编号Idc的cell的第I个点；
	return Vertex_triple(vertex(ch, Triangulation_utils_3::vertex_triple_index(i, 0)),
		vertex(ch, Triangulation_utils_3::vertex_triple_index(i, 1)),
		vertex(ch, Triangulation_utils_3::vertex_triple_index(i, 2)));
}

template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
insert_first_finite_cell()
{
	const double *p0=point(1);
	const double *p1=point(2);
	const double *p2=point(3);
	const double *p3=point(4);

	if(GeometricTraits<T>::orientation(p0, p1, p2, p3) == NEGATIVE)
	{
		double pTemp[3];
		pTemp[0]=p0[0];
		pTemp[1]=p0[1];
		pTemp[2]=p0[2];
		set_point(1,p1);
		set_point(2,pTemp);
	}
	set_dimension(3);
	int v0=1,v1=2,v2=3,v3=4;
	//create_cell四个形参是顶点编号，依次对应四面体相对位置v0,v1,v2,v3
	Cell_Idtype c0123 = create_cell(v0, v1, v2, v3);         //符合右手定则
	Cell_Idtype ci012 = create_cell(Infinite, v0, v1, v2);   //符合右手定则
	Cell_Idtype ci103 = create_cell(Infinite, v1, v0, v3);   //符合右手定则
	Cell_Idtype ci023 = create_cell(Infinite, v0, v2, v3);   //符合右手定则
	Cell_Idtype ci132 = create_cell(Infinite, v1, v3, v2);   //符合右手定则


	//将点Infinite的cell设定为ci012
	set_cell(Infinite,ci012);
	set_cell(v0,c0123);   //将点v0的cell设定为c0123
	set_cell(v1,c0123);   //将点v1的cell设定为c0123
	set_cell(v2,c0123);   //将点v2的cell设定为c0123
	set_cell(v3,c0123);   //将点v3的cell设定为c0123


	//设定相邻关系，c0123的第0个领域是ci132；ci132的第0个邻域是c0123
	set_adjacency(c0123, 0, ci132, 0);  
	set_adjacency(c0123, 1, ci023, 0);    //设定相邻关系，c0123的第1个领域是ci023；ci023的第0个邻域是c0123
	set_adjacency(c0123, 2, ci103, 0);    //设定相邻关系，c0123的第2个领域是ci103；ci103的第0个邻域是c0123
	set_adjacency(c0123, 3, ci012, 0);    //设定相邻关系，c0123的第3个领域是ci012；ci012的第0个邻域是c0123

	set_adjacency(ci012, 3, ci103, 3);
	set_adjacency(ci012, 2, ci023, 3);
	set_adjacency(ci012, 1, ci132, 2);
	set_adjacency(ci103, 1, ci023, 2);
	set_adjacency(ci023, 1, ci132, 1);
	set_adjacency(ci103, 2, ci132, 3);
	//HHHHHHHHHHHHHHHHHHH-------------------surface extractor-------------------HHHHHHHHHHHHHHHHHHHHHH
	if (this->Surface)
	{
		Facet_Idtype f012 = create_surface_facet(c0123, 3);
		Facet_Idtype f013 = create_surface_facet(c0123, 2);
		Facet_Idtype f023 = create_surface_facet(c0123, 1);
		Facet_Idtype f123 = create_surface_facet(c0123, 0);
		//利用make_vertex_triple中取点的方法
		//FacetNeighbor暂时不考虑
		Indextype i0,i1;
		

		//vertex_index_facet(Facet F,Vertex_Idtype IdV)求点IdV是facet F的第几个点
		i0=vertex_index_facet(Facet(c0123,3),3);
		i1=vertex_index_facet(Facet(c0123,2),4);
		//f012的i0个相邻面是f013,f013的第i1个相邻面是f012
		set_facet_adjacency(f012,i0,f013,i1);

		i0=vertex_index_facet(Facet(c0123,3),2);
		i1=vertex_index_facet(Facet(c0123,1),4);
		set_facet_adjacency(f012,i0,f023,i1);

		i0=vertex_index_facet(Facet(c0123,3),1);
		i1=vertex_index_facet(Facet(c0123,0),4);
		set_facet_adjacency(f012,i0,f123,i1);

		i0=vertex_index_facet(Facet(c0123,2),2);
		i1=vertex_index_facet(Facet(c0123,1),3);
		set_facet_adjacency(f013,i0,f023,i1);

		i0=vertex_index_facet(Facet(c0123,2),1);
		i1=vertex_index_facet(Facet(c0123,0),3);
		set_facet_adjacency(f013,i0,f123,i1);

		i0=vertex_index_facet(Facet(c0123,1),1);
		i1=vertex_index_facet(Facet(c0123,0),2);
		set_facet_adjacency(f023,i0,f123,i1);

		//标注
		label_cell_side(c0123, true);
		label_cell_side(ci012, false);
		label_cell_side(ci023, false);
		label_cell_side(ci103, false);
		label_cell_side(ci132, false);
	}

}

template<typename T, typename T_INFO>
Vertex_Idtype DataStructure<T, T_INFO>::
insert_first_finite_cell(Vertex_Idtype &v0, Vertex_Idtype &v1, Vertex_Idtype &v2, Vertex_Idtype &v3,
Vertex_Idtype Infinite)
{

	if (Infinite == -1)
		Infinite = create_vertex();
	//--------------------在data structure中加入infinite vertex-------------------//
	Infinite=Infinite;
	//============================加入end=========================================//

	set_dimension(3);

	v0 = create_vertex();
	v1 = create_vertex();
	v2 = create_vertex();
	v3 = create_vertex();

	Cell_Idtype c0123 = create_cell(v0, v1, v2, v3);
	Cell_Idtype ci012 = create_cell(Infinite, v0, v1, v2);
	Cell_Idtype ci103 = create_cell(Infinite, v1, v0, v3);
	Cell_Idtype ci023 = create_cell(Infinite, v0, v2, v3);
	Cell_Idtype ci132 = create_cell(Infinite, v1, v3, v2);

	set_cell(Infinite,ci012);
	set_cell(v0,c0123);
	set_cell(v1,c0123);
	set_cell(v2,c0123);
	set_cell(v3,c0123);

	set_adjacency(c0123, 0, ci132, 0);
	set_adjacency(c0123, 1, ci023, 0);
	set_adjacency(c0123, 2, ci103, 0);
	set_adjacency(c0123, 3, ci012, 0);

	set_adjacency(ci012, 3, ci103, 3);
	set_adjacency(ci012, 2, ci023, 3);
	set_adjacency(ci012, 1, ci132, 2);
	set_adjacency(ci103, 1, ci023, 2);
	set_adjacency(ci023, 1, ci132, 1);
	set_adjacency(ci103, 2, ci132, 3);
	//HHHHHHHHHHHHHHHHHHH-------------------surface extractor-------------------HHHHHHHHHHHHHHHHHHHHHH
	if (this->Surface)
	{
		Facet_Idtype f012 = create_surface_facet(c0123, 3);
		Facet_Idtype f013 = create_surface_facet(c0123, 2);
		Facet_Idtype f023 = create_surface_facet(c0123, 1);
		Facet_Idtype f123 = create_surface_facet(c0123, 0);
		//??????????????三角面的数据结构怎么存储（facet并没有存点，也不知道点的index）
		//利用make_vertex_triple中取点的方法
		//FacetNeighbor暂时不考虑
		Indextype i0,i1;
		
		i0=vertex_index_facet(Facet(c0123,3),3);
		i1=vertex_index_facet(Facet(c0123,2),4);
		set_facet_adjacency(f012,i0,f013,i1);

		i0=vertex_index_facet(Facet(c0123,3),2);
		i1=vertex_index_facet(Facet(c0123,1),4);
		set_facet_adjacency(f012,i0,f023,i1);

		i0=vertex_index_facet(Facet(c0123,3),1);
		i1=vertex_index_facet(Facet(c0123,0),4);
		set_facet_adjacency(f012,i0,f123,i1);

		i0=vertex_index_facet(Facet(c0123,2),2);
		i1=vertex_index_facet(Facet(c0123,1),3);
		set_facet_adjacency(f013,i0,f023,i1);

		i0=vertex_index_facet(Facet(c0123,2),1);
		i1=vertex_index_facet(Facet(c0123,0),3);
		set_facet_adjacency(f013,i0,f123,i1);

		i0=vertex_index_facet(Facet(c0123,1),1);
		i1=vertex_index_facet(Facet(c0123,0),2);
		set_facet_adjacency(f023,i0,f123,i1);

		//标注
		label_cell_side(c0123, true);
		label_cell_side(ci012, false);
		label_cell_side(ci023, false);
		label_cell_side(ci103, false);
		label_cell_side(ci132, false);
	}
	

	return Infinite;
}



template<typename T, typename T_INFO>
template<typename OutputIteratorCell>
void DataStructure<T, T_INFO>::
incident_cells(Vertex_Idtype v, OutputIteratorCell outputCell) 
{


	if (dimension() < 2)
		return;

	//Visitor visit(v, outputCell, this, f);

	std::vector<Cell_Idtype> tmp_cells;
	tmp_cells.reserve(64);
	std::vector<Facet> tmp_facets;
	tmp_facets.reserve(64);
	std::vector<Vertex_Idtype> temp_vertices;
	if (dimension() == 3)
		incident_cells_3(v, cell(v), std::make_pair(std::back_inserter(tmp_cells), std::back_inserter(tmp_facets)));
	/*else
		incident_cells_2(v, cell(v), std::back_inserter(tmp_cells));*/

	typename std::vector<Cell_Idtype>::iterator cit;
	for (cit = tmp_cells.begin();
		cit != tmp_cells.end();
		++cit)
	{
		clear(*cit);

		*outputCell++ = *cit;
		
	}
	
}

template<typename T, typename T_INFO>
template<typename OutputIteratorCell,typename OutputIteratorVertex,typename OutputIteratorFacet>
void DataStructure<T, T_INFO>::
incident_cells_mark_side(Vertex_Idtype v, OutputIteratorCell outputCell,OutputIteratorVertex outputVertex,
	OutputIteratorFacet outputFacet,bool& IsIso,bool& PInside) 
{
	std::vector<Cell_Idtype> tmpCells;
	tmpCells.reserve(32);
	std::stack<Cell_Idtype> cell_stack;
	Cell_Idtype d=cell(v);
	cell_stack.push(d);
	mark_in_conflict(d);
	*outputCell++ = d;
	tmpCells.push_back(d);
	IsIso=true;
	PInside=is_label_inside(d);
	do {
		Cell_Idtype c = cell_stack.top();
		cell_stack.pop();

		for (int i = 0; i<4; ++i) {
			if (vertex(c,i) == v)
				continue;
			Cell_Idtype next = neighbor(c,i);
			if (is_surface(Facet(c,i)))
			{
				IsIso=false;
				*outputFacet++=Facet(c,i);
				for(int j=0;j<4;j++)
				{
					if(j!=i)
					{
						Vertex_Idtype vEdge=vertex(c,j);
						if(vEdge!=v&&visited_for_vertex_extractor(vEdge)==0)
						{
							*outputVertex++=vEdge;
							mark_vertex_visited(vEdge);
						}
					}
				}

			}
			if (!is_clear(next))
				continue;
			cell_stack.push(next);
			//next->tds_data().mark_in_conflict();
			mark_in_conflict(next);
			*outputCell++ = next;
			tmpCells.push_back(next);
		}
	} while (!cell_stack.empty());

	typename std::vector<Cell_Idtype>::iterator cit;
	for (cit = tmpCells.begin();
		cit != tmpCells.end();
		++cit)
	{
		clear(*cit);
	
		//=================vertex extractor===============
		for (int j = 0; j <= dimension(); ++j) {
			Vertex_Idtype w = vertex(*cit,j);
			
			if (w != v){

				if (visited_for_vertex_extractor(w)==0){
					//w->visited_for_vertex_extractor = true;
					if(is_label_inside(*cit))
					{
						mark_vertex_visited(w);
						mark_vertex_visited(w);
						mark_vertex_visited(w);
					}
					else
					{
						mark_vertex_visited(w);
						mark_vertex_visited(w);
					}
					
					*outputVertex++ = w;
					//treat(c, v, j);
				}
			}
		}
		//----------------------vertex extractor end-------------------
	}
	


}

template<typename T, typename T_INFO>
template<typename OutputIteratorFacet, typename OutputIteratorVertex>
void DataStructure<T, T_INFO>::
incident_cells(Vertex_Idtype v,  
				OutputIteratorFacet outputFacets, OutputIteratorVertex outputVertex) 
{

	if (dimension() < 2)
		return ;
	std::vector<Cell_Idtype> tmp_cells;
	tmp_cells.reserve(64);

	std::vector<Vertex_Idtype> tmp_vertices;

	if (dimension() == 3)
		incident_cells_3(v, cell(v), std::make_pair(std::back_inserter(tmp_cells), outputFacets));


	typename std::vector<Cell_Idtype>::iterator cit;
	for (cit = tmp_cells.begin();
		cit != tmp_cells.end();
		++cit)
	{
		clear(*cit);

		for (int j = 0; j <= dimension(); ++j) {
			Vertex_Idtype w = vertex(*cit,j);
			
			if (w != v){

				if (visited_for_vertex_extractor(w)==0){
					mark_vertex_visited(w);
					tmp_vertices.push_back(w);
					*outputVertex++ = w;
				}
			}
		}
	}
	for (std::size_t i = 0; i < tmp_vertices.size(); ++i){
		clear_vertex_visited(tmp_vertices[i]);
	}

}



template<typename T, typename T_INFO>
template<typename IncidentCellIterator, typename IncidentFacetIterator>
void DataStructure<T, T_INFO>::
incident_cells_3(Vertex_Idtype v, Cell_Idtype d,
						std::pair<IncidentCellIterator,
						IncidentFacetIterator> it) 
{

	std::stack<Cell_Idtype> cell_stack;
	cell_stack.push(d);
	mark_in_conflict(d);
	*it.first++ = d;

	do {
		Cell_Idtype c = cell_stack.top();
		cell_stack.pop();

		for (int i = 0; i<4; ++i) {
			if (vertex(c,i) == v)
				continue;
			Cell_Idtype next = neighbor(c,i);
			if (c < next)
				*it.second++ = Facet(c, i); // Incident facet.
			if (!is_clear(next))
				continue;
			cell_stack.push(next);
			mark_in_conflict(next);
			*it.first++ = next;
		}
	} while (!cell_stack.empty());

}



template<typename T, typename T_INFO>
template<typename OutputIteratorVertex, typename OutputIteratorFacet>
void DataStructure<T,T_INFO>::
adjacent_vertices_facets(Vertex_Idtype v,
		OutputIteratorFacet outputFacets, OutputIteratorVertex vertices) 
{
	if (dimension() == -1)
		return ;

	if (dimension() == 0) {
		//Vertex_Idtype v1 = v->cell()->neighbor(0)->vertex(0);
		Vertex_Idtype v1 = vertex(neighbor(cell(v),0),0);
		 *vertices++ = v1;
		return ;
	}

	if (dimension() == 1) {
		Cell_Idtype n0 = cell(v);
		const int index_v_in_n0 = vertex_index(n0,v);
		Cell_Idtype n1 = neighbor(n0,(1 - index_v_in_n0));
		const int index_v_in_n1 = vertex_index(n1,v);
		Vertex_Idtype v1 = vertex(n0,(1 - index_v_in_n0));
		Vertex_Idtype v2 = vertex(n1,(1 - index_v_in_n1));
		 *vertices++ = v1;
		 *vertices++ = v2;
		return ;
	}
	incident_cells(v, outputFacets, vertices);

}



//求p到面f的三个顶点中，最长的边
template<typename T, typename T_INFO>
double DataStructure<T, T_INFO>::
longest_edge_to_facet(const T* p,Facet f)
{
	double longestLength = 0;
	Vertex_triple vertexFacet = make_vertex_triple(f);
	Vertex_Idtype vertexTemp = -1;
	for (int i = 0; i < 3; i++)
	{
		//vertexTemp = vertexFacet.get<i>();
		if (i == 0)
			vertexTemp = vertexFacet.first;
		else if (i == 1)
			vertexTemp = vertexFacet.second;
		else
			vertexTemp = vertexFacet.third;
			
		
		double distantTemp = distant_to_vertex(p,vertexTemp);
		if (distantTemp > longestLength)
			longestLength = distantTemp;
		
	}
	return longestLength;
}


template<typename T, typename T_INFO>
double DataStructure<T, T_INFO>::
shortest_edge_to_facet(const T* p,Facet f)
{
	double shortestLength = 0;
	Vertex_triple vertexFacet = make_vertex_triple(f);
	Vertex_Idtype vertexTemp = -1;
	for (int i = 0; i < 3; i++)
	{
		if (i == 0)
			vertexTemp = vertexFacet.first;
		else if (i == 1)
			vertexTemp = vertexFacet.second;
		else
			vertexTemp = vertexFacet.third;
			
		
		double distantTemp = distant_to_vertex(p,vertexTemp);
		if (distantTemp > shortestLength)
			shortestLength = distantTemp;
		
	}
	return shortestLength;
}

template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
neighbor_surface_facet(EdgeInFacet SurfaceEdge,Facet& NeighborFacet, bool& IsToDelete)
{
	Facet surfaceFacet = SurfaceEdge.first;
	int ii = SurfaceEdge.second;
	Cell_Idtype c= surfaceFacet.first;
	int li = surfaceFacet.second;
	Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
	Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
	Cell_Idtype cur = c;
	int zz = ii;
	Cell_Idtype n = neighbor(cur,zz);
	// turn around the oriented edge vj1 vj2
	while (is_label_inside(n)) {
		
		cur = n;
		zz = Triangulation_utils_3::next_around_edge(index(n, vj1), index(n, vj2));
		n = neighbor(cur,zz);
	}
	if (is_facet_to_delete(Facet(cur, zz)))
	{
		IsToDelete = true;
	}
	else
		IsToDelete = false;
	if (is_surface(Facet(cur, zz)))
	{
		NeighborFacet = Facet(cur, zz);
		return true;
	}
	else
		return false;
	
}


//在孤立点情况下,使用Gabriel的规则决定替换面的延展
//V孤立点;InitFacet为种子三角面片;PInside为孤立点V关于表面的位置
//VertexOnBoundary为边界线上的顶点，用于孤立点的检查;
//VertexOnConflictFacet为被替换面的三角面片中的顶点,用于孤立点的检查;Conflicts是被替换面;
//BoundEdges为返回值，被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换
//NewCreateSurfacets待新生成表面的集合
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
find_isolate_conflict_surfacet(Vertex_Idtype V, Facet InitFacet,bool PInside,std::set<Vertex_Idtype>& VertexOnBoundary,
							set<Vertex_Idtype>& VertexOnConflictFacet,vector<Facet>& Conflicts,list<Edge>& BoundEdges,vector<Facet>& NewCreateSurfacets,bool& quit)
{
	const T* P=point(V);
	list<Facet> illSurfacets;
	vector<Facet> tmpSurfacets;
	list<Facet> neighborSurfacets;
	neighborSurfacets.push_back(InitFacet);
	tmpSurfacets.push_back(InitFacet);
	mark_facet_visited(InitFacet);
	mark_facet_visited(InitFacet);

	vector<Facet_Idtype> idtmpSurfacets;
	Facet_Idtype fid=CellSurface.get_tuple_element(InitFacet.first,InitFacet.second);
	idtmpSurfacets.push_back(fid);
	//-----------------求InitFacet邻域的表面三角面片（即与InitFacet有公共点的）------------------//
	while(!neighborSurfacets.empty())
	{
		Facet tmpSurfacet=*neighborSurfacets.begin();
		//--------------for test------------//
		Vertex_triple vtriT0=make_vertex_triple(tmpSurfacet);
		//===============test end==========//
		neighborSurfacets.pop_front();
		Cell_Idtype c = tmpSurfacet.first;//HHHHHHHHHHHHHHHHHHHHHHHHHH-------c肯定被label为inside
		int li = tmpSurfacet.second;
		Facet nei;
		for(int i=0;i<3;i++)
		{
			nei=neighbor_surfacet(tmpSurfacet,i);
			//--------------对于孤立点的处理------------//
			Cell_Idtype cur=nei.first;
			Indextype zz=nei.second;
			Cell_Idtype cOpp=neighbor(cur,zz);
			Indextype indOpp=neighbor_index(cOpp,cur);
			//==============对于孤立点的处理===========//
			//--------------for test------------//
			Vertex_triple vtriT1=make_vertex_triple(nei);
			//===============test end==========//
			if(visited_for_facet_extractor(nei)==0)
			{
				//孤立点下不会出现
				if(is_facet_in_conflict(nei))
				{
					mark_facet_visited(nei);
					mark_facet_visited(nei);
					neighborSurfacets.push_back(nei);
					tmpSurfacets.push_back(nei);
				}
				//待选被替换面应该应该有一个cell,该面相对孤立点
				if((PInside&&vertex(cur,zz)==V)||(!PInside&&vertex(cOpp,indOpp)==V))
				//if((!IsIsolate&&((is_on_boundary(nei.first)&&is_in_conflict(neighbor(nei.first,nei.second)))||
					//(is_on_boundary(neighbor(nei.first,nei.second))&&is_in_conflict(nei.first))))||(IsIsolate&&((PInside&&vertex(cur,zz)==V)||(!PInside&&vertex(cOpp,indOpp)==V))))
				{
					T cosDihedral[3];
					//if(PInside)
					cos_dihedral_point_to_facet(P,nei,cosDihedral);
				
					
		
					//加入方向判断
					Vertex_triple vf=make_vertex_triple(nei);
					Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), P);
					if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
					{
						//判断加入nei是否为凸延展
						bool isConvexExtension=true;
						for(int i=0;i<3;i++)
						{
							Facet neiEdge=neighbor_surfacet(nei,i);
							//--------------------------for test-------------------------//
							Vertex_triple vtriT3=make_vertex_triple(neiEdge);
							//==========================test end=========================//
							//if(cosDihedral[i]<0&&neighbor_surfacet())
							if(cosDihedral[i]<0&&!is_in_conflict(neiEdge.first)&&!is_in_conflict(neighbor(neiEdge.first,neiEdge.second)))
							{
								isConvexExtension=false;
								break;
							}
						}
						//改变选择替换面的规则，以pInside为基础，以是否符合Gabriel规则为准则
						if (PInside)
						{
							bool toBeSurface=false;

							for (int j = 0; j < 3; j++)
							{
								Facet NNei=neighbor_surfacet(nei,i);
								if (visited_for_facet_extractor(NNei)==2)
								{
									Vertex_Idtype vertexOpp=vertex_facet(nei,i);
									Vertex_pair vertexOppEdge=make_vertex_pair(EdgeInFacet(nei,j));
									bool isInShpere=GeometricTraits<T>::side_of_bounded_sphereC3(point(vertexOppEdge.first),P,
										point(vertexOppEdge.second),point(vertexOpp))==ON_BOUNDED_SIDE;
									if (isInShpere)
									{
										toBeSurface=true;
									}
								}
							}
							if (!toBeSurface)
							{
								Vertex_triple ftrip=make_vertex_triple(nei);
								toBeSurface=GeometricTraits<T>::side_of_bounded_sphereC3(point(ftrip.first),
									point(ftrip.second),point(ftrip.third),P)==ON_BOUNDED_SIDE;
							}
							if(toBeSurface)//&&isConvexExtension)
							{
								mark_facet_visited(nei);
								mark_facet_visited(nei);
								neighborSurfacets.push_back(nei);
								tmpSurfacets.push_back(nei);

								Facet_Idtype fid=CellSurface.get_tuple_element(nei.first,nei.second);
								idtmpSurfacets.push_back(fid);
							}
							else
							{
								mark_facet_visited(nei);
								illSurfacets.push_back(nei);
							}
						}
						else
						{
							bool toBeSurface=true;
							//Facet facetOpp=mirror_facet(nei);
							for (int j = 0; j < 3; j++)
							{
								Facet NNei=neighbor_surfacet(nei,i);
								if (visited_for_facet_extractor(NNei)!=2)
								{
									Vertex_Idtype vertexOpp=vertex_facet(nei,i);
									Vertex_pair vertexOppEdge=make_vertex_pair(EdgeInFacet(nei,j));
									bool isInShpere=GeometricTraits<T>::side_of_bounded_sphereC3(point(vertexOppEdge.first),
										point(vertexOppEdge.second),P,point(vertexOpp))==ON_BOUNDED_SIDE;
									if (isInShpere)
									{
										toBeSurface=false;
									}
								}
							}
							if(toBeSurface)//&&isConvexExtension)
							{
								mark_facet_visited(nei);
								mark_facet_visited(nei);
								neighborSurfacets.push_back(nei);
								tmpSurfacets.push_back(nei);

								Facet_Idtype fid=CellSurface.get_tuple_element(nei.first,nei.second);
								idtmpSurfacets.push_back(fid);
							}
							else
							{
								mark_facet_visited(nei);
								illSurfacets.push_back(nei);
							}
						}
					}
				}
			}
		}
	}
	//对上一步中求得的交界面删除其中的可能造成折叠的surface facet
	while(!illSurfacets.empty())
	{
		Facet tmpIllSurfacet=*(illSurfacets.begin());
		clear_facet_visited(tmpIllSurfacet);
		illSurfacets.pop_front();
		for(int i=0;i<3;i++)
		{
			Facet nei=neighbor_surfacet(tmpIllSurfacet,i);
			//---------------for test------------//
			Vertex_triple vtriNei=make_vertex_triple(nei);
			//==============test end============//
			if(visited_for_facet_extractor(nei)==2&&!is_facet_in_conflict(nei)&&nei!=InitFacet)
			{
				Indextype iNei=neighbor_index_facet(nei,tmpIllSurfacet);
				if(cos_dihedral_point_to_Edge(P,nei,iNei)<0)
				{
					illSurfacets.push_back(nei);
					clear_facet_visited(nei);

					Facet_Idtype fid=CellSurface.get_tuple_element(nei.first,nei.second);
					idtmpSurfacets.erase(find(idtmpSurfacets.begin(),idtmpSurfacets.end(),fid));
				}
			}
		}
	}

 //   if(idtmpSurfacets.size()>0)
	//{

	//	for(auto itt=tmpSurfacets.begin();itt!=tmpSurfacets.end();itt++)
	//	{
	//		clear_facet_visited(*itt);
	//	}
	//	quit=true;
	//	return;
	//}

	//重新以种子三角面片为基础搜索符合条件的交界面三角面片及其边界边
	//重新以种子三角面片为基础搜索符合条件的交界面三角面片及其边界边
	//EdgeInFacet是(Facet,0/1/2)表示,其中Facet以(cell,0/1/2/3)表示
	//BoundEdges被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换
	//边界边输出是按着顺时针顺序，初始化也必须顺时针(从表面内看,表面外看则逆时针),所以按着相对顶点0->2->1加入初始化边界边
	mark_facet_visited(InitFacet);
	EdgeInFacet tmpBoundEdge=mirror_edge(EdgeInFacet(InitFacet,0));
	BoundEdges.push_back(turn_edgeinfacet_to_edge(tmpBoundEdge));
	tmpBoundEdge=mirror_edge(EdgeInFacet(InitFacet,2));
	BoundEdges.push_back(turn_edgeinfacet_to_edge(tmpBoundEdge));
	tmpBoundEdge=mirror_edge(EdgeInFacet(InitFacet,1));
	BoundEdges.push_back(turn_edgeinfacet_to_edge(tmpBoundEdge));
	list<Facet> tmpConflitFacets;
	tmpConflitFacets.push_back(InitFacet);
	/*cout<<"被替换面: "<<endl;*/
	while(!tmpConflitFacets.empty())
	{
		Facet tmp=*(tmpConflitFacets.begin());
		Conflicts.push_back(tmp);
		tmpConflitFacets.pop_front();
		Indextype li=tmp.second;
		Vertex_triple vtrif=make_vertex_triple(tmp);
		VertexOnConflictFacet.insert(vtrif.first);
		VertexOnConflictFacet.insert(vtrif.second);
		VertexOnConflictFacet.insert(vtrif.third);
		/*VertexOnConflictFacet[vtrif.first]++;
		VertexOnConflictFacet[vtrif.second]++;
		VertexOnConflictFacet[vtrif.third]++;*/
		//------------------for test--------------//
		Vertex_triple vtriTC0=make_vertex_triple(tmp);
		Cell_Idtype tmpc=tmp.first;
		Facet_Idtype fid=CellSurface.get_tuple_element(tmpc,tmp.second);
		//cout<<fid<<"("<<vtriTC0.first<<","<<vtriTC0.second<<","<<vtriTC0.third<<")"
		//	<<"/("<<tmp.first<<","<<tmp.second<<")["<<tmpc<<"("
		//	<<vertex(tmpc, 0)<<","<<vertex(tmpc, 1)<<","<<vertex(tmpc, 2)<<","<<vertex(tmpc, 3)<<")]"<<endl;
		//==================test end==============//
		
		for(int i=0;i<3;i++)
		{
			Facet nei=neighbor_surfacet(tmp,i);
			//---------------for test------------//
			Vertex_triple vtriT=make_vertex_triple(nei);
			//===============test end============//
			Indextype iNei=neighbor_index_facet(nei,tmp);
			if(visited_for_facet_extractor(nei)==2)
			{
				mark_facet_visited(nei);
				tmpConflitFacets.push_back(nei);
				//更新交界面的边界边队列
				vector<int> numNei;
				int numOfList=0;

				//--------------修改------------//
				auto itEdge=BoundEdges.begin();
				for(;itEdge!=BoundEdges.end();itEdge++)
				{
					numOfList++;
					Facet tmpBoundFacet0=get_facet_cell_link((*itEdge).first);
					//--------------for test----------//
					Vertex_triple vtriT0=make_vertex_triple(tmpBoundFacet0);
					//==============test end==========//
					if(tmpBoundFacet0==nei)
					{
						numNei.push_back((*itEdge).second);
						BoundEdges.erase(itEdge++);
						break;
					}
				}
				//当list链的第一个被删除时，应该检查最后一个是否也应该被删除
				if(numOfList==1)
				{
					auto endEdge=BoundEdges.back();
					Facet tmpBoundFacet00=get_facet_cell_link(endEdge.first);
					//--------------for test----------//
					Vertex_triple vtriT00=make_vertex_triple(tmpBoundFacet00);
					//==============test end=========//
					if(tmpBoundFacet00==nei)
					{
					
						numNei.push_back(endEdge.second);						
						BoundEdges.pop_back();
						endEdge=BoundEdges.back();
						Facet tmpBoundFacet01=get_facet_cell_link(endEdge.first);
						//--------------for test----------//
						Vertex_triple vtriT01=make_vertex_triple(tmpBoundFacet00);
						//==============test end=========//
						if(tmpBoundFacet01==nei)
						{
							numNei.push_back(endEdge.second);						
							BoundEdges.pop_back();
						}
					}
					
				}
				if(numNei.size()!=3)
				{
					//list 中间的插入与删除要处理好
					Facet tmpBoundFacet1=make_pair(-1,-1);
					if(itEdge!=BoundEdges.end())
					{
						tmpBoundFacet1=get_facet_cell_link((*itEdge).first);
						//----------------for test------------//
						Vertex_triple vtriT1=make_vertex_triple(tmpBoundFacet1);
					}
					//================test end===========//
					if(itEdge!=BoundEdges.end()&&tmpBoundFacet1==nei)
					{
						numNei.push_back((*itEdge).second);
						BoundEdges.erase(itEdge++);
						Facet tmpBoundFacet1=make_pair(-1,-1);
						if(itEdge!=BoundEdges.end())
						{
							tmpBoundFacet1=get_facet_cell_link((*itEdge).first);
							//----------------for test------------//
							Vertex_triple vtriT1=make_vertex_triple(tmpBoundFacet1);
						}
						if(itEdge!=BoundEdges.end()&&tmpBoundFacet1==nei)
						{
							numNei.push_back((*itEdge).second);
							BoundEdges.erase(itEdge++);
						}
					}
					
				}
				if(numNei.size()==1)
				{
					int i0=numNei[0];
					auto leftRight=make_edge_index(i0);
					EdgeInFacet oppEdge0=mirror_edge(EdgeInFacet(nei,leftRight.second));
					EdgeInFacet oppEdge1=mirror_edge(EdgeInFacet(nei,leftRight.first));
					//--------------for test-----------//
					Vertex_triple vtriT3=make_vertex_triple(oppEdge0.first);
					Vertex_triple vtriT4=make_vertex_triple(oppEdge1.first);
					//==============test end==========//
					BoundEdges.insert(itEdge,turn_edgeinfacet_to_edge(oppEdge0));
					BoundEdges.insert(itEdge,turn_edgeinfacet_to_edge(oppEdge1));
				}
				else if (numNei.size()==2)
				{
					for(int i2=0;i2<3;i2++)
					{
						if(i2!=numNei[0]&&i2!=numNei[1])
						{
							EdgeInFacet oppEdge=mirror_edge(EdgeInFacet(nei,i2));
							//--------------for test-----------//
							Vertex_triple vtriT2=make_vertex_triple(oppEdge.first);
							//==============test end==========//
							BoundEdges.insert(itEdge,turn_edgeinfacet_to_edge(oppEdge));
							break;
						}
					}
				}


			}
			else if(visited_for_facet_extractor(nei)==0)
			{
				Vertex_Idtype vOpp=vertex_facet(tmp,i);
				
				if(PInside)
				{
					Indextype ii=vertex_index(tmp.first,vOpp);
					//create_surface_facet(InitFacet.first,ii);
					Facet newFacet=mirror_facet(Facet(tmp.first,ii));
					NewCreateSurfacets.push_back(newFacet);
				}
				else
				{
					Cell_Idtype cOpp=neighbor(tmp.first,li);
					//Vertex_Idtype vOppEdge=vertex(tmp.first,ii);
					Indextype iiNewSur=vertex_index(cOpp,vOpp);
					//create_surface_facet(cOpp,iiNewSur);
					NewCreateSurfacets.push_back(Facet(cOpp,iiNewSur));
				}
				Vertex_pair vEdge=make_vertex_pair(EdgeInFacet(tmp,i));
				VertexOnBoundary.insert(vEdge.first);
				VertexOnBoundary.insert(vEdge.second);
			}
		}
	}
	//-------------------for test(打印 BoundEdges)------------------//
	//EdgeInFacet是(Facet,0/1/2)表示,其中Facet以(cell,0/1/2/3)表示
	//BoundEdges被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换
	//BoundEdge也是边界线集合，以端点(vertex_id,vertex_id)形式表示
	//cout<<"conflict boundary edge:"<<endl;
	//for(auto itBE=BoundEdges.begin();itBE!=BoundEdges.end();itBE++)
	//{
	//	Edge tmp=*itBE;
	//	Facet tmpF=get_facet_cell_link(tmp.first);
	//	//--------------for test----------//
	//	Cell_Idtype tmpc=tmpF.first;
	//	Vertex_triple vtrit = make_vertex_triple(tmpF);
	//	//==============test end=========//
	//	Vertex_pair vp=make_vertex_pair(EdgeInFacet(tmpF,tmp.second));
	//	cout<<"("<<vp.first<<","<<vp.second<<")"
	//		<<"/Edge("<<tmp.first<<","<<tmp.second<<") "<<tmp.first<<"["<<vtrit.first<<","<<vtrit.second<<","<<vtrit.third
	//		<<"]/EdgeInFacet(("<<tmpF.first<<","<<tmpF.second<<"),"<<tmp.second<<") "<<tmpF.first<<"["
	//		<<vertex(tmpc, 0)<<","<<vertex(tmpc, 1)<<","<<vertex(tmpc, 2)<<","<<vertex(tmpc, 3)<<"]"<<endl;
	//}
	//===================test end==================//
	
	quit=false;
	//边界线顶点只能关联两条边界边,否则被替换面不是正则的单纯2-D复形
	for(auto itv=VertexOnBoundary.begin();itv!=VertexOnBoundary.end();itv++)
	{
		int singular_point=0;
		for (auto itBE=BoundEdges.begin();itBE!=BoundEdges.end();itBE++)
		{
			Edge tmp=*itBE;
			Facet tmpF=get_facet_cell_link(tmp.first);
			Vertex_pair vp=make_vertex_pair(EdgeInFacet(tmpF,tmp.second));
			if (*itv==vp.first||*itv==vp.second)
			{
				singular_point++;
				if (singular_point>2)
				{
					quit=true;
					break;
				}
			}

		}
		if (quit)
		{
			break;
		}
	}

	for(auto itt=tmpSurfacets.begin();itt!=tmpSurfacets.end();itt++)
	{
		clear_facet_visited(*itt);
	}

}

//求交界面，使用Gabriel的规则决定替换面的延展（主要针对符合pseudo-concave的surface），凸延展
//以种子三角面片为初始，按BFS搜索邻域内可能被替换的表面三角面片
//p插入点；InitFacet为种子三角面片（插入点表面内,所在的cell在影响域中,插入点表面外,其mirror facet所在cell不属于影响域）；
//BoundEdges为返回值，被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换；BoundEdge也是边界线集合，以端点(vertex_id,vertex_id)形式表示
//Begin所在cell属于影响域,mirror面所在cell不属于影响域
//VertexOnBoundary为边界线上的顶点，用于孤立点的检查；
//VertexOnConflictFacet为(被替换面上的顶点id,与之关联的表面三角面片被删除的个数)，用于孤立点的检查
//PInside为p关于表面的位置
//SurfacetsInConflict为影响域内的表面三角面片(该面及其mirror facet所在的cell都是conflict)
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
find_conflict_surfacet(const T* P,Facet InitFacet,list<Edge>& BoundEdges,set_pair_int& BoundEdge,Facet& Begin,std::set<Vertex_Idtype>& VertexOnBoundary,
															map<Vertex_Idtype,size_t>& VertexOnConflictFacet,bool PInside,
															vector<Facet_Idtype>& SurfacetsInConflict,bool& quit)
{
	list<Facet> illSurfacets;
	vector<Facet> tmpSurfacets;
	list<Facet> neighborSurfacets;
	neighborSurfacets.push_back(InitFacet);
	tmpSurfacets.push_back(InitFacet);
	mark_facet_visited(InitFacet);
	mark_facet_visited(InitFacet);

	vector<Facet_Idtype> idtmpSurfacets;
	Facet_Idtype fid=CellSurface.get_tuple_element(InitFacet.first,InitFacet.second);
	idtmpSurfacets.push_back(fid);
	//-----------------求InitFacet邻域的表面三角面片（即与InitFacet有公共点的）------------------//
	while(!neighborSurfacets.empty())
	{
		Facet tmpSurfacet=*neighborSurfacets.begin();
	
		//--------------for test------------//
		Vertex_triple vtriT0=make_vertex_triple(tmpSurfacet);
		//===============test end==========//
		neighborSurfacets.pop_front();
		Cell_Idtype c = tmpSurfacet.first;//HHHHHHHHHHHHHHHHHHHHHHHHHH-------c肯定被label为inside
		int li = tmpSurfacet.second;
		Facet nei;
		for(int i=0;i<3;i++)
		{
			nei=neighbor_surfacet(tmpSurfacet,i);
	
			//--------------for test------------//
			Vertex_triple vtriT1=make_vertex_triple(nei);
			//===============test end==========//
			if(visited_for_facet_extractor(nei)==0)
			{
				//如果是broken surface facet
				if(is_facet_in_conflict(nei))//is_facet_in_conflict(nei)面所关联的两个cell都是conflict
				{
					//如果此面是表面则，标记为2；否则标记为1
					mark_facet_visited(nei);
					mark_facet_visited(nei);
					neighborSurfacets.push_back(nei);
					tmpSurfacets.push_back(nei);

					Facet_Idtype fid=CellSurface.get_tuple_element(nei.first,nei.second);
					idtmpSurfacets.push_back(fid);
				}
				//判断如果是boundary surface facet(该面对应的一个cell在影响域，另一个是boundary)
				else if((is_on_boundary(nei.first)&&is_in_conflict(neighbor(nei.first,nei.second)))||
					(is_on_boundary(neighbor(nei.first,nei.second))&&is_in_conflict(nei.first)))
				{
					T cosDihedral[3];//cosDihedral[0]存零号点相对与底面夹角的余弦
					//if(PInside)
					cos_dihedral_point_to_facet(P,nei,cosDihedral);
					/*else
						cos_dihedral_point_to_facet(P,mirror_facet(nei),cosDihedral);*/
					//------------------for test----------------//
					/*T cosT[3];
					cos_dihedral_point_to_facet(P,nei,cosT);*/
					//==================test end================//

					
					
					//加入方向判断
					Vertex_triple vf=make_vertex_triple(nei);
					Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), P);
					if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
					{
						//判断加入nei是否为凸延展
						bool isConvexExtension=true;
						for(int j=0;j<3;j++)
						{
							Facet neiEdge=neighbor_surfacet(nei,j);
							//--------------------------for test-------------------------//
							Vertex_triple vtriT3=make_vertex_triple(neiEdge);
							//==========================test end=========================//
							//if(cosDihedral[i]<0&&neighbor_surfacet())
							if(cosDihedral[j]<0&&!is_in_conflict(neiEdge.first)&&!is_in_conflict(neighbor(neiEdge.first,neiEdge.second)))
							{
								isConvexExtension=false;
								break;
							}
						}
						//改变选择替换面的规则，以pInside为基础，以是否符合Gabriel规则为准则
						if (PInside)
						{
							bool toBeSurface=false;

							for (int j = 0; j < 3; j++)
							{
								Facet NNei=neighbor_surfacet(nei,j);
								if (visited_for_facet_extractor(NNei)==2)
								{
									Vertex_Idtype vertexOpp=vertex_facet(nei,j);
									Vertex_pair vertexOppEdge=make_vertex_pair(EdgeInFacet(nei,j));
									bool isInShpere=GeometricTraits<T>::side_of_bounded_sphereC3(point(vertexOppEdge.first),P,
										point(vertexOppEdge.second),point(vertexOpp))==ON_BOUNDED_SIDE;
									if (isInShpere)
									{
										toBeSurface=true;
									}
								}
							}
							if (!toBeSurface)
							{
								Vertex_triple ftrip=make_vertex_triple(nei);
								toBeSurface=GeometricTraits<T>::side_of_bounded_sphereC3(point(ftrip.first),
									point(ftrip.second),point(ftrip.third),P)==ON_BOUNDED_SIDE;
							}
							if(toBeSurface)//&&isConvexExtension)
							{
								mark_facet_visited(nei);
								mark_facet_visited(nei);
								neighborSurfacets.push_back(nei);
								tmpSurfacets.push_back(nei);

								Facet_Idtype fid=CellSurface.get_tuple_element(nei.first,nei.second);
								idtmpSurfacets.push_back(fid);
							}
							else
							{
								mark_facet_visited(nei);
								illSurfacets.push_back(nei);
							}
						}
						else
						{
							bool toBeSurface=true;
							for (int j = 0; j < 3; j++)
							{
								Facet NNei=neighbor_surfacet(nei,j);
								if (visited_for_facet_extractor(NNei)!=2)
								{
									Vertex_Idtype vertexOpp=vertex_facet(nei,j);
									Vertex_pair vertexOppEdge=make_vertex_pair(EdgeInFacet(nei,j));
									bool isInShpere=GeometricTraits<T>::side_of_bounded_sphereC3(point(vertexOppEdge.first),
										point(vertexOppEdge.second),P,point(vertexOpp))==ON_BOUNDED_SIDE;
									if (isInShpere)
									{
										toBeSurface=false;
									}
								}
							}
							if(toBeSurface)//&&isConvexExtension)
							{
								mark_facet_visited(nei);
								mark_facet_visited(nei);
								neighborSurfacets.push_back(nei);
								tmpSurfacets.push_back(nei);

								Facet_Idtype fid=CellSurface.get_tuple_element(nei.first,nei.second);
								idtmpSurfacets.push_back(fid);
							}
							else
							{
								mark_facet_visited(nei);
								illSurfacets.push_back(nei);
							}
						}

					
					}//end of if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
				}// end of else if(...)
			
			}//end of if(visited_for_facet_extractor(nei)==0)
		}//end of for
	}//end of while

	//对上一步中求得的交界面删除其中的可能造成折叠的surface facet
	while(!illSurfacets.empty())
	{
		Facet tmpIllSurfacet=*(illSurfacets.begin());
		clear_facet_visited(tmpIllSurfacet);
		illSurfacets.pop_front();
		for(int i=0;i<3;i++)
		{
			Facet nei=neighbor_surfacet(tmpIllSurfacet,i);
			//---------------for test------------//
			Vertex_triple vtriNei=make_vertex_triple(nei);
			//==============test end============//
			if(visited_for_facet_extractor(nei)==2&&!is_facet_in_conflict(nei)&&nei!=InitFacet)
			{
				Indextype iNei=neighbor_index_facet(nei,tmpIllSurfacet);
				if(cos_dihedral_point_to_Edge(P,nei,iNei)<0)
				{
					illSurfacets.push_back(nei);
					clear_facet_visited(nei);

					Facet_Idtype fid=CellSurface.get_tuple_element(nei.first,nei.second);
					idtmpSurfacets.erase(find(idtmpSurfacets.begin(),idtmpSurfacets.end(),fid));
				}
			}
		}
	}



	quit=false;
	if (SurfacetsInConflict.size()>idtmpSurfacets.size())
	{

		quit=true;
	}

	for (auto its=SurfacetsInConflict.begin();its!=SurfacetsInConflict.end();its++)
	{
		if (find(idtmpSurfacets.begin(),idtmpSurfacets.end(),*its)==idtmpSurfacets.end())
		{
			quit=true;
			break;
		}
	}

	//if (idtmpSurfacets.size()>4)
	//{
	//	quit=true;
	//}

	if (quit)
	{
		for(auto itt=tmpSurfacets.begin();itt!=tmpSurfacets.end();itt++)
		{
			clear_facet_visited(*itt);
		}
		return;
	}



	//重新以种子三角面片为基础搜索符合条件的交界面三角面片及其边界边
	//EdgeInFacet是(Facet,0/1/2)表示,其中Facet以(cell,0/1/2/3)表示
	//BoundEdges被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换
	//BoundEdge也是边界线集合，以端点(vertex_id,vertex_id)形式表示
	//边界边输出是按着顺时针顺序，初始化也必须顺时针(从表面内看,表面外看则逆时针),所以按着相对顶点0->2->1加入初始化边界边
	mark_facet_visited(InitFacet);
	//---------------for test------------//
	Vertex_triple vtri_tmp = make_vertex_triple(InitFacet);
	//===============test end============//

	EdgeInFacet tmpBoundEdge=mirror_edge(EdgeInFacet(InitFacet,0));
	//turn_edgeinfacet_to_edge(tmpBoundEdge)把tmpBoundEdge从（（cell_id,0/1/2/3）,0/1/2）改成(面id,0/1/2)
	Edge bound0= turn_edgeinfacet_to_edge(tmpBoundEdge);
	BoundEdges.push_back(bound0);
	//---------------for test------------//
	Vertex_triple vtri_tmp1 = make_vertex_triple(tmpBoundEdge.first);
	//===============test end============//

	tmpBoundEdge=mirror_edge(EdgeInFacet(InitFacet,2));
	//turn_edgeinfacet_to_edge(tmpBoundEdge)把tmpBoundEdge从（（cell_id,0/1/2/3）,0/1/2）改成(面id,0/1/2)
	Edge bound2= turn_edgeinfacet_to_edge(tmpBoundEdge);
	BoundEdges.push_back(bound2);
	//---------------for test------------//
	Vertex_triple vtri_tmp2 = make_vertex_triple(tmpBoundEdge.first);
	//===============test end============//

	tmpBoundEdge=mirror_edge(EdgeInFacet(InitFacet,1));
	//turn_edgeinfacet_to_edge(tmpBoundEdge)把tmpBoundEdge从（（cell_id,0/1/2/3）,0/1/2）改成(面id,0/1/2)
	Edge bound1= turn_edgeinfacet_to_edge(tmpBoundEdge);
	BoundEdges.push_back(bound1);
	//---------------for test------------//
	Vertex_triple vtri_tmp3 = make_vertex_triple(tmpBoundEdge.first);
	//===============test end============//

	list<Facet> tmpConflitFacets;
	tmpConflitFacets.push_back(InitFacet);
	
	//-----------------for test-------------//
	/*cout<<"被替换面：     "<<endl;*/
	//=================test end==============//

	while(!tmpConflitFacets.empty())
	{
		Facet tmp=*(tmpConflitFacets.begin());
		tmpConflitFacets.pop_front();
		Vertex_triple vtrif=make_vertex_triple(tmp);
		/*VertexOnConflictFacet.insert(vtrif.first);
		VertexOnConflictFacet.insert(vtrif.second);
		VertexOnConflictFacet.insert(vtrif.third);*/
		VertexOnConflictFacet[vtrif.first]++;
		VertexOnConflictFacet[vtrif.second]++;
		VertexOnConflictFacet[vtrif.third]++;

		//------------------for test--------------//
		Vertex_triple vtriTC0=make_vertex_triple(tmp);
		Cell_Idtype tmpc=tmp.first;
		Facet_Idtype fid=CellSurface.get_tuple_element(tmpc,tmp.second);
		//cout<<fid<<"("<<vtriTC0.first<<","<<vtriTC0.second<<","<<vtriTC0.third<<")"
		//	<<"/("<<tmp.first<<","<<tmp.second<<")["<<tmpc<<"("
		//	<<vertex(tmpc, 0)<<","<<vertex(tmpc, 1)<<","<<vertex(tmpc, 2)<<","<<vertex(tmpc, 3)<<")]"<<endl;
		//==================test end==============//
		
		for (int i = 0; i<3; i++)
		{
			Facet nei = neighbor_surfacet(tmp, i);
			//---------------for test------------//
			Vertex_triple vtriT = make_vertex_triple(nei);
			//===============test end============//
			//求tmp是nei的第几个相邻的surface facet
			Indextype iNei = neighbor_index_facet(nei, tmp);
			if (visited_for_facet_extractor(nei) == 2)
			{
				
				mark_facet_visited(nei);
				tmpConflitFacets.push_back(nei);

				//更新交界面的边界边队列
				vector<int> numNei;
				int numOfList = 0;

				//--------------修改------------//
				auto itEdge = BoundEdges.begin();
				for (; itEdge != BoundEdges.end(); itEdge++)
				{
					numOfList++;
					Facet tmpBoundFacet0 = get_facet_cell_link((*itEdge).first);
					//--------------for test----------//
					Vertex_triple vtriT0 = make_vertex_triple(tmpBoundFacet0);
					//==============test end==========//
					if (tmpBoundFacet0 == nei)
					{
						numNei.push_back((*itEdge).second);
						BoundEdges.erase(itEdge++);
						break;
					}
				}
				//当list链的第一个被删除时，应该检查最后一个是否也应该被删除
				if (numOfList == 1)
				{
					auto endEdge = BoundEdges.back();
					Facet tmpBoundFacet00 = get_facet_cell_link(endEdge.first);
					//--------------for test----------//
					Vertex_triple vtriT00 = make_vertex_triple(tmpBoundFacet00);
					//==============test end=========//
					if (tmpBoundFacet00 == nei)
					{

						numNei.push_back(endEdge.second);
						BoundEdges.pop_back();
						endEdge = BoundEdges.back();
						Facet tmpBoundFacet01 = get_facet_cell_link(endEdge.first);
						//--------------for test----------//
						Vertex_triple vtriT01 = make_vertex_triple(tmpBoundFacet00);
						//==============test end=========//
						if (tmpBoundFacet01 == nei)
						{
							numNei.push_back(endEdge.second);
							BoundEdges.pop_back();
						}
					}

				}
				if (numNei.size() != 3)
				{
					//list 中间的插入与删除要处理好
					Facet tmpBoundFacet1 = make_pair(-1, -1);
					if (itEdge != BoundEdges.end())
					{
						tmpBoundFacet1 = get_facet_cell_link((*itEdge).first);
						//----------------for test------------//
						Vertex_triple vtriT1 = make_vertex_triple(tmpBoundFacet1);
						int ttttt=0;
						//================test end===========//
					}				
					if (itEdge != BoundEdges.end() && tmpBoundFacet1 == nei)
					{
						numNei.push_back((*itEdge).second);
						BoundEdges.erase(itEdge++);
						Facet tmpBoundFacet1 = make_pair(-1, -1);
						if (itEdge != BoundEdges.end())
						{
							tmpBoundFacet1 = get_facet_cell_link((*itEdge).first);
							//----------------for test------------//
							Vertex_triple vtriT1 = make_vertex_triple(tmpBoundFacet1);
							//int ttttt=0;
							//================test end===========//
						}
						if (itEdge != BoundEdges.end() && tmpBoundFacet1 == nei)
						{
							numNei.push_back((*itEdge).second);
							BoundEdges.erase(itEdge++);
						}
					}

				}
				if (numNei.size() == 1)  //新加两条边界边
				{
					int i0 = numNei[0];
					auto leftRight = make_edge_index(i0);
					EdgeInFacet oppEdge0 = mirror_edge(EdgeInFacet(nei, leftRight.second));
					EdgeInFacet oppEdge1 = mirror_edge(EdgeInFacet(nei, leftRight.first));
					//--------------for test-----------//
					Vertex_triple vtriT3 = make_vertex_triple(oppEdge0.first);
					Vertex_triple vtriT4 = make_vertex_triple(oppEdge1.first);
					//==============test end==========//
					BoundEdges.insert(itEdge, turn_edgeinfacet_to_edge(oppEdge0));
					BoundEdges.insert(itEdge, turn_edgeinfacet_to_edge(oppEdge1));
				}
				else if (numNei.size() == 2)  //新加一条边界边
				{
					for (int i2 = 0; i2<3; i2++)
					{
						if (i2 != numNei[0] && i2 != numNei[1])
						{
							EdgeInFacet oppEdge = mirror_edge(EdgeInFacet(nei, i2));
							//--------------for test-----------//
							Vertex_triple vtriT2 = make_vertex_triple(oppEdge.first);
							//==============test end==========//
							BoundEdges.insert(itEdge, turn_edgeinfacet_to_edge(oppEdge));
							break;
						}
					}
				}



			}
			else if (visited_for_facet_extractor(nei) == 0)
			{

				Vertex_pair vEdge = make_vertex_pair(EdgeInFacet(tmp, i));
				BoundEdge.insert(make_pair(vEdge.first, vEdge.second));
				BoundEdge.insert(make_pair(vEdge.second, vEdge.first));
				VertexOnBoundary.insert(vEdge.first);
				VertexOnBoundary.insert(vEdge.second);
			}
		}
	}


	//-------------------for test(打印 BoundEdges)------------------//
	//EdgeInFacet是(Facet,0/1/2)表示,其中Facet以(cell,0/1/2/3)表示
	//BoundEdges被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换
	//BoundEdge也是边界线集合，以端点(vertex_id,vertex_id)形式表示
	//cout<<"conflict boundary edge:"<<endl;
	//for(auto itBE=BoundEdges.begin();itBE!=BoundEdges.end();itBE++)
	//{
	//	Edge tmp=*itBE;
	//	Facet tmpF=get_facet_cell_link(tmp.first);
	//	//--------------for test----------//
	//	Cell_Idtype tmpc=tmpF.first;
	//	Vertex_triple vtrit = make_vertex_triple(tmpF);
	//	//==============test end=========//
	//	Vertex_pair vp=make_vertex_pair(EdgeInFacet(tmpF,tmp.second));
	//	cout<<"("<<vp.first<<","<<vp.second<<")"
	//		<<"/Edge("<<tmp.first<<","<<tmp.second<<") "<<tmp.first<<"["<<vtrit.first<<","<<vtrit.second<<","<<vtrit.third
	//		<<"]/EdgeInFacet(("<<tmpF.first<<","<<tmpF.second<<"),"<<tmp.second<<") "<<tmpF.first<<"["
	//		<<vertex(tmpc, 0)<<","<<vertex(tmpc, 1)<<","<<vertex(tmpc, 2)<<","<<vertex(tmpc, 3)<<"]"<<endl;
	//}
	//===================test end==================//


	quit=false;
	//EdgeInFacet是(Facet,0/1/2)表示,其中Facet以(cell,0/1/2/3)表示
	//BoundEdges被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换
	//BoundEdge也是边界线集合，以端点(vertex_id,vertex_id)形式表示
	//边界线顶点只能关联两条边界边,否则被替换面不是正则的单纯2-D复形
	for(auto itv=VertexOnBoundary.begin();itv!=VertexOnBoundary.end();itv++)
	{
		int singular_point=0;
		for (auto itBE=BoundEdges.begin();itBE!=BoundEdges.end();itBE++)
		{
			Edge tmp=*itBE;
			Facet tmpF=get_facet_cell_link(tmp.first);
			Vertex_pair vp=make_vertex_pair(EdgeInFacet(tmpF,tmp.second));
			if (*itv==vp.first||*itv==vp.second)
			{
				singular_point++;
				if (singular_point>2)
				{
					quit=true;
					break;
				}
			}
		}
		if (quit)
		{
			break;
		}
	}

	if (BoundEdges.empty())
	{
		quit=true;
	}

	if (quit)
	{
		for(auto itt=tmpSurfacets.begin();itt!=tmpSurfacets.end();itt++)
		{
			clear_facet_visited(*itt);
		}
		return;
	}


	//求Begin
	auto edgeBound=BoundEdges.begin();
	Facet facetEdge=get_facet_cell_link((*edgeBound).first);
	Facet conflictBoundary=neighbor_surfacet(facetEdge,(*edgeBound).second);
	//---------------------------for test-----------------------//
	Vertex_triple vtriT=make_vertex_triple(conflictBoundary);
	//===========================test end=======================//
	if(!PInside&&is_on_boundary(conflictBoundary.first))//neighbor(conflictBoundary.first,conflictBoundary.second)))
	{
		Begin=mirror_facet(conflictBoundary);
	}
	else
	{
		Cell_Idtype c=conflictBoundary.first;
		Indextype li=conflictBoundary.second;
		Indextype iif=neighbor_index_facet(conflictBoundary,facetEdge);
		Vertex_Idtype vii=vertex_facet(conflictBoundary,iif);
		Indextype ii=vertex_index(c,vii);

		Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
		Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
		Cell_Idtype cur1 = c;
		int zz1 = ii;
		Cell_Idtype n1 = neighbor(cur1, zz1);
		// turn around the oriented edge vj1 vj2
		//沿着有向边vj1->vj2依次搜索与vj1->vj2关联的四面体,搜索方向与vj1->vj2方向（右手拇指方向）符合右手螺旋定则
		//vj1,vj2为被替换面边界线顶点
		while (is_in_conflict(n1)) {

			cur1 = n1;
			zz1 = Triangulation_utils_3::next_around_edge(vertex_index(n1, vj1), vertex_index(n1, vj2));
			n1 = neighbor(cur1, zz1);
		}
		// Now n1 is outside region, cur1 is inside.
		//此时面Begin含边界线,并且cur1属于影响域,但neighbor(cur1,zz1)不属于影响域
		Begin = Facet(cur1,zz1);
	}
	for(auto itt=tmpSurfacets.begin();itt!=tmpSurfacets.end();itt++)
	{
		clear_facet_visited(*itt);
	}

}



template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::surface_facet_nearest_to_point(const T* p, Facet& NearestFacet,
	std::vector<Facet>& FacetToDelete, std::vector<Facet>& FacetUndelete)
{
	Vertex_Idtype nearestVertex = PointDataSet.nearest_inexact(p);
	
	incident_surface_facet(nearestVertex, std::back_inserter(FacetToDelete),std::back_inserter(FacetUndelete));
	if (FacetToDelete.empty())
	{
		NearestFacet = incident_facet_with_min_longest_edge(p,FacetUndelete.begin(),FacetUndelete.end(),nearestVertex);
		return false;
	}
	else
	{
		NearestFacet = incident_facet_with_min_longest_edge(p,FacetToDelete.begin(),FacetToDelete.end(),nearestVertex);
		return true;
	}
}

template<typename T, typename T_INFO>
template<class IteratorIncidentFacet>
Facet DataStructure<T, T_INFO>::
incident_facet_with_min_longest_edge(const T* p, IteratorIncidentFacet FacetBegin, 
	IteratorIncidentFacet FacetEnd, Vertex_Idtype VertexIncident)
{
	double minLongestLength=COMPUTATION_DOUBLE_MAX;
	Facet nearestFacet;
	IteratorIncidentFacet itFacet;
	for (itFacet = FacetBegin; itFacet != FacetEnd; itFacet++)
	{
		//-----------------itFacet必须为boundary facet或boundary facet对立面或conflict facet----------------//

		if (!is_in_conflict((*itFacet).first))
		{
			if (!is_in_conflict(neighbor((*itFacet).first, (*itFacet).second)))
				continue;
		}
		double longestLength = 0;
		Vertex_triple vertexFacet = make_vertex_triple(*itFacet);
		Vertex_Idtype vertexTemp = -1;
		for (int i = 0; i < 3; i++)
		{
			//vertexTemp = vertexFacet.get<i>();
			if (i == 0)
				vertexTemp = vertexFacet.first;
			else if (i == 1)
				vertexTemp = vertexFacet.second;
			else
				vertexTemp = vertexFacet.third;

			if (vertexTemp != VertexIncident)
			{
				double distantTemp = distant_to_vertex(p,vertexTemp);
				if (distantTemp > longestLength)
					longestLength = distantTemp;
			}
		}
		if (longestLength < minLongestLength)
		{
			minLongestLength = longestLength;
			nearestFacet = *itFacet;
		}
	}
	return nearestFacet;
}

template<typename T, typename T_INFO>
double DataStructure<T, T_INFO>::distant_to_vertex(const T* p, Vertex_Idtype IdV)
{
	const T* ptr = point(IdV);
	return NumericComputation<T>::SquareDistance(p,ptr);
}


template<typename T, typename T_INFO>
template<class OutputIteratorFacet>
void DataStructure<T, T_INFO>::
incident_surface_facet(Vertex_Idtype IdV, OutputIteratorFacet FacetToDelete, OutputIteratorFacet FacetUndelete)
{
	if (dimension() < 2)
		return;

	std::vector<Cell_Idtype> tmp_cells;
	tmp_cells.reserve(64);
	std::vector<Facet> tmp_facets;
	tmp_facets.reserve(64);
	std::vector<Vertex_Idtype> temp_vertices;
	if (dimension() == 3)
		incident_cells_3_withoutmark(IdV, cell(IdV), std::back_inserter(tmp_cells), FacetToDelete, FacetUndelete);


	typename std::vector<Cell_Idtype>::iterator cit;
	for (cit = tmp_cells.begin();
		cit != tmp_cells.end();
		++cit)
	{
		clear_cell_visited(*cit);
	}
}


template<typename T, typename T_INFO>
template<typename IncidentCellIterator, typename OutputIteratorFacet>
void DataStructure<T, T_INFO>::
incident_cells_3_withoutmark(Vertex_Idtype IdV, Cell_Idtype d, IncidentCellIterator IncidentCell,
	OutputIteratorFacet FacetToDelete, OutputIteratorFacet FacetUndelete)
{

	std::stack<Cell_Idtype> cell_stack;
	cell_stack.push(d);
	mark_cell_visited(d);
	*IncidentCell++ = d;

	do {
		Cell_Idtype c = cell_stack.top();
		cell_stack.pop();

		for (int i = 0; i < 4; ++i) {
			if (vertex(c, i) == IdV)
				continue;
			Cell_Idtype next = neighbor(c, i);

			if (is_surface(Facet(c, i)))
			{
				if (is_facet_in_conflict(Facet(c, i)))
				{
					*FacetToDelete++ = Facet(c, i);
				}
				else
				{
					*FacetUndelete++ = Facet(c, i);
				}
			}


			if (visited_for_cell_extractor(next))
				continue;
			cell_stack.push(next);
			mark_cell_visited(next);
			*IncidentCell++ = next;
		}
	} while (!cell_stack.empty());

}


//不对cell做conflict mark的前提下求incident_cell 3d, 不标记是否被删除
//找出顶点IdV关联的四面体和表面面片,其中d是IdV的关联四面体,
//IncidentCells是关联四面体集合的迭代器（back_inserter），Surfacets是关联的表面集合的迭代器（back_inserter）
template<typename T, typename T_INFO>
template<typename IteratorCell, class IteratorFacet>
void DataStructure<T, T_INFO>::
incident_cells_surfacet_3_withoutmark(Vertex_Idtype IdV, Cell_Idtype d,IteratorCell IncidentCells,
		IteratorFacet Surfacets)
{
	std::stack<Cell_Idtype> cell_stack;
	cell_stack.push(d);
	mark_cell_visited(d);
	*IncidentCells++ = d;

	do {
		Cell_Idtype c = cell_stack.top();
		cell_stack.pop();

		for (int i = 0; i<4; ++i) {
			if (vertex(c, i) == IdV)
				continue;
			Cell_Idtype next = neighbor(c, i);

		    if (is_surface(Facet(c,i)))
			{
			   *Surfacets++=Facet(c,i);
			}
				


			if (visited_for_cell_extractor(next))
				continue;
			cell_stack.push(next);
			mark_cell_visited(next);
			*IncidentCells++ = next;
		}
	} while (!cell_stack.empty());
}

//判断V是否为孤立点;SurfaceFacet存V关联四面体,并与V相对面,使得V到该面三个顶点中最长边最小的面
//若为孤立点,Inside为该孤立点在表面内还是表面外
template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
is_vertex_isolate(Vertex_Idtype V,Facet& SurfaceFacet,bool& Inside)
{
	std::stack<Cell_Idtype> cell_stack;
	std::vector<Cell_Idtype> cellsIncident;
	Cell_Idtype d=cell(V);
	cell_stack.push(d);
	mark_cell_visited(d);
	cellsIncident.push_back(d);
	Inside=is_label_inside(d);
	bool isIso=true;
	double dis=COMPUTATION_DOUBLE_MAX;
	//-----------------for test----------------//
	//cout<<"to decide if vertex"<<V<<"is isolate"<<endl;
	//=================test end================//
	do {
		Cell_Idtype c = cell_stack.top();
		cell_stack.pop();

		for (int i = 0; i<4; ++i) {
			if (vertex(c, i) == V)
				continue;
			Cell_Idtype next = neighbor(c, i);

			if(Inside!=is_label_inside(next))
			{
				isIso=false;
				break;
			}

			if (visited_for_cell_extractor(next))
				continue;
			cell_stack.push(next);
			mark_cell_visited(next);
			cellsIncident.push_back(next);
		}
		if(isIso==false)
		{
			break;
		}

	} while (!cell_stack.empty());

	for(auto ic=cellsIncident.begin();ic!=cellsIncident.end();ic++)
	{
		clear_cell_visited(*ic);
	}

	return isIso;

}

//插入孤立点V
template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
insert_vertex_into_surface(Vertex_Idtype V)
{
	//寻找最近邻点
	vector<Vertex_Idtype> vs;
	vector<Facet> incidentSurfacets;
	list<Vertex_Idtype> nearest(5,-1); //nearest初始化为五个-1
	bool pInside=is_label_inside(cell(V)); //pInside为孤立点在表面内还是表面外
	//---------------for test------------//
	//cout<<"插入孤立点："<<V<<endl;
	//cout<<"孤立点在表面内: "<<pInside<<endl;
	//===============test end============//

	
	const T* p=point(V);
	vs.reserve(32);
	vector<Facet> initFacet;
	incidentSurfacets.reserve(32);
	//求孤立点V关联顶点,vs是V关联顶点的集合,incidentSurfacets是关联的表面(umbrella)的集合
	incident_vertex_and_surface_facet(V,std::back_inserter(vs),std::back_inserter(incidentSurfacets));
	//更新孤立点p的最近邻点
	for(auto vsit=vs.begin();vsit!=vs.end();++vsit)
	{
		refresh_K_nearest(p,*vsit,nearest);			
	}
	//----------------------for test(输出最近邻点)------------------//
	//cout<<"最近邻点:  ";
	//for(auto iNV=nearest.begin();iNV!=nearest.end();iNV++)
	//{
	//	cout<<*iNV<<" ";
	//}
	//cout<<endl;
	//==========================test end============================//

	vector<Facet> incidentSurfacets_KNN;
	vector<Facet> incidentSurfacets_Sparse;
	vector<Facet> surfacetInterConflict;
	vector<Facet_Idtype> surfacetsInConflict;
	std::map<Vertex_Idtype,size_t> vertexWithDeletedSurfacet;
	Facet nearestFacet(-1,-1);
	Facet seedFacet(-1,-1);
	auto fit=nearest.begin();

	if (!incidentSurfacets.empty())
	{
		return true;
	}

	//找出顶点*fit关联的四面体和表面面片,incidentSurfacets是关联的表面集合
	incident_surface_facet(*fit,std::back_inserter(incidentSurfacets));
	
	//pInside是孤立点关于表面的位置,true为表面内部,false为表面外部;V为孤立点顶点编号;
	//nearest为孤立点的K个最近点,incidentSurfacets为与最邻近点关联的初始表面三角面片;
	//incidentSurfacets_KNN作为膨胀和收缩后与KNN关联的所有表面三角面片;
	//incidentSurfacets_Sparse为KNN领域内含有较长边的表面三角面片集合
	update_surfacets_KNN(pInside,V,nearest,incidentSurfacets,std::back_inserter(incidentSurfacets_KNN),std::back_inserter(incidentSurfacets_Sparse));

	//求孤立点的最近的三角面片
	//surfacetInterConflict是表面面片(该面是表面,该面所在的cell属于影响域,或其mirror facet所在cell属于影响域,满足其一即可),插入孤立点时候不需要
	//incidentSurfacets_KNN为与KNN关联的备选表面三角面片;incidentSurfacets_Sparse为KNN领域内含有较长边的备选表面三角面片集合
	//PInside表示孤立点关于表面的位置;nearestFacet为返回的最邻近表面三角面片;
	//surfacetsInConflict为影响域内的表面三角面片(该面及其mirror facet所在的cell都是conflict),插入孤立点时候不需要
	//vertexWithDeletedSurfacet为(删除的表面三角面片中的点id,与之关联的表面三角面片被删除的个数),用于新生成的孤立点的检测
	nearest_surfacet(p,surfacetInterConflict,incidentSurfacets_KNN,incidentSurfacets_Sparse,pInside,nearestFacet,surfacetsInConflict,
		vertexWithDeletedSurfacet);

	if (nearestFacet.first==-1)
	{
		return false;
	}
	//------------------------for test---------------------//
	//if (nearestFacet.first!=-1)
	//{
	//	Vertex_triple vtriF0=make_vertex_triple(nearestFacet);
	//	cout<<"初始化最近表面: ("<<vtriF0.first<<","<<vtriF0.second<<","<<vtriF0.third<<")"
	//		<<"/("<<nearestFacet.first<<","<<nearestFacet.second<<")"<<endl;
	//}
	//由孤立点V最邻近表面三角面片nearestFacet求种子三角面片seedFacet(种子面片关联的两个cell,有一个cell中,该种子面片相对的顶点为孤立点)
	//V为孤立点,pInside为孤立点关于表面的位置
	seed_facet_of_isolate_vertex(V,nearestFacet,pInside,seedFacet);

	//------------------------for test---------------------//
	//if (seedFacet.first!=-1)
	//{
	//	Vertex_triple vtriF1=make_vertex_triple(seedFacet);
	//	cout<<"种子面片: ("<<vtriF1.first<<","<<vtriF1.second<<","<<vtriF1.third<<")"
	//		<<"/("<<seedFacet.first<<","<<seedFacet.second<<")"<<endl;
	//}
	//=======================test end======================//
	Facet tmpF;
	if (!is_vertex_isolate(V,tmpF,pInside))
	{
		return true;
	}

	if (seedFacet.first!=-1)
	{
		//插入孤立点V
		//V孤立点,seedFacet种子表面三角面片,pInside为孤立点关于表面的位置
		link_isolate_to_surface(V,seedFacet,pInside);
		return true;
	}
	else
	{
		return false;
	}



}


template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
insert_vertex_into_surface(Vertex_Idtype V,Facet Surfacet,bool Inside)
{
	if(Inside)
	{
		Cell_Idtype c=Surfacet.first;
		label_cell_side(c,false);
		delete_surface_facet(Surfacet);
		int iv=vertex_index(c,V);
		for(int i=0;i<4;i++)
		{
			if(i!=iv)
			{
				Facet sf=mirror_facet(Facet(c,i));
				create_surface_facet(sf.first,sf.second);
			}
		}
	}
	else
	{
		Facet sf=mirror_facet(Surfacet);
		Cell_Idtype c=sf.first;
		label_cell_side(c,true);
		delete_surface_facet(Surfacet);
		int iv=vertex_index(c,V);
		for(int i=0;i<4;i++)
		{
			if(i!=iv)
			{
				create_surface_facet(c,i);
			}
		}
	}
	return true;
}


//处理孤立点
//VertexWithSurfacetDeleted为(被替换面上的顶点id,与之关联的表面三角面片被删除的个数)
//VertexWithSurfacetCreated为边界线上的顶点
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
isolate_vertices(std::map<Vertex_Idtype,size_t> VertexWithSurfacetDeleted,std::set<Vertex_Idtype> VertexWithSurfacetCreated)
{
	//检测是否有孤立点,同时将孤立点加入孤立点链表集合中,且不重复
	if(has_isolate_vertex(VertexWithSurfacetDeleted,VertexWithSurfacetCreated))
	{
		
			
		bool loop=false;
		Vertex_Idtype loopStart=*(IsolateVertices.begin());
		bool insertSuccess=false;
		bool handle=false;
		int t=IsolateVertices.size()+1;
		while((!IsolateVertices.empty()&&!loop)&&t--)
		{
			//bool inside;
			Facet surfacet=Facet(-1,-1);
			auto iiv=IsolateVertices.begin();
			Vertex_Idtype viiv=*iiv;
			if(viiv==loopStart&&handle)
			{
				loop=true;
				continue;
			}
			handle=true;	
			//---------------------for test------------------//
			/*cout<<endl<<"处理孤立点: "<<viiv<<endl;*/
			//====================test end===================//
			//is_vertex_isolate(*iiv,surfacet,inside);
			//插入孤立点viiv
			if(insert_vertex_into_surface(viiv))
			{
				insertSuccess=true;
				if (!IsolateVertices.empty())
				{
					IsolateVertices.pop_front();
				}	
				//---------------------for test------------------//
				/*cout<<"insert isolate point:"<<viiv<<" success"<<endl<<endl;*/
				//====================test end===================//
			}
			else
			{
				if(insertSuccess)
				{
					loopStart=viiv;
				}
				if (!IsolateVertices.empty())
				{
					IsolateVertices.pop_front();
				}				
				IsolateVertices.push_back(viiv);
				//insert_to_isolate_vertices(viiv);
				insertSuccess=false;
				//---------------------for test------------------//
				/*cout<<"insert isolate point:"<<viiv<<" failure"<<endl<<endl;*/
				//====================test end===================//
			}
			
		}
	}
}

//VertexWithSurfacetDeleted为(被替换面上的顶点id,与之关联的表面三角面片被删除的个数)
//VertexWithSurfacetCreated为边界线上的顶点
//检测是否有孤立点,同时将孤立点加入孤立点链表集合中,且不重复
template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
has_isolate_vertex(std::map<Vertex_Idtype,size_t> VertexWithSurfacetDeleted,std::set<Vertex_Idtype> VertexWithSurfacetCreated)
{
	
	auto ivd=VertexWithSurfacetDeleted.begin();
	while (ivd!=VertexWithSurfacetDeleted.end())
	{
		//被替换面中某顶点的umbrella(至少包含3个表面)全被包含了,产生孤立点,一定不在被替换面边界线上
		if(VertexWithSurfacetCreated.find(ivd->first)==VertexWithSurfacetCreated.end()&&
			ivd->second>2&&ivd->first!=0)
		{
			Vertex_Idtype vi=ivd->first;
			Facet fi=Facet(-1,-1);
			bool inside=false;
			bool is_isolate=is_vertex_isolate(vi,fi,inside);
			if(is_isolate)
			{
				//将孤立点V加入孤立点链表集合中,且不重复
				insert_to_isolate_vertices(vi);
			}
		}
		ivd++;
	}
	if(IsolateVertices.empty())
		return false;
	else
		return true;
}



template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
create_face()
{
	return create_cell();
}

template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
create_face(Vertex_Idtype v0, Vertex_Idtype v1,
	Vertex_Idtype v2)
{
	return create_cell(v0, v1, v2, -1);
}

// The following functions come from TDS_2.
template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
create_face(Cell_Idtype f0, int i0,
	Cell_Idtype f1, int i1,
	Cell_Idtype f2, int i2)
{
	Cell_Idtype newf = create_face(vertex(f0,Triangulation_cw_ccw_2::cw(i0)),
		vertex(f1,Triangulation_cw_ccw_2::cw(i1)),
		vertex(f2,Triangulation_cw_ccw_2::cw(i2)));
	set_adjacency(newf, 2, f0, i0);
	set_adjacency(newf, 0, f1, i1);
	set_adjacency(newf, 1, f2, i2);
	return newf;
}

template<typename T, typename T_INFO>
Cell_Idtype DataStructure<T, T_INFO>::
create_face(Cell_Idtype f0, int i0,
	Cell_Idtype f1, int i1)
{
	Cell_Idtype newf = create_face(vertex(f0,Triangulation_cw_ccw_2::cw(i0)),
		vertex(f1,Triangulation_cw_ccw_2::cw(i1)),
		vertex(f1,Triangulation_cw_ccw_2::ccw(i1)));
	set_adjacency(newf, 2, f0, i0);
	set_adjacency(newf, 0, f1, i1);
	return newf;
}




template<typename T, typename T_INFO>
Facet DataStructure<T, T_INFO>::
mirror_facet(Facet f) 
{
	Cell_Idtype neighbor_cell = neighbor(f.first,f.second);
	const int opposite_index = neighbor_index(neighbor_cell,f.first);
	return Facet(neighbor_cell, opposite_index);
}

template<typename T, typename T_INFO>
EdgeInFacet DataStructure<T, T_INFO>::
mirror_edge(EdgeInFacet E)
{
	Facet fE=E.first;
	
	Facet oppFacet=neighbor_surfacet(fE,E.second);
	
	Indextype ii=neighbor_index_facet(oppFacet,fE);
	return EdgeInFacet(oppFacet,ii);
}


//求某一点的incident vertex和在surface上的umbrella，不做标记
//求IdV关联顶点，IncidentVertices是IdV关联顶点的迭代器（back_inserter），IncidentFacets是关联的表面(为umbrella)的迭代器（back_inserter）
template<typename T, typename T_INFO>
template<class IteratorVertex,class IteratorFacet>
void DataStructure<T, T_INFO>::
incident_vertex_and_surface_facet(Vertex_Idtype IdV,IteratorVertex IncidentVertices,IteratorFacet IncidentFacets)
{
	if (dimension() < 2)
		return ;

	//Visitor visit(IdV, outputCell, this, f);

	std::vector<Cell_Idtype> tmp_cells;
	tmp_cells.reserve(64);
	/*std::vector<Facet> temp_facets;
	tmp_facets.reserve(64);*/
	std::vector<Vertex_Idtype> tmp_vertices;

	if (dimension() == 3)
		//找出顶点IdV关联的四面体和表面面片，其中cell(IdV)是IdV在CellIncident记录的关联四面体
		//tmp_cells是IdV关联四面体，IncidentFacets是关联的表面的迭代器（back_inserter）
		incident_cells_surfacet_3_withoutmark(IdV, cell(IdV),std::back_inserter(tmp_cells),IncidentFacets);


	typename std::vector<Cell_Idtype>::iterator cit;
	for (cit = tmp_cells.begin();
		cit != tmp_cells.end();
		++cit)
	{
		clear_cell_visited(*cit);
		

		//=================vertex extractor===============
		for (int j = 0; j <= dimension(); ++j) {
			Vertex_Idtype w = vertex(*cit,j);
			
			if (w != IdV){

				if (visited_for_vertex_extractor(w)==0){
					//w->visited_for_vertex_extractor = true;
					mark_vertex_visited(w);
					tmp_vertices.push_back(w);
					*IncidentVertices++ = w;
					//treat(c, IdV, j);
				}
			}
		}
		//----------------------vertex extractor end-------------------
	}
	for (std::size_t i = 0; i < tmp_vertices.size(); ++i){
		clear_vertex_visited(tmp_vertices[i]);
	}
}

template<typename T, typename T_INFO>
template<class IteratorFacet>
void DataStructure<T, T_INFO>::
incident_surface_facet(Vertex_Idtype IdV, IteratorFacet Surfacets)
{
	if (dimension() < 2)
		return ;

	std::vector<Cell_Idtype> tmp_cells;
	tmp_cells.reserve(64);

	std::vector<Vertex_Idtype> tmp_vertices;

	if (dimension() == 3)
    {
		//不对cell做conflict mark的前提下求incident_cell 3d, 不标记是否被删除
		//找出顶点IdV关联的四面体和表面面片,
		//tmp_cells是关联四面体集合,Surfacets是关联的表面集合
		incident_cells_surfacet_3_withoutmark(IdV, cell(IdV),std::back_inserter(tmp_cells),Surfacets);
	}


	typename std::vector<Cell_Idtype>::iterator cit;
	for (cit = tmp_cells.begin();
		cit != tmp_cells.end();
		++cit)
	{
		clear_cell_visited(*cit);
	}
}


//更新点p最近邻列表nearest,如果v比nearest里的点距离p更近，就更新nearest，否则不更新nearest
template<typename T, typename T_INFO>
void  DataStructure<T, T_INFO>::
refresh_K_nearest(const T* p,Vertex_Idtype v,list<Vertex_Idtype>& nearest)
{
	auto in=nearest.begin();
	for(;in!=nearest.end();in++)
	{
	//nearest_vertex(p,v,*in)比较一个点p到两个点的最近点
		if(*in==-1||(nearest_vertex(p,v,*in)==v))
		{
			if(*in!=v)
			{
				nearest.insert(in,v);
				nearest.pop_back();
			}
			break;
		}
	}
}

template<typename T, typename T_INFO>
Vertex_Idtype  DataStructure<T, T_INFO>::
nearest_vertex(const T* p,Vertex_Idtype v,Vertex_Idtype w)
{
	if(is_infinite(v))
		return w;
	if(is_infinite(w))
		return v;
	if(NumericComputation<T>::SquareDistance(p,point(v))<NumericComputation<T>::SquareDistance(p,point(w)))
		return v;
	else
		return w;
}

//是否要加入inflate/sculpture的处理
//处理不均匀的稀疏区域时用三角面片的最长边的长度做判断
//InitFacetConflict为在Delaunay conflict region内求得的最近邻面，第一个为facet，第二个为cosDihedral,第三个为当cosDihedral>-0.1时的distance
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
nearest_surface_facet(const T* p,list<Vertex_Idtype> NearestVertices,vector<Facet> InitFacet,
	Triple<Facet,double,double> InitFacetConflict, bool PInside,Facet& NearestFacet)
{

	//step1:求K-NN的incident surface facet中最有可能的nearest surface facet，并标记其中最大的三角面片
	//最小二面角的余弦值
	double cosMinDihedral=-1;
	//最小二面角对应的三角面片
	Facet facetMinDihedral(-1,-1);
	//当P与三角面片形成的二面角均为锐角时，P到此面的最小距离
	double minDistance=COMPUTATION_DOUBLE_MAX;
	//最小距离对应的三角面片
	Facet facetMinDistance(-1,-1);
	//K-NN的incident surfacet中的最大三角面片
	double maxLengthEdge=0;
	Facet facetLarge(-1,-1);

	//markForVertices用于标记第i个最近邻点是否处理，未处理未0，待处理为2,；已经处理且不是孤立点为1，是孤立点或者不包含可能为最近邻的表面为-1
	int markForVertices[5]={0};
	//用于存储K-NN的incident surfacet，且对于P可见
	Facet facetNVertices[5]={Facet(-1,-1)};
	//存储正在处理的vertex 
	Vertex_Idtype vertexHandled=*NearestVertices.begin();
	int rankVertexHandled=0;
	//用于存储vertexHandled的一个incident surface facet并且第一个最邻近点时此facet对于P可见  
	Facet facetHandledInit(-1,-1);
	if (!InitFacet.empty())
	{
		facetHandledInit=*InitFacet.begin();
	}
	//markForVertices[0]=2;
	//facetNVertices[0]=facetHandledInit;
	if (facetHandledInit.first!=-1)
	{
		markForVertices[0]=1;
		facetNVertices[0]=facetHandledInit;
	}
	else
	{
		markForVertices[0]=-1;
	}	
	//表征所有的K-NN是否处理完毕，是为false，否为true
	bool testForAll=true;

	while (testForAll)
	{
		double cosDihedral;
		
		double distance;
		Facet tmpFacet=facetHandledInit;
		bool markForCircleBegin=true;
		if (facetHandledInit.first!=-1)
		{
			while (facetHandledInit!=tmpFacet||markForCircleBegin)
			{
				//将tmpFacet绕vertexHandled旋转时的当前facet		
				Indextype idVertexOfFacet=vertex_index_facet(tmpFacet,vertexHandled);
				std::pair<Indextype,Indextype> idVertexOfEdge=make_edge_index(idVertexOfFacet);
				markForCircleBegin=false;
				Vertex_triple vf=make_vertex_triple(tmpFacet);
				Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), p);
				if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
				{
					if (markForVertices[rankVertexHandled]==2)
					{
						markForVertices[rankVertexHandled]=1;
					}
					if(!PInside)
					{
						cosDihedral=cos_min_dihedral_point_to_facet(p,mirror_facet(tmpFacet));
					}
					else
					{
						cosDihedral=cos_min_dihedral_point_to_facet(p,tmpFacet);
					}
					if(cosDihedral>-0.0)
					{
						distance =GeometricTraits<T>::distance_point_to_facet(p,point(vf.first),
															point(vf.second),point(vf.third));
						if(!PInside)
						{
							distance=-distance;
						}
						if (distance<minDistance)
						{
							minDistance=distance;
							facetMinDistance=tmpFacet;
						}
					}
					else if(cosDihedral>cosMinDihedral)
					{
						cosMinDihedral=cosDihedral;
						facetMinDihedral=tmpFacet;
					}
					//判断是否含有大三角面
					std::pair<Vertex_Idtype,Vertex_Idtype> verticesEdge=make_vertex_pair(EdgeInFacet(tmpFacet,idVertexOfFacet));
					double d0=NumericComputation<double>::SquareDistance(point(verticesEdge.first),point(verticesEdge.second));
					double d1=NumericComputation<double>::SquareDistance(point(vertexHandled),point(verticesEdge.second));
					if (d0>maxLengthEdge)
					{
						maxLengthEdge=d0;
						facetLarge=tmpFacet;
					}
					else if(d1>maxLengthEdge)
					{
						maxLengthEdge=d1;
						facetLarge=tmpFacet;
					}
				}
					
				
				Vertex_Idtype vertexLeftEdge=vertex_facet(tmpFacet,idVertexOfEdge.first);
				int iInNN=rank_in_list(vertexLeftEdge,NearestVertices.begin(),NearestVertices.end());
				if(iInNN!=-1&&markForVertices[iInNN]==0)
				{
					markForVertices[iInNN]=2;
					facetNVertices[iInNN]=tmpFacet;
				}
				tmpFacet=neighbor_surfacet(tmpFacet,idVertexOfEdge.first);
			}
			if (markForVertices[rankVertexHandled]==2)
			{
				markForVertices[rankVertexHandled]=-1;
			}
		
		}
		//检查下一个邻近点的incident surface facet
		testForAll=false;
		auto itN=NearestVertices.begin(); 
		for (int i = 0; i < K_NN; i++)
		{
			if (*itN<=0)
			{
				markForVertices[i]=-1;
				itN++;
				continue;
			}
			else if (markForVertices[i]==2)
			{
				testForAll=true;
				rankVertexHandled=i;
				vertexHandled=*itN;
				facetHandledInit=facetNVertices[i];
				break;
			}
			itN++;
		}
		itN=NearestVertices.begin();
		if (!testForAll)
		{
			for (int i = 0; i < K_NN; i++)
			{
				if (markForVertices[i]==0)
				{
					std::vector<Facet> incidentFacets;
					incident_surface_facet(*itN,std::back_inserter(incidentFacets));
					if (incidentFacets.empty())
					{
						markForVertices[i]=-1;
						continue;
					}
					else
					{
						testForAll=true;
						markForVertices[i]=2;
						vertexHandled=*itN;
						facetNVertices[i]=*incidentFacets.begin();
						facetHandledInit=*incidentFacets.begin();
					}
				}
				itN++;
			}
			
		}

	}
	//搜索大三角面片中非K-NN的点的incident surface facet
	std::vector<Vertex_Idtype> VertexNearSparse;
	//稀疏区域内最小二面角的余弦值
	double cosMinDihedralSparse=-1;
	//稀疏区域最小二面角对应的三角面片
	Facet facetMinDihedralSparse(-1,-1);
	//稀疏区域二面角均为锐角且到P距离最近的三角面片的距离
	double minDistanceSparse=0;
	//稀疏区域二面角均为锐角且到P距离最近的三角面片
	Facet facetMinDistanceSparse(-1,-1);

	for (int i = 0; i < 3; i++)
	{
		Vertex_Idtype vLargeFacet=vertex_facet(facetLarge,i);
		int rankVL=rank_in_list(vLargeFacet,NearestVertices.begin(),NearestVertices.end());
		if (rankVL==-1)
		{
			VertexNearSparse.push_back(vLargeFacet);
		}
	}
	for (auto itF=VertexNearSparse.begin();itF!=VertexNearSparse.end();itF++)
	{
		bool markForCircleBegin=true;
		Facet tmpFacet=facetLarge;
		vertexHandled=*itF;
		while (tmpFacet!=facetLarge||markForCircleBegin)
		{
			double cosDihedral;
			double distance;
			markForCircleBegin=false;
			//将tmpFacet绕vertexHandled旋转时的当前facet
			Vertex_triple vf=make_vertex_triple(tmpFacet);
			Indextype idVertexOfFacet=vertex_index_facet(tmpFacet,vertexHandled);
			std::pair<Indextype,Indextype> idVertexOfEdge=make_edge_index(idVertexOfFacet);

			
			Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), p);
			if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
			{
				
				if(!PInside)
				{
					cosDihedral=cos_min_dihedral_point_to_facet(p,mirror_facet(tmpFacet));
				}
				else
				{
					cosDihedral=cos_min_dihedral_point_to_facet(p,tmpFacet);
				}
				if(cosDihedral>-0.0)
				{
					distance =GeometricTraits<T>::distance_point_to_facet(p,point(vf.first),
														point(vf.second),point(vf.third));
					if(!PInside)
					{
						distance=-distance;
					}
					if (distance<minDistance)
					{
						minDistance=distance;
						facetMinDistanceSparse=tmpFacet;
					}
				}
				else if(cosDihedral>cosMinDihedral)
				{
					cosMinDihedral=cosDihedral;
					facetMinDihedralSparse=tmpFacet;
				}
			}			
			Vertex_Idtype vertexLeftEdge=vertex_facet(tmpFacet,idVertexOfEdge.first);
		
			tmpFacet=neighbor_surfacet(tmpFacet,idVertexOfEdge.first);
		}
	}

	if (facetMinDistance.first!=-1)
	{
		double mminDistance=minDistance;
		NearestFacet=facetMinDistance;
		if (facetMinDistanceSparse.first!=-1)
		{
			if (minDistanceSparse<mminDistance)
			{
				mminDistance=minDistanceSparse;
				NearestFacet=facetMinDistanceSparse;
			}
		}
		if (InitFacetConflict.second>-0.1)
		{
			if (InitFacetConflict.third<mminDistance)
			{
				mminDistance=InitFacetConflict.third;
				NearestFacet=InitFacetConflict.first;
			}
		}
	}
	else
	{
		double mminDistance=shortest_edge_to_facet(p,facetMinDihedral);
		NearestFacet=facetMinDihedral;
		if (facetMinDistanceSparse.first!=-1)
		{
			if (minDistanceSparse<mminDistance)
			{
				mminDistance=minDistanceSparse;
				NearestFacet=facetMinDistanceSparse;
			}
		}
		else if(facetMinDihedralSparse.first!=-1)
		{
			double tmpDis=shortest_edge_to_facet(p,facetMinDihedralSparse);
			if (tmpDis<mminDistance)
			{
				if (cosMinDihedralSparse>cosMinDihedral)
				{
					mminDistance=tmpDis;
					NearestFacet=facetMinDihedralSparse;
				}
			}
		}
		if (InitFacetConflict.second>-0.1)
		{
			if (InitFacetConflict.third<mminDistance)
			{
				mminDistance=InitFacetConflict.third;
				NearestFacet=InitFacetConflict.first;
			}
		}

	}
}

template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
nearest_surface_facet(const T* p,list<Vertex_Idtype> NearestVertices,vector<Facet> InitFacet,bool PInside,Facet& NearestFacet,bool& IsConcave)
{
	IsConcave=false;
	
	std::vector<Facet> tmpSurfacets;
	std::vector<bool> isConcaveNearestAcuteDihedral;
	std::vector<Facet> facetsOfNearestAcuteDihedral;

	Facet facetOfNearestMinDihedral(-1,-1);

	double cosNearestMinDiheral=-1;
	double minNearestDistance=COMPUTATION_DOUBLE_MAX;
	
	//初始化NearestFacet
	NearestFacet=make_pair(-1,-1);
	
	for(auto fit=InitFacet.begin();fit!=InitFacet.end();fit++)
	{
		Facet tmpInitFacet=*fit;	
		//------------------for test----------------//
		//Vertex_triple vtriT0=make_vertex_triple(tmpInitFacet);
		//==================test end================//

		if(visited_for_facet_extractor(tmpInitFacet)==1)
		{	
			std::vector<bool> isConcaveAcuteDihedral;
			std::vector<Facet> facetsOfAcuteDihedral;
								
			Facet facetOfMinDihedral=tmpInitFacet;
			bool isConcaveMinDihedral;
			Facet tmpNearestFacet=tmpInitFacet;
			list<Facet> neighborSurfacets;
			neighborSurfacets.push_back(tmpInitFacet);
			double cosMinDihedral=-1;
			double minDistance=COMPUTATION_DOUBLE_MAX;
			if(!PInside)
			{
				cosMinDihedral=cos_min_dihedral_point_to_facet(p,mirror_facet(tmpInitFacet));
			}
			else
			{
				cosMinDihedral=cos_min_dihedral_point_to_facet(p,tmpInitFacet);
			}
			if(cosMinDihedral>-0.1)
			{
				facetsOfAcuteDihedral.push_back(tmpInitFacet);
			}
			while(!neighborSurfacets.empty())
			{
				Facet tmpSurfacet=*neighborSurfacets.begin();
				tmpSurfacets.push_back(tmpSurfacet);
				mark_facet_visited(tmpSurfacet);
				neighborSurfacets.pop_front();
				//-----------------------for test--------------------------//
				//Vertex_triple vTriT0=make_vertex_triple(tmpSurfacet);
				//=======================test end==========================//
		
				Cell_Idtype c = tmpSurfacet.first;
				int li = tmpSurfacet.second;
				// Look for the other neighbors of c.
				for(int ii=0;ii<3;++ii)	{
			
					Facet fCurNei=neighbor_surfacet(tmpSurfacet,ii);
					Vertex_triple vftT0=make_vertex_triple(fCurNei);
					//==================for test================//
					//bool isVis=(visited_for_facet_extractor(fCurNei)!=0);
					//bool isNei=is_facet_in_neighbor(fCurNei,NearestVertices);
					//==================test end================//
					if((visited_for_facet_extractor(fCurNei)==0)&&is_facet_in_neighbor(fCurNei,NearestVertices))
					{
						neighborSurfacets.push_back(fCurNei);
						
						Vertex_triple vf=make_vertex_triple(fCurNei);						
						Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), p);
						if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
						{
							double cosDihedral;				
					

							if(!PInside)
							{						
								cosDihedral=cos_min_dihedral_point_to_facet(p,mirror_facet(fCurNei));
							}
							else
							{
								cosDihedral=cos_min_dihedral_point_to_facet(p,fCurNei);
							}
							if(cosDihedral>cosMinDihedral)
							{
								facetOfMinDihedral=fCurNei;
								cosMinDihedral=cosDihedral;
							}
							if(cosDihedral>-0.1)
							{
								facetsOfAcuteDihedral.push_back(fCurNei);
							}
					
						}

					}
	
				}

			}
			if(facetsOfAcuteDihedral.size()>1)
			{
				for(auto fit=facetsOfAcuteDihedral.begin();fit!=facetsOfAcuteDihedral.end();fit++)
				{
					Vertex_triple vf=make_vertex_triple(Facet(*fit));
					double distance =GeometricTraits<T>::distance_point_to_facet(p,point(vf.first),
																						point(vf.second),point(vf.third));
					if(!PInside)
					{
						distance=-distance;										
					}
					if(distance<minDistance)
					{
						minDistance=distance;
						tmpNearestFacet=*fit;
						facetsOfNearestAcuteDihedral.push_back(tmpNearestFacet);
						isConcaveNearestAcuteDihedral.push_back(true);
					}
				}
			
			}
			else
			{		
				if(facetOfMinDihedral.first!=-1)
				{
					if(cosMinDihedral>cosNearestMinDiheral)
					{
						tmpNearestFacet=facetOfMinDihedral;
						cosNearestMinDiheral=cosMinDihedral;
						facetOfNearestMinDihedral=tmpNearestFacet;
						if(cosMinDihedral>-0.1)
						{
							facetsOfNearestAcuteDihedral.push_back(tmpNearestFacet);
							isConcaveNearestAcuteDihedral.push_back(false);
						}
					}
				}
			}
		}
	}

	if(!facetsOfNearestAcuteDihedral.empty())
	{
		auto itConcave=isConcaveNearestAcuteDihedral.begin();
		for(auto fit=facetsOfNearestAcuteDihedral.begin();
			fit!=facetsOfNearestAcuteDihedral.end();fit++,itConcave++)
		{
			Vertex_triple vf=make_vertex_triple(Facet(*fit));
			double distance =GeometricTraits<T>::distance_point_to_facet(p,point(vf.first),
														point(vf.second),point(vf.third));
			if(!PInside)
			{
				distance=-distance;										
			}
			if(distance<minNearestDistance)
			{
				minNearestDistance=distance;
				NearestFacet=*fit;				
				IsConcave=*itConcave;
			}
		}
	}
	else
	{
		IsConcave=false;
		NearestFacet=facetOfNearestMinDihedral;
		//if(cosNearestMinDiheral<-0.1&&PInside)
		//	cout<<"此处插入点的投影不在表面上"<<endl;
	}
	
	for(auto fit=tmpSurfacets.begin();fit!=tmpSurfacets.end();fit++)
	{
		clear_facet_visited(*fit);
	}
	
	

}

//插入孤立点V
//V孤立点,SeedFacet种子表面三角面片,PInside为孤立点关于表面的位置
template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
link_isolate_to_surface(Vertex_Idtype V,Facet SeedFacet,bool PInside)
{

	if(SeedFacet.first!=-1)
	{
		std::vector<Facet> conflictFacets;
		std::set<Vertex_Idtype> vertexOnBoundary;
		std::set<Vertex_Idtype> vertexOnConflict;
		std::vector<Facet> newCreateSurfacets;
		std::list<Edge> conflictBoundEdges;
		std::vector<Facet> newCreateBoundFacets;//没用
		vector<Facet> newBoundarySurfacets;
		vector<Facet> newConflictSurfacets;

		//在孤立点情况下,使用Gabriel的规则决定替换面的延展
		//V孤立点;SeedFacet为种子三角面片;PInside为孤立点V关于表面的位置
		//VertexOnBoundary为边界线上的顶点，用于孤立点的检查;
		//VertexOnConflict为删除的表面三角面片中的顶点,用于孤立点的检查;conflictFacets是被替换面;
		//conflictBoundEdges为返回值，被替换表面区域的边界线集合，以Edge形式(面id,0/1/2)表示,并且该面不替换
		//newCreateSurfacets待新生成表面的集合
		bool quit=false;
		find_isolate_conflict_surfacet(V,SeedFacet,PInside,vertexOnBoundary,vertexOnConflict,conflictFacets,conflictBoundEdges,newCreateSurfacets,quit);
		if (quit)
		{
			/*cout<<"[iso] quit point"<<endl;*/
			return false;
		}
		for(auto fit=conflictFacets.begin();fit!=conflictFacets.end();fit++)
		{
			delete_surface_facet(*fit);
			if(PInside)
			{
				label_cell_side((*fit).first,false);
			}
			else
			{
				Cell_Idtype cOpp=neighbor((*fit).first,(*fit).second);
				label_cell_side(cOpp,true);
			}
		}
		/*cout<<"新生成表面: "<<endl;*/
		for(auto fit=newCreateSurfacets.begin();fit!=newCreateSurfacets.end();fit++)
		{
			//-------------------for test-----------------//
			//Vertex_triple vftriT0=make_vertex_triple(*fit);
			//cout<<"("<<vftriT0.first<<","<<vftriT0.second<<","<<vftriT0.third<<")"
			//	<<"/("<<(*fit).first<<","<<(*fit).second<<")"<<endl;
			//===================test end=================//

			create_surface_facet(*fit);			
		}

		update_surface_connection(V,conflictBoundEdges,newCreateBoundFacets,newBoundarySurfacets);

		for(auto itF=newCreateSurfacets.begin();itF!=newCreateSurfacets.end();itF++)
		{
			vector<Facet> newUpdateSur;
			Facet tmpF(-1,-1);
			//迭代膨胀,实际操作为向表面内部移入cell
			//*itF为初始膨胀的表面三角面片;第二个参数(形参VertexIso)为可能存在的孤立点,取-1;newUpdateSur为迭代膨胀过程中新形成的表面三角面片集合
			//第二个-1(形参IdV)为某一顶点;tmpF为新形成的且与IdV关联的表面三角面片,初始(-1,-1);true(形参MarkIso)不处理孤立点
			bool isFlated=iterative_inflate_with_mark(*itF,-1,std::back_inserter(newUpdateSur),-1,tmpF,true);
			//迭代雕刻,实际操作剔除表面内部cell
			//*itF为初始雕刻的表面三角面片;第二个参数(形参VertexIso)为可能存在的孤立点,取-1;newUpdateSur为迭代膨胀过程中新形成的表面三角面片集合
			//第二个-1(形参IdV)为某一顶点;tmpF为新形成的且与IdV关联的表面三角面片,初始(-1,-1);true(形参MarkIso)不处理孤立点
			bool isSculptured=iterative_sculpture_with_mark(*itF,-1,std::back_inserter(newUpdateSur),-1,tmpF,true);
		}

		return true;
	}
	else 
		return false;
}


template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
recursive_find_isolate_conflict_surfacet(Vertex_Idtype V, Facet SeedFacet,bool PInside,
	std::set<Vertex_Idtype>& VertexOnBoundary,set<Vertex_Idtype>& VertexOnConflictFacet,
	vector<Facet>& Conflicts,list<Edge>& ConflictBoundEdges,vector<Facet>& NewCreateSurfacets,int prev_ind3)
{
	const T* p=point(V);
	//----------for test----------//
	//BeginFacet=Facet(-1,-1);
	Facet fT=mirror_facet(SeedFacet);
	Vertex_triple ftrip=make_vertex_triple(SeedFacet);
	Cell_Idtype cT=fT.first;
	bool isHull=has_vertex(cT,0);
	bool isNoGabriel;

	if(PInside)
	{
		isNoGabriel=GeometricTraits<T>::side_of_bounded_sphereC3(point(ftrip.first),point(ftrip.second),point(ftrip.third),p)==ON_BOUNDED_SIDE;
	}
	else
	{
		isNoGabriel=GeometricTraits<T>::side_of_bounded_sphereC3(point(ftrip.second),point(ftrip.first),point(ftrip.third),p)==ON_BOUNDED_SIDE;
	}
	cout<<ftrip.first<<" "<<ftrip.second<<" "<<ftrip.third<<"    "<<"is or not hull_"<<isHull<<"      "<<"is Gabriel_"<<isNoGabriel<<endl;
	
	Conflicts.push_back(SeedFacet);
	mark_facet_visited(SeedFacet);
	Vertex_triple vertexOnFacet=make_vertex_triple(SeedFacet);
	VertexOnConflictFacet.insert(vertexOnFacet.first);
	VertexOnConflictFacet.insert(vertexOnFacet.second);
	VertexOnConflictFacet.insert(vertexOnFacet.third);

	Cell_Idtype c = SeedFacet.first;
	
	int li = SeedFacet.second;
	// Look for the other neighbors of c.
	for (int ii = 0; ii<4; ++ii) {
		//if (ii == prev_ind3 || neighbor(c,ii) != -1)
		if (ii==prev_ind3||ii==li)
			continue;
		//cnew->vertex(ii)->set_cell(cnew);

		// Indices of the vertices of cnew such that ii,vj1,vj2,li positive.
		Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
		Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
		Cell_Idtype cur = c;
		int zz = ii;
		Cell_Idtype n = neighbor(cur,zz);
		// turn around the oriented edge vj1 vj2
		while (is_label_inside(n)) {
		
			cur = n;
			zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
			n = neighbor(cur,zz);
		}
		// Now n is outside region, cur is inside.
		//-------------------for test Facet(cur,zz)---------------------//
		//bool isSurfceFacet=is_surface(Facet(cur,zz));
		//Vertex_triple vfnext=make_vertex_triple(Facet(cur,zz));
		//Vertex_Idtype vOpp=vertex(cur,zz);
		//======================test end================================//
		int jj1 = vertex_index(n, vj1);
		int jj2 = vertex_index(n, vj2);
		Vertex_Idtype vvv = vertex(n, Triangulation_utils_3::next_around_edge(jj1, jj2));
		Cell_Idtype nnn = neighbor(n, Triangulation_utils_3::next_around_edge(jj2, jj1));
		int zzz = vertex_index(nnn, vvv);

		bool isNoGabrielNew=false;
		if(PInside)
		{
			Vertex_triple ftrip=make_vertex_triple(Facet(cur,zz));
			isNoGabrielNew=GeometricTraits<T>::side_of_bounded_sphereC3(point(ftrip.first),
				point(ftrip.second),point(ftrip.third),p)==ON_BOUNDED_SIDE;
		}

		if(visited_for_facet_extractor(Facet(cur,zz))==0)
		{
			//判断facet是否在Delaunay影响域边界上或者内部
			/*if(((is_in_conflict(cur)&&PInside)||(is_in_conflict(neighbor(cur,zz))&&!PInside))
				&&((!PInside)||(PInside&&isNoGabrielNew)))*/

			//判断此facet的顶点是否为孤立点			
			Cell_Idtype cOpp=neighbor(cur,zz);
			int indOpp=neighbor_index(cOpp,cur);
			
			if((PInside&&vertex(cur,zz)==V)||(!PInside&&vertex(cOpp,indOpp)==V))
			//判断facet是否为boundary facet
			/*if(((is_in_conflict(cur)&&is_on_boundary(neighbor(cur,zz)))||(is_in_conflict(neighbor(cur,zz))&&is_on_boundary(cur)))
				&&((!PInside)||(PInside&&isNoGabrielNew)))*/
			{	
				//更新交界面的边界边队列
				Facet nei=make_pair(cur,zz);
				int i0=-1,i1=-1;
				int numOfList=0;
				bool handled=false;
				auto itEdge=ConflictBoundEdges.begin();
				for(;itEdge!=ConflictBoundEdges.end();itEdge++)
				{
					numOfList++;
					Facet fCurEdge=get_facet_cell_link((*itEdge).first);
					//------------for test-------------//
					Vertex_triple vTriT0=make_vertex_triple(fCurEdge);
					//============test end=============//
					if(fCurEdge==nei)
					{
						
						i0=(*itEdge).second;
						ConflictBoundEdges.erase(itEdge++);
						break;
					}
				}
					//当list链的第一个被删除时，应该检查最后一个是否也应该被删除
				if(numOfList==1)
				{
					auto endEdge=ConflictBoundEdges.back();
					Facet tmpBoundFacet00=get_facet_cell_link(endEdge.first);
					//--------------for test----------//
					Vertex_triple vtriT00=make_vertex_triple(tmpBoundFacet00);
					//==============test end=========//
					if(tmpBoundFacet00==nei)
					{
						handled=true;
						i1=endEdge.second;
						ConflictBoundEdges.pop_back();
						for(int i2=0;i2<3;i2++)
						{
							if(i2!=i0&&i2!=i1)
							{
								EdgeInFacet oppEdge=mirror_edge(EdgeInFacet(nei,i2));
								//--------------for test-----------//
								Vertex_triple vtriT2=make_vertex_triple(oppEdge.first);
								//==============test end==========//
								ConflictBoundEdges.insert(itEdge,turn_edgeinfacet_to_edge(oppEdge));
								break;
							}
						}
					}
					
				}
				//list 中间的插入与删除要处理好
				if(!handled)
				{
					if(itEdge!=ConflictBoundEdges.end()&&get_facet_cell_link((*itEdge).first)==nei)
					{

						i1=(*itEdge).second;
						ConflictBoundEdges.erase(itEdge++);
						for(int i2=0;i2<3;i2++)
						{
							if(i2!=i0&&i2!=i1)
							{
								EdgeInFacet oppEdge=mirror_edge(EdgeInFacet(nei,i2));
								ConflictBoundEdges.insert(itEdge,turn_edgeinfacet_to_edge(oppEdge));
								break;
							}
						}
					}
					else
					{
						auto leftRight=make_edge_index(i0);
						EdgeInFacet oppEdge0=mirror_edge(EdgeInFacet(nei,leftRight.second));
						EdgeInFacet oppEdge1=mirror_edge(EdgeInFacet(nei,leftRight.first));
						//--------------for test-----------//
						Vertex_triple vtriT3=make_vertex_triple(oppEdge0.first);
						Vertex_triple vtriT4=make_vertex_triple(oppEdge1.first);
						//=============test end============//
						ConflictBoundEdges.insert(itEdge,turn_edgeinfacet_to_edge(oppEdge0));
						ConflictBoundEdges.insert(itEdge,turn_edgeinfacet_to_edge(oppEdge1));
					}
				}
				recursive_find_isolate_conflict_surfacet(V,Facet(cur,zz),PInside,VertexOnBoundary,VertexOnConflictFacet,Conflicts,ConflictBoundEdges,NewCreateSurfacets,zzz);		
			}
			else
			{
				if(PInside)
				{
					//create_surface_facet(SeedFacet.first,ii);
					Facet newFacet=mirror_facet(Facet(SeedFacet.first,ii));
					NewCreateSurfacets.push_back(newFacet);
				}
				else
				{
					Cell_Idtype cOpp=neighbor(SeedFacet.first,li);
					Vertex_Idtype vOppEdge=vertex(SeedFacet.first,ii);
					Indextype iiNewSur=vertex_index(cOpp,vOppEdge);
					//create_surface_facet(cOpp,iiNewSur);
					NewCreateSurfacets.push_back(Facet(cOpp,iiNewSur));
				}
			}
		}
	}

}

//由孤立点V最邻近表面三角面片InitFacet求种子三角面片SeedFacet(种子面片关联的两个cell,有一个cell中,该种子面片相对的顶点为孤立点)
//V为孤立点,PInside为孤立点关于表面的位置
template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
seed_facet_of_isolate_vertex(Vertex_Idtype V,Facet InitFacet,bool PInside,Facet& SeedFacet)
{
	const T* p=point(V);
	SeedFacet=Facet(-1,-1);
	if(PInside)
	{
		if(vertex(InitFacet.first,InitFacet.second)==V)
		{
			SeedFacet=InitFacet;
			return true;
		}
		else  //使用sculpture（不可以产生一个外部孤立点）
		{
			//搜索与InitFacet邻域的surface facet，即与InitFacet有公共点的表面三角面片
			//那么这些三角面片依据什么排序？距离/二面角？即inflate/sculpture的顺序是什么？
			//inflate/sculpture有三种方式，即flip，产生一个孤立点，消除一个孤立点
			//InitFacet邻域的凹凸性质？即凹时以距离为约束，凸时以二面角为约束。或许此时可以用nearest_surface_facet（）函数中得到的凹凸性质
			//暂时不严格排序，利用广度优先的顺序，辅助nearest_surface_facet（）得到的凹凸性质
			std::list<Facet> neighborSurfacets;
			std::list<std::pair<Facet,double>> tmpSurfacetsDihedral;
			Vertex_triple vf=make_vertex_triple(InitFacet);
			std::list<Vertex_Idtype> nearestVertices;//(vf.first,vf.second,vf.third);
			nearestVertices.push_back(vf.first);
			nearestVertices.push_back(vf.second);
			nearestVertices.push_back(vf.third);
			double cosDihedral;								
			cosDihedral=cos_min_dihedral_point_to_facet(p,InitFacet);
	
			neighborSurfacets.push_back(InitFacet);
			tmpSurfacetsDihedral.push_back(make_pair(InitFacet,cosDihedral));
			mark_facet_visited(InitFacet);
			//-----------------求InitFacet邻域的表面三角面片（即与InitFacet有公共点的）------------------//
			while(!neighborSurfacets.empty())
			{
				Facet tmpSurfacet=*neighborSurfacets.begin();
				neighborSurfacets.pop_front();
		
				for (int ii = 0; ii<3; ++ii) {
					
					Facet surfacetNei=neighbor_surfacet(tmpSurfacet,ii);
					if(!visited_for_facet_extractor(surfacetNei)&&is_facet_in_neighbor(surfacetNei,nearestVertices))
					{
						Vertex_triple vf=make_vertex_triple(surfacetNei);						
						Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), p);
						//p表面内,方向与面片内法线一致;p表面外,方向与面片外法线一致
						if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
						{
							mark_facet_visited(surfacetNei);
						
							double cosDihedral;								
		
							cosDihedral=cos_min_dihedral_point_to_facet(p,surfacetNei);
							bool isInserted=false;	
							for(auto iNF=tmpSurfacetsDihedral.begin();iNF!=tmpSurfacetsDihedral.end();iNF++)
							{
								if(cosDihedral>=(*iNF).second)
								{
									tmpSurfacetsDihedral.insert(iNF,make_pair(surfacetNei,cosDihedral));
									neighborSurfacets.push_back(surfacetNei);
									isInserted=true;
									break;
								}
							}
							if (!isInserted)
							{
								tmpSurfacetsDihedral.push_back(make_pair(surfacetNei,cosDihedral));
								neighborSurfacets.push_back(surfacetNei);
							}
				
						}

					}
	
				}
			

			}
			for(auto fit=tmpSurfacetsDihedral.begin();fit!=tmpSurfacetsDihedral.end();fit++)
			{
				clear_facet_visited((*fit).first);
			}
			//======================求InitFacet邻域的表面三角面片（即与InitFacet有公共点的）end===================//


			
			bool isSculpture=false;
			//----------------------在InitFacet邻域中寻找种子三角面片，或sculpture-------------------------//
			for(auto fit=tmpSurfacetsDihedral.begin();fit!=tmpSurfacetsDihedral.end();fit++)
			{
				if (isSculpture)
				{
					fit=tmpSurfacetsDihedral.begin();
					isSculpture=false;
				}
				Facet tmpSurfacet=(*fit).first;
				if (!is_surface(tmpSurfacet))
				{
					continue;
				}
				Cell_Idtype c=(tmpSurfacet).first;
				int ii=(tmpSurfacet).second;
				vector<Facet> newSurfacets;
				if(is_surface(tmpSurfacet))
				{
					//-------------------------for test----------------------//
					Vertex_triple vft0=make_vertex_triple(tmpSurfacet);
					//=========================test end======================//
					
					if(vertex(c,ii)==V)
					{
						SeedFacet=tmpSurfacet;					
						return true;
					}
					else //sculpture
					{
						int fIndv[4];
						int numOfSurface=0;
						int numOfNSurface=0;
						for(int i=0;i<4;i++)
						{
							if(is_surface(Facet(c,i)))
							{
								
								fIndv[numOfSurface]=i;
								numOfSurface++;
							}
							else
							{
								numOfNSurface++;
								fIndv[4-numOfNSurface]=i;
							}
						}

						bool newEdgeShorter=true;
						if (numOfSurface==2)
						{
							double distanceNewEdge=NumericComputation<T>::SquareDistance(point(vertex(c,fIndv[0])),point(vertex(c,fIndv[1])));
							double distanceOldEdge=NumericComputation<T>::SquareDistance(point(vertex(c,fIndv[2])),point(vertex(c,fIndv[3])));
							if (distanceNewEdge>distanceOldEdge)
							{
								newEdgeShorter=false;
							}
						}

						Facet tmpNearFacet(-1,-1);
						bool isoPInside;
						//当numOfSurface=2时要判断inflate/sculpture是否会产生奇异边
						bool isSingularityEdge=false;
						if(numOfNSurface==2)
						{
							isSingularityEdge=is_edge_in_surface(c,fIndv[2],fIndv[3]);

						}

						if((numOfSurface==3)||(numOfSurface==2&&!isSingularityEdge&&newEdgeShorter))//可以优化,此处用flip或者消除一个除V以外的孤立点
						{
							
							sculpture(c,numOfSurface,fIndv,newSurfacets);
							for(auto itNS=newSurfacets.begin();itNS!=newSurfacets.end();itNS++)
							{
								Vertex_triple vf=make_vertex_triple(*itNS);						
								Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), p);
								if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
								{
									double cosDihedral;								
		
									cosDihedral=cos_min_dihedral_point_to_facet(p,*itNS);
									bool isInserted=false;
									for(auto iNF=tmpSurfacetsDihedral.begin();iNF!=tmpSurfacetsDihedral.end();iNF++)
									{
										if(cosDihedral>=(*iNF).second)
										{
											tmpSurfacetsDihedral.insert(iNF,make_pair(*itNS,cosDihedral));
											
											isInserted=true;
											break;
										}
									}
									if (!isInserted)
									{
										tmpSurfacetsDihedral.push_back(make_pair(*itNS,cosDihedral));
										neighborSurfacets.push_back(*itNS);
									}
								}
							}
						
							fit=tmpSurfacetsDihedral.begin();
							isSculpture=true;

						}
					}
				}
			}
			return false;
			//========================在InitFacet邻域中寻找种子三角面片，或sculpture===========================//
		}
		
	}
	else
	{
		Facet fOpp=mirror_facet(InitFacet);
		Cell_Idtype cOpp=fOpp.first;
		int iiOpp=fOpp.second;
		if(vertex(cOpp,iiOpp)==V)
		{
			SeedFacet=InitFacet;	
			return true;
		}
		else  //使用inflate（可能会加入一个内部孤立点）
		{
			//搜索与InitFacet邻域的surface facet，即与InitFacet有公共点的表面三角面片
			//那么这些三角面片依据什么排序？距离/二面角？即inflate/sculpture的顺序是什么？
			//InitFacet邻域的凹凸性质？即凹时以距离为约束，凸时以二面角为约束。或许此时可以用nearest_surface_facet（）函数中得到的凹凸性质
			//暂时不严格排序，利用广度优先的顺序，辅助nearest_surface_facet（）得到的凹凸性质

			std::list<Facet> neighborSurfacets;
			std::list<std::pair<Facet,double>> tmpSurfacetsDihedral;
			Vertex_triple vf=make_vertex_triple(InitFacet);
			std::list<Vertex_Idtype> nearestVertices;//(vf.first,vf.second,vf.third);
			nearestVertices.push_back(vf.first);
			nearestVertices.push_back(vf.second);
			nearestVertices.push_back(vf.third);
			double cosDihedral;								
			cosDihedral=cos_min_dihedral_point_to_facet(p,InitFacet);
	
			neighborSurfacets.push_back(InitFacet);
			tmpSurfacetsDihedral.push_back(make_pair(InitFacet,cosDihedral));
			mark_facet_visited(InitFacet);

			while(!neighborSurfacets.empty())
			{
				Facet tmpSurfacet=*neighborSurfacets.begin();
				neighborSurfacets.pop_front();
			
				for (int ii = 0; ii<3; ++ii) {
					
					Facet surfacetNei=neighbor_surfacet(tmpSurfacet,ii);

					if(!visited_for_facet_extractor(surfacetNei)&&is_facet_in_neighbor(surfacetNei,nearestVertices))
					{
						Vertex_triple vf=make_vertex_triple(surfacetNei);						
						Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), p);
						//p表面内,方向与面片内法线一致;p表面外,方向与面片外法线一致
						if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
						{
							mark_facet_visited(surfacetNei);
	
							double cosDihedral;								
						
							cosDihedral=cos_min_dihedral_point_to_facet(p,mirror_facet(surfacetNei));
							bool isInserted=false;	
							for(auto iNF=tmpSurfacetsDihedral.begin();iNF!=tmpSurfacetsDihedral.end();iNF++)
							{
								if(cosDihedral>=(*iNF).second)
								{
									tmpSurfacetsDihedral.insert(iNF,make_pair(surfacetNei,cosDihedral));
									neighborSurfacets.push_back(surfacetNei);
									isInserted=true;
									break;
								}
							}
							if (!isInserted)
							{
								tmpSurfacetsDihedral.push_back(make_pair(surfacetNei,cosDihedral));
								neighborSurfacets.push_back(surfacetNei);
							}	
					
					
						}

					}
	
				}
			

			}
			for(auto fit=tmpSurfacetsDihedral.begin();fit!=tmpSurfacetsDihedral.end();fit++)
			{
				clear_facet_visited((*fit).first);
			}

			bool isInflate=false;
			for(auto fit=tmpSurfacetsDihedral.begin();fit!=tmpSurfacetsDihedral.end();fit++)
			{
				if (isInflate)
				{
					fit=tmpSurfacetsDihedral.begin();
					isInflate=false;
				}
				Facet tmpSurfacet=(*fit).first;
				if (!is_surface(tmpSurfacet))
				{
					continue;
				}
				vector<Facet> newSurfacets;
				Vertex_triple vtriF=make_vertex_triple(tmpSurfacet);
				if(is_surface(tmpSurfacet))
				{
					Cell_Idtype c=(tmpSurfacet).first;
					int ii=(tmpSurfacet).second;
					Cell_Idtype cOpp=neighbor(c,ii);
					int iiOpp=neighbor_index(cOpp,c);
					if(vertex(cOpp,iiOpp)==V)
					{
						SeedFacet=tmpSurfacet;
						
						return true;
					}
					else //inflate
					{
						int fIndv[4];
						int numOfSurface=0;
						int numOfNSurface=0;
						for(int i=0;i<4;i++)
						{
							if(is_surface(mirror_facet(Facet(cOpp,i))))
							{
								
								fIndv[numOfSurface]=i;
								numOfSurface++;
							}
							else
							{
								numOfNSurface++;
								fIndv[4-numOfNSurface]=i;
							}
						}

						Facet tmpNearFacet(-1,-1);
						bool isoPInside;
						//当numOfSurface=2时要判断inflate/sculpture是否会产生奇异边
						bool isSingularityEdge=false;
						if(numOfNSurface==2)
						{
							isSingularityEdge=is_edge_in_surface(cOpp,fIndv[2],fIndv[3]);

						}
						//-------------------inflate之后不会有infinite点在表面上---------------------//
					    if((numOfSurface==3)||(numOfSurface==2&&!isSingularityEdge))
						{
							
							inflate(cOpp,numOfSurface,fIndv,newSurfacets);
							
							for(auto itNS=newSurfacets.begin();itNS!=newSurfacets.end();itNS++)
							{
								Vertex_triple vf=make_vertex_triple(*itNS);						
								Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), p);
								if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
								{
									double cosDihedral;								
		
									cosDihedral=cos_min_dihedral_point_to_facet(p,*itNS);
									bool isInserted=false;
									for(auto iNF=tmpSurfacetsDihedral.begin();iNF!=tmpSurfacetsDihedral.end();iNF++)
									{
										if(cosDihedral>=(*iNF).second)
										{
											tmpSurfacetsDihedral.insert(iNF,make_pair(*itNS,cosDihedral));
											
											isInserted=true;
											break;
										}
									}
									if (!isInserted)
									{
										tmpSurfacetsDihedral.push_back(make_pair(*itNS,cosDihedral));
										neighborSurfacets.push_back(*itNS);
									}
								}
							}
							fit=tmpSurfacetsDihedral.begin();
							isInflate=true;
						
						}
					}
				}
			}

			return false;
		}
	}
}


//V待插入顶点,ConflictBound为被替换面边界线，以Edge(面id,0/1/2)表示（该面不是替换面）
//NewBoundFacet为新形成的非替换面(即非umbrella上面)的表面三角面片;NewBoundarySurfacets(cell,0/1/2/3)表示与新形成表面共替换面(即umbrella上面)边界线的表面;
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
update_surface_connection(Vertex_Idtype V,list<Edge> ConflictBound,vector<Facet> NewBoundFacet,
	vector<Facet>& NewBonudarySurfacets)
{
	//----------------testing--------------------
	/*cout<<"starting update surface connection:"<<endl;*/
	//---------------testing end

	//更新交界面的邻近关系
	if(!ConflictBound.empty())
	{
		//Edge boundEdge=*(ConflictBound.begin());
		//ConflictBound.push_back(boundEdge);
		EdgeInFacet firstBoundEdge=make_pair(Facet(-1,-1),-1);
		EdgeInFacet former=make_pair(Facet(-1,-1),-1);
		for(auto itE=ConflictBound.begin();itE!=ConflictBound.end();itE++)
		{
			//先求出itE对偶边所在的surface facet即与插入点新形成的surface facet
			//绕itE指向的有向边旋转cell 和 facet
			Facet_Idtype idF=(*itE).first;
			
		    //把面id换成(cell,id)
			Facet boundF=get_facet_cell_link(idF);
			NewBonudarySurfacets.push_back(boundF);
			//求面boundF第(*itE).second号顶点实际编号,即该面中边界边相对的顶点
			Vertex_Idtype vii=vertex_facet(boundF,(*itE).second);
			//求vii是四面体boundF.first的第几号顶点
			Indextype ii=vertex_index(boundF.first,vii);
			//-----------------for test-------------//
			Vertex_triple vtriT0=make_vertex_triple(boundF);
			//=================test end============//
			Indextype li=boundF.second;
			Cell_Idtype c=boundF.first;
			//vj1,vj2为边界边
			Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
			Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
			Cell_Idtype cur = c;
			int zz = ii;
			Cell_Idtype n = neighbor(cur,zz);

			while (1)
			{
				if (!is_label_inside(n))
				{
					Indextype indvT=Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
					Vertex_Idtype vT=vertex(n,indvT);
					bool isTure=vT==V;
					if (isTure&&is_surface(Facet(cur,zz)))
					{
						break;
					}
				}
				cur = n;
				zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
				n = neighbor(cur,zz);

			}
			
			//求新插入点V是Facet(cur,zz)的第几个点
			Indextype iiNei=vertex_index_facet(Facet(cur,zz),V);

			EdgeInFacet curE=make_pair(Facet(cur,zz),iiNei);
			set_facet_adjacency(Facet(cur,zz),iiNei,boundF,(*itE).second);

			//----------------for test--------------//
			//Vertex_triple vtriT1=make_vertex_triple(Facet(cur,zz));
			//Vertex_triple vtriT1b=make_vertex_triple(boundF);
			//cout<<"connect conflict facet: ("<<cur<<","<<zz<<")/("<<vtriT1.first<<","<<vtriT1.second<<","<<vtriT1.third<<")-"<<iiNei
			//	<<"-("<<boundF.first<<","<<boundF.second<<")/("<<vtriT1b.first<<","<<vtriT1b.second<<","<<vtriT1b.third
			//	<<")-"<<(*itE).second<<endl;
			//================test end==============//
	

	

			if(former.second!=-1)
			{
				Vertex_pair vEdge0=make_vertex_pair(former);
				Indextype i0=vertex_index_facet(former.first,vEdge0.second);
				Vertex_pair vEdge1=make_vertex_pair(curE);
				Indextype i1=vertex_index_facet(Facet(cur,zz),vEdge1.first);
				set_facet_adjacency(former.first,i0,Facet(cur,zz),i1);
				//----------------------for test-----------------//
				Vertex_triple vtriT2=make_vertex_triple(former.first);
				Vertex_triple vtriT2b=make_vertex_triple(Facet(cur,zz));
				//cout<<"connect conflict facet: ("<<(former.first).first<<","<<(former.first).second<<")/("<<vtriT2.first<<","<<vtriT2.second<<","<<vtriT2.third
				//	<<")-"<<i0<<"-("<<cur<<","<<zz<<")/("<<vtriT2b.first<<","<<vtriT2b.second<<","<<vtriT2b.third
				//	<<")-"<<i1<<endl;
				//======================test end=================//
			}
			else
			{
				firstBoundEdge=curE;
			}
			former=curE;

		}
		//建立最后一个交界面表面与第一个表面的关系
		Vertex_pair vEdge0=make_vertex_pair(former);
		Indextype i0=vertex_index_facet(former.first,vEdge0.second);
		Vertex_pair vEdge1=make_vertex_pair(firstBoundEdge);
		Indextype i1=vertex_index_facet(firstBoundEdge.first,vEdge1.first);

		set_facet_adjacency(former.first,i0,firstBoundEdge.first,i1);

		//----------------------for test-----------------//
		//Vertex_triple vtriT3=make_vertex_triple(former.first);
        //Vertex_triple vtriT4=make_vertex_triple(firstBoundEdge.first);
		//cout<<"connect conflict facet: ("<<(former.first).first<<","<<(former.first).second<<")/("<<vtriT3.first<<","<<vtriT3.second<<","<<vtriT3.third
		//	<<")-"<<i0<<"-("<<(firstBoundEdge.first).first<<","<<(firstBoundEdge.first).second<<")/"<<vtriT4.first<<","<<vtriT4.second<<","<<vtriT4.third
		//	<<")-"<<i1<<endl;
		//======================test end=================//

	}
	

	//更新新形成的boundary surface facet的邻近关系
	for(auto itB=NewBoundFacet.begin();itB!=NewBoundFacet.end();itB++)
	{
		//---------for test---------//
		Vertex_triple vtriF=make_vertex_triple(*itB);
		//=========test end=========//
		//NewConflictSurfacets.push_back(*itB);
		for(int ii=0;ii<4;ii++)
		{
			
			Cell_Idtype c=(*itB).first;
			Indextype li=(*itB).second;

			if(ii==li)
				continue;
			Vertex_Idtype vTemp=vertex(c,ii);
			Indextype i1=vertex_index_facet((*itB),vTemp);
			Facet_Idtype idFB=get_facet_index(*itB);
			/*if(neighbor_surfacet(idFB,i1)!=-1)
				continue;*/
			Vertex_Idtype vj1 = vertex(c, Triangulation_utils_3::next_around_edge(ii, li));
			Vertex_Idtype vj2 = vertex(c, Triangulation_utils_3::next_around_edge(li, ii));
			Cell_Idtype cur = c;
			int zz = ii;
			Cell_Idtype n = neighbor(cur,zz);
			// turn around the oriented edge vj1 vj2
			while (is_label_inside(n)) {
		
				cur = n;
				zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
				n = neighbor(cur,zz);
			}
			
			Indextype indvOpp=Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
			Vertex_Idtype vOpp=vertex(n,indvOpp);
			Indextype i0=vertex_index_facet(Facet(cur,zz),vOpp);
			set_facet_adjacency(Facet(cur,zz),i0,(*itB),i1);

			//----------------------for test---------------------//
			Vertex_triple vtriFN=make_vertex_triple(Facet(cur,zz));
			Vertex_triple vtriFNB=make_vertex_triple((*itB));
			//cout<<"connect boundary facet: ("<<cur<<","<<zz<<")/("<<vtriFN.first<<","<<vtriFN.second<<","<<vtriFN.third<<")-"<<i0
			//	<<"("<<(*itB).first<<","<<(*itB).second<<")/"<<vtriFNB.first<<","<<vtriFNB.second<<","<<vtriFNB.third
			//	<<")-"<<i1<<endl;
			//======================test end=====================//
		}
	}


}

template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
neighbor_surfacet_around_outside(Facet F0,Indextype I0,Facet& F1,Indextype& I1)
{
	Cell_Idtype c=F0.first;
	Indextype li=F0.second;
	
	
	Vertex_Idtype vTemp=vertex_facet(F0,I0);
	Indextype ii=vertex_index(c,vTemp);

	Cell_Idtype cOpp=neighbor(c,li);
	Vertex_Idtype liOpp=neighbor_index(cOpp,c);
	Vertex_Idtype iiOpp=vertex_index(cOpp,vTemp);
	Vertex_Idtype vj1Opp = vertex(cOpp, Triangulation_utils_3::next_around_edge(iiOpp, liOpp));
	Vertex_Idtype vj2Opp = vertex(cOpp, Triangulation_utils_3::next_around_edge(liOpp, iiOpp));
	Cell_Idtype cur = cOpp;
	int zz = iiOpp;
	Cell_Idtype n = neighbor(cur,zz);
	// turn around the oriented edge vj1 vj2
	while (!is_label_inside(n)) {
		
		cur = n;
		zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1Opp), vertex_index(n, vj2Opp));
		n = neighbor(cur,zz);
	}
	Indextype indvN=Triangulation_utils_3::next_around_edge(vertex_index(n, vj1Opp), vertex_index(n, vj2Opp));
	Indextype indvF=Triangulation_utils_3::next_around_edge(vertex_index(n, vj2Opp), vertex_index(n, vj1Opp));
	Vertex_Idtype vN=vertex(n,indvN);
	Indextype i0=vertex_index_facet(Facet(n,indvF),vN);

	//set_facet_adjacency(Facet(n,indvF),i0,(*itB),i1);
	F1=Facet(n,indvF);
	I1=i0;
		
	
}

//雕刻,把表面内的C标记为表面外,NumOfSurface为C中表面个数
//FIndv为C中前数NumOfSurface个存储的为表面的相对索引,并且FIndv的元素个数一定为4
//NewCreateSurfacets雕刻后新生成表面
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
sculpture(Cell_Idtype C,int NumOfSurface,int* FIndv,vector<Facet>& NewCreateSurfacets)
{
	//存储原表面的相邻关系,
	//存储方式:删除面0号点的邻接面id,删除面相对其邻接面的索引;删除面1号点的邻接面id,删除面相对其邻接面的索引;删除面2号点的邻接面id,删除面相对其邻接面的索引...
	//一组六个,至多三组
	Facet_Idtype SurfacetNeighbor[18];

	//--------------------for test----------------//
	/*cout<<"sculpture cell "<<C<<":"<<is_label_inside(C)<<"->"<<false<<endl;*/
	//====================test end===============//

	label_cell_side(C,false);
	for(int i=0;i<4;i++)
	{
		if(i<NumOfSurface)
		{
			
			//保存原表面的相邻关系
			//存储方式:被删除面0号点的邻接面id,被删除面相对其邻接面的索引;
			//         被删除面1号点的邻接面id,被删除面相对其邻接面的索引;
			//         被删除面2号点的邻接面id,被删除面相对其邻接面的索引...
			//一组六个,至多三组
			Facet_Idtype idF=get_facet_index(Facet(C,FIndv[i]));
			SurfacetNeighbor[i*6]=neighbor_surfacet(idF,0);
			SurfacetNeighbor[i*6+1]=neighbor_index_facet(SurfacetNeighbor[i*6],idF);
			SurfacetNeighbor[i*6+2]=neighbor_surfacet(idF,1);
			SurfacetNeighbor[i*6+3]=neighbor_index_facet(SurfacetNeighbor[i*6+2],idF);
			SurfacetNeighbor[i*6+4]=neighbor_surfacet(idF,2);
			SurfacetNeighbor[i*6+5]=neighbor_index_facet(SurfacetNeighbor[i*6+4],idF);

			//-------------------for test-------------//
			//Vertex_triple vft1=make_vertex_triple(Facet(C,FIndv[i]));	
			//Facet_Idtype fid=CellSurface.get_tuple_element(C,FIndv[i]);
			//cout<<"sculpture delete facet:"<<fid<<"("<<vft1.first<<","<<vft1.second<<","<<vft1.third<<")"
			//	<<"/("<<C<<","<<FIndv[i]<<")["<<C<<"("
			//	<<vertex(C, 0)<<","<<vertex(C, 1)<<","<<vertex(C, 2)<<","<<vertex(C, 3)<<")]"<<endl;
			//==================test end==============//

			//删除C中原NumOfSurface个表面
			delete_surface_facet(Facet(C,FIndv[i]));

			if(NumOfSurface==1)
			{
				remove_from_isolate_vertices(vertex(C,FIndv[i]));
			}

		}
		else
		{
			Facet tmpF=mirror_facet(Facet(C,FIndv[i]));
			//把C中原来不是表面的mirror_facet变表面
			create_surface_facet(tmpF);
			NewCreateSurfacets.push_back(tmpF);
			if(NumOfSurface==3)
			{
				insert_to_isolate_vertices(vertex(C,FIndv[i]));
			}
			//-------------------for test-------------//
			//Vertex_triple vft1=make_vertex_triple(tmpF);
			//Cell_Idtype tmpc=tmpF.first;
			//Facet_Idtype fid=CellSurface.get_tuple_element(tmpc,tmpF.second);
			//cout<<"sculpture create facet:"<<fid<<"("<<vft1.first<<","<<vft1.second<<","<<vft1.third<<")"
			//	<<"/("<<tmpF.first<<","<<tmpF.second<<")["<<tmpc<<"("
			//	<<vertex(tmpc, 0)<<","<<vertex(tmpc, 1)<<","<<vertex(tmpc, 2)<<","<<vertex(tmpc, 3)<<")]"<<endl;
			//==================test end==============//
									
		}
	}

	//C中删一个表面,生成三个mirror face表面
	if(NumOfSurface==1)
	{
		//依次更新三个新成面与被删除面的邻接面的邻接关系
		for(int i=0;i<3;i++)
		{
			Indextype iindv=FIndv[1+i];
			//fCur是相对iindv新生成表面
			Facet fCur=mirror_facet(Facet(C,iindv));
			Facet_Idtype idFCur=get_facet_index(fCur);
			Facet_Idtype i0=vertex_index_facet(fCur,vertex(C,FIndv[0]));		
			Indextype iTemp=Triangulation_utils_3::cell_index_in_triple_index(iindv,FIndv[0]);
			Facet_Idtype idFNei=SurfacetNeighbor[2*iTemp];
			Indextype i1=SurfacetNeighbor[2*iTemp+1];

			//----------------------------for test----------------------//
			Vertex_Idtype vT0=vertex(C,iindv);
			Vertex_Idtype vT1=vertex(C,FIndv[0]);
			Vertex_triple vTriT0=make_vertex_triple(fCur);
			Vertex_triple vTriT1=make_vertex_triple(get_facet_cell_link(idFNei));
			//===========================test end======================//
			
			set_facet_adjacency(idFCur,i0,idFNei,i1);			
		}
		//依次更新三个新生成面彼此邻接关系
		for(int i=0;i<3;i++)
		{
			for(int j=i+1;j<3;j++)
			{
				Indextype iindv0=FIndv[i+1];
				Indextype iindv1=FIndv[j+1];
				Facet f1=mirror_facet(Facet(C,iindv0));
				Facet f2=mirror_facet(Facet(C,iindv1));
				Indextype i0=vertex_index_facet(f1,vertex(C,iindv1));
				Indextype i1=vertex_index_facet(f2,vertex(C,iindv0));
				set_facet_adjacency(f1,i0,f2,i1);
			}
		}
	}

	//C中删两个表面,生成两个mirror face表面
	if(NumOfSurface==2)
	{
		for(int i=0;i<2;i++)
		{
			//fCur是新生成表面
			Facet fCur=mirror_facet(Facet(C,FIndv[2+i]));
			Facet_Idtype idFCur=get_facet_index(fCur);
			for(int j=0;j<2;j++)
			{
				//fDel是被删除表面
				Facet fDel=make_pair(C,FIndv[j]);
				Indextype i0=vertex_index_facet(fCur,vertex(C,FIndv[j]));
				Indextype iTemp=Triangulation_utils_3::cell_index_in_triple_index(FIndv[2+i], FIndv[j]);
				Facet_Idtype idFNei=SurfacetNeighbor[6*j+2*iTemp];
				Indextype i1=SurfacetNeighbor[6*j+2*iTemp+1];
				set_facet_adjacency(idFCur,i0,idFNei,i1);

			}
		}
		//更新两个新生成表面邻接关系
		Facet fNew0=mirror_facet(Facet(C,FIndv[2]));
		Facet fNew1=mirror_facet(Facet(C,FIndv[3]));
		Indextype i0=vertex_index_facet(fNew0,vertex(C,FIndv[3]));
		Indextype i1=vertex_index_facet(fNew1,vertex(C,FIndv[2]));
		set_facet_adjacency(fNew0,i0,fNew1,i1);
	}

	//C中删三个表面,生成一个mirror face表面
	if(NumOfSurface==3)
	{
		Facet fNew=mirror_facet(Facet(C,FIndv[3]));
		Facet_Idtype idFNew=get_facet_index(fNew);
		//依次更新该新成面与被删除面的邻接面的邻接关系
		for(int i=0;i<3;i++)
		{
			Indextype i0=vertex_index_facet(fNew,vertex(C,FIndv[i]));
			int iTemp=Triangulation_utils_3::cell_index_in_triple_index(FIndv[3], FIndv[i]);
			Facet_Idtype idFNei=SurfacetNeighbor[6*i+2*iTemp];
			Indextype i1=SurfacetNeighbor[6*i+2*iTemp+1];
			set_facet_adjacency(idFNew,i0,idFNei,i1);
		}
	}

}


//膨胀,把表面外的C标记为表面内,NumOfSurface为C中mirror face是表面的个数
//FIndv为C中前数NumOfSurface个存储的为表面的相对索引,并且FIndv的元素个数一定为4
//NewCreateSurfacets膨胀后新生成表面
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
inflate(Cell_Idtype C,int NumOfSurface,int* FIndv,vector<Facet>& NewCreateSurfacets)
{
	//存储原表面的相邻关系,
	//存储方式:删除面(其在C中的mirror face 0号点)相对邻接面id,删除面相对其邻接面的索引;
	//         删除面(其在C中的mirror face 1号点)相对邻接面id,删除面相对其邻接面的索引;
	//         删除面(其在C中的mirror face 2号点)相对邻接面id,删除面相对其邻接面的索引...
	//一组六个,共NumOfSurface组
	Facet_Idtype* SurfacetNeighbor=new Facet_Idtype[NumOfSurface*6];
	//-------------for test----------//
	//Facet_Idtype SurfacetNeighbor[18];
	/*cout<<"inflate cell "<<C<<":"<<is_label_inside(C)<<"->"<<true<<endl;*/
	//============test end==========//
	label_cell_side(C,true);
	for(int i=0;i<4;i++)
	{
		if(i<NumOfSurface)
		{
			Facet tmpF=mirror_facet(Facet(C,FIndv[i]));
			
			//保存原表面的相邻关系
			//存储方式:删除面(其在C中的mirror face 0号点)相对邻接面id,删除面相对其邻接面的索引;
			//         删除面(其在C中的mirror face 1号点)相对邻接面id,删除面相对其邻接面的索引;
			//         删除面(其在C中的mirror face 2号点)相对邻接面id,删除面相对其邻接面的索引...
			//一组六个,共NumOfSurface组			
			for(int j=0;j<3;j++)
			{
				Facet_Idtype idF=get_facet_index(tmpF);
				Indextype iTri=vertex_index_facet(tmpF,vertex_facet(Facet(C,FIndv[i]),j));
				SurfacetNeighbor[i*6+2*j]=neighbor_surfacet(idF,iTri);
				SurfacetNeighbor[i*6+2*j+1]=neighbor_index_facet(SurfacetNeighbor[i*6+2*j],idF);
			}

			//-------------------for test-------------//
			//Vertex_triple vft1=make_vertex_triple(tmpF);
			//Cell_Idtype tmpc=tmpF.first;
			//Facet_Idtype fid=CellSurface.get_tuple_element(tmpc,tmpF.second);
			//cout<<"inflate delete facet:"<<fid<<"("<<vft1.first<<","<<vft1.second<<","<<vft1.third<<")"
			//	<<"/("<<tmpF.first<<","<<tmpF.second<<")["<<tmpc<<"("
			//	<<vertex(tmpc, 0)<<","<<vertex(tmpc, 1)<<","<<vertex(tmpc, 2)<<","<<vertex(tmpc, 3)<<")]"<<endl;
			//==================test end==============//

			delete_surface_facet(tmpF);
			if(NumOfSurface==1)
			{
				remove_from_isolate_vertices(vertex(C,FIndv[i]));
			}

		}
		else
		{
			//把C中原来不是表面的面片变表面
			create_surface_facet(Facet(C,FIndv[i]));
			NewCreateSurfacets.push_back(Facet(C,FIndv[i]));
			if(NumOfSurface==3)
			{
				insert_to_isolate_vertices(vertex(C,FIndv[i]));
			}
			//-------------------for test-------------//
			Vertex_triple vft1=make_vertex_triple(Facet(C,FIndv[i]));	
			Facet_Idtype fid=CellSurface.get_tuple_element(C,FIndv[i]);
			//cout<<"inflate create facet:"<<fid<<"("<<vft1.first<<","<<vft1.second<<","<<vft1.third<<")"
			//	<<"/("<<C<<","<<FIndv[i]<<")["<<C<<"("
			//	<<vertex(C, 0)<<","<<vertex(C, 1)<<","<<vertex(C, 2)<<","<<vertex(C, 3)<<")]"<<endl;
			//==================test end==============//
		}
	}

	//C中删一个mirror face表面,生成三个表面
	if(NumOfSurface==1)
	{
		
		for(int i=0;i<3;i++)
		{
			Indextype iindv=FIndv[1+i];
			Indextype i0= Triangulation_utils_3::cell_index_in_triple_index(FIndv[0], iindv);
				
			Indextype iTemp=Triangulation_utils_3::cell_index_in_triple_index(iindv,FIndv[0]);
			Facet_Idtype idFNei=SurfacetNeighbor[2*iTemp];
			Indextype i1=SurfacetNeighbor[2*iTemp+1];
			Facet_Idtype idFCur=get_facet_index(Facet(C,iindv));
			set_facet_adjacency(idFCur,i0,idFNei,i1);
				
		}
		for(int i=0;i<3;i++)
			for(int j=i+1;j<3;j++)
			{
				Indextype iindv0=FIndv[i+1];
				Indextype iindv1=FIndv[j+1];
				Indextype i0=Triangulation_utils_3::cell_index_in_triple_index(iindv1, iindv0);
				Indextype i1=Triangulation_utils_3::cell_index_in_triple_index(iindv0, iindv1);
				set_facet_adjacency(Facet(C,iindv0),i0,Facet(C,iindv1),i1);
			}
	}

	//C中删两个mirror face表面,生成两个表面
	if(NumOfSurface==2)
	{
		for(int i=0;i<2;i++)
		{
			//fCur是新生成表面
			Facet fCur=Facet(C,FIndv[2+i]);
			Facet_Idtype idFCur=get_facet_index(fCur);
			for(int j=0;j<2;j++)
			{
				//fDel是被删除表面
				Facet fDel=make_pair(C,FIndv[j]);
				Indextype i0=Triangulation_utils_3::cell_index_in_triple_index(FIndv[j], FIndv[2+i]);
				Indextype iTemp=Triangulation_utils_3::cell_index_in_triple_index(FIndv[2+i], FIndv[j]);
				Facet_Idtype idFNei=SurfacetNeighbor[6*j+2*iTemp];
				Indextype i1=SurfacetNeighbor[6*j+2*iTemp+1];
				set_facet_adjacency(idFCur,i0,idFNei,i1);

			}
		}
		//更新两个新生成表面邻接关系
		Facet fNew0=(Facet(C,FIndv[2]));
		Facet fNew1=(Facet(C,FIndv[3]));
		Indextype i0=Triangulation_utils_3::cell_index_in_triple_index(FIndv[3], FIndv[2]);
		Indextype i1=Triangulation_utils_3::cell_index_in_triple_index(FIndv[2], FIndv[3]);
		set_facet_adjacency(fNew0,i0,fNew1,i1);
	}

	//C中删三个mirror face表面,生成一个表面
	if(NumOfSurface==3)
	{
		Facet fNew=(Facet(C,FIndv[3]));
		Facet_Idtype idFNew=get_facet_index(fNew);
		//依次更新该新成面与被删除面的邻接面的邻接关系
		for(int i=0;i<3;i++)
		{
			Indextype i0=Triangulation_utils_3::cell_index_in_triple_index(FIndv[i], FIndv[3]);
			int iTemp=Triangulation_utils_3::cell_index_in_triple_index(FIndv[3], FIndv[i]);
			Facet_Idtype idFNei=SurfacetNeighbor[6*i+2*iTemp];
			Indextype i1=SurfacetNeighbor[6*i+2*iTemp+1];
			set_facet_adjacency(idFNew,i0,idFNei,i1);
		}
	}

	delete []SurfacetNeighbor;


}

//C中,I0-I1相对边是否为表面的边
template<typename T, typename T_INFO>
bool DataStructure<T, T_INFO>::
is_edge_in_surface(Cell_Idtype C,Indextype I0,Indextype I1)
{
			
	bool isCellInside=is_label_inside(C);
	Vertex_Idtype vj1 = vertex(C, Triangulation_utils_3::next_around_edge(I0, I1));
	Vertex_Idtype vj2 = vertex(C, Triangulation_utils_3::next_around_edge(I1, I0));
	Cell_Idtype cur = C;
	int zz = I0;
	Cell_Idtype n = neighbor(cur,zz);
	
	// turn around the oriented edge vj1 vj2
	while (n!=C) {

		if(is_label_inside(n)!=isCellInside)
		{
			return true;
		}
		cur = n;
		zz = Triangulation_utils_3::next_around_edge(vertex_index(n, vj1), vertex_index(n, vj2));
		n = neighbor(cur,zz);		
	}
	return false;
}


//迭代收缩，实际操作剔除表面内部cell
//F为初始收缩表面三角面片;VertexIso为可能存在的孤立点,初始为-1;ItSurf为新形成的表面三角面片集合;
//IdV为某一顶点,初始为-1;Fcur为新形成且与IdV关联的表面三角面片,初始(-1,-1);MarkIso为是否处理孤立点,初始true
template<typename T, typename T_INFO>
template<typename IteratorFacet>
bool DataStructure<T, T_INFO>:: 
iterative_sculpture_with_mark(Facet F,Vertex_Idtype VertexIso,IteratorFacet ItSurF,Vertex_Idtype IdV,Facet& FCur,bool MarkIso)
{
	//----------------for test-------------//
	bool isS=is_surface(F);
	if (!isS)
	{
		//cout<<"error!:iterative_sculpture_with_mark no surface"<<endl;
		return false;
	}
	//================test end============//	
	Cell_Idtype c=F.first;
	int ii=F.second;
	int fIndv[4];//fIndv为cOpp前数numOfSurface个存储的是mirror_facet为表面的相对索引，后数numOfNSurface是非表面的相对索引
	int numOfSurface=0;
	int numOfNSurface=0;//numOfSurface+numOfNSurface=4成立
	vector<Facet> newSurfacets;
	bool isInSphere=false;
	bool containFCur=false;

	for(int i=0;i<4;i++)
	{
		if(is_surface(Facet(c,i)))
		{
								
			fIndv[numOfSurface]=i;
			numOfSurface++;
			Vertex_triple vNSur=make_vertex_triple(Facet(c,i));

			//处理无穷远点
			if (is_infinite(vNSur.first)||is_infinite(vNSur.second)||
				is_infinite(vNSur.third)||is_infinite(vertex(c,i)))
			{
				return false;
			}

			bool isInSphereT=GeometricTraits<T>::side_of_bounded_sphereC3(point(vNSur.first),point(vNSur.second),
										point(vNSur.third),point(vertex(c,i)))==ON_BOUNDED_SIDE;
			if (isInSphereT)
			{
				isInSphere=true;
			}
			if (Facet(c,i)==FCur)
			{
				containFCur=true;
			}
		}
		else
		{
			numOfNSurface++;
			fIndv[4-numOfNSurface]=i;
		}
	}

	//当替换的边小于被替换的边时flip sculpture
	bool newEdgeShorter=true;
	if (numOfSurface==2)
	{
		double distanceNewEdge=NumericComputation<T>::SquareDistance(point(vertex(c,fIndv[0])),point(vertex(c,fIndv[1])));
		double distanceOldEdge=NumericComputation<T>::SquareDistance(point(vertex(c,fIndv[2])),point(vertex(c,fIndv[3])));
		if (distanceNewEdge>distanceOldEdge)
		{
			newEdgeShorter=false;
		}
	}

	Facet tmpNearFacet(-1,-1);
	bool isoPInside;
	//当numOfSurface=2时要判断inflate/sculpture是否会产生奇异边
	bool isSingularityEdge=false;
	if(numOfNSurface==2)
	{
		isSingularityEdge=is_edge_in_surface(c,fIndv[2],fIndv[3]);
	}

	if (containFCur&&IdV==-1)
	{
		return false;
	}
	if (MarkIso&&(numOfSurface==3))
	{
		if (IdV==-1)
		{
			return false;
		}
		else
		{
			FCur=F;
			return false;
		}
	}
	if ((vertex(c,ii)==VertexIso)&&(numOfSurface==1))
	{
		if (IdV==-1)
		{
			return false;
		}
		else
		{
			FCur=F;
			return false;
		}
	}

	if(isInSphere&&(numOfSurface==2&&!isSingularityEdge&&newEdgeShorter))//可以优化,此处用flip或者消除一个除V以外的孤立点
	{
		sculpture(c,numOfSurface,fIndv,newSurfacets);
		if (IdV!=-1)
		{
			int iFNext=-1;
			Facet* newSur=new Facet[numOfNSurface];
			int numOfIncidentSur=0;
			int numOfNIncidentSur=0;
			int iVF=-1;
			for(auto itNS=newSurfacets.begin();itNS!=newSurfacets.end();itNS++)
			{
				*ItSurF++=*itNS;
				int iVF_tmp=vertex_index_facet(*itNS,IdV);
				if (iVF_tmp>=0)
				{
					mark_facet_visited(*itNS);
					newSur[numOfIncidentSur]=*itNS;
					if(numOfIncidentSur==0)
					{
						iVF=iVF_tmp;
					}
					numOfIncidentSur++;
				}
				else
				{
					numOfNIncidentSur++;
					newSur[numOfNSurface-numOfNIncidentSur]=*itNS;
				}
			}
			
			//按逆时针方向查找某一点关联的表面三角面片
			if (numOfIncidentSur==1)
			{
				iFNext=0;
			}
			else if (numOfIncidentSur==2)
			{
				std::pair<int,int> idVEdge=make_edge_index(iVF);
				Facet FNext_tmp=neighbor_surfacet(newSur[0],idVEdge.first);
				iFNext=0;
				if (FNext_tmp==newSur[1])
				{
					iFNext=1;
				}
			}
			if (numOfSurface!=3)
			{
				
				for (int i = 0; i < numOfNSurface; i++)
				{
					Facet FCur_tmp=newSur[iFNext];
					if (i==iFNext)
					{
						continue;
					}
					else
					{
						if (is_surface(newSur[i]))
						{
							iterative_sculpture_with_mark(newSur[i],VertexIso,ItSurF,-1,FCur_tmp,MarkIso);
						}
						
					}
				}
				iterative_sculpture_with_mark(newSur[iFNext],VertexIso,ItSurF,IdV,FCur,MarkIso);
			}
			else
			{
				if (numOfIncidentSur=0)
				{
					FCur=Facet(-1,-1);
				}
				else
				{
					FCur=newSur[iFNext];
				}
			}
			delete [] newSur;
		}//end of if (IdV!=-1)
		else
		{
			if (numOfSurface!=3)
			{
				for(auto itNS=newSurfacets.begin();itNS!=newSurfacets.end();itNS++)
				{
					
					if (is_surface(*itNS))
					{
						*ItSurF++=*itNS;
						iterative_sculpture_with_mark(*itNS,VertexIso,ItSurF,-1,FCur,MarkIso);
					}
				
				}
			}
			
			
		}
		return true;
	}//end of if(isInSphere&&...
	else
	{
		if (IdV!=-1)
		{
			FCur=F;
		
		}
	
		return false;
	}
}


//迭代膨胀,实际操作为向表面内部移入cell
//F为初始膨胀的表面三角面片;VertexIso为可能存在的孤立点,初始-1;ItSurF为迭代膨胀过程中新形成的表面三角面片集合
//IdV为某一顶点,初始-1;FCur为新形成的且与IdV关联的表面三角面片,初始(-1,-1);MarkIso为是否处理孤立点标识,初始true
template<typename T, typename T_INFO>
template<typename IteratorFacet>
bool DataStructure<T, T_INFO>::
iterative_inflate_with_mark(Facet F,Vertex_Idtype VertexIso,IteratorFacet ItSurF,Vertex_Idtype IdV,Facet& FCur,bool MarkIso)
{
	//----------------for test-------------//
	bool isS=is_surface(F);
	if (!isS)
	{
		//--------------------for test--------------//
		//Vertex_triple vtriF=make_vertex_triple(F);
		//====================test end==============//
		//cout<<"error!:iterative_inflate_with_mark no surface"<<endl;
		return false;
	}
	//================test end============//
	Cell_Idtype c=F.first;
	int ii=F.second;
	//取c的第ii个相邻Cell
	Cell_Idtype cOpp=neighbor(c,ii);
	//求c是cOpp的第几个临近Cell
	int iiOpp=neighbor_index(cOpp,c);
	int fIndv[4];//fIndv为cOpp前数numOfSurface个存储的是mirror_facet为表面的相对索引，后数numOfNSurface是非表面的相对索引
	int numOfSurface=0;
	int numOfNSurface=0;//numOfSurface+numOfNSurface=4成立
	vector<Facet> newSurfacets;
	bool isInSphere=false;
	bool containFCur=false;

	int numInsphere=0;
	for(int i=0;i<4;i++)
	{
		if(is_surface(mirror_facet(Facet(cOpp,i))))
		{
								
			fIndv[numOfSurface]=i;
			numOfSurface++;
			if (mirror_facet(Facet(cOpp,i))==FCur)
			{
				containFCur=true;
			}
		}
		else
		{
			numOfNSurface++;
			fIndv[4-numOfNSurface]=i;
			Vertex_triple vNSur=make_vertex_triple(Facet(cOpp,i));

			//---------for test---------//
			//const double* p0=point(vNSur.first);
			//const double* p1=point(vNSur.second);
			//const double* p2=point(vNSur.third);
			//const double* p3=point(vertex(cOpp,i));
			//========test end==========//
			Vertex_Idtype vT = vertex(cOpp, i);

			//处理无穷远点
			if (is_infinite(vNSur.first)||is_infinite(vNSur.second)||
				is_infinite(vNSur.third)||is_infinite(vT))
			{
				return false;
			}

			//side_of_bounded_sphereC3判断点point(vertex(cOpp,i))在point(vNSur.first)，point(vNSur.second)和point(vNSur.third)最大外接球内/外/上
			//ON_BOUNDED_SIDE为表示在后三点的bounded sphere内部;ON_UNBOUNDED_SIDE表示在外部
			Bounded_side relative_position=GeometricTraits<T>::side_of_bounded_sphereC3(point(vNSur.first),point(vNSur.second),
				point(vNSur.third),point(vertex(cOpp,i)));
			bool isInSphereT=(relative_position==ON_BOUNDED_SIDE);
			if (isInSphereT)
			{
				numInsphere++;
			}
		}
	}
	
	if (numInsphere>0)
	{
		isInSphere=true;
	}


	Facet tmpNearFacet(-1,-1);
	bool isoPInside;
	//当numOfSurface=2时要判断inflate/sculpture是否会产生奇异边
	bool isSingularityEdge=false;
	if(numOfNSurface==2)
	{
		isSingularityEdge=is_edge_in_surface(cOpp,fIndv[2],fIndv[3]);

	}
	//-------------------inflate之后不会有infinite点在表面上---------------------//
	//加入Gabriel判断
	if (containFCur&&IdV==-1)
	{
		return false;
	}
	if (MarkIso&&(numOfSurface==3))
	{
		if (IdV==-1)
		{
			return false;
		}
		else
		{
			FCur=F;
			return false;
		}
	}
	if ((vertex(cOpp,iiOpp)==VertexIso)&&(numOfSurface==1))
	{
		if (IdV==-1)
		{
			return false;
		}
		else
		{
			FCur=F;
			return false;
		}
	}

	
	if(!isInSphere&&(numOfSurface==2&&!isSingularityEdge))
	{
		//膨胀,把表面外的cOpp标记为表面内,numOfSurface为cOpp中mirror face是表面的个数
		//fIndv为cOpp中前数numOfSurface个存储的为表面的相对索引,并且fIndv的元素个数一定为4
		//newSurfacets膨胀后新生成表面		
		inflate(cOpp,numOfSurface,fIndv,newSurfacets);
		
		if (IdV!=-1)
		{
			int iFNext=-1;
			Facet* newSur=new Facet[numOfNSurface];
			int numOfIncidentSur=0;
			int numOfNIncidentSur=0;
			int iVF=-1;
			for(auto itNS=newSurfacets.begin();itNS!=newSurfacets.end();itNS++)
			{
				*ItSurF++=*itNS;
				int iVF_tmp=vertex_index_facet(*itNS,IdV);
				if (iVF_tmp>=0)
				{
					mark_facet_visited(*itNS);
					newSur[numOfIncidentSur]=*itNS;
					if(numOfIncidentSur==0)
					{
						iVF=iVF_tmp;
					}
					numOfIncidentSur++;
				}
				else
				{
					numOfNIncidentSur++;
					newSur[numOfNSurface-numOfNIncidentSur]=*itNS;
				}
			}
			
		
			if (numOfIncidentSur==1)
			{
				iFNext=0;
			}
			else if (numOfIncidentSur==2)
			{
				std::pair<int,int> idVEdge=make_edge_index(iVF);
				Facet FNext_tmp=neighbor_surfacet(newSur[0],idVEdge.first);
				iFNext=0;
				if (FNext_tmp==newSur[1])
				{
					iFNext=1;
				}
			}
			if (numOfSurface!=3)
			{
				
				for (int i = 0; i < numOfNSurface; i++)
				{
					Facet FCur_tmp=newSur[iFNext];
					if (i==iFNext)
					{
						continue;
					}
					else
					{
						if (is_surface(newSur[i]))
						{
							iterative_inflate_with_mark(newSur[i],VertexIso,ItSurF,-1,FCur_tmp,MarkIso);
						}
						
					}
				}
				iterative_inflate_with_mark(newSur[iFNext],VertexIso,ItSurF,IdV,FCur,MarkIso);
			}
			else
			{
				if (numOfIncidentSur=0)
				{
					FCur=Facet(-1,-1);
				}
				else
				{
					FCur=newSur[iFNext];
				}
	
			}
			
		}//end of if(IdV!=-1)
		else
		{
			if (numOfSurface!=3)
			{
				for(auto itNS=newSurfacets.begin();itNS!=newSurfacets.end();itNS++)
				{
					if (is_surface(*itNS))
					{
						*ItSurF++=*itNS;
						iterative_inflate_with_mark(*itNS,VertexIso,ItSurF,-1,FCur,MarkIso);
					}
					
				}
			}
	
			
		}
		return true;	
	}//end of if(!isInSphere&&...
	else
	{
		if (IdV!=-1)
		{
			FCur=F;
		}
	
		return false;
	}
}


//插入点p领域内膨胀或收缩，使领域内的表面到采样点集的距离最小。
//pInside是p关于表面的位置，true为表面内部,false为表面外部;VertexIso处理孤立点(已剖分)时候使用,一般为-1;
//NearestVertices为p的K个最近点，InitFacets为与最邻近点关联的初始三角面片;
//ItSurF为输出,作为膨胀和收缩后与KNN关联的所有表面三角面片;
//ItSurSparse为KNN领域内含有较长边的表面三角面片集合，用来处理采样不均匀
template<typename T, typename T_INFO>
template<typename IteratorFacet>
void DataStructure<T, T_INFO>::
update_surfacets_KNN(bool PInside,Vertex_Idtype VertexIso,list<Vertex_Idtype> NearestVertices,vector<Facet> InitFacets,IteratorFacet ItSurF,IteratorFacet ItSurSparse)
{
	//markForVertices用于标记第i个最近邻点是否处理，未处理为0，待处理为2,；已经处理且不是孤立点为1，是孤立点或者不包含可能为最近邻的表面为-1
	int markForVertices[5]={0};
	//用于存储K-NN的incident surfacet，且对于P可见
	Facet facetNVertices[5]={Facet(-1,-1)};
	//存储正在处理的vertex 
	Vertex_Idtype vertexHandled=*NearestVertices.begin();
	int rankVertexHandled=0;
	//用于存储vertexHandled的一个incident surface facet
	Facet facetHandledInit(-1,-1);
	//incident surfacets的最长边的距离
	double maxLengthEdge=0;
	//最长边所在的表面三角面片
	Facet facetLarge(-1,-1);
	if (!InitFacets.empty())
	{
		facetHandledInit=*InitFacets.begin();
	}
	//markForVertices[0]=2;
	//facetNVertices[0]=facetHandledInit;
	if (facetHandledInit.first!=-1)
	{
		markForVertices[0]=1;
		facetNVertices[0]=facetHandledInit;
	}
	else
	{
		markForVertices[0]=-1;
	}
	//存储KNN的incident surfacets但是其中可能有重复
	vector<Facet> tmpSurIncident;
	//表征所有的K-NN是否处理完毕，是为false，否为true
	bool testForAll=true;

	while (testForAll)
	{
		Facet tmpFacet=facetHandledInit;
		
		if (facetHandledInit.first!=-1)
		{
			//存储与KNN中某一点关联的所有表面三角面片，可能有重复
			vector<Facet> tmpIncidentSurfacets;
			while (visited_for_facet_extractor(tmpFacet)==0)
			{
				//---------------for test-----------//
				//Vertex_triple vtriFT0=make_vertex_triple(tmpFacet);
				//==============test end============//
				mark_facet_visited(tmpFacet);
				tmpIncidentSurfacets.push_back(tmpFacet);

				//---------------for test-----------//
				//Vertex_triple vtriFN1=make_vertex_triple(tmpFacet);
				//==============test end============//
				//将tmpFacet绕vertexHandled旋转时的当前facet
				//求点vertexHandled是面tmpFacet的第几个点（0/1/2）
				Indextype idVertexOfFacet=vertex_index_facet(tmpFacet,vertexHandled);
				//idVertexOfEdge是相对于顶点idVertexOfFacet（0/1/2）的边
				std::pair<Indextype,Indextype> idVertexOfEdge=make_edge_index(idVertexOfFacet);
				//求面的相对顶点编号（0/1/2）idVertexOfEdge.first在面fNext中实际顶点编号，
				Vertex_Idtype vertexLeftEdge=vertex_facet(tmpFacet,idVertexOfEdge.first);
				int iInNN=rank_in_list(vertexLeftEdge,NearestVertices.begin(),NearestVertices.end());
				if(iInNN!=-1&&markForVertices[iInNN]==0)
				{
					markForVertices[iInNN]=2;
					facetNVertices[iInNN]=tmpFacet;
				}
			    //求表面tmpFacet的第idVertexOfEdge.first个相邻的表面
				tmpFacet=neighbor_surfacet(tmpFacet,idVertexOfEdge.first);
				//判断是否含有大三角面
				
				//make_vertex_pair把边以EdgeInFacet<面，0/1/2>形式换成<顶点实际编号，顶点实际编号>的形式
				std::pair<Vertex_Idtype,Vertex_Idtype> verticesEdge=make_vertex_pair(EdgeInFacet(tmpFacet,idVertexOfFacet));
				const double d0=NumericComputation<double>::SquareDistance(point(verticesEdge.first),point(verticesEdge.second));
				const double d1=NumericComputation<double>::SquareDistance(point(vertexHandled),point(verticesEdge.first));
				if (d0>maxLengthEdge)
				{
					maxLengthEdge=d0;
					facetLarge=tmpFacet;
				}
				else if(d1>maxLengthEdge)
				{
					maxLengthEdge=d1;
					facetLarge=tmpFacet;
				}
			}
			for(auto itIST=tmpIncidentSurfacets.begin();itIST!=tmpIncidentSurfacets.end();itIST++)
			{
				if (is_surface(*itIST))
				{
					//--------------------for test-------------//
					//Vertex_triple vtriF=make_vertex_triple(*itIST);
					//====================test end=============//
					if (visited_for_facet_extractor(*itIST)!=0)
					{
						clear_facet_visited(*itIST);					
					}
					tmpSurIncident.push_back(*itIST);
				}
				
			}
			markForVertices[rankVertexHandled]=1;
		}
		else
		{
			markForVertices[rankVertexHandled]=-1;
		}
		//检查下一个邻近点的incident surface facet
		testForAll=false;
		auto itN=NearestVertices.begin(); 
		for (int i = 0; i < K_NN; i++)
		{
			if (*itN<=0)
			{
				markForVertices[i]=-1;
				itN++;
				continue;
			}
			else if (markForVertices[i]==2)
			{
				testForAll=true;
				rankVertexHandled=i;
				vertexHandled=*itN;
				facetHandledInit=facetNVertices[i];
				break;
			}
			itN++;
		}
		itN=NearestVertices.begin();
		if (!testForAll)
		{
			for (int i = 0; i < K_NN; i++)
			{
				if (markForVertices[i]==0)
				{
					std::vector<Facet> incidentFacets;
					incident_surface_facet(*itN,std::back_inserter(incidentFacets));
					if (incidentFacets.empty())
					{
						markForVertices[i]=-1;
						itN++;
						continue;
					}
					else
					{
						testForAll=true;
						markForVertices[i]=2;
						vertexHandled=*itN;
						facetNVertices[i]=*incidentFacets.begin();
						facetHandledInit=*incidentFacets.begin();
					}
				}
				itN++;
			}
			
		}

	}
	vector<Facet> tmpSurfacetsIncident;//tmpSurfacetsIncident为tmpSurIncident去除重复后的关连KNN顶点的面片
	for (auto iTS=tmpSurIncident.begin();iTS!=tmpSurIncident.end();iTS++)
	{
		if (is_surface(*iTS))
		{
			if (visited_for_facet_extractor(*iTS)==0)
			{
				//--------------------for test-------------//
				//Vertex_triple vtriF=make_vertex_triple(*iTS);
				//====================test end=============//
				tmpSurfacetsIncident.push_back(*iTS);
				mark_facet_visited(*iTS);
			}
		}
	}
	for(auto iTS=tmpSurfacetsIncident.begin();iTS!=tmpSurfacetsIncident.end();iTS++)
	{
		*ItSurF++=*iTS;
		clear_facet_visited(*iTS);
	}

	if (tmpSurIncident.empty()&&tmpSurfacetsIncident.empty())
	{
		return;
	}
	//搜索大三角面片中非K-NN的点的incident surface facet
	std::vector<Vertex_Idtype> VertexNearSparse;
	
	if (is_surface(facetLarge))
	{
		for (int i = 0; i < 3; i++)
		{
			Vertex_Idtype vLargeFacet=vertex_facet(facetLarge,i);
			int rankVL=rank_in_list(vLargeFacet,NearestVertices.begin(),NearestVertices.end());
			if (rankVL==-1)
			{
				VertexNearSparse.push_back(vLargeFacet);
			}
		}
	}
	
	vector<Facet> tmpSurS;
	for (auto itF=VertexNearSparse.begin();itF!=VertexNearSparse.end();itF++)
	{
		
		vertexHandled=*itF;
		Facet tmpFacet=facetLarge;
		vector<Facet> tmpIncidentSparse;
		while (visited_for_facet_extractor(tmpFacet)==0)
		{
		
			mark_facet_visited(tmpFacet);
			tmpIncidentSparse.push_back(tmpFacet);

			//将tmpFacet绕vertexHandled旋转时的当前facet
			Vertex_triple vf=make_vertex_triple(tmpFacet);
			Indextype idVertexOfFacet=vertex_index_facet(tmpFacet,vertexHandled);
			std::pair<Indextype,Indextype> idVertexOfEdge=make_edge_index(idVertexOfFacet);
				
			Vertex_Idtype vertexLeftEdge=vertex_facet(tmpFacet,idVertexOfEdge.first);
		
			tmpFacet=neighbor_surfacet(tmpFacet,idVertexOfEdge.first);
		}
		for(auto itIST=tmpIncidentSparse.begin();itIST!=tmpIncidentSparse.end();itIST++)
		{
			if (is_surface(*itIST))
			{
				clear_facet_visited(*itIST);
				tmpSurS.push_back(*itIST);
			}
				
		}
	}
	vector<Facet> tmpSurfacetsIncidentSparse;
	for (auto iTS=tmpSurS.begin();iTS!=tmpSurS.end();iTS++)
	{
		if (visited_for_facet_extractor(*iTS)==0)
		{
			tmpSurfacetsIncidentSparse.push_back(*iTS);
			mark_facet_visited(*iTS);
		}
	}
	for(auto iTS=tmpSurfacetsIncidentSparse.begin();iTS!=tmpSurfacetsIncidentSparse.end();iTS++)
	{
		*ItSurSparse++=*iTS;
		clear_facet_visited(*iTS);
	}
}


//求插入点p的最近的三角面片
//SurfacetsInternalConflict是表面面片（该面是表面,该面所在的cell属于影响域,或其mirror facet所在cell属于影响域,满足其一即可）
//SurfacetsIncidentKNN为与KNN关联的备选表面三角面片；SurfacetsIncidentSparse为KNN领域内含有较长边的备选表面三角面片集合
//PInside表示p关于表面的位置；NearestFacet为返回的最邻近表面三角面片；
//SurfacetsInConflict为影响域内的表面三角面片(该面及其mirror facet所在的cell都是conflict)
// VerticesDeletedSurfacet(删除的表面三角面片中的点id,与之关联的表面三角面片被删除的个数)，用于孤立点的检测
template<typename T, typename T_INFO>
void DataStructure<T, T_INFO>::
nearest_surfacet(const T* P,vector<Facet> SurfacetsInternalConflict,vector<Facet> SurfacetsIncidentKNN,vector<Facet> SurfacetsIncidentSparse,
	bool PInside,Facet& NearestFacet,vector<Facet_Idtype>& SurfacetsInConflict,std::map<Vertex_Idtype,size_t>& VerticesDeletedSurfacet)
{
	//step1:求K-NN的incident surface facet中最有可能的nearest surface facet，并标记其中最大的三角面片
	//最小二面角的余弦值
	double cosMinDihedral_KNN=-1;
	//最小二面角对应的三角面片
	Facet facetMinDihedral_KNN(-1,-1);
	//当P与三角面片形成的二面角均为锐角时，P到此面的最小距离
	double minDistance_KNN=COMPUTATION_DOUBLE_MAX;
	//最小距离对应的三角面片
	Facet facetMinDistance_KNN(-1,-1);

	for(auto itFacetKNN=SurfacetsIncidentKNN.begin();itFacetKNN!=SurfacetsIncidentKNN.end();itFacetKNN++)
	{
		double cosDihedral=0;
		double distance=0;
		Facet tmpFacet=*itFacetKNN;
		if (tmpFacet.first==-1||!is_surface(tmpFacet))
		{
			continue;
		}
		Vertex_triple vf=make_vertex_triple(tmpFacet);
		if (is_facet_in_conflict(tmpFacet)) //判断该面及其mirror facet所在的cell都是conflict
		{
			if (visited_for_facet_extractor(tmpFacet)==0)
			{
				Facet_Idtype idF=get_facet_index(tmpFacet);
				SurfacetsInConflict.push_back(idF);
				VerticesDeletedSurfacet[vf.first]++;
				VerticesDeletedSurfacet[vf.second]++;
				VerticesDeletedSurfacet[vf.third]++;
				mark_facet_visited(tmpFacet);
			}
			
		}
		Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), P);
		if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
		{
			
			cosDihedral=cos_min_dihedral_point_to_facet(P,mirror_facet(tmpFacet));
		
			if(cosDihedral>-0.0)
			{
				distance =GeometricTraits<T>::distance_point_to_facet(P,point(vf.first),
													point(vf.second),point(vf.third));
				if(!PInside)
				{
					distance=-distance;
				}
				if (distance<minDistance_KNN)
				{
					minDistance_KNN=distance;
					facetMinDistance_KNN=tmpFacet;
				}
			}
			else if(cosDihedral>cosMinDihedral_KNN)
			{
				cosMinDihedral_KNN=cosDihedral;
				facetMinDihedral_KNN=tmpFacet;
			}
		
		}
	}


	//step2:求sparse vertices的incident surface facet中最有可能的nearest surface facet
	//最小二面角的余弦值
	double cosMinDihedral_Sparse=-1;
	//最小二面角对应的三角面片
	Facet facetMinDihedral_Sparse(-1,-1);
	//当P与三角面片形成的二面角均为锐角时，P到此面的最小距离
	double minDistance_Sparse=COMPUTATION_DOUBLE_MAX;
	//最小距离对应的三角面片
	Facet facetMinDistance_Sparse(-1,-1);

	for(auto itFacetS=SurfacetsIncidentSparse.begin();itFacetS!=SurfacetsIncidentSparse.end();itFacetS++)
	{
		double cosDihedral=0;
		double distance=0;
		Facet tmpFacet=*itFacetS;
		if (tmpFacet.first==-1||!is_surface(tmpFacet))
		{
			continue;
		}
		Vertex_triple vf=make_vertex_triple(tmpFacet);
		if (is_facet_in_conflict(tmpFacet))  //判断该面及其mirror facet所在的cell都是conflict
		{
			if (visited_for_facet_extractor(tmpFacet)==0)
			{
				Facet_Idtype idF=get_facet_index(tmpFacet);
				SurfacetsInConflict.push_back(idF);
				VerticesDeletedSurfacet[vf.first]++;
				VerticesDeletedSurfacet[vf.second]++;
				VerticesDeletedSurfacet[vf.third]++;
				mark_facet_visited(tmpFacet);
			}
			
		}
		Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), P);
		if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
		{
			
			cosDihedral=cos_min_dihedral_point_to_facet(P,mirror_facet(tmpFacet));
		
			if(cosDihedral>-0.0)
			{
				distance =GeometricTraits<T>::distance_point_to_facet(P,point(vf.first),
													point(vf.second),point(vf.third));
				if(!PInside)
				{
					distance=-distance;
				}
				if (distance<minDistance_Sparse)
				{
					minDistance_Sparse=distance;
					facetMinDistance_Sparse=tmpFacet;
				}
			}
			else if(cosDihedral>cosMinDihedral_Sparse)
			{
				cosMinDihedral_Sparse=cosDihedral;
				facetMinDihedral_Sparse=tmpFacet;
			}
		
		}
	}

	//step3:求conflict region中最有可能的nearest surface facet,包括内部和边界上的表面三角面片
	//最小二面角的余弦值
	double cosMinDihedral_Conflict=-1;
	//最小二面角对应的三角面片
	Facet facetMinDihedral_Conflict(-1,-1);
	//当P与三角面片形成的二面角均为锐角时，P到此面的最小距离
	double minDistance_Conflict=COMPUTATION_DOUBLE_MAX;
	//最小距离对应的三角面片
	Facet facetMinDistance_Conflict(-1,-1);

	for(auto itFacetConflict=SurfacetsInternalConflict.begin();itFacetConflict!=SurfacetsInternalConflict.end();itFacetConflict++)
	{
		double cosDihedral=0;
		double distance=0;
		Facet tmpFacet=*itFacetConflict;
		if (tmpFacet.first==-1||!is_surface(tmpFacet))
		{
			continue;
		}
		
		Vertex_triple vf=make_vertex_triple(tmpFacet);
		if (is_facet_in_conflict(tmpFacet))  //判断该面及其mirror facet所在的cell都是conflict
		{
			if (visited_for_facet_extractor(tmpFacet)==0)
			{
				Facet_Idtype idF=get_facet_index(tmpFacet);
				SurfacetsInConflict.push_back(idF);
				VerticesDeletedSurfacet[vf.first]++;
				VerticesDeletedSurfacet[vf.second]++;
				VerticesDeletedSurfacet[vf.third]++;
				mark_facet_visited(tmpFacet);
			}
			
		}
		Orientation o = GeometricTraits<T>::orientation(point(vf.first), point(vf.second), point(vf.third), P);
		if ((o!=NEGATIVE&&PInside)||(o!=POSITIVE&&!PInside))
		{
			
			cosDihedral=cos_min_dihedral_point_to_facet(P,mirror_facet(tmpFacet));
		
			if(cosDihedral>-0.0)
			{
				distance =GeometricTraits<T>::distance_point_to_facet(P,point(vf.first),
													point(vf.second),point(vf.third));
				if(!PInside)
				{
					distance=-distance;
				}
				if (distance<minDistance_Conflict)
				{
					minDistance_Conflict=distance;
					facetMinDistance_Conflict=tmpFacet;
				}
			}
			else if(cosDihedral>cosMinDihedral_Conflict)
			{
				cosMinDihedral_Conflict=cosDihedral;
				facetMinDihedral_Conflict=tmpFacet;
			}
		
		}
	}

	for (auto itS=SurfacetsInConflict.begin();itS!=SurfacetsInConflict.end();itS++)
	{
		clear_facet_visited(*itS);
	}

	if (facetMinDistance_KNN.first!=-1)
	{
		double mminDistance=minDistance_KNN;
		NearestFacet=facetMinDistance_KNN;
		if (facetMinDistance_Sparse.first!=-1)
		{
			if (minDistance_Sparse<mminDistance)
			{
				mminDistance=minDistance_Sparse;
				NearestFacet=facetMinDistance_Sparse;
			}
		}
		if (facetMinDistance_Conflict.first!=-1)
		{
			if (minDistance_Conflict<mminDistance)
			{
				mminDistance=minDistance_Conflict;
				NearestFacet=facetMinDistance_Conflict;
			}
		}
	}
	else if (facetMinDihedral_KNN.first!=-1)
	{

		double mminDistance=shortest_edge_to_facet(P,facetMinDihedral_KNN);
		NearestFacet=facetMinDihedral_KNN;
		if (facetMinDistance_Sparse.first!=-1)
		{
			if (minDistance_Sparse<mminDistance)
			{
				mminDistance=minDistance_Sparse;
				NearestFacet=facetMinDistance_Sparse;
			}
		}
		else if(facetMinDihedral_Sparse.first!=-1)
		{
			double tmpDis=shortest_edge_to_facet(P,facetMinDihedral_Sparse);
			if (tmpDis<mminDistance)
			{
				if (cosMinDihedral_Sparse>cosMinDihedral_KNN)
				{
					mminDistance=tmpDis;
					NearestFacet=facetMinDihedral_Sparse;
				}
			}
		}
		if (facetMinDistance_Conflict.first!=-1)
		{
			if (minDistance_Conflict<mminDistance)
			{
				mminDistance=minDistance_Conflict;
				NearestFacet=facetMinDistance_Conflict;
			}
		}

	}
	else if (facetMinDistance_Sparse.first!=-1)
	{
		double mminDistance=shortest_edge_to_facet(P,facetMinDistance_Sparse);
		NearestFacet=facetMinDistance_Sparse;

		if(facetMinDihedral_Sparse.first!=-1)
		{
			double tmpDis=shortest_edge_to_facet(P,facetMinDihedral_Sparse);
			if (tmpDis<mminDistance)
			{
				if (cosMinDihedral_Sparse>cosMinDihedral_KNN)
				{
					mminDistance=tmpDis;
					NearestFacet=facetMinDihedral_Sparse;
				}
			}
		}
		if (facetMinDistance_Conflict.first!=-1)
		{
			if (minDistance_Conflict<mminDistance)
			{
				mminDistance=minDistance_Conflict;
				NearestFacet=facetMinDistance_Conflict;
			}
		}
	}
	else if (facetMinDihedral_Sparse.first!=-1)
	{
		double mminDistance=shortest_edge_to_facet(P,facetMinDihedral_Sparse);
		NearestFacet=facetMinDihedral_Sparse;

		if (facetMinDistance_Conflict.first!=-1)
		{
			if (minDistance_Conflict<mminDistance)
			{
				mminDistance=minDistance_Conflict;
				NearestFacet=facetMinDistance_Conflict;
			}
		}
	}
	else
	{
		NearestFacet=facetMinDistance_Conflict;
	}


}
#endif // !DATASTRUCTURE_H

