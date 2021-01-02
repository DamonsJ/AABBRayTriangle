///////////////////////////////////////////////////////////////////
/*

///--------------------------------------------------------
/// Ray Triangle intersect data structure is designed to make intersect
/// between ray and triangle.
/// we use AABB tree as main data structure and check if the
/// ray is intersected with box and triangle
/// The AABB tree is written learned from OPCODE (http://www.codercorner.com/Opcode.htm) 
/// which is an open source collision detect.
/// File:			AABBTree.h
/// Programmer:		Damons(shlkl99@163.com)
/// Description:    Standard AABB-tree intersect
/// Last modified:	2016/04/07 (Version 1.0)
///--------------------------------------------------------

*/


#ifndef     _AABBTREE_H_
#define		_AABBTREE_H_

#include <iostream>
#include <vector>
#include <map>
#include <assert.h>

#define  MAX_FLOAT FLT_MAX
#define  MIN_FLOAT (-FLT_MAX)

#define  MAX_VAL2(a,b)  (( (a) > (b))?(a):(b))
#define  MAX_VAL3(a,b,c)  ( (MAX_VAL2(a,b)) > (b))?(MAX_VAL2(a,b)):(b)

typedef  unsigned int  DUInt;
typedef  float		   DType;
/*
///@class:AABBTreeNode
///@breif:树的节点
/// 每个节点存储左右子树，以及当前包围盒
*/
class AABBTreeNode
{
public:
	AABBTreeNode()
	{
		left  = NULL;
		right = NULL;
	}
	~AABBTreeNode()
	{
		if (left)
		{
			//delete left;
			left = NULL;
		}
		if (right)
		{
			//delete right;
			right = NULL;
		}
	}
	//判断是否为叶子节点（如果是则进行三角形与射线的求交）
	bool isLeaf(){
		return (NULL == right);
	}

	DType box[6];//当前节点的包围盒
	DUInt *m_NodePrimitives;//当前节点的包围盒内所包含的三角形的索引
	DUInt  mNbPrimitives;//当前节点的包围盒内所包含的三角形的数目

	AABBTreeNode *left;
	AABBTreeNode *right;

};
/*
///@class:AABBTreeMeshInterface
///@breif:三角网格接口
/// 存储待划分的所有的三角形的信息，包括点与索引
*/
class AABBTreeMeshInterface
{
public:
	AABBTreeMeshInterface();
	~AABBTreeMeshInterface();
public:
	void SetPointList(DType* ptl,DUInt ptsize);//设置所有的点 每个点为xyz坐标顺序排放，深拷贝
	void SetTriangleIndex(DUInt* ids,DUInt trisize);//设置三角形点的索引 每个三角形三个点，深拷贝
public:
	DUInt GetTriangleSize() const;//获取三角形面片的总个数
	DUInt GetPointListSize() const;//获取三角形顶点的总个数
	void   GetBoundBoxByTriangleIndex(DUInt triID,DType *bdbox);//通过三角形的索引，获取包围盒
	void   GetMeshBoundBox(DType *bdbox);//整个网格的包围盒
	DType  GetMeanValueByAxis(DUInt triID, DUInt axis);//获取第triID个三角形的 第axis维坐标的平均值
	void   GetTriangle(DUInt triID, DType *p1, DType *p2, DType *p3);//获取三角形的点
private:
	//传入的为列表的指针，外部析构，内部不析构，仅置空
	DType*    ptlist;//点的列表
	DUInt*    triindex;//索引的列表
	DUInt     ptSize;
	DUInt     triSize;
};

/*
///@class:AABBTree
///@breif:树的入口类
/// 存储所有树的节点以及进行射线与三角形求交功能
*/

class AABBTree
{
public:
	AABBTree();
	~AABBTree();
public:
	void SetMesh(AABBTreeMeshInterface *_mesh);//设置网格接口
	void Build(); //构建AABB树
	//射线求交
	void RayHit(DType orgx, DType orgy, DType orgz, DType dirx, DType diry, DType dirz);
	//获取三角形与射线的相交个数
	DUInt  GetNbIntersects() const { return resInfo.size(); }
	//获取相交三角形的index
	DUInt GetIntersectTriangleIndex(DUInt idx) const;
	//获取相交点距原点的距离
	DType GetIntersectDistance(DUInt idx) const;
	//设置是否第一个相交之后就返回，交到一次就返回
	void SetFirstContact(bool isfirst){ isFirstContact = isfirst; }
	//获取三角形的三个点
	void GetTriangle(DUInt triID, DType *p1, DType *p2, DType *p3) const
	{
		mesh->GetTriangle(triID, p1, p2, p3);
	}
	bool IsValid(){
		return (NULL != root);
	}
protected:
	//初始化树的参数
	void InitParameters();
	//获取当前节点的分裂轴 选择最长轴
	inline DUInt GetSplitAxis(AABBTreeNode *node) const;
    //在最长轴划分节点
	DUInt Split(DUInt axis, AABBTreeNode *node) const;
	//计算当前节点包围盒
	void ComputeGlobeBox(AABBTreeNode *node) const;
	//递归构建树
	void BuildHierarchy(AABBTreeNode *node);
	//进行射线与节点相交判断
	void _RayStab(AABBTreeNode *node) ;
	//是否相交
	bool ContactFound();
	//进行射线与三角形相交判断
	void RayPrim(AABBTreeNode *node);

protected:
	//进行射线与包围盒相交算法
	bool RayAABBOverlap(DType *center, DType *extents) const;
	//进行射线与三角形相交算法
	bool TriangleAABBOverlap(DType *p1, DType * p2, DType * p3, DType &mDistance);
	//简单的数学函数
	inline void    cross(DType* v1, DType* v2, DType* n) const
	{
		n[0] = v1[1] * v2[2] - v2[1] * v1[2];
		n[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
		n[2] = v1[0] * v2[1] - v2[0] * v1[1];
	}

	inline DType   dot(DType* v1, DType* v2) const
	{
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	}

	inline void    sub(DType* v1, DType* v2, DType *n) const
	{
		n[0] = v1[0] - v2[0];
		n[1] = v1[1] - v2[1];
		n[2] = v1[2] - v2[2];
	}

private:
	DType *boundbox;//临时变量 为了计算包围盒
	DUInt *mNodePrimitives;//所有三角形的索引，不是三角行点的索引
	DUInt mNbPrimitives;//三角形个数
	AABBTreeMeshInterface *mesh;//网格接口
	AABBTreeNode *root;//所有节点的数组
	int cnt;//以添加节点的个数
	DType ray[6];//射线
	bool  isFirstContact;
	bool  isContact;
	std::map<DType, DUInt> resInfo;//结果信息
};
#endif