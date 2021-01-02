#include "AABBTree.h"

AABBTreeMeshInterface::AABBTreeMeshInterface()
{
	ptlist = NULL;
	triindex = NULL;
	ptSize = 0;
	triSize = 0;
}
AABBTreeMeshInterface::~AABBTreeMeshInterface()
{
	ptlist = NULL;
	triindex = NULL;
	ptSize = 0;
	triSize = 0;
}
void AABBTreeMeshInterface::SetPointList(DType* ptl, DUInt ptsize)
{
	ptlist = ptl;
	ptSize = ptsize;
}
void AABBTreeMeshInterface::SetTriangleIndex(DUInt* ids, DUInt trisize)
{
	triindex = ids;
	triSize = trisize;
}

DUInt AABBTreeMeshInterface::GetTriangleSize() const{

	return triSize;
}
DUInt AABBTreeMeshInterface::GetPointListSize() const{
	return ptSize;
}

void   AABBTreeMeshInterface::GetBoundBoxByTriangleIndex(DUInt triID, DType *bdbox)
{
	DType xmin, ymin, zmin, xmax, ymax, zmax;
	xmin = MAX_FLOAT; ymin = MAX_FLOAT; zmin = MAX_FLOAT;
	xmax = MIN_FLOAT; ymax = MIN_FLOAT; zmax = MIN_FLOAT;

	for (int i = 0; i < 3;++i)
	{
		DUInt id  = triindex[3 * triID + i];
		DType x = ptlist[3 * id];
		DType y = ptlist[3 * id + 1];
		DType z = ptlist[3 * id + 2];

		xmin = xmin > x ? x : xmin;
		ymin = ymin > y ? y : ymin;
		zmin = zmin > z ? z : zmin;

		xmax = xmax < x ? x : xmax;
		ymax = ymax < y ? y : ymax;
		zmax = zmax < z ? z : zmax;
	}

	bdbox[0] = xmin; bdbox[1] = ymin; bdbox[2] = zmin;
	bdbox[3] = xmax; bdbox[4] = ymax; bdbox[5] = zmax;
}
void   AABBTreeMeshInterface::GetMeshBoundBox(DType *bdbox)
{
	DType xmin, ymin, zmin, xmax, ymax, zmax;
	xmin = MAX_FLOAT; ymin = MAX_FLOAT; zmin = MAX_FLOAT;
	xmax = MIN_FLOAT; ymax = MIN_FLOAT; zmax = MIN_FLOAT;
 
	for (unsigned int i = 0; i < ptSize;++i)
	{
		DType x = ptlist[3 * i];
		DType y = ptlist[3 * i + 1];
		DType z = ptlist[3 * i + 2];

		xmin = xmin > x ? x : xmin;
		ymin = ymin > y ? y : ymin;
		zmin = zmin > z ? z : zmin;

		xmax = xmax < x ? x : xmax;
		ymax = ymax < y ? y : ymax;
		zmax = zmax < z ? z : zmax;
	}

	bdbox[0] = xmin; bdbox[1] = ymin; bdbox[2] = zmin;
	bdbox[3] = xmax; bdbox[4] = ymax; bdbox[5] = zmax;
}
DType   AABBTreeMeshInterface::GetMeanValueByAxis(DUInt triID, DUInt axis)
{
	DType val = 0.0f;
	for (int i = 0; i < 3; ++i)
	{
		DUInt id = triindex[3 * triID + i];
		val     += ptlist[3 * id + axis];
	}

	return val / 3.0f;
}

void   AABBTreeMeshInterface::GetTriangle(DUInt triID, DType *p1, DType *p2, DType *p3)
{
	DUInt p = triindex[3 * triID];
	DUInt q = triindex[3 * triID + 1];
	DUInt r = triindex[3 * triID + 2];

	p1[0] = ptlist[3 * p]; p1[1] = ptlist[3 * p + 1]; p1[2] = ptlist[3 * p + 2];
	p2[0] = ptlist[3 * q]; p2[1] = ptlist[3 * q + 1]; p2[2] = ptlist[3 * q + 2];
	p3[0] = ptlist[3 * r]; p3[1] = ptlist[3 * r + 1]; p3[2] = ptlist[3 * r + 2];
}

AABBTree::AABBTree()
{
	boundbox = NULL;
	mNodePrimitives = NULL;
	mesh = NULL;
	root = NULL;
	cnt = 0;
	isFirstContact = false;
	isContact = false;
}
AABBTree::~AABBTree()
{
	if (mNodePrimitives)
	{
		delete mNodePrimitives;
		mNodePrimitives = NULL;
	}
	if (mesh)
	{
		delete mesh;
		mesh = NULL;
	}
	if (root)
	{
		delete[] root;
		root = NULL;
	}
	if (boundbox)
	{
		delete[] boundbox;
		boundbox = NULL;
	}
}
 
void AABBTree::SetMesh(AABBTreeMeshInterface *_mesh)
{
	mesh = _mesh;
}
 
void AABBTree::InitParameters()
{
	mNbPrimitives = mesh->GetTriangleSize();
	mNodePrimitives = new DUInt[mNbPrimitives];
	for (unsigned int i = 0; i < mNbPrimitives; ++i)
	{
		mNodePrimitives[i] = i;
	}
	boundbox = new DType[6];
	root = new AABBTreeNode[2 * mNbPrimitives-1];
}
void AABBTree::Build()
{
	if (mesh && (mesh->GetTriangleSize() > 1))
	{
		InitParameters();
		//root = new AABBTreeNode();
		AABBTreeNode *node = &(root[cnt]);
		node->mNbPrimitives = mNbPrimitives;
		node->m_NodePrimitives = mNodePrimitives;
		ComputeGlobeBox(node);
		BuildHierarchy(node);
		//printf("all the node number is :%d\n", cnt);
	}
}

DUInt AABBTree::Split(DUInt axis, AABBTreeNode *node) const
{
	DType SplitValue = 0.0f;
	DUInt nb_prims = node->mNbPrimitives;
	for (unsigned int i = 0; i < nb_prims; i++)
	{
		SplitValue += mesh->GetMeanValueByAxis(node->m_NodePrimitives[i],axis);
	}
	SplitValue = SplitValue / float(nb_prims);

	DUInt NbPos = 0;
	// Loop through all node-related primitives. Their indices range from mNodePrimitives[0] to mNodePrimitives[mNbPrimitives-1].
	// Those indices map the global list in the tree builder.
	for (unsigned int i = 0; i<nb_prims; i++)
	{
		// Test against the splitting value. The primitive value is tested against the enclosing-box center.
		// [We only need an approximate partition of the enclosing box here.]
		DType PrimitiveValue = mesh->GetMeanValueByAxis(node->m_NodePrimitives[i], axis);

		// Reorganize the list of indices in this order: positive - negative.
		if (PrimitiveValue > SplitValue)
		{
			// Swap entries
			DUInt Tmp = node->m_NodePrimitives[i];
			node->m_NodePrimitives[i] = node->m_NodePrimitives[NbPos];
			node->m_NodePrimitives[NbPos] = Tmp;
			// Count primitives assigned to positive space
			NbPos++;
		}
	}

	if (!NbPos || NbPos == node->mNbPrimitives)
		NbPos = node->mNbPrimitives / 2;

	return NbPos;
}
DUInt AABBTree::GetSplitAxis(AABBTreeNode *node) const
{
	DType ext_x = node->box[3] - node->box[0];
	DType ext_y = node->box[4] - node->box[1];
	DType ext_z = node->box[5] - node->box[2];

	DType max_ext = MAX_VAL3(ext_x, ext_y, ext_z);

	if (fabs(ext_x - max_ext) < 1.0e-4)
		return 0;
	else if (fabs(ext_y - max_ext) < 1.0e-4)
		return 1;
	else
		return 2;
}
void AABBTree::ComputeGlobeBox(AABBTreeNode *node) const
{
	DUInt nb_prims = node->mNbPrimitives;
	DType xmin, ymin, zmin, xmax, ymax, zmax;
	xmin = MAX_FLOAT; ymin = MAX_FLOAT; zmin = MAX_FLOAT;
	xmax = MIN_FLOAT; ymax = MIN_FLOAT; zmax = MIN_FLOAT;

	for (unsigned int i = 0; i < nb_prims; i++)
	{
		mesh->GetBoundBoxByTriangleIndex(node->m_NodePrimitives[i], boundbox);
		xmin = xmin > boundbox[0] ? boundbox[0] : xmin;
		ymin = ymin > boundbox[1] ? boundbox[1] : ymin;
		zmin = zmin > boundbox[2] ? boundbox[2] : zmin;

		xmax = xmax < boundbox[3] ? boundbox[3] : xmax;
		ymax = ymax < boundbox[4] ? boundbox[4] : ymax;
		zmax = zmax < boundbox[5] ? boundbox[5] : zmax;
	}

	node->box[0] = xmin;node->box[1] = ymin;node->box[2] = zmin;
	node->box[3] = xmax;node->box[4] = ymax;node->box[5] = zmax;
}
void AABBTree::BuildHierarchy(AABBTreeNode *node)
{
    /// get axis to split
	DUInt axis = GetSplitAxis(node);
	/// split value
	DUInt nbPos = Split(axis, node);
	
	//AABBTreeNode *leftnd  = new AABBTreeNode();
	AABBTreeNode *leftnd = &(root[cnt+1]);
	leftnd->m_NodePrimitives = &node->m_NodePrimitives[0];
	leftnd->mNbPrimitives = nbPos;
	ComputeGlobeBox(leftnd);

	//AABBTreeNode *rightnd = new AABBTreeNode();
	AABBTreeNode *rightnd = &(root[cnt+2]);
	rightnd->m_NodePrimitives = &node->m_NodePrimitives[nbPos];
	rightnd->mNbPrimitives = node->mNbPrimitives - nbPos;
	ComputeGlobeBox(rightnd);
	
	cnt += 2;

	node->left = leftnd;
	node->right = rightnd;

	if (nbPos > 1)
		BuildHierarchy(leftnd);
	if (rightnd->mNbPrimitives > 1)
		BuildHierarchy(rightnd);
}
void AABBTree::_RayStab(AABBTreeNode *node) 
{
	DType center[3];
	DType extents[3];

	center[0] = (node->box[3] + node->box[0])*0.5f;
	center[1] = (node->box[4] + node->box[1])*0.5f;
	center[2] = (node->box[5] + node->box[2])*0.5f;

	extents[0] = (node->box[3] - node->box[0])*0.5f;
	extents[1] = (node->box[4] - node->box[1])*0.5f;
	extents[2] = (node->box[5] - node->box[2])*0.5f;

	// Perform Ray-AABB overlap test
	if (!RayAABBOverlap(center, extents))	return;

	if (node->isLeaf())
	{
		RayPrim(node);
	}
	else
	{
		_RayStab(node->left);

		if (ContactFound()) return;

		_RayStab(node->right);
	}
}
bool AABBTree::ContactFound()
{
	if (isFirstContact && isContact)
		return true;
	else
		return false;
}
void AABBTree::RayHit(DType orgx, DType orgy, DType orgz, DType dirx, DType diry, DType dirz)
{
	ray[0] = orgx; ray[1] = orgy; ray[2] = orgz;
	ray[3] = dirx; ray[4] = diry; ray[5] = dirz;

	std::map<DType, DUInt>().swap(resInfo);
	isContact = false;

	_RayStab(&root[0]);
}

bool AABBTree::RayAABBOverlap(DType *center, DType *extents) const
{
	DType Dx = ray[0] - center[0];	if ((fabs(Dx) > extents[0]) && Dx*ray[3] >= 0.0f)	return false;
	DType Dy = ray[1] - center[1];	if ((fabs(Dy) > extents[1]) && Dy*ray[4] >= 0.0f)	return false;
	DType Dz = ray[2] - center[2];	if ((fabs(Dz) > extents[2]) && Dz*ray[5] >= 0.0f)	return false;

	DType f;
	f = ray[4] * Dz - ray[5] * Dy;		if (fabs(f) > extents[1] * fabs(ray[5]) + extents[2] * fabs(ray[4]))	return false;
	f = ray[5] * Dx - ray[3] * Dz;		if (fabs(f) > extents[0] * fabs(ray[5]) + extents[2] * fabs(ray[3]))	return false;
	f = ray[3] * Dy - ray[4] * Dx;		if (fabs(f) > extents[0] * fabs(ray[4]) + extents[1] * fabs(ray[3]))	return false;

	return true;
}

void AABBTree::RayPrim(AABBTreeNode *node)
{
	DUInt nbPrim = node->mNbPrimitives;
	assert(1==nbPrim);
	DType p1[3];DType p2[3];DType p3[3];

	DUInt id = node->m_NodePrimitives[0];
	mesh->GetTriangle(id, p1, p2, p3);
	DType distance = 0.0f;

	if (TriangleAABBOverlap(p1, p2, p3, distance))
	{
		isContact = true;
		resInfo.insert(std::make_pair(distance, id));
	}
}

bool AABBTree::TriangleAABBOverlap(DType *p1, DType * p2, DType * p3, DType &mDistance)
{
	DType edge1[3], edge2[3];
	sub(p2, p1, edge1);
	sub(p3, p1, edge2);

	DType p[3];
	cross(&ray[3], edge2, p);
	DType det = dot(edge1, p);;

	if (det < 1.0e-6 && det > 1.0e-6)
		return false;

	DType OneOverDet = 1.0f / det;
	DType t[3];
	sub(&ray[0], p1, t);

	DType u = OneOverDet*dot(t, p);
	if (u < 0.0 || u > 1.0)
		return false;

	DType q[3];
	cross(t, edge1, q);
	DType v = OneOverDet*dot(&ray[3], q);
	if (v < 0.0 || v > 1.0)
		return false;

	mDistance = dot(edge2, q) * OneOverDet;
	if (mDistance < 0.0)
		return false;

	return true;
}
DUInt AABBTree::GetIntersectTriangleIndex(DUInt idx) const
{
	auto itor = resInfo.begin();
	int id = 0;
	for (; itor != resInfo.end(); ++itor,++id)
	{
		if (id == idx)
			return itor->second;
	}
	
	return -1;
}
DType AABBTree::GetIntersectDistance(DUInt idx) const
{
	auto itor = resInfo.begin();
	int id = 0;
	for (; itor != resInfo.end(); ++itor, ++id)
	{
		if (id == idx)
			return itor->first;
	}

	return -1;
}