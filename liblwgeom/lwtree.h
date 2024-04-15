/**********************************************************************
 *
 * PostGIS - Spatial Types for PostgreSQL
 * http://postgis.net
 *
 * PostGIS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * PostGIS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PostGIS.  If not, see <http://www.gnu.org/licenses/>.
 *
 **********************************************************************
 *
 * Copyright (C) 2009-2012 Paul Ramsey <pramsey@cleverelephant.ca>
 *
 **********************************************************************/

#define RECT_NODE_SIZE 8

#define RECT_TREE_QUEUE_SIZE 64

typedef enum
{
	RECT_NODE_INTERNAL_TYPE,
	RECT_NODE_LEAF_TYPE
} RECT_NODE_TYPE;

typedef enum
{
	RECT_NODE_RING_NONE = 0,
	RECT_NODE_RING_EXTERIOR,
	RECT_NODE_RING_INTERIOR
} RECT_NODE_RING_TYPE;

typedef enum
{
	RECT_NODE_SEG_UNKNOWN = 0,
	RECT_NODE_SEG_POINT,
	RECT_NODE_SEG_LINEAR,
	RECT_NODE_SEG_CIRCULAR
} RECT_NODE_SEG_TYPE;

typedef struct
{
	const POINTARRAY *pa;
	RECT_NODE_SEG_TYPE seg_type;
	int seg_num;
} RECT_NODE_LEAF;

struct rect_node;

typedef struct
{
	unsigned int num_nodes;
	RECT_NODE_RING_TYPE ring_type;
	struct rect_node *nodes[RECT_NODE_SIZE];
	int sorted;
} RECT_NODE_INTERNAL;

typedef struct rect_node
{
	RECT_NODE_TYPE type;
	unsigned char geom_type;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double d;
	union {
		RECT_NODE_INTERNAL i;
		RECT_NODE_LEAF l;
	};
} RECT_NODE;

/*
* Allocate a flat array of RECT_NODE for
* locality reasons, and then parcel out references
* into that array in the *_new() functions.
*/
typedef struct rect_tree
{
	RECT_NODE *root;
	RECT_NODE *buffer;
	uint32_t next;
	uint32_t capacity;
} RECT_TREE;


typedef struct rect_tree_query_state
{
	double threshold;
	double min_dist;
	double max_dist;
	POINT2D p1;
	POINT2D p2;
} RECT_TREE_QUERY_STATE;

/**
* Create a tree index on top an LWGEOM. Do not free the LWGEOM
* until the tree is freed.
*/
RECT_TREE * rect_tree_from_lwgeom(const LWGEOM *geom);

/**
* Test if two RECT_NODE trees intersect one another.
*/
int rect_tree_intersects_tree(const RECT_TREE *tree1, const RECT_TREE *tree2);

/**
* Return the distance between two RECT_NODE trees.
*/
double rect_tree_distance_tree(const RECT_TREE *t1, const RECT_TREE *t2, double threshold);

/**
* Free the rect-tree memory
*/
void rect_tree_free(RECT_TREE *node);

int rect_tree_contains_point(const RECT_TREE *tree, const POINT2D *pt);
RECT_NODE * rect_node_from_ptarray(const POINTARRAY *pa, int geom_type, RECT_TREE *tree);
LWGEOM * rect_node_to_lwgeom(const RECT_NODE *tree);
char * rect_node_to_wkt(const RECT_NODE *node);
void rect_node_printf(const RECT_NODE *node, int depth);
uint32_t rect_node_count_from_lwgeom(const LWGEOM *lwgeom);

