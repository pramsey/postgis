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

#include "liblwgeom_internal.h"
#include "lwgeom_log.h"
#include "lwtree.h"
#include "measures.h"
#include "pqueue.h"

static inline int
rect_node_is_leaf(const RECT_NODE *node)
{
	return node->type == RECT_NODE_LEAF_TYPE;
}

/*
* Support qsort of nodes for collection/multi types so nodes
* are in "spatial adjacent" order prior to merging.
*/
static int
rect_node_cmp(const void *pn1, const void *pn2)
{
	GBOX b1, b2;
	RECT_NODE *n1 = *((RECT_NODE**)pn1);
	RECT_NODE *n2 = *((RECT_NODE**)pn2);
	uint64_t h1, h2;
	b1.flags = 0;
	b1.xmin = n1->xmin;
	b1.xmax = n1->xmax;
	b1.ymin = n1->ymin;
	b1.ymax = n1->ymax;

	b2.flags = 0;
	b2.xmin = n2->xmin;
	b2.xmax = n2->xmax;
	b2.ymin = n2->ymin;
	b2.ymax = n2->ymax;

	h1 = gbox_get_sortable_hash(&b1, 0);
	h2 = gbox_get_sortable_hash(&b2, 0);
	return h1 < h2 ? -1 : (h1 > h2 ? 1 : 0);
}

/**
* Recurse from top of node tree and free all children.
* does not free underlying point array.
*/
void
rect_tree_free(RECT_TREE *tree)
{
	if (!tree) return;
	/* All the nodes are in the buffer */
	if (tree->buffer)
	{
		lwfree(tree->buffer);
	}
	lwfree(tree);
}

static RECT_NODE *
rect_tree_get_node(RECT_TREE *tree)
{
	if (tree->next >= tree->capacity)
		lwerror("%s: tree capacity (%d) exceeded", __func__, tree->capacity);
	RECT_NODE *node = tree->buffer + tree->next;
	tree->next++;
	return node;
}


static int
rect_leaf_node_intersects(const RECT_NODE_LEAF *n1, const RECT_NODE_LEAF *n2)
{
	const POINT2D *p1, *p2, *p3, *q1, *q2, *q3;
	DISTPTS dl;
	lw_dist2d_distpts_init(&dl, 1);
	switch (n1->seg_type)
	{
		case RECT_NODE_SEG_POINT:
		{
			p1 = getPoint2d_cp(n1->pa, n1->seg_num);

			switch (n2->seg_type)
			{
				case RECT_NODE_SEG_POINT:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					lw_dist2d_pt_pt(q1, p1, &dl);
					return dl.distance == 0.0;

				case RECT_NODE_SEG_LINEAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num+1);
					lw_dist2d_pt_seg(p1, q1, q2, &dl);
					return dl.distance == 0.0;

				case RECT_NODE_SEG_CIRCULAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num*2);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num*2+1);
					q3 = getPoint2d_cp(n2->pa, n2->seg_num*2+2);
					lw_dist2d_pt_arc(p1, q1, q2, q3, &dl);
					return dl.distance == 0.0;

				default:
					lwerror("%s: unsupported segment type", __func__);
					break;
			}

			break;
		}

		case RECT_NODE_SEG_LINEAR:
		{
			p1 = getPoint2d_cp(n1->pa, n1->seg_num);
			p2 = getPoint2d_cp(n1->pa, n1->seg_num+1);

			switch (n2->seg_type)
			{
				case RECT_NODE_SEG_POINT:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					lw_dist2d_pt_seg(q1, p1, p2, &dl);
					return dl.distance == 0.0;

				case RECT_NODE_SEG_LINEAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num+1);
					return lw_segment_intersects(p1, p2, q1, q2) > 0;

				case RECT_NODE_SEG_CIRCULAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num*2);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num*2+1);
					q3 = getPoint2d_cp(n2->pa, n2->seg_num*2+2);
					lw_dist2d_seg_arc(p1, p2, q1, q2, q3, &dl);
					return dl.distance == 0.0;

				default:
					lwerror("%s: unsupported segment type", __func__);
					break;
			}

			break;
		}
		case RECT_NODE_SEG_CIRCULAR:
		{
			p1 = getPoint2d_cp(n1->pa, n1->seg_num*2);
			p2 = getPoint2d_cp(n1->pa, n1->seg_num*2+1);
			p3 = getPoint2d_cp(n1->pa, n1->seg_num*2+2);

			switch (n2->seg_type)
			{
				case RECT_NODE_SEG_POINT:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					lw_dist2d_pt_arc(q1, p1, p2, p3, &dl);
					return dl.distance == 0.0;

				case RECT_NODE_SEG_LINEAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num+1);
					lw_dist2d_seg_arc(q1, q2, p1, p2, p3, &dl);
					return dl.distance == 0.0;

				case RECT_NODE_SEG_CIRCULAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num*2);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num*2+1);
					q3 = getPoint2d_cp(n2->pa, n2->seg_num*2+2);
					lw_dist2d_arc_arc(p1, p2, p3, q1, q2, q3, &dl);
					return dl.distance == 0.0;

				default:
					lwerror("%s: unsupported segment type", __func__);
					break;
			}

			break;
		}
		default:
			return LW_FALSE;
	}
	return LW_FALSE;
}


/*
* Returns 1 if segment is to the right of point.
*/
static inline int
rect_leaf_node_segment_side(const RECT_NODE_LEAF *node, const POINT2D *q, int *on_boundary)
{
	const POINT2D *p1, *p2, *p3;
	switch (node->seg_type)
	{
		case RECT_NODE_SEG_LINEAR:
		{
			int side;
			p1 = getPoint2d_cp(node->pa, node->seg_num);
			p2 = getPoint2d_cp(node->pa, node->seg_num+1);

			side = lw_segment_side(p1, p2, q);

			/* Always note case where we're on boundary */
			if (side == 0 && lw_pt_in_seg(q, p1, p2))
			{
				*on_boundary = LW_TRUE;
				return 0;
			}

			/* Segment points up and point is on left */
			if (p1->y < p2->y && side == -1 && q->y != p2->y)
			{
				return 1;
			}

			/* Segment points down and point is on right */
			if (p1->y > p2->y && side == 1 && q->y != p2->y)
			{
				return 1;
			}

			/* Segment is horizontal, do we cross first point? */
			if (p1->y == p2->y && q->x < p1->x)
			{
				return 1;
			}

			return 0;
		}
		case RECT_NODE_SEG_CIRCULAR:
		{
			int arc_side, seg_side;

			p1 = getPoint2d_cp(node->pa, node->seg_num*2);
			p2 = getPoint2d_cp(node->pa, node->seg_num*2+1);
			p3 = getPoint2d_cp(node->pa, node->seg_num*2+2);

			/* Always note case where we're on boundary */
			arc_side = lw_arc_side(p1, p2, p3, q);
			if (arc_side == 0)
			{
				*on_boundary = LW_TRUE;
				return 0;
			}

			seg_side = lw_segment_side(p1, p3, q);
			if (seg_side == arc_side)
			{
				/* Segment points up and point is on left */
				if (p1->y < p3->y && seg_side == -1 && q->y != p3->y)
				{
					return 1;
				}

				/* Segment points down and point is on right */
				if (p1->y > p3->y && seg_side == 1 && q->y != p3->y)
				{
					return 1;
				}
			}
			else
			{
				/* Segment points up and point is on left */
				if (p1->y < p3->y && seg_side == 1 && q->y != p3->y)
				{
					return 1;
				}

				/* Segment points down and point is on right */
				if (p1->y > p3->y && seg_side == -1 && q->y != p3->y)
				{
					return 1;
				}

				/* Segment is horizontal, do we cross first point? */
				if (p1->y == p3->y)
				{
					return 1;
				}
			}

			return 0;

		}
		default:
		{
			lwerror("%s: unsupported seg_type - %d", __func__, node->seg_type);
			return 0;
		}
	}

	return 0;
}

/*
* Only pass the head of a ring. Either a LinearRing from a polygon
* or a CompoundCurve from a CurvePolygon.
* Takes a horizontal line through the ring, and adds up the
* crossing directions. An "inside" point will be on the same side
* ("both left" or "both right") of the edges the line crosses,
* so it will have a non-zero return sum. An "outside" point will
* be on both sides, and will have a zero return sum.
*/
static int
rect_node_ring_contains_point(const RECT_NODE *node, const POINT2D *pt, int *on_boundary)
{
	/* Only test nodes that straddle our stabline vertically */
	/* and might be to the right horizontally */
	if (node->ymin <= pt->y && pt->y <= node->ymax && pt->x <= node->xmax)
	{
		if (rect_node_is_leaf(node))
		{
			return rect_leaf_node_segment_side(&node->l, pt, on_boundary);
		}
		else
		{
			uint32_t i, r = 0;
			for (i = 0; i < node->i.num_nodes; i++)
			{
				r += rect_node_ring_contains_point(node->i.nodes[i], pt, on_boundary);
			}
			return r;
		}
	}
	return 0;
}


/*
* Only pass in the head of an "area" type. Polygon or CurvePolygon.
* Sums up containment of exterior (+1) and interior (-1) rings, so
* that zero is uncontained, +1 is contained and negative is an error
* (multiply contained by interior rings?)
*/
static int
rect_node_area_contains_point(const RECT_NODE *node, const POINT2D *pt)
{
	/* Can't do anything with a leaf node, makes no sense */
	if (rect_node_is_leaf(node))
		return 0;

	/* Iterate into area until we find ring heads */
	if (node->i.ring_type == RECT_NODE_RING_NONE)
	{
		uint32_t i, sum = 0;
		for (i = 0; i < node->i.num_nodes; i++)
			sum += rect_node_area_contains_point(node->i.nodes[i], pt);
		return sum;
	}
	/* See if the ring encloses the point */
	else
	{
		int on_boundary = 0;
		int edge_crossing_count = rect_node_ring_contains_point(node, pt, &on_boundary);
		/* Odd number of crossings => contained */
		int contained = (edge_crossing_count % 2 == 1);
		/* External rings return positive containment, interior ones negative, */
		/* so that a point-in-hole case nets out to zero (contained by both */
		/* interior and exterior rings. */
		if (node->i.ring_type == RECT_NODE_RING_INTERIOR)
		{
			return on_boundary ? 0 : -1 * contained;
		}
		else
		{
			return contained || on_boundary;
		}

	}
}

/*
* Simple containment test for node/point inputs
*/
static inline int
rect_node_bounds_point(const RECT_NODE *node, const POINT2D *pt)
{
	if (pt->y < node->ymin || pt->y > node->ymax ||
		pt->x < node->xmin || pt->x > node->xmax)
		return LW_FALSE;
	else
		return LW_TRUE;
}

/*
* Pass in arbitrary tree, get back true if point is contained or on boundary,
* and false otherwise.
*/
static int
rect_node_contains_point(const RECT_NODE *node, const POINT2D *pt)
{
	uint32_t i;
	int c;

	/* Empty/invalid things cannot contain or be contained */
	if (!(node && pt)) return LW_FALSE;

	/* Object cannot contain point if bounds don't */
	if (!rect_node_bounds_point(node, pt))
		return LW_FALSE;

	switch (node->geom_type)
	{
		case POLYGONTYPE:
		case CURVEPOLYTYPE:
			return rect_node_area_contains_point(node, pt) > 0;

		case MULTIPOLYGONTYPE:
		case MULTISURFACETYPE:
		case COLLECTIONTYPE:
		{
			for (i = 0; i < node->i.num_nodes; i++)
			{
				c = rect_node_area_contains_point(node->i.nodes[i], pt);
				if (c) return LW_TRUE;
			}
			return LW_FALSE;
		}

		default:
			return LW_FALSE;
	}
}

int rect_tree_contains_point(const RECT_TREE *tree, const POINT2D *pt)
{
    return rect_node_contains_point(tree->root, pt);
}


/*
* For area types, doing intersects and distance, we will
* need to do a point-in-poly test first to find the full-contained
* case where an intersection exists without any edges actually
* intersecting.
*/
static int
rect_node_is_area(const RECT_NODE *node)
{
	switch (node->geom_type)
	{
		case POLYGONTYPE:
		case CURVEPOLYTYPE:
		case MULTISURFACETYPE:
			return LW_TRUE;

		case COLLECTIONTYPE:
		{
			if (rect_node_is_leaf(node))
				return LW_FALSE;
			else
			{
				uint32_t i;
				for (i = 0; i < node->i.num_nodes; i++)
				{
					if (rect_node_is_area(node->i.nodes[i]))
						return LW_TRUE;
				}
				return LW_FALSE;
			}
		}

		default:
			return LW_FALSE;
	}
}

static RECT_NODE_SEG_TYPE lwgeomTypeArc[] =
{
	RECT_NODE_SEG_UNKNOWN,  // "Unknown"
	RECT_NODE_SEG_POINT,    // "Point"
	RECT_NODE_SEG_LINEAR,   // "LineString"
	RECT_NODE_SEG_LINEAR,   // "Polygon"
	RECT_NODE_SEG_UNKNOWN,  // "MultiPoint"
	RECT_NODE_SEG_LINEAR,   // "MultiLineString"
	RECT_NODE_SEG_LINEAR,   // "MultiPolygon"
	RECT_NODE_SEG_UNKNOWN,  // "GeometryCollection"
	RECT_NODE_SEG_CIRCULAR, // "CircularString"
	RECT_NODE_SEG_UNKNOWN,  // "CompoundCurve"
	RECT_NODE_SEG_UNKNOWN,  // "CurvePolygon"
	RECT_NODE_SEG_UNKNOWN,  // "MultiCurve"
	RECT_NODE_SEG_UNKNOWN,  // "MultiSurface"
	RECT_NODE_SEG_LINEAR,   // "PolyhedralSurface"
	RECT_NODE_SEG_LINEAR,   // "Triangle"
	RECT_NODE_SEG_LINEAR    // "Tin"
};

/*
* Create a new leaf node.
*/
static RECT_NODE *
rect_node_leaf_new(const POINTARRAY *pa, int seg_num, int geom_type, RECT_TREE *tree)
{
	const POINT2D *p1, *p2, *p3;
	RECT_NODE *node;
	GBOX gbox;
	RECT_NODE_SEG_TYPE seg_type = lwgeomTypeArc[geom_type];

	switch (seg_type)
	{
		case RECT_NODE_SEG_POINT:
		{
			p1 = getPoint2d_cp(pa, seg_num);
			gbox.xmin = gbox.xmax = p1->x;
			gbox.ymin = gbox.ymax = p1->y;
			break;
		}

		case RECT_NODE_SEG_LINEAR:
		{
			p1 = getPoint2d_cp(pa, seg_num);
			p2 = getPoint2d_cp(pa, seg_num+1);
			/* Zero length edge, doesn't get a node */
			if ((p1->x == p2->x) && (p1->y == p2->y))
				return NULL;
			gbox.xmin = FP_MIN(p1->x, p2->x);
			gbox.xmax = FP_MAX(p1->x, p2->x);
			gbox.ymin = FP_MIN(p1->y, p2->y);
			gbox.ymax = FP_MAX(p1->y, p2->y);
			break;
		}

		case RECT_NODE_SEG_CIRCULAR:
		{
			p1 = getPoint2d_cp(pa, 2*seg_num);
			p2 = getPoint2d_cp(pa, 2*seg_num+1);
			p3 = getPoint2d_cp(pa, 2*seg_num+2);
			/* Zero length edge, doesn't get a node */
			if ((p1->x == p2->x) && (p2->x == p3->x) &&
				(p1->y == p2->y) && (p2->y == p3->y))
				return NULL;
			lw_arc_calculate_gbox_cartesian_2d(p1, p2, p3, &gbox);
			break;
		}

		default:
		{
			lwerror("%s: unsupported seg_type - %d", __func__, seg_type);
			return NULL;
		}
	}

	node = rect_tree_get_node(tree);
	node->type = RECT_NODE_LEAF_TYPE;
	node->geom_type = geom_type;
	node->xmin = gbox.xmin;
	node->xmax = gbox.xmax;
	node->ymin = gbox.ymin;
	node->ymax = gbox.ymax;
	node->l.seg_num = seg_num;
	node->l.seg_type = seg_type;
	node->l.pa = pa;
	return node;
}


static void
rect_node_internal_add_node(RECT_NODE *node, RECT_NODE *add)
{
	if (rect_node_is_leaf(node))
		lwerror("%s: call on leaf node", __func__);
	node->xmin = FP_MIN(node->xmin, add->xmin);
	node->xmax = FP_MAX(node->xmax, add->xmax);
	node->ymin = FP_MIN(node->ymin, add->ymin);
	node->ymax = FP_MAX(node->ymax, add->ymax);
	node->i.nodes[node->i.num_nodes++] = add;
	return;
}


static RECT_NODE *
rect_node_internal_new(const RECT_NODE *seed, RECT_TREE *tree)
{
	RECT_NODE *node = rect_tree_get_node(tree);
	node->xmin = seed->xmin;
	node->xmax = seed->xmax;
	node->ymin = seed->ymin;
	node->ymax = seed->ymax;
	node->geom_type = seed->geom_type;
	node->type = RECT_NODE_INTERNAL_TYPE;
	node->i.num_nodes = 0;
	node->i.ring_type = RECT_NODE_RING_NONE;
	node->i.sorted = 0;
	return node;
}

/*
* We expect the incoming nodes to be in a spatially coherent
* order. For incoming nodes derived from point arrays,
* the very fact that they are
* a vertex list implies a reasonable ordering: points nearby in
* the list will be nearby in space. For incoming nodes from higher
* level structures (collections, etc) the caller should sort the
* nodes using a z-order first, so that this merge step results in a
* spatially coherent structure.
*/
static RECT_NODE *
rect_nodes_merge(RECT_NODE ** nodes, uint32_t num_nodes, RECT_TREE *tree)
{
	if (num_nodes < 1)
	{
		return NULL;
	}

	while (num_nodes > 1)
	{
		uint32_t i, k = 0;
		RECT_NODE *node = NULL;
		for (i = 0; i < num_nodes; i++)
		{
			if (!node)
				node = rect_node_internal_new(nodes[i], tree);

			rect_node_internal_add_node(node, nodes[i]);

			if (node->i.num_nodes == RECT_NODE_SIZE)
			{
				nodes[k++] = node;
				node = NULL;
			}
		}
		if (node)
			nodes[k++] = node;
		num_nodes = k;
	}

	return nodes[0];
}

/*
* Build a tree of nodes from a point array, one node per edge.
*/
RECT_NODE *
rect_node_from_ptarray(const POINTARRAY *pa, int geom_type, RECT_TREE *tree)
{
	int num_edges = 0, i = 0, j = 0;
	RECT_NODE_SEG_TYPE seg_type = lwgeomTypeArc[geom_type];
	RECT_NODE **nodes = NULL;
	RECT_NODE *root = NULL;

	/* No-op on empty ring/line/pt */
	if ( pa->npoints < 1 )
		return NULL;

	/* For arcs, 3 points per edge, for lines, 2 per edge */
	switch(seg_type)
	{
		case RECT_NODE_SEG_POINT:
			return rect_node_leaf_new(pa, 0, geom_type, tree);
			break;
		case RECT_NODE_SEG_LINEAR:
			num_edges = pa->npoints - 1;
			break;
		case RECT_NODE_SEG_CIRCULAR:
			num_edges = (pa->npoints - 1)/2;
			break;
		default:
			lwerror("%s: unsupported seg_type - %d", __func__, seg_type);
	}

	/* First create a flat list of nodes, one per edge. */
	nodes = lwalloc(sizeof(RECT_NODE*) * num_edges);
	for (i = 0; i < num_edges; i++)
	{
		RECT_NODE *node = rect_node_leaf_new(pa, i, geom_type, tree);
		if (node) /* Not zero length? */
			nodes[j++] = node;
	}

	/* Merge the list into a tree */
	root = rect_nodes_merge(nodes, j, tree);

	/* Free the old list structure, leaving the tree in place */
	lwfree(nodes);

	/* Return top of tree */
	return root;
}

LWGEOM *
rect_node_to_lwgeom(const RECT_NODE *node)
{
	LWGEOM *poly = (LWGEOM*)lwpoly_construct_envelope(0, node->xmin, node->ymin, node->xmax, node->ymax);
	if (rect_node_is_leaf(node))
	{
		return poly;
	}
	else
	{
		uint32_t i;
		LWCOLLECTION *col = lwcollection_construct_empty(COLLECTIONTYPE, 0, 0, 0);
		lwcollection_add_lwgeom(col, poly);
		for (i = 0; i < node->i.num_nodes; i++)
		{
			lwcollection_add_lwgeom(col, rect_node_to_lwgeom(node->i.nodes[i]));
		}
		return (LWGEOM*)col;
	}
}

char *
rect_node_to_wkt(const RECT_NODE *node)
{
	LWGEOM *geom = rect_node_to_lwgeom(node);
	char *wkt = lwgeom_to_wkt(geom, WKT_ISO, 12, 0);
	lwgeom_free(geom);
	return wkt;
}

void
rect_node_printf(const RECT_NODE *node, int depth)
{
	printf("%*s----\n", depth, "");
	printf("%*stype: %d\n", depth, "", node->type);
	printf("%*sgeom_type: %d\n", depth, "", node->geom_type);
	printf("%*sbox: %g %g, %g %g\n", depth, "", node->xmin, node->ymin, node->xmax, node->ymax);
	if (node->type == RECT_NODE_LEAF_TYPE)
	{
		printf("%*sseg_type: %d\n", depth, "", node->l.seg_type);
		printf("%*sseg_num: %d\n", depth, "", node->l.seg_num);
	}
	else
	{
		uint32_t i;
		for (i = 0; i < node->i.num_nodes; i++)
		{
			rect_node_printf(node->i.nodes[i], depth+2);
		}
	}
}

static RECT_NODE *
rect_node_from_lwpoint(const LWGEOM *lwgeom, RECT_TREE *tree)
{
	const LWPOINT *lwpt = (const LWPOINT*)lwgeom;
	return rect_node_from_ptarray(lwpt->point, lwgeom->type, tree);
}

static RECT_NODE *
rect_node_from_lwline(const LWGEOM *lwgeom, RECT_TREE *tree)
{
	const LWLINE *lwline = (const LWLINE*)lwgeom;
	return rect_node_from_ptarray(lwline->points, lwgeom->type, tree);
}

static RECT_NODE *
rect_node_from_lwpoly(const LWGEOM *lwgeom, RECT_TREE *tree)
{
	RECT_NODE **nodes;
	RECT_NODE *root;
	uint32_t i, j = 0;
	const LWPOLY *lwpoly = (const LWPOLY*)lwgeom;

	if (lwpoly->nrings < 1)
		return NULL;

	nodes = lwalloc(sizeof(RECT_NODE*) * lwpoly->nrings);
	for (i = 0; i < lwpoly->nrings; i++)
	{
		RECT_NODE *node = rect_node_from_ptarray(lwpoly->rings[i], lwgeom->type, tree);
		if (node)
		{
			node->i.ring_type = i ? RECT_NODE_RING_INTERIOR : RECT_NODE_RING_EXTERIOR;
			nodes[j++] = node;
		}
	}
	/* Put the top nodes in a z-order curve for a spatially coherent */
	/* tree after node merge */
	qsort(nodes, j, sizeof(RECT_NODE*), rect_node_cmp);
	root = rect_nodes_merge(nodes, j, tree);
	root->geom_type = lwgeom->type;
	lwfree(nodes);
	return root;
}

/* Prototype */
static RECT_NODE * rect_node_from_lwgeom(const LWGEOM *lwgeom, RECT_TREE *tree);

static RECT_NODE *
rect_node_from_lwcurvepoly(const LWGEOM *lwgeom, RECT_TREE *tree)
{
	RECT_NODE **nodes;
	RECT_NODE *root;
	uint32_t i, j = 0;
	const LWCURVEPOLY *lwcol = (const LWCURVEPOLY*)lwgeom;

	if (lwcol->nrings < 1)
		return NULL;

	nodes = lwalloc(sizeof(RECT_NODE*) * lwcol->nrings);
	for (i = 0; i < lwcol->nrings; i++)
	{
		RECT_NODE *node = rect_node_from_lwgeom(lwcol->rings[i], tree);
		if (node)
		{
			/*
			* In the case of arc circle, it's possible for a ring to consist
			* of a single closed edge. That will arrive as a leaf node. We
			* need to wrap that node in an internal node with an appropriate
			* ring type so all the other code can try and make sense of it.
			*/
			if (node->type == RECT_NODE_LEAF_TYPE)
			{
				RECT_NODE *internal = rect_node_internal_new(node, tree);
				rect_node_internal_add_node(internal, node);
				node = internal;
			}
			/* Each subcomponent is a ring */
			node->i.ring_type = i ? RECT_NODE_RING_INTERIOR : RECT_NODE_RING_EXTERIOR;
			nodes[j++] = node;
		}
	}
	/* Put the top nodes in a z-order curve for a spatially coherent */
	/* tree after node merge */
	qsort(nodes, j, sizeof(RECT_NODE*), rect_node_cmp);

	root = rect_nodes_merge(nodes, j, tree);

	root->geom_type = lwgeom->type;
	lwfree(nodes);
	return root;
}

static RECT_NODE *
rect_node_from_lwcollection(const LWGEOM *lwgeom, RECT_TREE *tree)
{
	RECT_NODE **nodes;
	RECT_NODE *root;
	uint32_t i, j = 0;
	const LWCOLLECTION *lwcol = (const LWCOLLECTION*)lwgeom;

	if (lwcol->ngeoms < 1)
		return NULL;

	/* Build one tree for each sub-geometry, then below */
	/* we merge the root notes of those trees to get a single */
	/* top node for the collection */
	nodes = lwalloc(sizeof(RECT_NODE*) * lwcol->ngeoms);
	for (i = 0; i < lwcol->ngeoms; i++)
	{
		RECT_NODE *node = rect_node_from_lwgeom(lwcol->geoms[i], tree);
		if (node)
		{
			/* Curvepolygons are collections where the sub-geometries */
			/* are the rings, and will need to doint point-in-poly */
			/* tests in order to do intersects and distance calculations */
			/* correctly */
			if (lwgeom->type == CURVEPOLYTYPE)
				node->i.ring_type = i ? RECT_NODE_RING_INTERIOR : RECT_NODE_RING_EXTERIOR;
			nodes[j++] = node;
		}
	}
	/* Sort the nodes using a z-order curve, so that merging the nodes */
	/* gives a spatially coherent tree (near things are in near nodes) */
	/* Note: CompoundCurve has edges already spatially organized, no */
	/* sorting needed */
	if (lwgeom->type != COMPOUNDTYPE)
		qsort(nodes, j, sizeof(RECT_NODE*), rect_node_cmp);

	root = rect_nodes_merge(nodes, j, tree);

	root->geom_type = lwgeom->type;
	lwfree(nodes);
	return root;
}

static RECT_NODE *
rect_node_from_lwgeom(const LWGEOM *lwgeom, RECT_TREE *tree)
{
	switch(lwgeom->type)
	{
		case POINTTYPE:
			return rect_node_from_lwpoint(lwgeom, tree);
		case TRIANGLETYPE:
		case CIRCSTRINGTYPE:
		case LINETYPE:
			return rect_node_from_lwline(lwgeom, tree);
		case POLYGONTYPE:
			return rect_node_from_lwpoly(lwgeom, tree);
		case CURVEPOLYTYPE:
			return rect_node_from_lwcurvepoly(lwgeom, tree);
		case COMPOUNDTYPE:
		case MULTICURVETYPE:
		case MULTISURFACETYPE:
		case MULTIPOINTTYPE:
		case MULTILINETYPE:
		case MULTIPOLYGONTYPE:
		case POLYHEDRALSURFACETYPE:
		case TINTYPE:
		case COLLECTIONTYPE:
			return rect_node_from_lwcollection(lwgeom, tree);
		default:
			lwerror("%s: Unknown geometry type: %s", __func__, lwtype_name(lwgeom->type));
			return NULL;
	}
	return NULL;
}


RECT_TREE *
rect_tree_from_lwgeom(const LWGEOM *lwgeom)
{
	RECT_TREE *tree = lwalloc0(sizeof(RECT_TREE));
	uint32_t num_nodes = rect_node_count_from_lwgeom(lwgeom);
	if (num_nodes > 0)
	{
		/* Flat array of nodes, nice and localized */
		tree->buffer = lwalloc0(num_nodes * sizeof(RECT_NODE));
	}
	tree->capacity = num_nodes;
	tree->next = 0;
	tree->root = rect_node_from_lwgeom(lwgeom, tree);
	return tree;
}


static uint32_t
rect_node_count_nodes(uint32_t level_nodes)
{
	uint32_t num_nodes = 0;
	while(level_nodes > 1)
	{
		uint32_t next_level_nodes = level_nodes / RECT_NODE_SIZE;
		/* For now, always slightly over-allocate space */
		/* next_level_nodes += ((level_nodes % RECT_NODE_SIZE) > 0); */
		next_level_nodes += 1;
		num_nodes += level_nodes;
		level_nodes = next_level_nodes;
	}
	/* and the root node */
	return num_nodes + 1;
}

static uint32_t
rect_node_count_from_lwline(const LWGEOM *lwgeom)
{
	const LWLINE* line = (const LWLINE*)lwgeom;
	return rect_node_count_nodes(line->points->npoints - 1);
}

static uint32_t
rect_node_count_from_lwpoly(const LWGEOM *lwgeom)
{
	const LWPOLY* poly = (const LWPOLY*)lwgeom;
	uint32_t num_nodes = rect_node_count_nodes(poly->nrings);
	for (uint32_t i = 0; i < poly->nrings; i++)
	{
		num_nodes += rect_node_count_nodes(poly->rings[i]->npoints - 1);
	}
	return num_nodes;
}

static uint32_t
rect_node_count_from_lwcollection(const LWGEOM *lwgeom)
{
	const LWCOLLECTION* col = (const LWCOLLECTION*)lwgeom;
	uint32_t num_nodes = rect_node_count_nodes(col->ngeoms);
	for (uint32_t i = 0; i < col->ngeoms; i++)
	{
		num_nodes += rect_node_count_from_lwgeom(col->geoms[i]);
	}
	return num_nodes;
}

uint32_t
rect_node_count_from_lwgeom(const LWGEOM *lwgeom)
{
	switch(lwgeom->type)
	{
		case POINTTYPE:
			return 1;
		case TRIANGLETYPE:
		case CIRCSTRINGTYPE:
		case LINETYPE:
			return rect_node_count_from_lwline(lwgeom);
		case POLYGONTYPE:
			return rect_node_count_from_lwpoly(lwgeom);
		case CURVEPOLYTYPE:
		case COMPOUNDTYPE:
		case MULTICURVETYPE:
		case MULTISURFACETYPE:
		case MULTIPOINTTYPE:
		case MULTILINETYPE:
		case MULTIPOLYGONTYPE:
		case POLYHEDRALSURFACETYPE:
		case TINTYPE:
		case COLLECTIONTYPE:
			return rect_node_count_from_lwcollection(lwgeom);
		default:
			lwerror("%s: Unknown geometry type: %s", __func__, lwtype_name(lwgeom->type));
			return 0;
	}
	return 0;
}

/*
* Get an actual coordinate point from a tree to use
* for point-in-polygon testing.
*/
static const POINT2D *
rect_node_get_point(const RECT_NODE *node)
{
	if (!node) return NULL;
	if (rect_node_is_leaf(node))
		return getPoint2d_cp(node->l.pa, 0);
	else
		return rect_node_get_point(node->i.nodes[0]);
}


static inline int
rect_node_intersects(const RECT_NODE *n1, const RECT_NODE *n2)
{
	if (n1->xmin > n2->xmax || n2->xmin > n1->xmax ||
		n1->ymin > n2->ymax || n2->ymin > n1->ymax)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}


#if POSTGIS_DEBUG_LEVEL >= 4
static char *
rect_node_to_str(const RECT_NODE *n)
{
	char *buf = lwalloc(256);
	snprintf(buf, 256, "(%.9g %.9g,%.9g %.9g)",
		n->xmin, n->ymin, n->xmax, n->ymax);
   	return buf;
}
#endif

/*
* Work down to leaf nodes, until we find a pair of leaf nodes
* that intersect. Prune branches that do not intersect.
*/
static int
rect_node_intersects_node_recursive(const RECT_NODE *n1, const RECT_NODE *n2)
{
	uint32_t i, j;
#if POSTGIS_DEBUG_LEVEL >= 4
	char *n1_str = rect_node_to_str(n1);
	char *n2_str = rect_node_to_str(n2);
	LWDEBUGF(4,"n1 %s  n2 %s", n1, n2);
	lwfree(n1_str);
	lwfree(n2_str);
#endif
	/* There can only be an edge intersection if the rectangles overlap */
	if (rect_node_intersects(n1, n2))
	{
		LWDEBUG(4," interaction found");
		/* We can only test for a true intersection if the nodes are both leaf nodes */
		if (rect_node_is_leaf(n1) && rect_node_is_leaf(n2))
		{
			LWDEBUG(4,"  leaf node test");
			/* Check for true intersection */
			return rect_leaf_node_intersects(&n1->l, &n2->l);
		}
		else if (rect_node_is_leaf(n2) && !rect_node_is_leaf(n1))
		{
			for (i = 0; i < n1->i.num_nodes; i++)
			{
				if (rect_node_intersects_node_recursive(n1->i.nodes[i], n2))
					return LW_TRUE;
			}
		}
		else if (rect_node_is_leaf(n1) && !rect_node_is_leaf(n2))
		{
			for (i = 0; i < n2->i.num_nodes; i++)
			{
				if (rect_node_intersects_node_recursive(n2->i.nodes[i], n1))
					return LW_TRUE;
			}
		}
		else
		{
			for (j = 0; j < n1->i.num_nodes; j++)
			{
				for (i = 0; i < n2->i.num_nodes; i++)
				{
					if (rect_node_intersects_node_recursive(n2->i.nodes[i], n1->i.nodes[j]))
						return LW_TRUE;
				}
			}
		}
	}
	return LW_FALSE;
}


int
rect_tree_intersects_tree(const RECT_TREE *t1, const RECT_TREE *t2)
{
	const RECT_NODE *n1 = t1->root;
	const RECT_NODE *n2 = t2->root;

	/*
	* It is possible for an area to intersect another object
	* without any edges intersecting, if the object is fully contained.
	* If that is so, then any point in the object will be contained,
	* so we do a quick point-in-poly test first for those cases
	*/
	if (rect_node_is_area(n1) &&
		rect_node_contains_point(n1, rect_node_get_point(n2)))
	{
		return LW_TRUE;
	}

	if (rect_node_is_area(n2) &&
		rect_node_contains_point(n2, rect_node_get_point(n1)))
	{
		return LW_TRUE;
	}

	/*
	* Not contained, so intersection can only happen if
	* edges actually intersect.
	*/
	return rect_node_intersects_node_recursive(n1, n2);
}

static inline double
distance(double x1, double y1, double x2, double y2)
{
	double dx = x1-x2;
	double dy = y1-y2;
	return sqrt(dx*dx + dy*dy);
}

/*
* The closest any two objects in two nodes can be is the smallest
* distance between the nodes themselves.
*/
static inline double
rect_node_min_distance(const RECT_NODE *n1, const RECT_NODE *n2)
{
	int   left = n1->xmin > n2->xmax;
	int  right = n1->xmax < n2->xmin;
	int bottom = n1->ymin > n2->ymax;
	int    top = n1->ymax < n2->ymin;

	//lwnotice("rect_node_min_distance");

	if (top && left)
		return distance(n1->xmin, n1->ymax, n2->xmax, n2->ymin);
	else if (top && right)
		return distance(n1->xmax, n1->ymax, n2->xmin, n2->ymin);
	else if (bottom && left)
		return distance(n1->xmin, n1->ymin, n2->xmax, n2->ymax);
	else if (bottom && right)
		return distance(n1->xmax, n1->ymin, n2->xmin, n2->ymax);
	else if (left)
		return n1->xmin - n2->xmax;
	else if (right)
		return n2->xmin - n1->xmax;
	else if (bottom)
		return n1->ymin - n2->ymax;
	else if (top)
		return n2->ymin - n1->ymax;
	else
		return 0.0;

	return 0.0;
}

/*
* The furthest apart the objects in two nodes can be is if they
* are at opposite corners of the bbox that contains both nodes
* TODO, this is not quite right. For the fully contained case, the
* maximum distance is actually smaller than the largest box.
*/
// static inline double
// rect_node_max_distance(const RECT_NODE *n1, const RECT_NODE *n2)
// {
// 	double d11 =
// 	double xmin = FP_MIN(n1->xmin, n2->xmin);
// 	double ymin = FP_MIN(n1->ymin, n2->ymin);
// 	double xmax = FP_MAX(n1->xmax, n2->xmax);
// 	double ymax = FP_MAX(n1->ymax, n2->ymax);
// 	double dx = xmax - xmin;
// 	double dy = ymax - ymin;
// 	//lwnotice("rect_node_max_distance");
// 	return sqrt(dx*dx + dy*dy);
// }

/*
* The maximum distance between opposing corners seems to
* capture the "maximum distance a child node could be
* from another child node" better than the old function.
*/
static inline double
rect_node_max_distance(const RECT_NODE *n1, const RECT_NODE *n2)
{
	double d00 = distance(n1->xmin, n1->ymin, n2->xmax, n2->ymax);
	double d11 = distance(n1->xmax, n1->ymax, n2->xmin, n2->ymin);
	double d01 = distance(n1->xmin, n1->ymax, n2->xmax, n2->ymin);
	double d10 = distance(n1->xmax, n1->ymin, n2->xmin, n2->ymax);
	return FP_MAX(FP_MAX(d00, d11), FP_MAX(d10, d01));
}

/*
* Leaf nodes represent individual edges from the original shape.
* As such, they can be either points (if original was a (multi)point)
* two-point straight edges (for linear types), or
* three-point circular arcs (for curvilinear types).
* The distance calculations for each possible combination of primitive
* edges is different, so there's a big case switch in here to match
* up the right combination of inputs to the right distance calculation.
*/
static double
rect_leaf_node_distance(const RECT_NODE_LEAF *n1, const RECT_NODE_LEAF *n2, RECT_TREE_QUERY_STATE *state)
{
	const POINT2D *p1, *p2, *p3, *q1, *q2, *q3;
	DISTPTS dl;

	//lwnotice("rect_leaf_node_distance, %d<->%d", n1->seg_num, n2->seg_num);

	lw_dist2d_distpts_init(&dl, DIST_MIN);

	switch (n1->seg_type)
	{
		case RECT_NODE_SEG_POINT:
		{
			p1 = getPoint2d_cp(n1->pa, n1->seg_num);

			switch (n2->seg_type)
			{
				case RECT_NODE_SEG_POINT:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					lw_dist2d_pt_pt(q1, p1, &dl);
					break;

				case RECT_NODE_SEG_LINEAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num+1);
					lw_dist2d_pt_seg(p1, q1, q2, &dl);
					break;

				case RECT_NODE_SEG_CIRCULAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num*2);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num*2+1);
					q3 = getPoint2d_cp(n2->pa, n2->seg_num*2+2);
					lw_dist2d_pt_arc(p1, q1, q2, q3, &dl);
					break;

				default:
					lwerror("%s: unsupported segment type", __func__);
			}
			break;
		}

		case RECT_NODE_SEG_LINEAR:
		{
			p1 = getPoint2d_cp(n1->pa, n1->seg_num);
			p2 = getPoint2d_cp(n1->pa, n1->seg_num+1);

			switch (n2->seg_type)
			{
				case RECT_NODE_SEG_POINT:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					lw_dist2d_pt_seg(q1, p1, p2, &dl);
					break;

				case RECT_NODE_SEG_LINEAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num+1);
					lw_dist2d_seg_seg(q1, q2, p1, p2, &dl);
					// lwnotice(
					// 	"%d\tLINESTRING(%g %g,%g %g)\t%d\tLINESTRING(%g %g,%g %g)\t%g\t%g\t%g",
					// 	n1->seg_num,
					// 	p1->x, p1->y, p2->x, p2->y,
					// 	n2->seg_num,
					// 	q1->x, q1->y, q2->x, q2->y,
					// 	dl.distance, state->min_dist, state->max_dist);
					break;

				case RECT_NODE_SEG_CIRCULAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num*2);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num*2+1);
					q3 = getPoint2d_cp(n2->pa, n2->seg_num*2+2);
					lw_dist2d_seg_arc(p1, p2, q1, q2, q3, &dl);
					break;

				default:
					lwerror("%s: unsupported segment type", __func__);
			}
			break;
		}
		case RECT_NODE_SEG_CIRCULAR:
		{
			p1 = getPoint2d_cp(n1->pa, n1->seg_num*2);
			p2 = getPoint2d_cp(n1->pa, n1->seg_num*2+1);
			p3 = getPoint2d_cp(n1->pa, n1->seg_num*2+2);

			switch (n2->seg_type)
			{
				case RECT_NODE_SEG_POINT:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					lw_dist2d_pt_arc(q1, p1, p2, p3, &dl);
					break;

				case RECT_NODE_SEG_LINEAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num+1);
					lw_dist2d_seg_arc(q1, q2, p1, p2, p3, &dl);
					break;

				case RECT_NODE_SEG_CIRCULAR:
					q1 = getPoint2d_cp(n2->pa, n2->seg_num*2);
					q2 = getPoint2d_cp(n2->pa, n2->seg_num*2+1);
					q3 = getPoint2d_cp(n2->pa, n2->seg_num*2+2);
					lw_dist2d_arc_arc(p1, p2, p3, q1, q2, q3, &dl);
					break;

				default:
					lwerror("%s: unsupported segment type", __func__);
			}
			break;
		}
		default:
			lwerror("%s: unsupported segment type", __func__);
	}

	/* If this is a new global minima, save it */
	if (dl.distance < state->min_dist)
	{
		state->min_dist = dl.distance;
		state->p1 = dl.p1;
		state->p2 = dl.p2;
	}

	return dl.distance;
}


static double
rect_node_distance_node_recursive(const RECT_NODE *n1, const RECT_NODE *n2, RECT_TREE_QUERY_STATE *state)
{
	double min, max;

	/* Short circuit if we've already hit the minimum */
	if (state->min_dist < state->threshold || state->min_dist == 0.0)
		return state->min_dist;

	/* If your minimum is greater than anyone's maximum, you can't hold the winner */
	min = rect_node_min_distance(n1, n2);
	if (min > state->max_dist)
	{
		//lwnotice("pruning pair %p, %p", n1, n2);
		LWDEBUGF(4, "pruning pair %p, %p", n1, n2);
		return FLT_MAX;
	}

	/* If your maximum is a new low, we'll use that as our new global tolerance */
	max = rect_node_max_distance(n1, n2);
	if (max < state->max_dist)
		state->max_dist = max;

	/* Both leaf nodes, do a real distance calculation */
	if (rect_node_is_leaf(n1) && rect_node_is_leaf(n2))
	{
		return rect_leaf_node_distance(&n1->l, &n2->l, state);
	}
	/* Recurse into nodes */
	else
	{
		uint32_t i, j;
		double d_min = FLT_MAX;
		if (rect_node_is_leaf(n1) && !rect_node_is_leaf(n2))
		{
			for (i = 0; i < n2->i.num_nodes; i++)
			{
				min = rect_node_distance_node_recursive(n1, n2->i.nodes[i], state);
				d_min = FP_MIN(d_min, min);
			}
		}
		else if (rect_node_is_leaf(n2) && !rect_node_is_leaf(n1))
		{
			for (i = 0; i < n1->i.num_nodes; i++)
			{
				min = rect_node_distance_node_recursive(n1->i.nodes[i], n2, state);
				d_min = FP_MIN(d_min, min);
			}
		}
		else
		{
			for (i = 0; i < n1->i.num_nodes; i++)
			{
				for (j = 0; j < n2->i.num_nodes; j++)
				{
					min = rect_node_distance_node_recursive(n1->i.nodes[i], n2->i.nodes[j], state);
					d_min = FP_MIN(d_min, min);
				}
			}
		}
		return d_min;
	}
}



//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

typedef struct node_pair_t
{
	pqueue_pri_t pri;
	size_t pos;
	double min_dist;
	double max_dist;
	const RECT_NODE *n1;
	const RECT_NODE *n2;
} node_pair_t;

typedef struct node_pair_list_t
{
	size_t size;
	size_t capacity;
	node_pair_t *items;
} node_pair_list_t;

node_pair_list_t node_pair_list;

static size_t
get_pos(void *node_pair)
{
	return ((node_pair_t*)node_pair)->pos;
}

static void
set_pos(void *node_pair, size_t pos)
{
	((node_pair_t*)node_pair)->pos = pos;
}

static pqueue_pri_t
get_pri(void *node_pair)
{
	return ((node_pair_t*)node_pair)->pri;
}

static void
set_pri(void *node_pair, pqueue_pri_t pri)
{
	((node_pair_t*)node_pair)->pri = pri;
}

static int
cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
	/* Smallest has most priority */
	return next > curr;
}

static node_pair_t *
rect_node_distance_node_pair(const RECT_NODE *n1, const RECT_NODE *n2)
{
	node_pair_t *node_pair;
	if (node_pair_list.size >= node_pair_list.capacity)
	{
		node_pair_list.capacity *= 2;
		node_pair_list.items = lwrealloc(node_pair_list.items, node_pair_list.capacity);
	}
	node_pair = &(node_pair_list.items[node_pair_list.size++]);
	memset(node_pair, 0, sizeof(node_pair_t));

	node_pair->min_dist = rect_node_min_distance(n1, n2);
	node_pair->max_dist = rect_node_max_distance(n1, n2);
	node_pair->n1 = n1;
	node_pair->n2 = n2;
	node_pair->pri = node_pair->min_dist;
	return node_pair;
}

static double
rect_node_distance_node_queue(const RECT_NODE *root1, const RECT_NODE *root2, RECT_TREE_QUERY_STATE *state)
{
	const node_pair_t *nodePair;
	size_t i, j;

	pqueue_t *queue = pqueue_init(RECT_TREE_QUEUE_SIZE, cmp_pri, get_pri, set_pri, get_pos, set_pos);

	/* Add the first node pair to prime the queue */
	pqueue_insert(queue, rect_node_distance_node_pair(root1, root2));

	while ((nodePair = pqueue_pop(queue)))
	{
		const RECT_NODE *n1 = nodePair->n1;
		const RECT_NODE *n2 = nodePair->n2;

		bool n1_is_leaf = rect_node_is_leaf(n1);
		bool n2_is_leaf = rect_node_is_leaf(n2);

		/* Short circuit if we've already hit the minimum */
		if (state->min_dist < state->threshold || state->min_dist == 0.0)
		{
			break;
		}

		/* If your minimum distance is greater than anyone's, */
		/* maximum distance you cannot hold the nearest pair of leaves */
		if (nodePair->min_dist > state->max_dist)
		{
			//lwnotice("pruning pair %p, %p", n1, n2);
			LWDEBUGF(4, "pruning pair %p, %p", n1, n2);
			continue;
		}

		/* If your maximum distance is a new low, */
		/* we'll use that as our new global tolerance */
		if (nodePair->max_dist < state->max_dist)
			state->max_dist = nodePair->max_dist;

		/* If both nodes are leaves, do a real distance calculation */
		if (n1_is_leaf && n2_is_leaf)
		{
			/* This updates the state nearest points and min_dist */
			/* if a new minimum is found */
			rect_leaf_node_distance(&n1->l, &n2->l, state);
			continue;
		}
		/* If one is a leaf, add node pair combos to the queue */
		else if (n2_is_leaf && !n1_is_leaf)
		{
			for (i = 0; i < n1->i.num_nodes; i++)
				pqueue_insert(queue, rect_node_distance_node_pair(n1->i.nodes[i], n2));
		}
		/* If one is a leaf, add node pair combos to the queue */
		else if (n1_is_leaf && !n2_is_leaf)
		{
			for (i = 0; i < n2->i.num_nodes; i++)
				pqueue_insert(queue, rect_node_distance_node_pair(n1, n2->i.nodes[i]));
		}
		/* Both nodes are internal, add all child combos to the queue */
		else
		{
			for (i = 0; i < n1->i.num_nodes; i++)
				for (j = 0; j < n2->i.num_nodes; j++)
					pqueue_insert(queue, rect_node_distance_node_pair(n1->i.nodes[i], n2->i.nodes[j]));
		}

		/* Go back and pop another off the top */
	}

	return state->min_dist;
}



double rect_tree_distance_tree(const RECT_TREE *t1, const RECT_TREE *t2, double threshold)
{
	double distance;

	const RECT_NODE *n1 = t1->root;
	const RECT_NODE *n2 = t2->root;

	RECT_TREE_QUERY_STATE state;

	/*
	* It is possible for an area to intersect another object
	* without any edges intersecting, if the object is fully contained.
	* If that is so, then any point in the object will be contained,
	* so we do a quick point-in-poly test first for those cases
	*/
	if (rect_node_is_area(n1) &&
	    rect_node_contains_point(n1, rect_node_get_point(n2)))
	{
		return 0.0;
	}

	if (rect_node_is_area(n2) &&
	    rect_node_contains_point(n2, rect_node_get_point(n1)))
	{
		return 0.0;
	}

	/* Initialize state */
	state.threshold = threshold;
	state.min_dist = FLT_MAX;
	state.max_dist = FLT_MAX;

	/* Initialize state */
	node_pair_list.capacity = RECT_TREE_QUEUE_SIZE;
	node_pair_list.size = 0;
	node_pair_list.items = lwalloc0(node_pair_list.capacity * sizeof(node_pair_t));

	distance = rect_node_distance_node_queue(n1, n2, &state);
	// distance = rect_node_distance_node_recursive(n1, n2, &state);

	lwfree(node_pair_list.items);

	// *p1 = state.p1;
	// *p2 = state.p2;
	return distance;
}




