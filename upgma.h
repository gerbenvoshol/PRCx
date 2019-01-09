/* This a mix of tree functions, data-structures and including 
 * the fast UPGMA algorithm (O(N^2)) as implemented in Bob Edgar's
 * Muscle (version 3.7) ported to pure C. This version is adapted 
 * from the C version included with clustal-Omega.
 *
 * Muscle's code is public domain and so is this code here.

 * From http://www.drive5.com/muscle/license.htm:
 * """
 * MUSCLE is public domain software
 *
 * The MUSCLE software, including object and source code and
 * documentation, is hereby donated to the public domain.
 *
 * Disclaimer of warranty
 * THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * """
 *
 */

#ifndef __UPGMA_H__
#define __UPGMA_H__

#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <stdbool.h>

#ifndef uint
/* limit use of uint (see coding_style_guideline.txt) */
typedef unsigned int uint;
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/* type boolean and false and true defined in stdbool.h */
#ifndef TRUE
#define TRUE true
#endif
#ifndef FALSE
#define FALSE false
#endif

static const uint NULL_NEIGHBOR = UINT_MAX;

/**
 * @brief guide-tree structure
 *     
 * @note We kept the original variable names here, to make it easy to
 * search through Muscle's source code.
 * From phy.cpp:
 *    Node has 0 to 3 neighbors:
 *    0 neighbors: singleton root
 *    1 neighbor:  leaf, neighbor is parent
 *    2 neigbors:  non-singleton root
 *    3 neighbors: internal node (other than root)
 *    
 *    Minimal rooted tree is single node.
 *    Minimal unrooted tree is single edge.
 *    Leaf node always has nulls in neighbors 2 and 3, neighbor 1 is parent.
 *    When tree is rooted, neighbor 1=parent, 2=left, 3=right.
 *
 */
typedef struct {
    uint m_uNodeCount;/**< number of nodes */
    uint m_uCacheCount;/**< reserved memory */
    
    uint *m_uNeighbor1;/**< parent node */
    uint *m_uNeighbor2;/**< left node */
    uint *m_uNeighbor3;/**< right node */
    
    /* do we have edge lengths info stored (m_dEdgeLength[123]) */
    bool *m_bhas_edge_length1;
    bool *m_bhas_edge_length2;
    bool *m_bhas_edge_length3;

    double *m_dEdgeLength1;
    double *m_dEdgeLength2;
    double *m_dEdgeLength3;
    
#if USE_HEIGHT
    /* unused in our version of the code. we might need it at some
     * stage so keep it in here, but disable via USE_HEIGHT throughout
     * the code */
    double *m_dHeight;
    bool *m_bHasHeight;
#endif

    /**
     * leaf labels.
     * index range: 0 -- (m_uNodeCount+1)/2
     */
    char **m_ptrName;
    
    /**
     * node id.
     * index range: 0 -- m_uNodeCount
     */
    uint *m_Ids;
    
    bool m_bRooted; /**< tree is rooted */
    uint m_uRootNodeIndex;
} tree_t;

enum linkage_e {
    LINKAGE_MIN = 0,
    LINKAGE_AVG,
    LINKAGE_MAX,
};
typedef enum linkage_e linkage_t;

typedef struct {
    int nrows; /* number of rows    */
    int ncols; /* number of columns */
    double **data; 
} matrix_t;

extern void
tree_create(tree_t *tree, uint uLeafCount, uint uRoot, const uint *Left,
           const uint  *Right, const double *LeftLength, const double* RightLength,
           const uint *LeafIds, char **LeafNames);

extern void
save_tree(FILE *fp, tree_t *tree);

extern void
save_distmat(FILE *fp, matrix_t *distmat, char **labels);

extern int
load_tree(tree_t *tree, char *ftree);

extern void
free_tree(tree_t *tree);

extern void
log_tree(tree_t *tree, FILE *fp);

extern bool
is_rooted(tree_t *tree);

extern uint
get_node_count(tree_t *tree);

extern uint
get_leaf_count(tree_t *tree);
        
extern uint
first_depth_first_node(tree_t *tree);

extern uint
next_depth_first_node(uint nodeindex, tree_t *tree);

extern bool
is_leaf(uint nodeindex, tree_t *tree);

extern void
set_leaf_ID(tree_t *tree, uint uNodeIndex, uint uId);
    
extern uint
get_leaf_ID(uint nodeindex, tree_t *tree);

extern char *
get_leaf_name(unsigned uNodeIndex, tree_t *tree);

extern uint
get_left(uint nodeindex, tree_t *tree);

extern uint
get_right(uint nodeindex, tree_t *tree);

extern uint
get_root_node_index(tree_t *tree);

extern bool
is_root(uint uNodeIndex, tree_t *tree);

extern uint
get_parent(unsigned uNodeIndex, tree_t *tree);

extern double
get_edge_length(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree);

extern uint
leaf_index_to_node_index(uint uLeafIndex, tree_t *prTree);

extern void
append_tree(tree_t *prDstTree,
          uint uDstTreeNodeIndex, tree_t *prSrcTree);

extern void
validate_tree(tree_t *tree);

extern void prc_upgma(tree_t *tree, matrix_t *distmat,
                   linkage_t linkage, char **names);

/* Allocates a matrix of rowsize*colsize with a single call to malloc, so
 * only a single call to free is necessary
 */
extern double **prc_allocmat(int rowsize, int colsize);

extern int prc_newmat(matrix_t **distmat, int nrows, int ncols);

extern int prc_freemat(matrix_t *distmat);

/* Reads an entire file into an array of strings, needs only a single call to free */
extern char **prc_stringfile(char *filename, size_t *number_of_lines);

#endif
