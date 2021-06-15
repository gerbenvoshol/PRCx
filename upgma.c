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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <ctype.h>

#include "upgma.h"

/*const double dInsane = VERY_NEGATIVE_DOUBLE;*/
static const double dInsane = -9e29;
static const unsigned uInsane = 8888888;

typedef enum 
{
    NTT_Unknown,

    /* Returned from Tree::get_token: */
    NTT_Lparen,
    NTT_Rparen,
    NTT_Colon,
    NTT_Comma,
    NTT_Semicolon,
    NTT_String,
    
    /* Following are never returned from Tree::get_token: */
    NTT_SingleQuotedString,
    NTT_DoubleQuotedString,
    NTT_Comment
} NEWICK_TOKEN_TYPE;


static void
init_cache(uint uCacheCount, tree_t *tree);
static void
zero_tree(tree_t *tree);
static uint
get_neighbor_count(unsigned uNodeIndex, tree_t *tree);
static bool
is_edge(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree);
static bool
has_edge_length(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree);
static void
tree_to_file_node_rooted(tree_t *tree, uint m_uRootNodeIndex, FILE *fp);
static void
validate_node(uint uNodeIndex, tree_t *tree);
static void
assert_are_neighbors(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree);
static void
expand_cache(tree_t *tree);
static void
tree_create_rooted(tree_t *tree);
static bool
get_group_from_file(FILE *fp, uint uNodeIndex, double *ptrdEdgeLength, tree_t *tree);
static NEWICK_TOKEN_TYPE
get_token(FILE *fp, char szToken[], uint uBytes);
/* stuff from textfile.cpp */
static void
file_skip_white(FILE *fp);
static bool
file_skip_whitex(FILE *fp);
static void
set_leaf_name(uint uNodeIndex, const char *ptrName, tree_t *tree);
uint
append_branch(tree_t *tree, uint uExistingLeafIndex);
static void
set_edge_length(uint uNodeIndex1, uint uNodeIndex2,
              double dLength, tree_t *tree);
static uint
unroot_from_file(tree_t *tree);
uint
get_neighbor(uint uNodeIndex, uint uNeighborSubscript, tree_t *tree);
static void
init_node(tree_t *prTree, uint uNodeIndex);


/* Allocates a matrix of rowsize*colsize with a single call to malloc, so
 * only a single call to free is necessary
 */
double **prc_allocmat(int rowsize, int colsize)
{
    double **matrix = NULL;
    int item_size = sizeof(double);
    int datasize = colsize * rowsize * item_size;
    int ptrsize = rowsize;
    double *rowptr;

    int totalsize = ptrsize * sizeof(double *) + datasize;

    if((matrix = (double **) calloc(totalsize, sizeof(double))) == NULL) {
        fprintf(stderr, "Can't allocate %d bytes\n.", totalsize);
        return NULL;
    }

    rowptr = (double *)(matrix + ptrsize);

    /* set pointers */
    for (int i = 0; i < rowsize; i++) {
        (matrix)[i] = rowptr + i * colsize * item_size;
    }

    return matrix;
}

int prc_newmat(matrix_t **distmat, int nrows, int ncols)
{   
    /* only symmetrical matrixes are allowed */
    assert(nrows > 0);
    assert(ncols == nrows);

    (*distmat) = (matrix_t *) malloc(1 * sizeof(matrix_t));
    (*distmat)->nrows = nrows;
    (*distmat)->ncols = ncols;
    (*distmat)->data = prc_allocmat(nrows, ncols);
    
    return 0;
}

int prc_freemat(matrix_t *distmat)
{
    free(distmat->data);
    free(distmat);

    return 0;
}

/* Reads an entire file into an array of strings, needs only a single call to free */
char **prc_stringfile(char *filename, size_t *number_of_lines)
{
    size_t count = 1;
    char **sfile = NULL;
    char *p;

    FILE *f = NULL;
    f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Unable to open: %s\n", filename);
        return NULL;
    }

    /* Determine the file size */
    fseek(f, 0, SEEK_END);
    size_t fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    /* Read the file into a temporary buffer */
    char *buffer = calloc(fsize + 1, sizeof(char));
    fread(buffer, fsize, 1, f);
    if (!buffer) {
        fprintf(stderr, "Failed to read %s into memory\n", filename);
        return NULL;
    }
    buffer[fsize] = 0;
    
    /* Close the file */
    fclose(f);

    /* Count the number of new lines */
    p = buffer;
    size_t i = 0;
    while (p[i]) {
        if (p[i] == '\n') {
            if ( p[i+1] == '\r') {
                count++;
                i++;
            } else {
                count++;
            }
        } else if (*p == '\r') {
            count++;
        }
        i++;
    }

    if (number_of_lines) {
        *number_of_lines = count;
    }

    /* Allocate space to keep the entire file */
    sfile = (char **) calloc(sizeof(char *) * (count + 1) + fsize + 1, 1);
    if (!sfile) {
        fprintf(stderr, "Could not copy the data\n");
        return NULL;
    }
    sfile[count] = NULL;
    /* Copy in the original data */
    memcpy(&sfile[count + 1], buffer, fsize + 1);
    
    free(buffer);
    buffer = (char *) &sfile[count + 1];

    /* Go over everything again and set the pointers */
    p = buffer;
    i = 0;
    count = 0;
    sfile[count] = &p[i];
    while (p[i]) {
        if (p[i] == '\n') {
            if ( p[i+1] == '\r') {
                p[i] = '\0';
                p[i+1] = '\0';
                count++;
                i++;
                if (p[i+1]) {
                    sfile[count] = &p[i+1];
                }
            } else {
                p[i] = '\0';
                count++;
                if (p[i+1]) {
                    sfile[count] = &p[i+1];
                }
            }
        } else if (*p == '\r') {
            p[i] = '\0';
            count++;
            if (p[i+1]) {
                sfile[count] = &p[i+1];
            }
        }
        i++;
    }

    return sfile;
}

/**
 * @param[in] uNodeIndex
 * node to examine
 * @param[in] tree
 * tree to examine
 * @return id of left node
 * @note called get_right in Muscle3.7
 */
uint
get_left(uint uNodeIndex, tree_t *tree) 
{
    assert(NULL != tree);
    assert(tree->m_bRooted && uNodeIndex < tree->m_uNodeCount);
    return tree->m_uNeighbor2[uNodeIndex];
}
/***   end: get_left   ***/



/**
 * @param[in] uNodeIndex
 * node to examine
 * @param[in] tree
 * tree to examine
 * @return id of right node
 * @note called get_right in Muscle3.7
 */
uint
get_right(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(tree->m_bRooted && uNodeIndex < tree->m_uNodeCount);
    return tree->m_uNeighbor3[uNodeIndex];
}
/***   end: get_right   ***/



/**
 * @param[in] uNodeIndex
 * node to examine
 * @param[in] tree
 * tree to examine
 * @return leaf id of current node
 */
uint
get_leaf_ID(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(is_leaf(uNodeIndex, tree));
    
    return tree->m_Ids[uNodeIndex];
}
/***   end: get_leaf_ID   ***/



/**
 * @note originally called get_leaf_name
 *
 */
char *
get_leaf_name(unsigned uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(is_leaf(uNodeIndex, tree));
    return tree->m_ptrName[uNodeIndex];
}
/***   end: get_leaf_name   ***/



/**
 * @brief returns first leaf node for a depth-first traversal of tree
 *
 * @param[in] tree
 * tree to traverse
 *
 * @return node index of first leaf node for depth-first traversal
 *
 * @note called first_depth_first_node in Muscle3.7
 *
 */
uint
first_depth_first_node(tree_t *tree)
{
    uint uNodeIndex;
        
    assert(NULL != tree);
    assert(is_rooted(tree));    
    
    /* Descend via left branches until we hit a leaf */
    uNodeIndex =  tree->m_uRootNodeIndex;
    while (!is_leaf(uNodeIndex, tree))
        uNodeIndex = get_left(uNodeIndex, tree);
    return uNodeIndex;
}
/***   end: first_depth_first_node   ***/


/**
 * @brief returns next leaf node index for depth-first traversal of
 * tree
 *
 * @param[in] tree
 * tree to traverse
 * @param[in] uNodeIndex
 * Current node index
 *
 * @return node index of next leaf node during depth-first traversal
 *
 * @note called next_depth_first_node in Muscle3.7
 */
uint
next_depth_first_node(uint uNodeIndex, tree_t *tree)
{
    uint uParent;
    
    assert(NULL != tree);
    assert(is_rooted(tree));
    assert(uNodeIndex < tree->m_uNodeCount);
    
    if (is_root(uNodeIndex, tree))
    {
        /* Node uNodeIndex is root, end of traversal */
        return NULL_NEIGHBOR;
    }
    
    uParent = get_parent(uNodeIndex, tree);
    if (get_right(uParent, tree) == uNodeIndex) {
        /* Is right branch, return parent=uParent */
        return uParent;
    }
    
    uNodeIndex = get_right(uParent, tree);
    /* Descend left from right sibling uNodeIndex */
    while (!is_leaf(uNodeIndex, tree)) {
        uNodeIndex = get_left(uNodeIndex, tree);
    }
    
    /*  bottom out at leaf uNodeIndex */
    return uNodeIndex;
}
/***   end: next_depth_first_node   ***/



/**
 * @brief check if tree is a rooted tree
 * @param[in] tree
 * tree to check
 * @return TRUE if given tree is rooted, FALSE otherwise
 *
 */
bool
is_rooted(tree_t *tree) 
{
    assert(NULL != tree);
    return tree->m_bRooted;
}
/***   end: is_rooted   ***/


/**
 *
 */
void
free_tree(tree_t *tree)
{
    uint i;

    assert(tree!=NULL);
    

    /* FIXME use GetLeafNodeIndex? or
       for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
       {
       if (tree.is_leaf(uNodeIndex))
       {
       const char *ptrName =
       tree.get_leaf_name(uNodeIndex);
    */
    /* is_leaf needs m_uNodeCount and all m_uNeighbor's
     * so free first
     */
    for (i=0; i<tree->m_uNodeCount; i++) {
        /* is_leaf needs neighbour count, so free those guys later */
        if (is_leaf(i, tree)) {
            free(tree->m_ptrName[i]);
        }
    }
    free(tree->m_ptrName);

    free(tree->m_uNeighbor1);
    free(tree->m_uNeighbor2);
    free(tree->m_uNeighbor3);
    
    free(tree->m_Ids);
    
    free(tree->m_dEdgeLength1);
    free(tree->m_dEdgeLength2);
    free(tree->m_dEdgeLength3);
#if USE_HEIGHT
    free(tree->m_dHeight);
    free(tree->m_bHasHeight);
#endif 
    free(tree->m_bhas_edge_length1);
    free(tree->m_bhas_edge_length2);
    free(tree->m_bhas_edge_length3);
    
    zero_tree(tree);

    //free(tree);
}
/***   end: free_tree   ***/



/**
 * @brief create a muscle tree
 *
 * @note Original comment in Muscle code: "Create rooted tree from a
 * vector description. Node indexes are 0..N-1 for leaves, N..2N-2 for
 * internal nodes. Vector subscripts are i-N and have values for
 * internal nodes only, but those values are node indexes 0..2N-2. So
 * e.g. if N=6 and Left[2]=1, this means that the third internal node
 * (node index 8) has the second leaf (node index 1) as its left
 * child. uRoot gives the vector subscript of the root, so add N to
 * get the node index."
 *
 * @param[out] tree
 * newly created tree
 * @param[in] uLeafCount
 * number of leaf nodes
 * @param[in] uRoot
 * Internal node index of root node
 * @param[in] Left
 * Node index of left sibling of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] Right
 * Node index of right sibling of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] LeftLength
 * Branch length of left branch of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] RightLength
 * Branch length of right branch of an internal node.
 * Index range: 0--uLeafCount-1
 * @param[in] LeafIds
 * Leaf id. Index range: 0--uLeafCount
 * @param[in] LeafNames
 * Leaf label. Index range: 0--uLeafCount
 *
 */
void
tree_create(tree_t *tree,
                  uint uLeafCount, uint uRoot,
                  const uint *Left, const uint *Right,
                  const double *LeftLength, const double* RightLength,
                  const uint *LeafIds, char **LeafNames)
{
    uint uNodeIndex;

	zero_tree(tree);
    tree->m_uNodeCount = 2*uLeafCount - 1;
	init_cache(tree->m_uNodeCount, tree);
    
	for (uNodeIndex = 0; uNodeIndex < uLeafCount; ++uNodeIndex) {
		tree->m_Ids[uNodeIndex] = LeafIds[uNodeIndex];
		tree->m_ptrName[uNodeIndex] = strdup(LeafNames[uNodeIndex]);
    }

	for (uNodeIndex = uLeafCount; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
		uint v = uNodeIndex - uLeafCount;
		uint uLeft = Left[v];
		uint uRight = Right[v];
		double fLeft = LeftLength[v];
		double fRight = RightLength[v];
        
		tree->m_uNeighbor2[uNodeIndex] = uLeft;
		tree->m_uNeighbor3[uNodeIndex] = uRight;
        
		tree->m_bhas_edge_length2[uNodeIndex] = TRUE;
		tree->m_bhas_edge_length3[uNodeIndex] = TRUE;

		tree->m_dEdgeLength2[uNodeIndex] = fLeft;
		tree->m_dEdgeLength3[uNodeIndex] = fRight;

		tree->m_uNeighbor1[uLeft] = uNodeIndex;
		tree->m_uNeighbor1[uRight] = uNodeIndex;

		tree->m_dEdgeLength1[uLeft] = fLeft;
		tree->m_dEdgeLength1[uRight] = fRight;

		tree->m_bhas_edge_length1[uLeft] = TRUE;
		tree->m_bhas_edge_length1[uRight] = TRUE;
    }

	tree->m_bRooted = TRUE;
	tree->m_uRootNodeIndex = uRoot + uLeafCount;
#ifndef NDEBUG
	validate_tree(tree);
#endif
}
/***   end: tree_create   ***/



/**
 * @param[in] tree
 * tree to write
 * @param[out] fp
 * filepointer to write to2
 * 
 * @brief write a muscle tree to a file in newick format (distances
 * and all names)
 *
 * 
 */
void
save_tree(FILE *fp, tree_t *tree)
{
    assert(NULL != tree);
    if (is_rooted(tree)) {
        tree_to_file_node_rooted(tree, tree->m_uRootNodeIndex, fp);
        fprintf(fp, ";\n");
        return;
    } else {
        printf("FIXME: output of unrooted muscle trees not implemented");
    }
}
/***   end: save_tree   ***/

extern void
save_distmat(FILE *fp, matrix_t *distmat, char **labels)
{
    int i, n;

    assert(NULL != distmat);
    assert(NULL != labels);
    fprintf(fp, "%d\n", distmat->nrows);

    for(i = 0; i <distmat->nrows; i++) {
        fprintf(fp, "%s", labels[i]);
        for(n = 0; n < distmat->ncols; n++) {
            fprintf(fp, "\t%f", distmat->data[i][n]);
        }
        fprintf(fp, "\n");
    }
}

/**
 * @brief check if given node is a leaf node
 *
 * @param[in] uNodeIndex
 * the node index
 * @param tree
 * the tree
 *
 * @return TRUE if given node is a leaf, FALSE otherwise
 *
 * @note called is_leaf in Muscle3.7. See tree.h in original code
 */
bool
is_leaf(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex < tree->m_uNodeCount);
    if (1 == tree->m_uNodeCount)
        return TRUE;
    return 1 == get_neighbor_count(uNodeIndex, tree);
}
/***   end: is_leaf   ***/



/**
 */
bool
is_root(uint uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    return is_rooted(tree) && tree->m_uRootNodeIndex == uNodeIndex;
}
/***   end: is_root   ***/



/**
 */
uint
get_neighbor_count(uint uNodeIndex, tree_t *tree)
{
    uint n1, n2, n3;
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(NULL != tree);
    assert(NULL != tree->m_uNeighbor1);
    assert(NULL != tree->m_uNeighbor2);
    assert(NULL != tree->m_uNeighbor3);
    n1 = tree->m_uNeighbor1[uNodeIndex];
    n2 = tree->m_uNeighbor2[uNodeIndex];
    n3 = tree->m_uNeighbor3[uNodeIndex];
    return (NULL_NEIGHBOR != n1) + (NULL_NEIGHBOR != n2) + (NULL_NEIGHBOR != n3);
}
/***   end: get_neighbor_count   ***/


/**
 */
uint
get_parent(unsigned uNodeIndex, tree_t *tree)
{
    assert(NULL != tree);
    assert(tree->m_bRooted && uNodeIndex < tree->m_uNodeCount);
    return tree->m_uNeighbor1[uNodeIndex];
}
/***   end: get_parent   ***/



/**
 */
bool
is_edge(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree)
{
    assert(uNodeIndex1 < tree->m_uNodeCount && uNodeIndex2 < tree->m_uNodeCount);
    assert(NULL != tree);
    
    return tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2 ||
        tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2 ||
        tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2;
}
/***   end: is_edge   ***/



/**
 */
bool
has_edge_length(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex1 < tree->m_uNodeCount);
    assert(uNodeIndex2 < tree->m_uNodeCount);
    assert(is_edge(uNodeIndex1, uNodeIndex2, tree));
    
    if (tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
        return tree->m_bhas_edge_length1[uNodeIndex1];
    else if (tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
        return tree->m_bhas_edge_length2[uNodeIndex1];
    assert(tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
    return tree->m_bhas_edge_length3[uNodeIndex1];
}
/***   end:   ***/


/**
 */
double
get_edge_length(uint uNodeIndex1, uint uNodeIndex2, tree_t *tree)
{
    assert(NULL != tree);
    assert(uNodeIndex1 < tree->m_uNodeCount && uNodeIndex2 < tree->m_uNodeCount);
    if (!has_edge_length(uNodeIndex1, uNodeIndex2, tree))
    {
        printf("Missing edge length in tree %u-%u", uNodeIndex1, uNodeIndex2);
    }
    
    if (tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
        return tree->m_dEdgeLength1[uNodeIndex1];
    else if (tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
        return tree->m_dEdgeLength2[uNodeIndex1];
    assert(tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
    return tree->m_dEdgeLength3[uNodeIndex1];
}
/***   end: get_edge_length   ***/


/**
 *
 */
void
init_node(tree_t *prTree, uint uNodeIndex)
{
    prTree->m_uNeighbor1[uNodeIndex] = NULL_NEIGHBOR;
    prTree->m_uNeighbor2[uNodeIndex] = NULL_NEIGHBOR;
    prTree->m_uNeighbor3[uNodeIndex] = NULL_NEIGHBOR;
    prTree->m_bhas_edge_length1[uNodeIndex] = FALSE;
    prTree->m_bhas_edge_length2[uNodeIndex] = FALSE;
    prTree->m_bhas_edge_length3[uNodeIndex] = FALSE;
#if USE_HEIGHT
    prTree->m_bHasHeight[uNodeIndex] = FALSE;
    prTree->m_dHeight[uNodeIndex] = dInsane;
#endif
    prTree->m_dEdgeLength1[uNodeIndex] = dInsane;
    prTree->m_dEdgeLength2[uNodeIndex] = dInsane;
    prTree->m_dEdgeLength3[uNodeIndex] = dInsane;
    
    prTree->m_ptrName[uNodeIndex] = NULL;
    prTree->m_Ids[uNodeIndex] = uInsane;
}
/***   end: init_node   ***/


/**
 *
 */
void
init_cache(uint uCacheCount, tree_t *tree)
{
    uint uNodeIndex;
    
    tree->m_uCacheCount = uCacheCount;
    
    tree->m_uNeighbor1 = (uint *) malloc(sizeof(uint) * tree->m_uCacheCount);
    tree->m_uNeighbor2 = (uint *) malloc(sizeof(uint) * tree->m_uCacheCount);
    tree->m_uNeighbor3 = (uint *) malloc(sizeof(uint) * tree->m_uCacheCount);
    
    tree->m_Ids = (uint *) malloc(sizeof(uint) * tree->m_uCacheCount);
    
    tree->m_dEdgeLength1 = (double *) malloc(sizeof(double) * tree->m_uCacheCount);
    tree->m_dEdgeLength2 = (double *) malloc(sizeof(double) * tree->m_uCacheCount);
    tree->m_dEdgeLength3 = (double *) malloc(sizeof(double) * tree->m_uCacheCount);
#if USE_HEIGHT
    tree->m_dHeight = (double *) malloc(sizeof(double) * tree->m_uCacheCount);
    tree->m_bHasHeight = (bool *) malloc(sizeof(bool) * tree->m_uCacheCount);
#endif    
    tree->m_bhas_edge_length1 = (bool *) malloc(sizeof(bool) * tree->m_uCacheCount);
    tree->m_bhas_edge_length2 = (bool *) malloc(sizeof(bool) * tree->m_uCacheCount);
    tree->m_bhas_edge_length3 = (bool *) malloc(sizeof(bool) * tree->m_uCacheCount);

    tree->m_ptrName = (char **) malloc(sizeof(char *) * tree->m_uCacheCount);
    
    for (uNodeIndex = 0; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
        init_node(tree, uNodeIndex);
    }
}
/***   end: init_cache   ***/




/**
 *
 * @note Replacing Tree::Clear but no freeing of memory! Just setting
 * everything to 0/NULL
 */
void
zero_tree(tree_t *tree)
{
    assert(NULL != tree);
    tree->m_uNodeCount = 0;
    tree->m_uCacheCount = 0;

    tree->m_uNeighbor1 = NULL;
    tree->m_uNeighbor2 = NULL;
    tree->m_uNeighbor3 = NULL;
    
    tree->m_dEdgeLength1 = NULL;
    tree->m_dEdgeLength2 = NULL;
    tree->m_dEdgeLength3 = NULL;
    
#if USE_HEIGHT
    tree->m_dHeight = NULL;
    tree->m_bHasHeight = NULL;
#endif
    tree->m_bhas_edge_length1 = NULL;
    tree->m_bhas_edge_length2 = NULL;
    tree->m_bhas_edge_length3 = NULL;

    tree->m_ptrName = NULL;
    tree->m_Ids = NULL;
    
    tree->m_bRooted = FALSE;
    tree->m_uRootNodeIndex = 0;
}
/* end: zero_tree */



/**
 *
 */
void
tree_to_file_node_rooted(tree_t *tree, uint uNodeIndex, FILE *fp)
{
    bool bGroup;
    
    assert(is_rooted(tree));
    assert(NULL != tree);
    bGroup = !is_leaf(uNodeIndex, tree) || is_root(uNodeIndex, tree);
    
    if (bGroup) 
        fprintf(fp, "(\n");
    

    if (is_leaf(uNodeIndex, tree)) {
        /* File.PutString(GetName(uNodeIndex)); */
        fprintf(fp, "%s", tree->m_ptrName[uNodeIndex]);
    } else {
        tree_to_file_node_rooted(tree, get_left(uNodeIndex, tree), fp);
        fprintf(fp, ",\n");
        tree_to_file_node_rooted(tree, get_right(uNodeIndex, tree), fp);
    }

    if (bGroup)
        fprintf(fp, ")");

    if (!is_root(uNodeIndex, tree)) {
        uint uParent = get_parent(uNodeIndex, tree);
        if (has_edge_length(uNodeIndex, uParent, tree))
            fprintf(fp, ":%g", get_edge_length(uNodeIndex, uParent, tree));
    }
    fprintf(fp, "\n");
}
/***   end: tree_to_file_node_rooted   ***/


/**
 *
 *
 */
void
validate_tree(tree_t *tree)
{
    uint uNodeIndex;
    assert(NULL != tree);
    for (uNodeIndex = 0; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
        validate_node(uNodeIndex, tree);
    }
    /* FIXME: maybe set negative length to zero? What impact would
     * that have? */
}
/***   end validate_tree   ***/



/**
 *
 *
 */
void
validate_node(uint uNodeIndex, tree_t *tree)
{
    uint uNeighborCount;
    uint n1, n2, n3;
    assert(NULL != tree);
    
    if (uNodeIndex >= tree->m_uNodeCount)
        printf("validate_node(%u), %u nodes", uNodeIndex, tree->m_uNodeCount);

    uNeighborCount = get_neighbor_count(uNodeIndex, tree);
    

    if (2 == uNeighborCount) {
        if (!tree->m_bRooted) {
            printf("Tree::validate_node: Node %u has two neighbors, tree is not rooted",
                  uNodeIndex);
        }
        if (uNodeIndex != tree->m_uRootNodeIndex) {
            printf("Tree::validate_node: Node %u has two neighbors, but not root node=%u",
                  uNodeIndex, tree->m_uRootNodeIndex);
        }
    }

    n1 = tree->m_uNeighbor1[uNodeIndex];
    n2 = tree->m_uNeighbor2[uNodeIndex];
    n3 = tree->m_uNeighbor3[uNodeIndex];
    
    if (NULL_NEIGHBOR == n2 && NULL_NEIGHBOR != n3) {
        printf("Tree::validate_node, n2=null, n3!=null, uNodeIndex=%u", uNodeIndex);
    }
    if (NULL_NEIGHBOR == n3 && NULL_NEIGHBOR != n2) {
        printf("Tree::validate_node, n3=null, n2!=null, uNodeIndex=%u", uNodeIndex);
    }
    
    if (n1 != NULL_NEIGHBOR)
        assert_are_neighbors(uNodeIndex, n1, tree);
    if (n2 != NULL_NEIGHBOR)
        assert_are_neighbors(uNodeIndex, n2, tree);
    if (n3 != NULL_NEIGHBOR)
        assert_are_neighbors(uNodeIndex, n3, tree);


    
    if (n1 != NULL_NEIGHBOR && (n1 == n2 || n1 == n3)) {
        printf("Tree::validate_node, duplicate neighbors in node %u", uNodeIndex);
    }
    if (n2 != NULL_NEIGHBOR && (n2 == n1 || n2 == n3)) {
        printf("Tree::validate_node, duplicate neighbors in node %u", uNodeIndex);
    }
    if (n3 != NULL_NEIGHBOR && (n3 == n1 || n3 == n2)) {
        printf("Tree::validate_node, duplicate neighbors in node %u", uNodeIndex);
    }


    if (is_rooted(tree)) {
        if (NULL_NEIGHBOR == get_parent(uNodeIndex, tree)) {
            if (uNodeIndex != tree->m_uRootNodeIndex) {
                printf("Tree::ValiateNode(%u), no parent", uNodeIndex);
            }
        } else if (get_left(get_parent(uNodeIndex, tree), tree) != uNodeIndex &&
                   get_right(get_parent(uNodeIndex, tree), tree) != uNodeIndex) {
            printf("Tree::validate_node(%u), parent / child mismatch", uNodeIndex);
        }
    }
}
/***   end: validate_node   ***/


/**
 *
 *
 */
void
assert_are_neighbors(unsigned uNodeIndex1, unsigned uNodeIndex2, tree_t *tree)
{
    bool Has12, Has21;
    assert(NULL != tree);

    if (uNodeIndex1 >= tree->m_uNodeCount || uNodeIndex2 >= tree->m_uNodeCount)
        printf("assert_are_neighbors(%u,%u), are %u nodes",
              uNodeIndex1, uNodeIndex2, tree->m_uNodeCount);

    if (tree->m_uNeighbor1[uNodeIndex1] != uNodeIndex2 &&
        tree->m_uNeighbor2[uNodeIndex1] != uNodeIndex2 &&
        tree->m_uNeighbor3[uNodeIndex1] != uNodeIndex2) {
        printf("assert_are_neighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
    }

    
    if (tree->m_uNeighbor1[uNodeIndex2] != uNodeIndex1 &&
        tree->m_uNeighbor2[uNodeIndex2] != uNodeIndex1 &&
        tree->m_uNeighbor3[uNodeIndex2] != uNodeIndex1) {
        printf("assert_are_neighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
    }


    Has12 = has_edge_length(uNodeIndex1, uNodeIndex2, tree);
    Has21 = has_edge_length(uNodeIndex2, uNodeIndex1, tree);
    if (Has12 != Has21) {
        has_edge_length(uNodeIndex1, uNodeIndex2, tree);
        has_edge_length(uNodeIndex2, uNodeIndex1, tree);
        printf("has_edge_length(%u, %u)=%c has_edge_length(%u, %u)=%c\n",
              uNodeIndex1,
              uNodeIndex2,
              Has12 ? 'T' : 'F',
              uNodeIndex2,
              uNodeIndex1,
              Has21 ? 'T' : 'F');
        printf("Tree::assert_are_neighbors, has_edge_length not symmetric");
    }

    
    if (Has12) {
        double d12 = get_edge_length(uNodeIndex1, uNodeIndex2, tree);
        double d21 = get_edge_length(uNodeIndex2, uNodeIndex1, tree);
        if (d12 != d21) {
            printf("Tree::assert_are_neighbors, Edge length disagrees %u-%u=%.3g, %u-%u=%.3g",
                  uNodeIndex1, uNodeIndex2, d12,
                  uNodeIndex2, uNodeIndex1, d21);
        }
    }
}
/***   end: assert_are_neighbors   ***/



/**
 *
 * @note see phyfromfile.cpp in Muscle3.7. Tree has to be a pointer to
 * an already allocated tree structure.
 *
 * return non-Zero on failure
 *
 * leafids will not be set here (FIXME:CHECK if true)
 */
int
load_tree(tree_t *tree, char *ftree)
{
    double dEdgeLength;
    bool bEdgeLength;
    char szToken[16];
    NEWICK_TOKEN_TYPE NTT;
    unsigned uThirdNode;
    FILE *fp = NULL;

    assert(tree!=NULL);
    assert(ftree!=NULL);

    if (NULL == (fp=fopen(ftree, "r"))) {
        printf("Couldn't open tree-file '%s' for reading. Skipping", ftree);
        return -1;
    }
    
    /* Assume rooted.
     * If we discover that it is unrooted, will convert on the fly.
     */
    tree_create_rooted(tree);

    bEdgeLength = get_group_from_file(fp, 0, &dEdgeLength, tree);


    /* Next token should be either ';' for rooted tree or ',' for
     * unrooted.
     */
    NTT = get_token(fp, szToken, sizeof(szToken));

    /* If rooted, all done. */
    if (NTT_Semicolon == NTT) {
        if (bEdgeLength)
            printf(" *** Warning *** edge length on root group in Newick file %s\n", ftree);
        validate_tree(tree);
        fclose(fp);
        return 0;
    }

    if (NTT_Comma != NTT)
        printf("Tree::FromFile, expected ';' or ',', got '%s'", szToken);

    uThirdNode = unroot_from_file(tree);
    bEdgeLength = get_group_from_file(fp, uThirdNode, &dEdgeLength, tree);
    if (bEdgeLength)
        set_edge_length(0, uThirdNode, dEdgeLength, tree);
    validate_tree(tree);

    fclose(fp);
    return 0;
}
/***   end load_tree   ***/



/**
 *
 */
void
expand_cache(tree_t *tree)
{
    const uint uNodeCount = 100;
    uint uNewCacheCount;
    uint *uNewNeighbor1, *uNewNeighbor2, *uNewNeighbor3;
    uint *uNewIds;
    double *dNewEdgeLength1, *dNewEdgeLength2, *dNewEdgeLength3;
#if USE_HEIGHT
    double *dNewHeight;
    bool *bNewHasHeight;
#endif
    bool *bNewhas_edge_length1, *bNewhas_edge_length2, *bNewhas_edge_length3;
    char **ptrNewName;

    assert(NULL != tree);
    uNewCacheCount = tree->m_uCacheCount + uNodeCount;
    uNewNeighbor1 = (uint *) malloc(
        uNewCacheCount * sizeof(uint));
    uNewNeighbor2 = (uint *) malloc(
        uNewCacheCount * sizeof(uint));
    uNewNeighbor3 = (uint *) malloc(
        uNewCacheCount * sizeof(uint));

    uNewIds = (uint *) calloc(
        uNewCacheCount, sizeof(uint));
    
    dNewEdgeLength1 = (double *) malloc(
        uNewCacheCount * sizeof(double));
    dNewEdgeLength2 =  (double *) malloc(
        uNewCacheCount * sizeof(double));
    dNewEdgeLength3 = (double *) malloc(
        uNewCacheCount * sizeof(double));
#if USE_HEIGHT
    dNewHeight = (double *) malloc(
        uNewCacheCount * sizeof(double));
    bNewHasHeight = (bool *) malloc(
        uNewCacheCount * sizeof(bool));
#endif
     
    bNewhas_edge_length1 = (bool *) malloc(
        uNewCacheCount * sizeof(bool));
    bNewhas_edge_length2 = (bool *) malloc(
        uNewCacheCount * sizeof(bool));
    bNewhas_edge_length3 = (bool *) malloc(
        uNewCacheCount * sizeof(bool));
    ptrNewName = (char **) calloc(uNewCacheCount, sizeof(char*));

    if (tree->m_uCacheCount > 0) {
        uint uUnsignedBytes, uEdgeBytes;
        uint uBoolBytes, uNameBytes;

        uUnsignedBytes = tree->m_uCacheCount*sizeof(uint);
        
        memcpy(uNewNeighbor1, tree->m_uNeighbor1, uUnsignedBytes);
        memcpy(uNewNeighbor2, tree->m_uNeighbor2, uUnsignedBytes);
        memcpy(uNewNeighbor3, tree->m_uNeighbor3, uUnsignedBytes);

        memcpy(uNewIds, tree->m_Ids, uUnsignedBytes);

        uEdgeBytes = tree->m_uCacheCount*sizeof(double);
        memcpy(dNewEdgeLength1, tree->m_dEdgeLength1, uEdgeBytes);
        memcpy(dNewEdgeLength2, tree->m_dEdgeLength2, uEdgeBytes);
        memcpy(dNewEdgeLength3, tree->m_dEdgeLength3, uEdgeBytes);
#if USE_HEIGHT
        memcpy(dNewHeight, tree->m_dHeight, uEdgeBytes);
#endif        
          
        uBoolBytes = tree->m_uCacheCount*sizeof(bool);
        memcpy(bNewhas_edge_length1, tree->m_bhas_edge_length1, uBoolBytes);
        memcpy(bNewhas_edge_length2, tree->m_bhas_edge_length2, uBoolBytes);
        memcpy(bNewhas_edge_length3, tree->m_bhas_edge_length3, uBoolBytes);
#if USE_HEIGHT
        memcpy(bNewHasHeight, tree->m_bHasHeight, uBoolBytes);
#endif
        uNameBytes = tree->m_uCacheCount*sizeof(char *);
        memcpy(ptrNewName, tree->m_ptrName, uNameBytes);

        /* similiar to free_tree
         */
        
        /* is_leaf needs m_uNodeCount and all m_uNeighbor's
         * so free first
         */
#if 0
        for (i=0; i<tree->m_uNodeCount; i++) {
            if (is_leaf(i, tree)) {
#ifndef NDEBUG
                if (NULL==tree->m_ptrName[i]) {
                    printf("FIXME tree->m_ptrName[%d] is already NULL", i);
                }
#endif
                free(tree->m_ptrName[i]);
            }
        }
#endif
        free(tree->m_ptrName);

        free(tree->m_uNeighbor1);
        free(tree->m_uNeighbor2);
        free(tree->m_uNeighbor3);

        free(tree->m_Ids);

        free(tree->m_dEdgeLength1);
        free(tree->m_dEdgeLength2);
        free(tree->m_dEdgeLength3);

        free(tree->m_bhas_edge_length1);
        free(tree->m_bhas_edge_length2);
        free(tree->m_bhas_edge_length3);
#if USE_HEIGHT
        free(tree->m_bHasHeight);
        free(tree->m_dHeight);
#endif
    }
    
    tree->m_uCacheCount = uNewCacheCount;
    tree->m_uNeighbor1 = uNewNeighbor1;
    tree->m_uNeighbor2 = uNewNeighbor2;
    tree->m_uNeighbor3 = uNewNeighbor3;
    tree->m_Ids = uNewIds;
    tree->m_dEdgeLength1 = dNewEdgeLength1;
    tree->m_dEdgeLength2 = dNewEdgeLength2;
    tree->m_dEdgeLength3 = dNewEdgeLength3;
        
#ifdef USE_HEIGHT
    tree->m_dHeight = dNewHeight;
    tree->m_bHasHeight = bNewHasHeight;
#endif
    tree->m_bhas_edge_length1 = bNewhas_edge_length1;
    tree->m_bhas_edge_length2 = bNewhas_edge_length2;
    tree->m_bhas_edge_length3 = bNewhas_edge_length3;
        
    tree->m_ptrName = ptrNewName;

}
/***   end: expand_cache   ***/



/**
 *
 * Tree must be pointer to an already allocated tree structure
 *
 */
void
tree_create_rooted(tree_t *tree)
{
    zero_tree(tree);
    expand_cache(tree);
    tree->m_uNodeCount = 1;

    tree->m_uNeighbor1[0] = NULL_NEIGHBOR;
    tree->m_uNeighbor2[0] = NULL_NEIGHBOR;
    tree->m_uNeighbor3[0] = NULL_NEIGHBOR;
    
    tree->m_bhas_edge_length1[0] = FALSE;
    tree->m_bhas_edge_length2[0] = FALSE;
    tree->m_bhas_edge_length3[0] = FALSE;
#if USE_HEIGHT
    tree->m_bHasHeight[0] = FALSE;
#endif
    
    tree->m_uRootNodeIndex = 0;
    tree->m_bRooted = TRUE;
    
#ifndef NDEBUG
    validate_tree(tree);
#endif
}
/***   end: tree_create_rooted   ***/


/**
 *
 */
uint
unroot_from_file(tree_t *tree)
{
    uint uThirdNode;

    
    if (!tree->m_bRooted)
        printf("Tree::Unroot, not rooted");
    
    /* Convention: root node is always node zero */
    assert(is_root(0, tree));
    assert(NULL_NEIGHBOR == tree->m_uNeighbor1[0]);

    uThirdNode = tree->m_uNodeCount++;
    
    tree->m_uNeighbor1[0] = uThirdNode;
    tree->m_uNeighbor1[uThirdNode] = 0;
    
    tree->m_uNeighbor2[uThirdNode] = NULL_NEIGHBOR;
    tree->m_uNeighbor3[uThirdNode] = NULL_NEIGHBOR;
    
    tree->m_dEdgeLength1[0] = 0;
    tree->m_dEdgeLength1[uThirdNode] = 0;
    tree->m_bhas_edge_length1[uThirdNode] = TRUE;
    
    tree->m_bRooted = FALSE;
    


    return uThirdNode;
}
/***   end: unroot_from_file   ***/



/**
 *
 *
 */
bool
get_group_from_file(FILE *fp, uint uNodeIndex, double *ptrdEdgeLength, tree_t *tree)
{
    char szToken[1024];
    NEWICK_TOKEN_TYPE NTT = get_token(fp, szToken, sizeof(szToken));
    bool bRightLength;
    bool bEof;
    char c;

    /* Group is either leaf name or (left, right). */
    if (NTT_String == NTT) {
        set_leaf_name(uNodeIndex, szToken, tree);

    } else if (NTT_Lparen == NTT) {
        const unsigned uLeft = append_branch(tree, uNodeIndex);
        const unsigned uRight = uLeft + 1;
        double dEdgeLength;
        bool bLeftLength = get_group_from_file(fp, uLeft, &dEdgeLength, tree);
        
        /* Left sub-group...
         */
#if TRACE
        printf("Got '(', group is compound, expect left sub-group");

        if (bLeftLength) {
            printf("Edge length for left sub-group: %.3g", dEdgeLength);
        } else {
            printf("No edge length for left sub-group");
        }
#endif
        if (bLeftLength)
            set_edge_length(uNodeIndex, uLeft, dEdgeLength, tree);

        /* ... then comma ...
         */

        NTT = get_token(fp, szToken, sizeof(szToken));
        if (NTT_Comma != NTT)
            printf("Tree::get_group_from_file, expected ',', got '%s'", szToken);
        
        /* ...then right sub-group...
         */

        bRightLength = get_group_from_file(fp, uRight, &dEdgeLength, tree);
        if (bRightLength)
            set_edge_length(uNodeIndex, uRight, dEdgeLength, tree);
        


        /* ... then closing parenthesis.
         */

        NTT = get_token(fp, szToken, sizeof(szToken));
        if (NTT_Rparen == NTT)
            ;
        else if (NTT_Comma == NTT) {
            if (ungetc(',', fp)==EOF)
                printf("%s", "ungetc failed");
            return FALSE;
        } else
            printf("Tree::get_group_from_file, expected ')' or ',', got '%s'", szToken);
    } else {
        printf("Tree::get_group_from_file, expected '(' or leaf name, got '%s'",
              szToken);
    }
    
    /* Group may optionally be followed by edge length.
     */
    bEof = file_skip_whitex(fp);
    if (bEof)
        return FALSE;
    if ((c = fgetc(fp))==EOF) /* GetCharX */
        printf("fgetc reached end of file");

    if (':' == c) {
        NTT = get_token(fp, szToken, sizeof(szToken));
        if (NTT_String != NTT)
            printf("Tree::get_group_from_file, expected edge length, got '%s'", szToken);
        *ptrdEdgeLength = atof(szToken);
        return TRUE;
    }
    if (ungetc(c, fp)==EOF)
        printf("%s", "ungetc failed");
    
    return FALSE;
}
/***   end: get_group_from_file   ***/




/**
 *
 */
void
file_skip_white(FILE *fp)
{
    bool bEof = file_skip_whitex(fp);
    if (bEof)
        printf("End-of-file skipping white space");
}
/***   end: file_skip_white   ***/




/**
 *
 */
bool
file_skip_whitex(FILE *fp)
{
    for (;;) {
        int c;
        bool bEof;
        
        /* GetChar */
        if ((c = fgetc(fp))==EOF) {
            bEof = TRUE;
        } else {
            bEof = FALSE;
        }

        if (bEof)
            return TRUE;
        if (!isspace(c)) {
            if (ungetc(c, fp)==EOF)
                printf("%s", "ungetc failed");
            break;
        }
    }
    return FALSE;
}
/***   end: file_skip_whitex   ***/




/**
 *
 */
NEWICK_TOKEN_TYPE
get_token(FILE *fp, char szToken[], uint uBytes)
{
    char c;
    unsigned uBytesCopied = 0;
    NEWICK_TOKEN_TYPE TT;

    /* Skip leading white space */
    file_skip_white(fp);

    if ((c = fgetc(fp))==EOF) /* GetCharX */
        printf("fgetc reached end of file");
    
    /* In case a single-character token */
    szToken[0] = c;
    szToken[1] = 0;

    switch (c) {

    case '(':
        return NTT_Lparen;
        
    case ')':
        return NTT_Rparen;
        
    case ':':
        return NTT_Colon;
        
    case ';':
        return NTT_Semicolon;
        
    case ',':
        return NTT_Comma;
        
    case '\'':
        TT = NTT_SingleQuotedString;
        if ((c = fgetc(fp))==EOF) /* GetCharX */
            printf("fgetc reached end of file");
        break;
        
    case '"':
        TT = NTT_DoubleQuotedString;
        if ((c = fgetc(fp))==EOF) /* GetCharX */
            printf("fgetc reached end of file");
        break;
        
    case '[':
        TT = NTT_Comment;
        break;
        
    default:
        TT = NTT_String;
        break;
    }
    
    for (;;)
    {
        bool bEof;
        if (TT != NTT_Comment) {
            if (uBytesCopied < uBytes - 2)  {
                szToken[uBytesCopied++] = c;
                szToken[uBytesCopied] = 0;
            } else {
                printf("Tree::get_token: input buffer too small, token so far='%s'", szToken);
            }
        }
        c = fgetc(fp); /* GetChar */
        bEof = (c==EOF ? TRUE : FALSE);
        if (bEof)
            return TT;
        
        switch (TT) {

        case NTT_String:
            if (0 != strchr("():;,", c)) {
                if (ungetc(c, fp)==EOF)
                    printf("%s", "ungetc failed");
                return NTT_String;
            }
            if (isspace(c))
                return NTT_String;
            break;
            
        case NTT_SingleQuotedString:
            if ('\'' == c)
                return NTT_String;
            break;
            
        case NTT_DoubleQuotedString:
            if ('"' == c)
                return NTT_String;
            break;
            
        case NTT_Comment:
            if (']' == c)
                return get_token(fp, szToken, uBytes);
            break;
            
        default:
            printf("Tree::get_token, invalid TT=%u", TT);
        }
    }
}
/***   end: get_token   ***/



/***   set_leaf_name
 *
 */
void
set_leaf_name(unsigned uNodeIndex, const char *ptrName, tree_t *tree)
{
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(is_leaf(uNodeIndex, tree));
    free(tree->m_ptrName[uNodeIndex]);
    /*LOG_DEBUG("Setting tree->m_ptrName[uNodeIndex=%d] to %s", uNodeIndex, ptrName);*/
    tree->m_ptrName[uNodeIndex] = strdup(ptrName);
}
/***   end: set_leaf_name   ***/




/**
 *
 * Append a new branch. This adds two new nodes and joins them to an
 * existing leaf node. Return value is k, new nodes have indexes k and
 * k+1 respectively.
 *
 */
uint
append_branch(tree_t *tree, uint uExistingLeafIndex)
{
    uint uNewLeaf1;
    uint uNewLeaf2;
    
    assert(tree!=NULL);
    if (0 == tree->m_uNodeCount) {
        printf("tree has not been created");
    }
    assert(NULL_NEIGHBOR == tree->m_uNeighbor2[uExistingLeafIndex]);
    assert(NULL_NEIGHBOR == tree->m_uNeighbor3[uExistingLeafIndex]);
    assert(uExistingLeafIndex < tree->m_uNodeCount);
#ifndef NDEBUG
    if (!is_leaf(uExistingLeafIndex, tree)) {
        printf("append_branch(%u): not leaf", uExistingLeafIndex);
    }
#endif

    if (tree->m_uNodeCount >= tree->m_uCacheCount - 2) {
        expand_cache(tree);
    }
    uNewLeaf1 = tree->m_uNodeCount;
    uNewLeaf2 = tree->m_uNodeCount + 1;

    tree->m_uNodeCount += 2;
    
    tree->m_uNeighbor2[uExistingLeafIndex] = uNewLeaf1;
    tree->m_uNeighbor3[uExistingLeafIndex] = uNewLeaf2;
    
    tree->m_uNeighbor1[uNewLeaf1] = uExistingLeafIndex;
    tree->m_uNeighbor1[uNewLeaf2] = uExistingLeafIndex;
    
    tree->m_uNeighbor2[uNewLeaf1] = NULL_NEIGHBOR;
    tree->m_uNeighbor2[uNewLeaf2] = NULL_NEIGHBOR;
    
    tree->m_uNeighbor3[uNewLeaf1] = NULL_NEIGHBOR;
    tree->m_uNeighbor3[uNewLeaf2] = NULL_NEIGHBOR;
    
    tree->m_dEdgeLength2[uExistingLeafIndex] = 0;
    tree->m_dEdgeLength3[uExistingLeafIndex] = 0;
    
    tree->m_dEdgeLength1[uNewLeaf1] = 0;
    tree->m_dEdgeLength2[uNewLeaf1] = 0;
    tree->m_dEdgeLength3[uNewLeaf1] = 0;
    
    tree->m_dEdgeLength1[uNewLeaf2] = 0;
    tree->m_dEdgeLength2[uNewLeaf2] = 0;
    tree->m_dEdgeLength3[uNewLeaf2] = 0;
    
    tree->m_bhas_edge_length1[uNewLeaf1] = FALSE;
    tree->m_bhas_edge_length2[uNewLeaf1] = FALSE;
    tree->m_bhas_edge_length3[uNewLeaf1] = FALSE;
    
    tree->m_bhas_edge_length1[uNewLeaf2] = FALSE;
    tree->m_bhas_edge_length2[uNewLeaf2] = FALSE;
    tree->m_bhas_edge_length3[uNewLeaf2] = FALSE;
    
#if USE_HEIGHT
    tree->m_bHasHeight[uNewLeaf1] = FALSE;
    tree->m_bHasHeight[uNewLeaf2] = FALSE;
#endif 
    tree->m_Ids[uNewLeaf1] = uInsane;
    tree->m_Ids[uNewLeaf2] = uInsane;
    
    return uNewLeaf1;
}
/***   end: append_branch   ***/


/**
 *
 *
 */
void
set_edge_length(uint uNodeIndex1, uint uNodeIndex2,
              double dLength, tree_t *tree)
{
    assert(uNodeIndex1 < tree->m_uNodeCount && uNodeIndex2 < tree->m_uNodeCount);
    assert(is_edge(uNodeIndex1, uNodeIndex2, tree));
    
    if (tree->m_uNeighbor1[uNodeIndex1] == uNodeIndex2) {
        tree->m_dEdgeLength1[uNodeIndex1] = dLength;
        tree->m_bhas_edge_length1[uNodeIndex1] = TRUE;
    } else if (tree->m_uNeighbor2[uNodeIndex1] == uNodeIndex2) {
        tree->m_dEdgeLength2[uNodeIndex1] = dLength;
        tree->m_bhas_edge_length2[uNodeIndex1] = TRUE;
        
    } else {
        assert(tree->m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
        tree->m_dEdgeLength3[uNodeIndex1] = dLength;
        tree->m_bhas_edge_length3[uNodeIndex1] = TRUE;
    }
    
    if (tree->m_uNeighbor1[uNodeIndex2] == uNodeIndex1) {
        tree->m_dEdgeLength1[uNodeIndex2] = dLength;
        tree->m_bhas_edge_length1[uNodeIndex2] = TRUE;
    } else if (tree->m_uNeighbor2[uNodeIndex2] == uNodeIndex1) {
        tree->m_dEdgeLength2[uNodeIndex2] = dLength;
        tree->m_bhas_edge_length2[uNodeIndex2] = TRUE;
    } else {
        assert(tree->m_uNeighbor3[uNodeIndex2] == uNodeIndex1);
        tree->m_dEdgeLength3[uNodeIndex2] = dLength;
        tree->m_bhas_edge_length3[uNodeIndex2] = TRUE;
    }
}
/***   end: set_edge_length   ***/



/**
 *
 * Debug output
 *
 * LogMe in phy.cpp
 *
 */
void
log_tree(tree_t *tree, FILE *fp)
{
    uint uNodeIndex;
    uint n1, n2, n3;
    char *ptrName;
    
    fprintf(fp, "This is a tree with %u nodes, which is ", tree->m_uNodeCount);
    
    if (is_rooted(tree)) {
        fprintf(fp, "rooted:\n");
        fprintf(fp, "Index  Parnt  LengthP  Left   LengthL  Right  LengthR     Id  Name\n");
        fprintf(fp, "-----  -----  -------  ----   -------  -----  -------  -----  ----\n");

    } else {
        fprintf(fp, "unrooted;\n");
        fprintf(fp, "Index  Nbr_1  Length1  Nbr_2  Length2  Nbr_3  Length3     Id  Name\n");
        fprintf(fp, "-----  -----  -------  -----  -------  -----  -------  -----  ----\n");
    }
    
    for (uNodeIndex = 0; uNodeIndex < tree->m_uNodeCount; ++uNodeIndex) {
        fprintf(fp, "%5u  ", uNodeIndex);
        n1 = tree->m_uNeighbor1[uNodeIndex];
        n2 = tree->m_uNeighbor2[uNodeIndex];
        n3 = tree->m_uNeighbor3[uNodeIndex];
        
        if (NULL_NEIGHBOR != n1) {
            fprintf(fp, "%5u  ", n1);
            if (tree->m_bhas_edge_length1[uNodeIndex])
                fprintf(fp, "%7.3g  ", tree->m_dEdgeLength1[uNodeIndex]);
            else
                fprintf(fp, "      *  ");
        } else {
            fprintf(fp, "                ");
        }
        
        if (NULL_NEIGHBOR != n2) {
            fprintf(fp, "%5u  ", n2);
            if (tree->m_bhas_edge_length2[uNodeIndex])
                fprintf(fp, "%7.3g  ", tree->m_dEdgeLength2[uNodeIndex]);
            else
                fprintf(fp, "      *  ");
        } else {
            fprintf(fp, "                ");
        }
        
        if (NULL_NEIGHBOR != n3) {
            fprintf(fp, "%5u  ", n3);
            if (tree->m_bhas_edge_length3[uNodeIndex])
                fprintf(fp, "%7.3g  ", tree->m_dEdgeLength3[uNodeIndex]);
            else
                fprintf(fp, "      *  ");
        } else {
            fprintf(fp, "                ");
        }
        
        if (tree->m_Ids != 0 && is_leaf(uNodeIndex, tree)) {
            unsigned uId = tree->m_Ids[uNodeIndex];
            if (uId == uInsane)
                fprintf(fp, "    *");
            else
                fprintf(fp, "%5u", uId);
        } else {
            fprintf(fp, "     ");
        }
        
        if (tree->m_bRooted && uNodeIndex == tree->m_uRootNodeIndex)
            fprintf(fp, "  [ROOT] ");
        ptrName = tree->m_ptrName[uNodeIndex];
        if (ptrName != 0)
            fprintf(fp, "  %s", ptrName);
        
        fprintf(fp, "\n");
    }     
}
/***   end: log_tree   ***/



/**
 *
 * replaces m_uLeafCount
 */
uint
get_leaf_count(tree_t *tree)
{
    assert(tree!=NULL);
    
    return (tree->m_uNodeCount+1)/2;
}
/***   get_leaf_count   ***/



/**
 *
 */
uint
get_node_count(tree_t *tree)
{
    assert(tree!=NULL);
    
    return 2*(get_leaf_count(tree)) - 1;
}
/***   end: get_node_count   ***/


/**
 *
 */
uint
get_neighbor(uint uNodeIndex, uint uNeighborSubscript, tree_t *prTree)
{
    assert(uNodeIndex < prTree->m_uNodeCount);
    switch (uNeighborSubscript)
    {
    case 0:
        return prTree->m_uNeighbor1[uNodeIndex];
    case 1:
        return prTree->m_uNeighbor2[uNodeIndex];
    case 2:
        return prTree->m_uNeighbor3[uNodeIndex];
    }
    printf("Internal error in %s: sub=%u", __FUNCTION__, uNeighborSubscript);
    return NULL_NEIGHBOR;
}
/***   end: get_neighbor   ***/





/**
 *
 */
void
set_leaf_ID(tree_t *tree, uint uNodeIndex, uint uId)
{
    assert(uNodeIndex < tree->m_uNodeCount);
    assert(is_leaf(uNodeIndex, tree));
    tree->m_Ids[uNodeIndex] = uId;
}
/***   end: set_leaf_ID    ***/


/**
 *
 */
uint
get_root_node_index(tree_t *tree)
{
    assert(NULL!=tree);
    return tree->m_uRootNodeIndex;
}
/***   end: get_root_node_index   ***/



/**
 * @note avoid usage if you want to iterate over all indices, because
 * this will be slow
 *
 */
uint
leaf_index_to_node_index(uint uLeafIndex, tree_t *prTree) {
    uint uLeafCount = 0;
    unsigned uNodeCount = get_node_count(prTree);
    uint uNodeIndex;
    
    for (uNodeIndex = 0; uNodeIndex < uNodeCount; uNodeIndex++) {
        if (is_leaf(uNodeIndex, prTree)) {
            if (uLeafCount == uLeafIndex) {
                return uNodeIndex;
            } else {
                uLeafCount++;
            }
        }
    }
    printf("Internal error: node index out of range");
    return 0;
}
/***   end: leaf_index_to_node_index   ***/




/**
 * @brief Append a (source) tree to a (dest) tree to a given node
 * which will be replaced. All other nodes in that tree will stay the
 * same.
 *
 * @param[out] prDstTree
 * The tree to append to
 * @param[in] uDstTreeReplaceNodeIndex
 * Dest tree node to which source tree will be appended
 * @param[in] prSrcTree
 * The tree to append
 * 
 * @note No nodes inside prDstTree will change except
 * uDstTreeReplaceNodeIndex
 *
 *
 * @warning: Function won't check or touch the m_Ids/leaf-indices!
 * That means if you want to join two trees with leaf indices 1-10 and
 * 1-10 your m_Ids/leaf-indices won't be unique anymore and the
 * association between your sequences and the tree are broken. Make
 * sure m_Ids are unique before calling me.
 *
 * The easiest would have been to do this by recursively calling
 * append_branch() (after adding uSrcTreeNodeIndex as extra argument to
 * this function). But recursion is evil. Yet another version would be
 * to setup all the data and call tree_create() to create a third
 * tree, which seems complicated and wasteful. The approach taken here
 * is the following: increase dest tree memory, copy over each src
 * tree node data and update the indices and counters. This is tricky
 * and has a lot of potential for bugs if tree interface is changed
 * (and even if it isn't).
 *
 */
void
append_tree(tree_t *prDstTree, uint uDstTreeReplaceNodeIndex, tree_t *prSrcTree)
{
    uint uSrcTreeNodeIndex;
    uint uOrgDstTreeNodeCount;

    assert(NULL!=prDstTree);
    assert(NULL!=prSrcTree);
    assert(uDstTreeReplaceNodeIndex < prDstTree->m_uNodeCount);
    assert(is_leaf(uDstTreeReplaceNodeIndex, prDstTree));
    assert(is_rooted(prDstTree) && is_rooted(prSrcTree));

    
    uOrgDstTreeNodeCount = prDstTree->m_uNodeCount;

    
    /* increase dest tree memory
     */
    while (prDstTree->m_uCacheCount
           <
           (get_node_count(prDstTree) + get_node_count(prSrcTree))) {
        expand_cache(prDstTree);
    }


    /* copy all src tree nodes
     *
     */
    for (uSrcTreeNodeIndex=0;
         uSrcTreeNodeIndex<get_node_count(prSrcTree); uSrcTreeNodeIndex++) {
        uint uNewDstNodeIndex = prDstTree->m_uNodeCount;
        
        /* distinguish 4 cases for src nodes to copy:
         *
         * 1. src node is the only node, i.e. root & leaf
         *
         * 2. src node is root: set only left & right, but not parent
         * and just replace the given dest index
         *
         * 3. src node is leaf: set only parent
         *
         * 4. src node is internal node: update all three neighbours
         *
         * FIXME: this is messy. Is there a cleaner way to do this by
         * merging all cases?
         *
         */
        if (is_root(uSrcTreeNodeIndex, prSrcTree) && is_leaf(uSrcTreeNodeIndex, prSrcTree)) {
            /* special case: if this is the only member in
             * tree, i.e. it's root and leaf. Copy leaf name and leaf
             * id. No neighbours to update
             */

            /* free dst node name if set */
            if (NULL != prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]) {
                free(prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]);
            }

            prDstTree->m_ptrName[uDstTreeReplaceNodeIndex] =
                strdup(get_leaf_name(uSrcTreeNodeIndex, prSrcTree));

            prDstTree->m_Ids[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_Ids[uSrcTreeNodeIndex];

            /* no updated of uNodeCount, because we used the replace node */


            
        } else if (is_root(uSrcTreeNodeIndex, prSrcTree)) {
            /* src node is root: replace uDstTreeReplaceNodeIndex
             * (not uNewDstNodeIndex) with src root, i.e. the
             * uDstTreeReplaceNodeIndex becomes an internal node now.
             *
             * We only have two neighbours 2 & 3 (no parent). Keep old
             * parent info (neighbor 1).
             */

            /* free dst node name if set */
            if (NULL != prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]) {
                free(prDstTree->m_ptrName[uDstTreeReplaceNodeIndex]);
            }
            
            prDstTree->m_uNeighbor2[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_uNeighbor2[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            prDstTree->m_uNeighbor3[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_uNeighbor3[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;          
            
            prDstTree->m_bhas_edge_length2[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_bhas_edge_length2[uSrcTreeNodeIndex];
            prDstTree->m_bhas_edge_length3[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_bhas_edge_length3[uSrcTreeNodeIndex];            

            prDstTree->m_dEdgeLength2[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_dEdgeLength2[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength3[uDstTreeReplaceNodeIndex] =
                prSrcTree->m_dEdgeLength3[uSrcTreeNodeIndex];

            /* make Id invalid */
            prDstTree->m_Ids[uDstTreeReplaceNodeIndex] = uInsane;

            /* no updated of uNodeCount, because we used the replace node */

#if TRACE
            printf("Updated dst rpl node %d with the src root node: (untouched) parent=%d (%f) left=%d (%f) right=%d (%f)",
                      uDstTreeReplaceNodeIndex,
                      prDstTree->m_uNeighbor1[uDstTreeReplaceNodeIndex], prDstTree->m_dEdgeLength1[uDstTreeReplaceNodeIndex],
                      prDstTree->m_uNeighbor2[uDstTreeReplaceNodeIndex], prDstTree->m_dEdgeLength2[uDstTreeReplaceNodeIndex],
                      prDstTree->m_uNeighbor3[uDstTreeReplaceNodeIndex], prDstTree->m_dEdgeLength3[uDstTreeReplaceNodeIndex]);
#endif
            
        } else if (is_leaf(uSrcTreeNodeIndex, prSrcTree)) {
            /* src node is a leaf, which means we only have one
             * neighbour, and that is its parent, i.e. n1
             *
             */

            /* initialise/zero new node to default values
             */
            init_node(prDstTree, uNewDstNodeIndex);

        
            /*  update m_ptrName/leaf name
             */
            prDstTree->m_ptrName[uNewDstNodeIndex] =
                strdup(get_leaf_name(uSrcTreeNodeIndex, prSrcTree));

            /* update parent node (beware of special case: parent was
               src tree root */
            if (is_root(prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex], prSrcTree)) {
                    prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                        uDstTreeReplaceNodeIndex;
            } else {
                prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                    prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            }

            /* update edge length info to parent
             */
            prDstTree->m_bhas_edge_length1[uNewDstNodeIndex] =
                prSrcTree->m_bhas_edge_length1[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength1[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength1[uSrcTreeNodeIndex];

            /* update sequence/object id
             */
            prDstTree->m_Ids[uNewDstNodeIndex] =
                prSrcTree->m_Ids[uSrcTreeNodeIndex];            

            /* we used a new node so increase their count */
            prDstTree->m_uNodeCount += 1;
            

            
        } else  {
            /* src node is not root neither leaf, means we have an
             * internal node. Update all neighbour info
             * 
             */

            /* initialise/zero node values to default values
             */
            init_node(prDstTree, uNewDstNodeIndex);
          
            /* update neigbours
             */
            /* parent: special case if parent was src tree root */
            if (is_root(prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex], prSrcTree)) {
                    prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                        uDstTreeReplaceNodeIndex;
            } else {
                prDstTree->m_uNeighbor1[uNewDstNodeIndex] =
                    prSrcTree->m_uNeighbor1[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            }
            /* left */
            prDstTree->m_uNeighbor2[uNewDstNodeIndex] =
                prSrcTree->m_uNeighbor2[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;
            /* right */
            prDstTree->m_uNeighbor3[uNewDstNodeIndex] =
                prSrcTree->m_uNeighbor3[uSrcTreeNodeIndex] + uOrgDstTreeNodeCount;

            /* update edge length info
             */
            /* parent */
            prDstTree->m_bhas_edge_length1[uNewDstNodeIndex] =
                prSrcTree->m_bhas_edge_length1[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength1[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength1[uSrcTreeNodeIndex];
            /* left */
            prDstTree->m_bhas_edge_length2[uNewDstNodeIndex] =
                prSrcTree->m_bhas_edge_length2[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength2[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength2[uSrcTreeNodeIndex];
            /* right */
            prDstTree->m_bhas_edge_length3[uNewDstNodeIndex] =
                prSrcTree->m_bhas_edge_length3[uSrcTreeNodeIndex];
            prDstTree->m_dEdgeLength3[uNewDstNodeIndex] =
                prSrcTree->m_dEdgeLength3[uSrcTreeNodeIndex];

            /* we used a new node so increase their count */
            prDstTree->m_uNodeCount += 1;
            
#if TRACE
            printf("Updated dst node %d with an internal src node: parent=%d (%f) left=%d (%f) right=%d (%f)",
                      uNewDstNodeIndex,
                      prDstTree->m_uNeighbor1[uNewDstNodeIndex], prDstTree->m_dEdgeLength1[uNewDstNodeIndex],
                      prDstTree->m_uNeighbor2[uNewDstNodeIndex], prDstTree->m_dEdgeLength2[uNewDstNodeIndex],
                      prDstTree->m_uNeighbor3[uNewDstNodeIndex], prDstTree->m_dEdgeLength3[uNewDstNodeIndex]);
#endif
        }

    }
    /* end for each src tree node */

    
    /*
     * m_uRootNodeIndex stays the same.
     *
     * No need to touch m_uCacheCount.
     *
     */    
#if USE_HEIGHT
    printf("Internal error: Height usage not implemented in %s", __FUNCTION__);
#endif 

    
#ifndef NDEBUG
    validate_tree(prDstTree);
#endif
        
    return;
}
/***   end: append_tree()   ***/

/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

/*
 * Notes:
 * ------
 * LINKAGE become linkage_t here
 *
 * Replaced the the following member functions for DistCalc DC:
 * DC.GetId = sequence id as int
 * DC.GetName = sequence name
 * DC.GetCount = matrix dim
 * DC.DistRange = vector / matrix row for object i with index j<i
 *
 * Log() has been replaced with Clustal's Info(), Quiet() with Log(&rLog, LOG_FATAL)
 *
 * Made TriangleSubscript() and g_ulTriangleSize ulong to prevent overflow for many sequences
 */

#ifndef ulint
/* limit use of unsigned vars (see coding_style_guideline.txt) */
typedef unsigned long int ulong;
#endif



#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* from distcalc.h */
typedef double dist_t;
static const dist_t BIG_DIST = (dist_t) 1e29;

/*static inline*/
ulong TriangleSubscript(uint uIndex1, uint uIndex2);




#define TRACE   0

#ifndef MIN
#define MIN(x, y)   ((x) < (y) ? (x) : (y))
#endif
#ifndef MIN
#define MAX(x, y)   ((x) > (y) ? (x) : (y))
#endif
#define AVG(x, y)   (((x) + (y))/2)

static uint g_uLeafCount;
static ulong g_ulTriangleSize;
static uint g_uInternalNodeCount;
static uint g_uInternalNodeIndex;

/* Triangular distance matrix is g_Dist, which is allocated
 * as a one-dimensional vector of length g_ulTriangleSize.
 * TriangleSubscript(i,j) maps row,column=i,j to the subscript
 * into this vector.
 * Row / column coordinates are a bit messy.
 * Initially they are leaf indexes 0..N-1.
 * But each time we create a new node (=new cluster, new subtree),
 * we re-use one of the two rows that become available (the children
 * of the new node). This saves memory.
 * We keep track of this through the g_uNodeIndex vector.
 */
static dist_t *g_Dist;

/* Distance to nearest neighbor in row i of distance matrix.
 * Subscript is distance matrix row.
 */
static dist_t *g_MinDist;

/* Nearest neighbor to row i of distance matrix.
 * Subscript is distance matrix row.
 */
static uint *g_uNearestNeighbor;

/* Node index of row i in distance matrix.
 * Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
 * Subscript is distance matrix row.
 */
static uint *g_uNodeIndex;

/* The following vectors are defined on internal nodes,
 * subscripts are internal node index 0..N-2.
 * For g_uLeft/Right, value is the node index 0 .. 2N-2
 * because a child can be internal or leaf.
 */
static uint *g_uLeft;
static uint *g_uRight;
static dist_t *g_Height;
static dist_t *g_LeftLength;
static dist_t *g_RightLength;


/***   CalcDistRange
 *
 * Imitation of DistCalc.DistRange
 *
 * Sets values of row (vector / matrix row) to distances for object i with index j<i
 *
 * row must be preallocated
 */
void CalcDistRange(matrix_t *distmat, uint i, dist_t *row)
{
    uint j;
    for (j = 0; j < i; ++j) {
        row[j] = distmat->data[i][j];
    }
}
/* end of CalcDistRange */



/*static inline*/
ulong
TriangleSubscript(uint uIndex1, uint uIndex2)
{
    ulong v;
#ifndef NDEBUG
    if (uIndex1 >= g_uLeafCount || uIndex2 >= g_uLeafCount)
        printf("TriangleSubscript(%u,%u) %u", uIndex1, uIndex2, g_uLeafCount);
#endif
    if (uIndex1 >= uIndex2)
        v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
    else
        v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
    assert(v < (g_uLeafCount*(g_uLeafCount - 1))/2);
    return v;
}

#ifdef UNUSED
static void ListState()
{
    uint i, j;
    Info("Dist matrix\n");
    Info("     ");
    for (i = 0; i < g_uLeafCount; ++i)
    {
        if (uInsane == g_uNodeIndex[i])
            continue;
        Info("  %5u", g_uNodeIndex[i]);
    }
    Info("\n");

    for (i = 0; i < g_uLeafCount; ++i)
    {
        if (uInsane == g_uNodeIndex[i])
            continue;
        Info("%5u  ", g_uNodeIndex[i]);
        for (j = 0; j < g_uLeafCount; ++j)
        {
            if (uInsane == g_uNodeIndex[j])
                continue;
            if (i == j)
                Info("       ");
            else
            {
                ulong v = TriangleSubscript(i, j);
                Info("%5.2g  ", g_Dist[v]);
            }
        }
        Info("\n");
    }

    Info("\n");
    Info("    i   Node   NrNb      Dist\n");
    Info("-----  -----  -----  --------\n");
    for (i = 0; i < g_uLeafCount; ++i)
    {
        if (uInsane == g_uNodeIndex[i])
            continue;
        Info("%5u  %5u  %5u  %8.3f\n",
             i,
             g_uNodeIndex[i],
             g_uNearestNeighbor[i],
             g_MinDist[i]);
    }

    Info("\n");
    Info(" Node      L      R  Height  LLength  RLength\n");
    Info("-----  -----  -----  ------  -------  -------\n");
    for (i = 0; i <= g_uInternalNodeIndex; ++i)
        Info("%5u  %5u  %5u  %6.2g  %6.2g  %6.2g\n",
             i,
             g_uLeft[i],
             g_uRight[i],
             g_Height[i],
             g_LeftLength[i],
             g_RightLength[i]);
}
#endif
/* ifdef UNUSED */

/**
 * @brief Creates a UPGMA in O(N^2) tree from given distmat
 *
 * @param[out] tree
 * newly created rooted UPGMA tree
 * @param[in] distmat
 * distance matrix to be clustered
 * @param[in] linkage
 * linkage type
 * @param[in] names
 * leaf names, will be copied
 *
 * @note called UPGMA2() in Muscle3.7.
 * caller has to free with free_tree()
 *
 * @see free_tree()
 */
void
prc_upgma(tree_t *tree, matrix_t *distmat, linkage_t linkage, char **names)
{
    int i, j;
    uint *Ids;

    /* only works on full symmetric matrices */
    assert (distmat->nrows==distmat->ncols);
   
    g_uLeafCount = distmat->ncols;    
    g_ulTriangleSize = (g_uLeafCount*(g_uLeafCount - 1))/2;
    g_uInternalNodeCount = g_uLeafCount - 1;

    g_Dist = (dist_t *) malloc(g_ulTriangleSize * sizeof(dist_t));

    g_uNodeIndex = (uint*) malloc(sizeof(uint) * g_uLeafCount);
    g_uNearestNeighbor = (uint*) malloc(sizeof(uint) * g_uLeafCount);
    g_MinDist = (dist_t *) malloc(sizeof(dist_t) * g_uLeafCount);
    Ids = (uint*) malloc(sizeof(uint) * g_uLeafCount);
    /* NOTE: we replaced Names with argument names */

    /**
     * left and right node indices, as well as left and right
     * branch-length and height for for internal nodes
     */
    g_uLeft =  (uint*) malloc(sizeof(uint) * g_uInternalNodeCount);
    g_uRight =  (uint*) malloc(sizeof(uint) * g_uInternalNodeCount);
    g_Height =  (dist_t*) malloc(sizeof(dist_t) * g_uInternalNodeCount);
    g_LeftLength =  (dist_t*) malloc(sizeof(dist_t) * g_uInternalNodeCount);
    g_RightLength =  (dist_t*) malloc(sizeof(dist_t) * g_uInternalNodeCount);
    
    for (i = 0; i < g_uLeafCount; ++i) {
        g_MinDist[i] = BIG_DIST;
        g_uNodeIndex[i] = i;
        g_uNearestNeighbor[i] = uInsane;
        Ids[i] = i;
    }
    
    for (i = 0; i < g_uInternalNodeCount; ++i) {
        g_uLeft[i] = uInsane;
        g_uRight[i] = uInsane;
        g_LeftLength[i] = BIG_DIST;
        g_RightLength[i] = BIG_DIST;
        g_Height[i] = BIG_DIST;
    }
    
/* Compute initial NxN triangular distance matrix.
 * Store minimum distance for each full (not triangular) row.
 * Loop from 1, not 0, because "row" is 0, 1 ... i-1,
 * so nothing to do when i=0.
 */
    for (i = 1; i < g_uLeafCount; ++i) {
        dist_t *Row = g_Dist + TriangleSubscript(i, 0);
        CalcDistRange(distmat, i, Row);
        for (j = 0; j < i; ++j) {
            const dist_t d = Row[j];
            if (d < g_MinDist[i]) {
                g_MinDist[i] = d;
                g_uNearestNeighbor[i] = j;
            }
            if (d < g_MinDist[j]) {
                g_MinDist[j] = d;
                g_uNearestNeighbor[j] = i;
            }
        }
    }



    for (g_uInternalNodeIndex = 0;
         g_uInternalNodeIndex < g_uLeafCount - 1;
         ++g_uInternalNodeIndex) {

        dist_t dtNewMinDist = BIG_DIST;
        uint uNewNearestNeighbor = uInsane;



        /* Find nearest neighbors */
        uint Lmin = uInsane;
        uint Rmin = uInsane;
        dist_t dtMinDist = BIG_DIST;
        for (j = 0; j < g_uLeafCount; ++j) {
            dist_t d;
            if (uInsane == g_uNodeIndex[j])
                continue;

            d = g_MinDist[j];
            if (d < dtMinDist) {
                dtMinDist = d;
                Lmin = j;
                Rmin = g_uNearestNeighbor[j];
                assert(uInsane != Rmin);
                assert(uInsane != g_uNodeIndex[Rmin]);
            }
        }

        assert(Lmin != uInsane);
        assert(Rmin != uInsane);
        assert(dtMinDist != BIG_DIST);

#if TRACE
        Info("Nearest neighbors Lmin %u[=%u] Rmin %u[=%u] dist %.3g\n",
             Lmin,
             g_uNodeIndex[Lmin],
             Rmin,
             g_uNodeIndex[Rmin],
             dtMinDist);
#endif

        /* Compute distances to new node
         * New node overwrites row currently assigned to Lmin
         */
        for ( j = 0; j < g_uLeafCount; ++j) {
            ulong vL, vR;
            dist_t dL, dR;
            dist_t dtNewDist;
            
            if (j == Lmin || j == Rmin)
                continue;
            if (uInsane == g_uNodeIndex[j])
                continue;

            vL = TriangleSubscript(Lmin, j);
            vR = TriangleSubscript(Rmin, j);
            dL = g_Dist[vL];
            dR = g_Dist[vR];
            dtNewDist = 0.0;

            switch (linkage) {
            case LINKAGE_AVG:
                dtNewDist = AVG(dL, dR);
                break;

            case LINKAGE_MIN:
                dtNewDist = MIN(dL, dR);
                break;

            case LINKAGE_MAX:
                dtNewDist = MAX(dL, dR);
                break;
/* couldn't be arsed to figure out proper usage of g_dSUEFF */
#if 0
            case LINKAGE_BIASED:
                dtNewDist = g_dSUEFF*AVG(dL, dR) + (1 - g_dSUEFF)*MIN(dL, dR);
                break;
#endif
            default:
                printf("UPGMA2: Invalid LINKAGE_%u", linkage);
            }

            /* Nasty special case.
             * If nearest neighbor of j is Lmin or Rmin, then make the new
             * node (which overwrites the row currently occupied by Lmin)
             * the nearest neighbor. This situation can occur when there are
             * equal distances in the matrix. If we don't make this fix,
             * the nearest neighbor pointer for j would become invalid.
             * (We don't need to test for == Lmin, because in that case
             * the net change needed is zero due to the change in row
             * numbering).
             */
            if (g_uNearestNeighbor[j] == Rmin)
                g_uNearestNeighbor[j] = Lmin;


            g_Dist[vL] = dtNewDist;
            if (dtNewDist < dtNewMinDist) {
                dtNewMinDist = dtNewDist;
                uNewNearestNeighbor = j;
            }
        }

        assert(g_uInternalNodeIndex < g_uLeafCount - 1 || BIG_DIST != dtNewMinDist);
        assert(g_uInternalNodeIndex < g_uLeafCount - 1 || uInsane != uNewNearestNeighbor);

        const ulong v = TriangleSubscript(Lmin, Rmin);
        const dist_t dLR = g_Dist[v];
        const dist_t dHeightNew = dLR/2;
        const uint uLeft = g_uNodeIndex[Lmin];
        const uint uRight = g_uNodeIndex[Rmin];
        const dist_t HeightLeft =
            uLeft < g_uLeafCount ? 0 : g_Height[uLeft - g_uLeafCount];
        const dist_t HeightRight =
            uRight < g_uLeafCount ? 0 : g_Height[uRight - g_uLeafCount];

        g_uLeft[g_uInternalNodeIndex] = uLeft;
        g_uRight[g_uInternalNodeIndex] = uRight;
        g_LeftLength[g_uInternalNodeIndex] = dHeightNew - HeightLeft;
        g_RightLength[g_uInternalNodeIndex] = dHeightNew - HeightRight;
        g_Height[g_uInternalNodeIndex] = dHeightNew;

        /* Row for left child overwritten by row for new node */
        g_uNodeIndex[Lmin] = g_uLeafCount + g_uInternalNodeIndex;
        g_uNearestNeighbor[Lmin] = uNewNearestNeighbor;
        g_MinDist[Lmin] = dtNewMinDist;

        /* Delete row for right child */
        g_uNodeIndex[Rmin] = uInsane;


    }

    uint uRoot = g_uLeafCount - 2;

    // printf("uRoot=%d g_uLeafCount=%d g_uInternalNodeCount=%d\n", uRoot, g_uLeafCount, g_uInternalNodeCount);
    // for (i=0; i<g_uInternalNodeCount; i++) {
    //     printf("internal node=%d:  g_uLeft=%d g_uRight=%d g_LeftLength=%f g_RightLength=%f g_Height=%f\n",
    //       i, g_uLeft[i], g_uRight[i],
    //       g_LeftLength[i], g_RightLength[i],
    //       g_Height[i]);
    // }
    // for (i=0; i<g_uLeafCount; i++) {
    //     printf("leaf node=%d:  Ids=%d names=%s\n",
    //               i, Ids[i], names[i]);
    // }

#if TRACE

    printf("uRoot=%d g_uLeafCount=%d g_uInternalNodeCount=%d", uRoot, g_uLeafCount, g_uInternalNodeCount);
    for (i=0; i<g_uInternalNodeCount; i++) {
        printf("internal node=%d:  g_uLeft=%d g_uRight=%d g_LeftLength=%f g_RightLength=%f g_Height=%f",
                  i, g_uLeft[i], g_uRight[i],
                  g_LeftLength[i], g_RightLength[i],
                  g_Height[i]);
    }
    for (i=0; i<g_uLeafCount; i++) {
        printf("leaf node=%d:  Ids=%d names=%s",
                  i, Ids[i], names[i]);
    }
#endif

    tree_create(tree, g_uLeafCount, uRoot,
                      g_uLeft, g_uRight,
                      g_LeftLength, g_RightLength,
                      Ids, names);



    free(g_Dist);

    free(g_uNodeIndex);
    free(g_uNearestNeighbor);
    free(g_MinDist);
    free(g_Height);

    free(g_uLeft);
    free(g_uRight);
    free(g_LeftLength);
    free(g_RightLength);

    /* NOTE: Muscle's "Names" variable is here the argument "names" */
    free(Ids);
}
/***   end of UPGMA2   ***/
