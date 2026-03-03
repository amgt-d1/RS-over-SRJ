#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <random>
#include <algorithm>
#include "pcg_random.hpp"


pcg32 mt;

// definition of point (geo-spatial)
struct point
{
    unsigned int idx = 0;
    float x = 0;
    float y = 0;
};

// definition of kd-tree node
struct Node
{    
    point location;
    Node *leftChild = NULL;
    Node *rightChild = NULL;
    unsigned int subtree_size = 0;

    unsigned int left_idx = 0;
    unsigned int right_idx = 0;
    
    float leftmost = 0;
    float rightmost = 0;
    float upmost = 0;
    float downmost = 0;
};

bool compare_x(const point &a, const point &b){ return  a.x < b.x; }
bool compare_y(const point &a, const point &b){ return  a.y < b.y; }

static void find_median(point *pts, unsigned int size, bool use_x)
{
    if (size == 0) return;
    unsigned int n = size;
    unsigned int m = n / 2;  // 0-based median
    if (use_x)
    {
        std::nth_element(pts, pts + m, pts + n, compare_x);
    }
    else
    {
        std::nth_element(pts, pts + m, pts + n, compare_y);
    }
}

Node *kdtree(point *pointlist, int right, int depth, unsigned int left_idx, unsigned int right_idx)
{
    if (right < 0)  return NULL;

    unsigned int n = (unsigned int)right + 1;
    unsigned int span = right_idx - left_idx + 1;

    if (n != span)
    {
        // debugging
        // fprintf(stderr, "kdtree: n != span (%u vs %u)\n", n, span);
        // exit(1);
    }

    // leaf node
    if (n == 1)
    {
        Node *node = (Node*)malloc(sizeof(Node));
        node->location = pointlist[0];
        node->leftChild = NULL;
        node->rightChild= NULL;

        node->left_idx = left_idx;
        node->right_idx = right_idx;
        node->subtree_size = 1;

        node->leftmost = pointlist[0].x;
        node->rightmost = pointlist[0].x;
        node->downmost = pointlist[0].y;
        node->upmost = pointlist[0].y;

        return node;
    }

    // intermediate node
    int axis = depth % 2;

    // MBR
    float lx = pointlist[0].x, rx = pointlist[0].x;
    float dy = pointlist[0].y, uy = pointlist[0].y;
    for (unsigned int i = 1; i < n; ++i)
    {
        if (lx > pointlist[i].x) lx = pointlist[i].x;
        if (rx < pointlist[i].x) rx = pointlist[i].x;
        if (dy > pointlist[i].y) dy = pointlist[i].y;
        if (uy < pointlist[i].y) uy = pointlist[i].y;
    }

    find_median(pointlist, n, axis == 0);
    unsigned int m = n / 2; 
    unsigned int mid_idx = left_idx + m;

    Node *node = (Node*)malloc(sizeof(Node));
    node->location  = pointlist[m];
    node->left_idx  = left_idx;
    node->right_idx = right_idx;
    node->subtree_size = span;

    node->leftmost  = lx;
    node->rightmost = rx;
    node->downmost  = dy;
    node->upmost    = uy;

    // left child
    unsigned int left_span = m;
    if (left_span > 0)
    {
        node->leftChild = kdtree(pointlist, (int)m - 1, depth + 1, left_idx, mid_idx - 1);
    }
    else
    {
        node->leftChild = NULL;
    }

    // right child
    unsigned int right_span = n - m - 1;
    if (right_span > 0)
    {
        node->rightChild = kdtree(pointlist + m + 1, (int)right_span - 1, depth + 1, mid_idx + 1, right_idx);
    }
    else
    {
        node->rightChild = NULL;
    }

    return node;
}

// check fully covered
int is_contained(const Node *region, const float sx, const float tx, const float sy, const float ty)
{    
    if (region->leftmost < sx || region->rightmost > tx || region->upmost > ty || region->downmost < sy) return 0;    
    return 1;
}

// orthogonal range search
void SearchKDtree(Node *v, std::vector<std::pair<Node*, bool>> &result, const float sx, const float tx, const float sy, const float ty)
{
    // disjoint case    
    if (v->rightmost < sx || v->leftmost > tx || v->downmost > ty || v->upmost < sy) return;
    
    // leaf node case
    if (v->leftChild == NULL && v->rightChild == NULL)
    {
        result.push_back({v, 0});
        return;
    }
    
    // intermediate node case
    if (sx <= v->location.x && v->location.x <= tx && sy <= v->location.y && v->location.y <= ty) result.push_back({v, 0});
    if (v->leftChild != NULL)
    {        
        if (is_contained(v->leftChild, sx, tx, sy, ty))
        {
            result.push_back({v->leftChild, 1});
        }
        else
        {
            SearchKDtree(v->leftChild, result, sx, tx, sy, ty);
        }
    }
    
    if (v->rightChild != NULL)
    {        
        if (is_contained(v->rightChild, sx, tx, sy, ty))
        {
            result.push_back({v->rightChild, 1});
        }
        else
        {
            SearchKDtree(v->rightChild, result, sx, tx, sy, ty);
        }
    }
    return;
}

// orthogonal range search
void range_search(Node *v, std::vector<unsigned int> &result, const float sx, const float tx, const float sy, const float ty)
{
    // disjoint case    
    if (v->rightmost < sx || v->leftmost > tx || v->downmost > ty || v->upmost < sy) return;
    
    // leaf node case
    if (v->leftChild == NULL && v->rightChild == NULL)
    {
        result.push_back(v->location.idx);
        return;
    }
    
    // intermediate node case
    if (sx <= v->location.x && v->location.x <= tx && sy <= v->location.y && v->location.y <= ty) result.push_back(v->location.idx);
    if (v->leftChild != NULL) range_search(v->leftChild, result, sx, tx, sy, ty);    
    if (v->rightChild != NULL) range_search(v->rightChild, result, sx, tx, sy, ty);
    
    return;
}

// range counting
void CountKDtree(const Node *v, unsigned int &count, const float sx, const float tx, const float sy, const float ty)
{
    if (!v) return;

    // disjoint case    
    if (v->rightmost < sx || v->leftmost > tx || v->downmost > ty || v->upmost < sy) return;
    
    // leaf node case
    if (v->leftChild == NULL && v->rightChild == NULL)
    {        
        if (sx <= v->location.x && sy <= v->location.y && tx >= v->location.x && ty >= v->location.y)
        {
            ++count;
        }
        return;
    }
    
    // intermediate node case
    if (sx <= v->location.x && sy <= v->location.y && tx >= v->location.x && ty >= v->location.y) ++count;

    if (v->leftChild != NULL)
    {        
        if (is_contained(v->leftChild, sx, tx, sy, ty))
        {
            count += v->leftChild->subtree_size;
        }
        else
        {
            CountKDtree(v->leftChild, count, sx, tx, sy, ty);
        }
    }
    
    if (v->rightChild != NULL)
    {        
        if (is_contained(v->rightChild, sx, tx, sy, ty))
        {
            count += v->rightChild->subtree_size;
        }
        else
        {
            CountKDtree(v->rightChild, count, sx, tx, sy, ty);
        }
    }
    return;
}
