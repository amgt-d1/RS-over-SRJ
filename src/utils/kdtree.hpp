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

// int min(const float a, const float b)
// {
//     if (a > b) return b;
//     return a;
// }

// int max(const float a, const float b)
// {
//     if (a < b) return b;
//     return a;
// }

// get median
void find_median(point *pointlist, unsigned int size, bool x)
{
    if (size % 2 == 0) --size;
    if (x)
    {
        std::nth_element(pointlist, pointlist + size / 2, pointlist + (size - 1), compare_x); 
        // return pointlist[size / 2].x;
    }
    else
    {
        std::nth_element(pointlist, pointlist + size / 2, pointlist + (size - 1), compare_y);    
        // return pointlist[size / 2].y;
    }
} 

// build a kd-tree
Node *kdtree(point *pointlist, int right, int depth, const unsigned int left_idx, const unsigned int right_idx)
{
    if (right < 0)
    {
        return NULL;
    }
    
    if (right == 0)
    {        
        Node *new_node = (Node*)malloc(sizeof(Node));        
        new_node->location = pointlist[right];        
        new_node->leftChild = NULL;
        new_node->rightChild = NULL;
        new_node->subtree_size = right + 1;
        new_node->leftmost = new_node->location.x;
        new_node->rightmost = new_node->location.x;
        new_node->upmost = new_node->location.y;
        new_node->downmost = new_node->location.y;

        new_node->left_idx = left_idx;
        new_node->right_idx = left_idx;
        
        return new_node;        
    }
    
        // determine divide dimension (x or y)
    int axis = depth % 2;
    
    // make a new node
    Node *new_node = (Node*)malloc(sizeof(Node));

    // set left & right indices
    new_node->left_idx = left_idx;
    new_node->right_idx = right_idx;

    // set coordinate of the space
    new_node->leftmost = pointlist[0].x;
    new_node->rightmost = pointlist[0].x;
    new_node->downmost = pointlist[0].y;
    new_node->upmost = pointlist[0].y;
    for (unsigned int i = 1; i < right + 1; ++i)
    {
        if (new_node->leftmost > pointlist[i].x) new_node->leftmost = pointlist[i].x;
        if (new_node->rightmost < pointlist[i].x) new_node->rightmost = pointlist[i].x;
        if (new_node->downmost > pointlist[i].y) new_node->downmost = pointlist[i].y;
        if (new_node->upmost < pointlist[i].y) new_node->upmost = pointlist[i].y;
    }

    // get median
    if (axis == 0)
    {
        find_median(pointlist, right + 1, 1);
    }
    else if (axis == 1)
    {
        find_median(pointlist, right + 1, 0);
    }    
    int median = right / 2;
    unsigned int median_idx = (left_idx + right_idx + 1) / 2;    
    new_node->location = pointlist[median];
    new_node->subtree_size = right + 1;

    // make child nodes
    new_node->rightChild = kdtree(pointlist + median + 1, right - (median + 1), depth + 1, median_idx + 1, right_idx);
    new_node->leftChild = kdtree(pointlist, median - 1, depth + 1, left_idx, median_idx - 1);
    
    // if (axis == 0)
    // {        
    //     if (new_node->rightChild != NULL && new_node->leftChild != NULL)
    //     {
    //         new_node->upmost = max(new_node->rightChild->upmost, new_node->leftChild->upmost);
    //         new_node->downmost = min(new_node->rightChild->downmost, new_node->leftChild->downmost);            
    //     }
    //     else if (new_node->rightChild != NULL)
    //     {
    //         new_node->upmost = new_node->rightChild->upmost;
    //         new_node->downmost = new_node->rightChild->downmost;            
    //     }
    //     else if (new_node->leftChild != NULL)
    //     {
    //         new_node->upmost = new_node->leftChild->upmost;
    //         new_node->downmost = new_node->leftChild->downmost;
    //     }
    //     else
    //     {
    //         new_node->upmost = new_node->location.y;
    //         new_node->downmost = new_node->location.y;
    //     }        
    // }
    // else
    // {        
    //     if (new_node->rightChild != NULL && new_node->leftChild != NULL)
    //     {
    //         new_node->rightmost = max(new_node->rightChild->rightmost,new_node->leftChild->rightmost);
    //         new_node->leftmost = min(new_node->rightChild->leftmost,new_node->leftChild->leftmost);            
    //     }
    //     else if (new_node->rightChild != NULL)
    //     {
    //         new_node->rightmost = new_node->rightChild->rightmost;
    //         new_node->leftmost = new_node->rightChild->leftmost;            
    //     }
    //     else if (new_node->leftChild != NULL)
    //     {
    //         new_node->rightmost = new_node->leftChild->rightmost;
    //         new_node->leftmost = new_node->leftChild->leftmost;
    //     }
    //     else
    //     {
    //         new_node->rightmost = new_node->location.x;
    //         new_node->leftmost = new_node->location.x;
    //     }        
    // }
    
    return new_node;
}

// kd-tree-based partitioning
void partition(point *pointlist, int right, int depth, const unsigned int left_idx, const unsigned int right_idx)
{
    if (right <= 0)
    {
        return;
    }
    
    // determine divide dimension (x or y)
    int axis = depth % 2;

    // get median
    if (axis == 0)
    {
        find_median(pointlist, right + 1, 1);
    }
    else if (axis == 1)
    {
        find_median(pointlist, right + 1, 0);
    }    
    int median = right / 2;
    unsigned int median_idx = (left_idx + right_idx + 1) / 2;    

    // make child nodes
    partition(pointlist + median + 1, right - (median + 1), depth + 1, median_idx + 1, right_idx);
    partition(pointlist, median - 1, depth + 1, left_idx, median_idx - 1);
}

// build a binary tree
Node *binary_tree(point *pointlist, int right, const unsigned int left_idx, const unsigned int right_idx)
{
    if (right < 0)
    {
        return NULL;
    }
    
    if (right == 0)
    {        
        Node *new_node = (Node*)malloc(sizeof(Node));        
        new_node->location = pointlist[right];        
        new_node->leftChild = NULL;
        new_node->rightChild = NULL;
        new_node->subtree_size = right + 1;
        new_node->leftmost = new_node->location.x;
        new_node->rightmost = new_node->location.x;
        new_node->upmost = new_node->location.y;
        new_node->downmost = new_node->location.y;

        new_node->left_idx = left_idx;
        new_node->right_idx = left_idx;
        
        return new_node;        
    }
    
    Node *new_node = (Node*)malloc(sizeof(Node));
    new_node->left_idx = left_idx;
    new_node->right_idx = right_idx;

    new_node->leftmost = pointlist[0].x;
    new_node->rightmost = pointlist[0].x;
    new_node->downmost = pointlist[0].y;
    new_node->upmost = pointlist[0].y;
    for (unsigned int i = 1; i < right + 1; ++i)
    {
        if (new_node->leftmost > pointlist[i].x) new_node->leftmost = pointlist[i].x;
        if (new_node->rightmost < pointlist[i].x) new_node->rightmost = pointlist[i].x;
        if (new_node->downmost > pointlist[i].y) new_node->downmost = pointlist[i].y;
        if (new_node->upmost < pointlist[i].y) new_node->upmost = pointlist[i].y;
    }
    
    int median = right / 2;
    unsigned int median_idx = (left_idx + right_idx + 1) / 2;
    
    new_node->location = pointlist[median];
    new_node->subtree_size = right + 1;

    new_node->rightChild = binary_tree(pointlist + median + 1, right - (median + 1), median_idx + 1, right_idx);
    new_node->leftChild = binary_tree(pointlist, median - 1, left_idx, median_idx - 1);
    
    return new_node;
}

// check fully covered
int is_contained(Node *region, const float sx, const float tx, const float sy, const float ty)
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

// range counting
void CountKDtree(Node *v, unsigned int &count, const float sx, const float tx, const float sy, const float ty)
{
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

// sample a point
int SamplePoint(point* points, Node *v, const float sx, const float tx, const float sy, const float ty)
{
    bool flag = 0;
    
    if (v->leftChild == NULL && v->rightChild == NULL)  // leaf node case
    {
        flag = 1;
    }
    else
    {
        if (v->leftChild != NULL)
        {
            if (is_contained(v->leftChild, sx, tx, sy, ty)) flag = 1;
        }
        if (flag == 0 && v->rightChild != NULL)
        {
            if (is_contained(v->rightChild, sx, tx, sy, ty)) flag = 1;
        }
    }

    if (flag)  // sample from this node (this node is a leaf or has a child covered by the query range)
    {
        std::uniform_int_distribution<> rnd_idx(v->left_idx, v->right_idx);
        int idx = rnd_idx(mt);
        if (sx <= points[idx].x && sy <= points[idx].y && tx >= points[idx].x && ty >= points[idx].y) return points[idx].idx;
    }
    else
    {
        bool left_exist = 1;
        if (v->leftChild != NULL)
        {
            if (v->leftChild->rightmost < sx || v->leftChild->leftmost > tx || v->leftChild->downmost > ty || v->leftChild->upmost < sy) left_exist = 0;
        }
        else
        {
            left_exist = 0;
        }

        bool right_exist = 1;
        if (v->rightChild != NULL)
        {
            if (v->rightChild->rightmost < sx || v->rightChild->leftmost > tx || v->rightChild->downmost > ty || v->rightChild->upmost < sy) right_exist = 0;
        }
        else
        {
            right_exist = 0;
        }

        std::uniform_int_distribution<> rnd_prob(0, v->subtree_size);

        if (left_exist == 1 && right_exist == 1)
        {
            unsigned int prob = rnd_prob(mt);
            if (prob <= v->leftChild->subtree_size)
            {
                SamplePoint(points, v->leftChild, sx, tx, sy, ty);
            }
            else
            {
                SamplePoint(points, v->rightChild, sx, tx, sy, ty);
            }
        }
        else if (left_exist == 1)
        {
            unsigned int prob = rnd_prob(mt);
            if (prob <= v->leftChild->subtree_size) SamplePoint(points, v->leftChild, sx, tx, sy, ty);
        }
        else if (right_exist == 1)
        {
            unsigned int prob = rnd_prob(mt);
            if (prob <= v->rightChild->subtree_size) SamplePoint(points, v->rightChild, sx, tx, sy, ty);
        }
    }
    
    return -1;
}