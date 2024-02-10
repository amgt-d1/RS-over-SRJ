#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <random>
#include <algorithm>
#include "../utils/utils.hpp"
#include <unordered_map>


// definition of point (geo-spatial)
struct point
{
    unsigned int idx = 0;
    float x = 0;
    float y = 0;
};

bool compare_x(const point &a, const point &b){ return  a.x < b.x; }
bool compare_y(const point &a, const point &b){ return  a.y < b.y; }

// definition of bucket
struct bucket
{
    float x_min = FLT_MAX;
    float x_max = 0;
    float y_min = FLT_MAX;
    float y_max = 0;
    std::vector<unsigned int> idx_set;
    // std::vector<std::pair<float, float>*> idx_set;

    // inner parameter initialization
    void init()
    {
        x_min = FLT_MAX;
        x_max = 0;
        y_min = FLT_MAX;
        y_max = 0;
    }
};

bool compare_y_min_bucket(const bucket &a, const bucket &b){ return  a.y_min < b.y_min; }
bool compare_y_max_bucket(const bucket &a, const bucket &b){ return  a.y_max < b.y_max; }

// definition of node of BBST
struct node
{
    unsigned int node_id = 0;
    // bucket* b = NULL;
    float x = 0;
    node* left = NULL;
    node* right = NULL;
    unsigned int height = 0;
};

// data structure (like range-tree)
struct tree
{
    node* root = NULL;
    std::vector<node> bbst;
    unsigned int height_max = 0;

    // buckets having the same x-coordinate
    // std::unordered_map<unsigned int, std::vector<bucket*>> bucket_map;
    std::unordered_map<unsigned int, std::vector<std::pair<float, bucket*>>> bucket_map_min, bucket_map_max;

    // y-coordinate of the buckets in sub-tree
    std::unordered_map<unsigned int, std::vector<std::pair<float, bucket*>>> range_map_min, range_map_max;
};

bool compare_y_(const std::pair<float, bucket*> &a, const std::pair<float, bucket*> &b) { return a.first < b.first; }

class bucket_bst
{
    // a set of buckets
    std::vector<bucket> bucket_set;

    // bucket size = log(m)
    unsigned int bucket_size = 0;

    // trees
    tree bbst_min, bbst_max;

    pcg32 _mt;


    // make a set of buckets on a cell
    void make_buckets(std::vector<unsigned int> &cell)
    {
        const unsigned int size = dataset_b.size();
        const unsigned int cell_size = cell.size();
        bucket_size = (unsigned int)std::log2f(size);
        // std::vector<std::pair<float, float>*> idx_set;
        std::vector<unsigned int> idx_set;
        bucket_set.reserve(cell.size() / bucket_size);
        bucket b;

        for (unsigned int i = 0; i < cell_size; ++i)
        {
            std::pair<float, float>* p = &dataset_b[cell[i]];
            idx_set.push_back(cell[i]);
            if (b.x_min > p->first) b.x_min = p->first;
            b.x_max = p->first;
            if (b.y_min > p->second) b.y_min = p->second;
            if (b.y_max < p->second) b.y_max = p->second;
            
            // make a bucket
            if (idx_set.size() == bucket_size || i == cell_size - 1)
            {
                b.idx_set = idx_set;
                bucket_set.push_back(b);

                idx_set.clear();
                b.init();
            }
        }
    }
    
    // function for making a new node
    node* make_node(std::vector<unsigned int> &bucket_idx_set, std::vector<std::pair<float, bucket*>> &y_range_min, std::vector<std::pair<float, bucket*>> &y_range_max, unsigned int height, bool flag)
    {
        const unsigned int b_size = bucket_idx_set.size();
        if (b_size == 0) return NULL;

        /*******************/
        /* make a new node */
        /*******************/
        node Node;
        node* new_node = NULL;
        if (flag)
        {
            bbst_min.bbst.push_back(Node);
            new_node = &bbst_min.bbst[bbst_min.bbst.size() - 1];
            new_node->node_id = bbst_min.bbst.size() - 1;
        }
        else
        {
            bbst_max.bbst.push_back(Node);
            new_node = &bbst_max.bbst[bbst_max.bbst.size() - 1];
            new_node->node_id = bbst_max.bbst.size() - 1;
        }

        // get median
        const unsigned int median_idx = bucket_idx_set[b_size / 2];
        new_node->x = bucket_set[median_idx].x_min;
        if (!flag) new_node->x = bucket_set[median_idx].x_max;

        // get height
        if (flag)
        {
            new_node->height = height;
            if (bbst_min.height_max < new_node->height) bbst_min.height_max = new_node->height;
        }
        else
        {
            new_node->height = height;
            if (bbst_max.height_max < new_node->height) bbst_max.height_max = new_node->height;
        }

        // split
        std::vector<unsigned int> bucket_idx_left, bucket_idx_right;
        for (unsigned int i = 0; i < b_size; ++i)
        {
            const unsigned int idx = bucket_idx_set[i];

            if (flag)
            {
                if (bucket_set[idx].x_min < new_node->x)
                {
                    bucket_idx_left.push_back(bucket_idx_set[i]);
                }
                else if (bucket_set[idx].x_min > new_node->x)
                {
                    bucket_idx_right.push_back(bucket_idx_set[i]);
                }
            }
            else
            {
                if (bucket_set[idx].x_max < new_node->x)
                {
                    bucket_idx_left.push_back(bucket_idx_set[i]);
                }
                else if (bucket_set[idx].x_max > new_node->x)
                {
                    bucket_idx_right.push_back(bucket_idx_set[i]);
                }
            }
        }

        // store the data structure for y-coordinate
        if (new_node->node_id > 0)
        {
            if (flag)
            {
                bbst_min.range_map_min.insert({new_node->node_id, y_range_min});
                bbst_min.range_map_max.insert({new_node->node_id, y_range_max});
            }
            else
            {
                bbst_max.range_map_min.insert({new_node->node_id, y_range_min});
                bbst_max.range_map_max.insert({new_node->node_id, y_range_max});
            }
        }

        // if (b_size == 1) return new_node;

        // split y-values
        std::vector<std::pair<float, bucket*>> y_range_min_left, y_range_min_right, y_range_max_left, y_range_max_right;
        const unsigned int y_size = y_range_min.size();
        for (unsigned int i = 0; i < y_size; ++i)
        {
            if (flag)
            {
                if (y_range_min[i].second->x_min < new_node->x)
                {
                    y_range_min_left.push_back(y_range_min[i]);
                }
                else if (y_range_min[i].second->x_min > new_node->x)
                {
                    y_range_min_right.push_back(y_range_min[i]);
                }
                else if (y_range_min[i].second->x_min == new_node->x)
                {
                    bbst_min.bucket_map_min[new_node->node_id].push_back(y_range_min[i]);
                }

                if (y_range_max[i].second->x_min < new_node->x)
                {
                    y_range_max_left.push_back(y_range_max[i]);
                }
                else if (y_range_max[i].second->x_min > new_node->x)
                {
                    y_range_max_right.push_back(y_range_max[i]);
                }
                else if (y_range_max[i].second->x_min == new_node->x)
                {
                    bbst_min.bucket_map_max[new_node->node_id].push_back(y_range_max[i]);
                }
            }
            else
            {
                if (y_range_min[i].second->x_max < new_node->x)
                {
                    y_range_min_left.push_back(y_range_min[i]);
                }
                else if (y_range_min[i].second->x_max > new_node->x)
                {
                    y_range_min_right.push_back(y_range_min[i]);
                }
                else if (y_range_min[i].second->x_max == new_node->x)
                {
                    bbst_max.bucket_map_min[new_node->node_id].push_back(y_range_min[i]);
                }

                if (y_range_max[i].second->x_max < new_node->x)
                {
                    y_range_max_left.push_back(y_range_max[i]);
                }
                else if (y_range_max[i].second->x_max > new_node->x)
                {
                    y_range_max_right.push_back(y_range_max[i]);
                }
                else if (y_range_max[i].second->x_max == new_node->x)
                {
                    bbst_max.bucket_map_max[new_node->node_id].push_back(y_range_max[i]);
                }
            }

            if (b_size == 1) return new_node;
        }

        // make child nodes
        new_node->left = make_node(bucket_idx_left, y_range_min_left, y_range_max_left, new_node->height + 1, flag);
        new_node->right = make_node(bucket_idx_right, y_range_min_right, y_range_max_right, new_node->height + 1, flag);

        return new_node;
    }

    // build bbsts on a cell
    void build_bbst()
    {
        // reserve spaces
        const unsigned int size = bucket_set.size();
        bbst_min.bbst.reserve(size);
        bbst_max.bbst.reserve(size);

        // get bucket idx
        std::vector<unsigned int> bucket_idx_set;
        for (unsigned int i = 0; i < bucket_set.size(); ++i) bucket_idx_set.push_back(i);
        const unsigned int b_size = bucket_idx_set.size();

        // prepare data structure for y-coordinate
        std::vector<std::pair<float, bucket*>> y_range_min, y_range_max;
        for (unsigned int i = 0; i < b_size; ++i)
        {
            y_range_min.push_back({bucket_set[i].y_min, &bucket_set[i]});
            y_range_max.push_back({bucket_set[i].y_max, &bucket_set[i]});
        }
        std::sort(y_range_min.begin(), y_range_min.end());
        std::sort(y_range_max.begin(), y_range_max.end());

        unsigned int height = 0;
        bbst_min.root = make_node(bucket_idx_set, y_range_min, y_range_max, height, 1);
        bbst_max.root = make_node(bucket_idx_set, y_range_min, y_range_max, height, 0);
    }

    // approx. range counting for case 1
    void approximate_range_count_1(const node* n, unsigned int &count, const float x, const float y)
    {
        if (n->x >= x)
        {
            // check right-sub-tree
            if (n->right != NULL)
            {
                auto itr = bbst_max.range_map_max.find(n->right->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count += (itr->second.end() - itr_bs) * bucket_size;
            }

            auto itr = bbst_max.bucket_map_max.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count += (itr->second.end() - itr_bs) * bucket_size;

            if (n->x == x) return;

            // go to left child
            if (n->left != NULL) approximate_range_count_1(n->left, count, x, y);
        }
        else
        {
            // go to right child
            if (n->right != NULL) approximate_range_count_1(n->right, count, x, y);
        }
    }

    // approx. range counting for case 3
    void approximate_range_count_3(const node* n, unsigned int &count, const float x, const float y)
    {
        if (n->x >= x)
        {
            // check right-sub-tree
            if (n->right != NULL)
            {
                auto itr = bbst_max.range_map_min.find(n->right->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count += (itr_bs - itr->second.begin()) * bucket_size;
            }

            auto itr = bbst_max.bucket_map_min.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count += (itr_bs - itr->second.begin()) * bucket_size;

            if (n->x == x) return;

            // go to left child
            if (n->left != NULL) approximate_range_count_3(n->left, count, x, y);
        }
        else
        {
            // go to right child
            if (n->right != NULL) approximate_range_count_3(n->right, count, x, y);
        }
    }

    // approx. range counting for case 7
    void approximate_range_count_7(const node* n, unsigned int &count, const float x, const float y)
    {
        if (n->x <= x)
        {
            // check left sub-tree
            if (n->left != NULL)
            {
                auto itr = bbst_min.range_map_max.find(n->left->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count += (itr->second.end() - itr_bs) * bucket_size;
            }

            auto itr = bbst_min.bucket_map_max.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count += (itr->second.end() - itr_bs) * bucket_size;

            if (n->x == x) return;

            // go to right child
            if (n->right != NULL) approximate_range_count_7(n->right, count, x, y);
        }
        else if (x < n->x)
        {
            // go to left child
            if (n->left != NULL) approximate_range_count_7(n->left, count, x, y);
        }
    }

    // approx. range counting for case 9
    void approximate_range_count_9(const node* n, unsigned int &count, const float x, const float y)
    {
        if (n->x <= x)
        {
            // check left sub-tree
            if (n->left != NULL)
            {
                auto itr = bbst_min.range_map_min.find(n->left->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count += (itr_bs - itr->second.begin()) * bucket_size;
            }

            auto itr = bbst_min.bucket_map_min.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count += (itr_bs - itr->second.begin()) * bucket_size;

            if (n->x == x) return;

            // go to right child
            if (n->right != NULL) approximate_range_count_9(n->right, count, x, y);
        }
        else if (x < n->x)
        {
            // go to left child
            if (n->left != NULL) approximate_range_count_9(n->left, count, x, y);
        }
    }

    // check in range
    bool check_within(const float x, const float y, const std::pair<float, float> &p)
    {
        bool f = 0;
        if (x - range <= p.first && p.first <= x + range && y - range <= p.second && p.second <= y + range)
        {
            f = 1;
        }

        return f;
    }

    // random sampling for case 1
    void random_sampling_1(std::pair<bool, unsigned int> &res, const node* n, const unsigned int count, const float x, const float y)
    {
        if (n->x >= x)
        {
            // get a random integer
            std::uniform_int_distribution<> rnd_weight(1, count);
            unsigned int rnd = rnd_weight(_mt);

            unsigned int count_right = 0;
            unsigned int count_center = 0;

            // check buckets of this node (n->x == x or leaf node)
            auto itr = bbst_max.bucket_map_max.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count_center = (itr->second.end() - itr_bs) * bucket_size;

            // sample from this node
            if (rnd <= count_center)
            {
                std::uniform_int_distribution<> rnd_b_idx(itr_bs - itr->second.begin(), itr->second.size() - 1);
                const unsigned int r_idx = rnd_b_idx(_mt);
                bucket* b = itr->second[r_idx].second;

                std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                const unsigned int _idx = rnd_p_idx(_mt);
                bool f = check_within(x + range, y + range, dataset_b[b->idx_set[_idx]]);
                if (f) res = {1, b->idx_set[_idx]};

                return;
            }

            if (n->x == x) return;

            // check right-sub-tree
            if (n->right != NULL)
            {
                std::uniform_int_distribution<> rnd_weight_2(1, count - count_center);
                rnd = rnd_weight_2(_mt);

                auto itr = bbst_max.range_map_max.find(n->right->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count_right = (itr->second.end() - itr_bs) * bucket_size;

                // sample from this node
                if (rnd <= count_right)
                {
                    // sample a bucket
                    std::uniform_int_distribution<> rnd_b_idx(itr->second.size() - (count_right / bucket_size), itr->second.size() - 1);
                    const unsigned int r_idx = rnd_b_idx(_mt);
                    bucket* b = itr->second[r_idx].second;

                    std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                    const unsigned int _idx = rnd_p_idx(_mt);
                    bool f = check_within(x + range, y + range, dataset_b[b->idx_set[_idx]]);
                    if (f) res = {1, b->idx_set[_idx]};

                    return;
                }
            }

            // go to left child
            if (n->left != NULL) random_sampling_1(res, n->left, count - count_center - count_right, x, y);
        }
        else
        {
            // go to right child
            if (n->right != NULL) random_sampling_1(res, n->right, count, x, y);
        }
    }

    // random sampling for case 3
    void random_sampling_3(std::pair<bool, unsigned int> &res, const node* n, const unsigned int count, const float x, const float y)
    {
        // get a random integer
        std::uniform_int_distribution<> rnd_weight(1, count);
        unsigned int rnd = rnd_weight(_mt);

        unsigned int count_right = 0;
        unsigned int count_center = 0;

        if (n->x >= x)
        {
            // check right-sub-tree
            if (n->right != NULL)
            {
                auto itr = bbst_max.range_map_min.find(n->right->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count_right = (itr_bs - itr->second.begin()) * bucket_size;

                // sample from this node
                if (rnd <= count_right)
                {
                    // sample a bucket
                    std::uniform_int_distribution<> rnd_b_idx(0, count_right / bucket_size - 1);
                    const unsigned int r_idx = rnd_b_idx(_mt);
                    bucket* b = itr->second[r_idx].second;

                    std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                    const unsigned int _idx = rnd_p_idx(_mt);
                    bool f = check_within(x + range, y - range, dataset_b[b->idx_set[_idx]]);
                    if (f) res = {1, b->idx_set[_idx]};

                    return;
                }
            }

            // check bucket_map
            auto itr = bbst_max.bucket_map_min.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count_center = (itr_bs - itr->second.begin()) * bucket_size;

            std::uniform_int_distribution<> rnd_weight_2(1, count - count_right);
            rnd = rnd_weight_2(_mt);

            // sample from this node
            if (rnd <= count_center)
            {
                std::uniform_int_distribution<> rnd_b_idx(0, itr_bs - itr->second.begin() - 1);
                const unsigned int r_idx = rnd_b_idx(_mt);
                bucket* b = itr->second[r_idx].second;

                std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                const unsigned int _idx = rnd_p_idx(_mt);
                bool f = check_within(x + range, y - range, dataset_b[b->idx_set[_idx]]);
                if (f) res = {1, b->idx_set[_idx]};

                return;
            }

            if (n->x == x) return;

            // go to left child
            if (n->left != NULL) random_sampling_3(res, n->left, count - count_center - count_right, x, y);
        }
        else
        {
            // go to right child
            if (n->right != NULL) random_sampling_3(res, n->right, count, x, y);
        }
    }

    // random sampling for case 7
    void random_sampling_7(std::pair<bool, unsigned int> &res, const node* n, const unsigned int count, const float x, const float y)
    {
        // get a random integer
        std::uniform_int_distribution<> rnd_weight(1, count);
        unsigned int rnd = rnd_weight(_mt);

        unsigned int count_left = 0;
        unsigned int count_center = 0;

        if (n->x <= x)
        {
            // check left sub-tree
            if (n->left != NULL)
            {
                auto itr = bbst_min.range_map_max.find(n->left->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count_left = (itr->second.end() - itr_bs) * bucket_size;

                // sample from this node
                if (rnd <= count_left)
                {
                    // sample a bucket
                    std::uniform_int_distribution<> rnd_b_idx(itr->second.size() - (count_left / bucket_size), itr->second.size() - 1);
                    const unsigned int r_idx = rnd_b_idx(_mt);
                    bucket* b = itr->second[r_idx].second;

                    std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                    const unsigned int _idx = rnd_p_idx(_mt);
                    bool f = check_within(x - range, y + range, dataset_b[b->idx_set[_idx]]);
                    if (f) res = {1, b->idx_set[_idx]};

                    return;
                }
            }

            auto itr = bbst_min.bucket_map_max.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count_center = (itr->second.end() - itr_bs) * bucket_size;

            std::uniform_int_distribution<> rnd_weight_2(1, count - count_left);
            rnd = rnd_weight_2(_mt);

            // sample from this node
            if (rnd <= count_center)
            {
                std::uniform_int_distribution<> rnd_b_idx(itr_bs - itr->second.begin(), itr->second.size() - 1);
                const unsigned int r_idx = rnd_b_idx(_mt);
                bucket* b = itr->second[r_idx].second;

                std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                const unsigned int _idx = rnd_p_idx(_mt);
                bool f = check_within(x - range, y + range, dataset_b[b->idx_set[_idx]]);
                if (f) res = {1, b->idx_set[_idx]};

                return;
            }

            if (n->x == x) return;

            // go to right child
            if (n->right != NULL) random_sampling_7(res, n->right, count - count_center - count_left, x, y);
        }
        else if (x < n->x)
        {
            // go to left child
            if (n->left != NULL) random_sampling_7(res, n->left, count, x, y);
        }
    }

    // random sampling for case 9
    void random_sampling_9(std::pair<bool, unsigned int> &res, const node* n, const unsigned int count, const float x, const float y)
    {
        // get a random integer
        std::uniform_int_distribution<> rnd_weight(1, count);
        unsigned int rnd = rnd_weight(_mt);

        unsigned int count_left = 0;
        unsigned int count_center = 0;

        if (n->x <= x)
        {
            // check left sub-tree
            if (n->left != NULL)
            {
                auto itr = bbst_min.range_map_min.find(n->left->node_id);
                bucket* b = NULL;
                std::pair p = {y,b};
                auto itr_bs = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
                count_left = (itr_bs - itr->second.begin()) * bucket_size;

                // sample from this node
                if (rnd <= count_left)
                {
                    // sample a bucket
                    std::uniform_int_distribution<> rnd_b_idx(0, count_left / bucket_size - 1);
                    const unsigned int r_idx = rnd_b_idx(_mt);
                    bucket* b = itr->second[r_idx].second;

                    std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                    const unsigned int _idx = rnd_p_idx(_mt);
                    bool f = check_within(x - range, y - range, dataset_b[b->idx_set[_idx]]);
                    if (f) res = {1, b->idx_set[_idx]};

                    return;
                }
            }

            auto itr = bbst_min.bucket_map_min.find(n->node_id);
            bucket* b = NULL;
            std::pair p = {y,b};
            auto itr_bs = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y_);
            count_center = (itr_bs - itr->second.begin()) * bucket_size;

            std::uniform_int_distribution<> rnd_weight_2(1, count - count_left);
            rnd = rnd_weight_2(_mt);

            // sample from this node
            if (rnd <= count_center)
            {
                std::uniform_int_distribution<> rnd_b_idx(0, itr_bs - itr->second.begin() - 1);
                const unsigned int r_idx = rnd_b_idx(_mt);
                bucket* b = itr->second[r_idx].second;

                std::uniform_int_distribution<> rnd_p_idx(0, b->idx_set.size() - 1);
                const unsigned int _idx = rnd_p_idx(_mt);
                bool f = check_within(x - range, y - range, dataset_b[b->idx_set[_idx]]);
                if (f) res = {1, b->idx_set[_idx]};

                return;
            }

            if (n->x == x) return;

            // go to right child
            if (n->right != NULL) random_sampling_9(res, n->right, count - count_center - count_left, x, y);
        }
        else if (x < n->x)
        {
            // go to left child
            if (n->left != NULL) random_sampling_9(res, n->left, count, x, y);
        }
    }

public:

    // constructor (do nothing)
    bucket_bst() {}

    // build two bbsts (on a cell)
    void build(std::vector<unsigned int> &cell)
    {
        // make a set of buckets
        make_buckets(cell);

        // build BSTs on buckets
        build_bbst();
    }

    // approx. range counting
    unsigned int approximate_range_counting(const unsigned int case_number, const float x, const float y)
    {
        unsigned int count = 0;

        if (case_number == 1)
        {
            approximate_range_count_1(bbst_max.root, count, x, y);
        }
        else if (case_number == 3)
        {
            approximate_range_count_3(bbst_max.root, count, x, y);
        }
        else if (case_number == 7)
        {
            approximate_range_count_7(bbst_min.root, count, x, y);
        }
        else if (case_number == 9)
        {
            approximate_range_count_9(bbst_min.root, count, x, y);
        }

        return count;
    }

    // random sampling from tree
    std::pair<bool, unsigned int> random_sampling(const unsigned int case_number, const unsigned int count, const float x, const float y)
    {
        std::pair<bool, unsigned int> res = {0,0};
        unsigned int r_count = 0;

        if (case_number == 1)
        {
            random_sampling_1(res, bbst_max.root, count, x, y);
        }
        else if (case_number == 3)
        {
            random_sampling_3(res, bbst_max.root, count, x, y);
        }
        else if (case_number == 7)
        {
            random_sampling_7(res, bbst_min.root, count, x, y);
        }
        else if (case_number == 9)
        {
            random_sampling_9(res, bbst_min.root, count, x, y);
        }
        // std::cout << " count: " << r_count << "\n";

        return res;
    }
};