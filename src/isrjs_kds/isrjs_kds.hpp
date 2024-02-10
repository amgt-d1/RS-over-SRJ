#include "../utils/utils.hpp"
#include "../utils/kdtree.hpp"
#include "../utils/weighted_sampling.hpp"
#include <unordered_map>


class isrjs_kds
{
    // data format for kd-tree
    point* datapoints = NULL;

    // root node of kd-tree
    Node* root = NULL;

    // arrays for maintaining range counts
    std::vector<unsigned int> count_array;

    // alias structure
    alias Alias;
    pcg32 mt;


    // range counting
    unsigned int range_counting(const float x_min, const float x_max, const float y_min, const float y_max)
    {
        unsigned int count = 0;
        CountKDtree(root, count, x_min, x_max, y_min, y_max);
        return count;
    }

    // range counting for all
    void all_range_counting()
    {
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        const unsigned int size_a = dataset_a.size();

        // init count array
        count_array.resize(size_a);

        // counting
        for (unsigned int i = 0; i < size_a; ++i)
        {
            if (i > 0 && dataset_a[i] == dataset_a[i-1])
            {
                count_array[i] = count_array[i - 1];
            }
            else
            {
                count_array[i] = range_counting(dataset_a[i].first - range, dataset_a[i].first + range, dataset_a[i].second - range, dataset_a[i].second + range);
            }
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_upperbounding = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " range counting time: " << time_upperbounding << "[msec]\t";
        time_total += time_upperbounding;

        unsigned long long all = 0;
        for (unsigned int i = 0; i < size_a; ++i) all += count_array[i];
        std::cout << " count sum: " << all << "\n";
    }

    // build alias structure for dataset_a
    void build_alias()
    {
        double _mem = process_mem_usage();

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // building an alias structure
        std::vector<unsigned int> arr;

        const unsigned int size_a = dataset_a.size();
        for (unsigned int i = 0; i < size_a; ++i) arr.push_back(count_array[i]);
        Alias.build(arr);

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_alias = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " alias building time: " << time_alias << "[msec]\t";
        time_total += time_alias;

        _mem = process_mem_usage() - _mem;
        mem += _mem;
        std::cout << " Memory: " << mem << "[MB]\n";
    }

    // range sampling
    unsigned int range_sampling(const unsigned idx_a)
    {
        // set of overlapping nodes
        std::vector<std::pair<Node*, bool>> result;
        SearchKDtree(root, result, dataset_a[idx_a].first - range, dataset_a[idx_a].first + range, dataset_a[idx_a].second - range, dataset_a[idx_a].second + range);      

        if (result.size() > 0)
        {
            // build an alias
            std::vector<unsigned int> arr;
            for (unsigned int i = 0; i < result.size(); ++i)
            {
                if (result[i].second == 0)
                {
                    arr.push_back(1);
                }
                else
                {
                    arr.push_back(result[i].first->right_idx - result[i].first->left_idx + 1);
                }
            }
            alias b(arr);

            // sample a node
            const unsigned int idx = b.get_index();
            Node* v = result[idx].first;

            if (result[idx].second == 1)
            {
                if (v->leftChild == NULL && v->rightChild == NULL)
                {
                    return v->location.idx;
                }
                else
                {
                    std::uniform_int_distribution<> rnd(result[idx].first->left_idx, result[idx].first->right_idx);
                    const unsigned int _idx = rnd(mt);
                    return datapoints[_idx].idx;
                }
            }
            else
            {
                return v->location.idx;
            }
        }

        return 0;
    }


public:

    // constructor
    isrjs_kds()
    {
        const unsigned int size_b = dataset_b.size();

        // prepare new format
        datapoints = new point[size_b];
        for (unsigned int i = 0; i < size_b; ++i)
        {
            datapoints[i].idx = i;
            datapoints[i].x = dataset_b[i].first;
            datapoints[i].y = dataset_b[i].second;
        }

        mem = process_mem_usage();
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // sort dataset_a (for duplication check)
        std::sort(dataset_a.begin(), dataset_a.end());

        // build a kd-tree fpr dataset_b
        root = kdtree(datapoints, size_b - 1, 0, 0, size_b - 1);

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_preprocess = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " Pre-processing time: " << time_preprocess << "[msec]\n";

        mem = process_mem_usage() - mem;
        std::cout << " Memory: " << mem << "[MB]\n";
        std::cout << " -----------------------------\n\n";
    }


    // join sampling
    void join_sampling()
    {
        // all range counting
        all_range_counting();

        // build alias based on upper-bounds
        build_alias();

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // sampling
        std::vector<std::pair<unsigned int, unsigned int>> samples;
        while (samples.size() < sample_size)
        {
            // get a sample from dataset_a
            const unsigned int idx_a = Alias.get_index();

            // spatial range sampling
            const unsigned int idx_b = range_sampling(idx_a);

            // store join pair
            samples.push_back({idx_a, idx_b});

            ++iteration_count;
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_sampling = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " join sampling time: " << time_sampling << "[msec]\t time per iteration: " << (time_sampling / sample_size) * 1000 << "[microsec]\n";
        time_total += time_sampling;
        std::cout << " total time: " << time_total << "[msec]\n\n";

        // result output
        output_result();
    }

};