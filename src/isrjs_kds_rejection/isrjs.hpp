#include "../utils/utils.hpp"
#include "../utils/kdtree.hpp"
#include "../utils/weighted_sampling.hpp"
#include <unordered_map>


class isrjs
{
    // data format for kd-tree
    point* datapoints = NULL;

    // root node of kd-tree
    Node* root = NULL;

    // grid (each cell has count, i.e., #points in the cell)
    std::unordered_map<std::string, std::vector<unsigned int>> grid;

    // arrays for maintaining (upper-bound of) range counts
    std::vector<unsigned int> upper_bound_array;

    // alias structure
    alias Alias;
    pcg32 mt;


    // grid mapping
    void map_grid()
    {
        double _mem = process_mem_usage();
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        const unsigned int size_b = dataset_b.size();
        for (unsigned int i = 0; i < size_b; ++i)
        {
            // compute key
            std::string key = std::to_string((unsigned int)(dataset_b[i].first / range)) + "+" + std::to_string((unsigned int)(dataset_b[i].second / range));
            grid[key].push_back(i);
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_mapping = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " grid mapping time: " << time_mapping << "[msec]\t";
        time_total += time_mapping;

        _mem = process_mem_usage() - _mem;
        mem += _mem;
        std::cout << " Memory: " << mem << "[MB]\t grid size: " << grid.size() << "\n";
    }

    // add count in the cell
    unsigned int add_count(const std::string &key)
    {
        auto itr = grid.find(key);
        if (itr != grid.end()) return itr->second.size();

        return 0;
    }

    // compute upper bound of count
    void upper_bounding()
    {
        double _mem = process_mem_usage();

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        const unsigned int size_a = dataset_a.size();
        upper_bound_array.resize(size_a);
        std::unordered_map<std::string, unsigned int> grid_count;

        for (unsigned int i = 0; i < size_a; ++i)
        {
            // compute key
            const unsigned int x_key = (unsigned int)(dataset_a[i].first / range);
            const unsigned int y_key = (unsigned int)(dataset_a[i].second / range);
            const std::string key = std::to_string(x_key) + "+" + std::to_string(y_key);

            auto it = grid_count.find(key);
            if (it == grid_count.end())
            {
                const std::string key_1 = std::to_string(x_key - 1) + "+" + std::to_string(y_key - 1);
                const std::string key_2 = std::to_string(x_key - 1) + "+" + std::to_string(y_key);
                const std::string key_3 = std::to_string(x_key - 1) + "+" + std::to_string(y_key + 1);
                const std::string key_4 = std::to_string(x_key) + "+" + std::to_string(y_key - 1);
                const std::string key_5 = std::to_string(x_key) + "+" + std::to_string(y_key);
                const std::string key_6 = std::to_string(x_key) + "+" + std::to_string(y_key + 1);
                const std::string key_7 = std::to_string(x_key + 1) + "+" + std::to_string(y_key - 1);
                const std::string key_8 = std::to_string(x_key + 1) + "+" + std::to_string(y_key);
                const std::string key_9 = std::to_string(x_key + 1) + "+" + std::to_string(y_key + 1);

                unsigned int cnt = add_count(key_1);
                cnt += add_count(key_2);
                cnt += add_count(key_3);
                cnt += add_count(key_4);
                cnt += add_count(key_5);
                cnt += add_count(key_6);
                cnt += add_count(key_7);
                cnt += add_count(key_8);
                cnt += add_count(key_9);
                upper_bound_array[i] = cnt;

                grid_count.insert({key, cnt});
            }
            else
            {
                upper_bound_array[i] = it->second;
            }
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_upperbounding = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " upper-bounding time: " << time_upperbounding << "[msec]\t";
        time_total += time_upperbounding;

        unsigned long long all = 0;
        for (unsigned int i = 0; i < size_a; ++i) all += upper_bound_array[i];
        std::cout << " upper-bound sum: " << all << "\t";

        _mem = process_mem_usage() - _mem;
        mem += _mem;
        std::cout << " Memory: " << mem << "[MB]\n";
    }

    // build alias structure for dataset_a
    void build_alias()
    {
        double _mem = process_mem_usage();

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // building an alias structure
        std::vector<unsigned int> arr;

        const unsigned int size_a = dataset_a.size();
        for (unsigned int i = 0; i < size_a; ++i) arr.push_back(upper_bound_array[i]);
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
    void range_sampling()
    {
        // unsigned int iteration_count = 0;
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // sampling
        std::uniform_real_distribution<> rnd_prob(0, 1.0);
        std::vector<std::pair<unsigned int, unsigned int>> samples;

        while (samples.size() < sample_size)
        {
            ++iteration_count;

            // get a sample from dataset_a
            const unsigned int idx_a = Alias.get_index();

            // range counting
            unsigned int count = 0;
            std::vector<std::pair<Node*, bool>> result;
            SearchKDtree(root, result, dataset_a[idx_a].first - range, dataset_a[idx_a].first + range, dataset_a[idx_a].second - range, dataset_a[idx_a].second + range);      

            if (result.size() > 0)
            {
                double range_count = 0;
                std::vector<unsigned int> arr;

                for (unsigned int i = 0; i < result.size(); ++i)
                {
                    if (result[i].second == 0)
                    {
                        ++range_count;
                        arr.push_back(1);
                    }
                    else
                    {
                        range_count += result[i].first->subtree_size;
                        arr.push_back(result[i].first->right_idx - result[i].first->left_idx + 1);
                    }
                }

                if (rnd_prob(mt) <= range_count / upper_bound_array[idx_a])
                {
                    // build an alias
                    alias b(arr);

                    // sample a node
                    const unsigned int idx = b.get_index();
                    Node* v = result[idx].first;

                    if (result[idx].second == 1)
                    {
                        if (v->leftChild == NULL && v->rightChild == NULL)
                        {
                            samples.push_back({idx_a, v->location.idx});
                        }
                        else
                        {
                            std::uniform_int_distribution<> rnd(result[idx].first->left_idx, result[idx].first->right_idx);
                            const unsigned int _idx = rnd(mt);
                            samples.push_back({idx_a, datapoints[_idx].idx});
                        }
                    }
                    else
                    {
                        samples.push_back({idx_a, v->location.idx});
                    }
                }
            }
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_sampling = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " join sampling time: " << time_sampling << "[msec]\t";
        std::cout << " iteration count: " << iteration_count << "\t";
        std::cout << " success ratio: " << (double)sample_size / iteration_count << "\t";
        std::cout << " time per iteration: " << (time_sampling / iteration_count) * 1000 << "[microsec]\n";
        time_total += time_sampling;
    }

    // grid sampling
    int grid_sampling(const unsigned int seed, const unsigned idx_a)
    {
        // compute key
        const unsigned int x_key = (unsigned int)(dataset_a[idx_a].first / range);
        const unsigned int y_key = (unsigned int)(dataset_a[idx_a].second / range);
        const std::string key_1 = std::to_string(x_key - 1) + "+" + std::to_string(y_key - 1);
        const std::string key_2 = std::to_string(x_key - 1) + "+" + std::to_string(y_key);
        const std::string key_3 = std::to_string(x_key - 1) + "+" + std::to_string(y_key + 1);
        const std::string key_4 = std::to_string(x_key) + "+" + std::to_string(y_key - 1);
        const std::string key_5 = std::to_string(x_key) + "+" + std::to_string(y_key);
        const std::string key_6 = std::to_string(x_key) + "+" + std::to_string(y_key + 1);
        const std::string key_7 = std::to_string(x_key + 1) + "+" + std::to_string(y_key - 1);
        const std::string key_8 = std::to_string(x_key + 1) + "+" + std::to_string(y_key);
        const std::string key_9 = std::to_string(x_key + 1) + "+" + std::to_string(y_key + 1);

        // cell sampling
        std::vector<unsigned int> arr;
        arr.push_back(add_count(key_1));
        arr.push_back(add_count(key_2));
        arr.push_back(add_count(key_3));
        arr.push_back(add_count(key_4));
        arr.push_back(add_count(key_5));
        arr.push_back(add_count(key_6));
        arr.push_back(add_count(key_7));
        arr.push_back(add_count(key_8));
        arr.push_back(add_count(key_9));

        alias a(arr);
        auto itr = grid.find(key_1);

        const unsigned int _idx = a.get_index() + 1;
        if (_idx == 2) itr = grid.find(key_2);
        else if (_idx == 3) itr = grid.find(key_3);
        else if (_idx == 4) itr = grid.find(key_4);
        else if (_idx == 5) itr = grid.find(key_5);
        else if (_idx == 6) itr = grid.find(key_6);
        else if (_idx == 7) itr = grid.find(key_7);
        else if (_idx == 8) itr = grid.find(key_8);
        else if (_idx == 9) itr = grid.find(key_9);

        // point sampling
        std::uniform_int_distribution<> rnd_idx(0, itr->second.size() - 1);
        const unsigned int idx = itr->second[rnd_idx(mt)];

        if (dataset_a[idx_a].first - range <= dataset_b[idx].first && dataset_a[idx_a].second - range <= dataset_b[idx].second && dataset_a[idx_a].first + range >= dataset_b[idx].first && dataset_a[idx_a].second + range >= dataset_b[idx].second)
        {
            return idx;
        }
        
        return -1;
    }

public:

    // constructor
    isrjs()
    {
        // prepare new format
        const unsigned int size_b = dataset_b.size();
        datapoints = new point[size_b];
        for (unsigned int i = 0; i < size_b; ++i)
        {
            datapoints[i].idx = i;
            datapoints[i].x = dataset_b[i].first;
            datapoints[i].y = dataset_b[i].second;
        }

        // kd-tree building
        mem = process_mem_usage();
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

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
        // grid mapping
        map_grid();

        // upper-bounding
        upper_bounding();

        // build alias based on upper-bounds
        build_alias();

        // range sampling
        range_sampling();
        
        std::cout << " total time: " << time_total << "[msec]\n\n";

        // output result
        output_result();
    }

    // grid-based  join sampling (beta ver.)
    void grid_join_sampling()
    {
        // grid mapping
        map_grid();

        // upper-bounding
        upper_bounding();

        // build alias based on upper-bounds
        build_alias();

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // sampling
        std::mt19937 mt(1);
        std::uniform_real_distribution<> rnd_prob(0, 1.0);
        std::vector<std::pair<unsigned int, unsigned int>> samples;
        samples.reserve(sample_size);
        unsigned int trial_number = 0;

        while (samples.size() < sample_size)
        {
            ++trial_number;

            // get a sample from dataset_a
            const unsigned int idx_a = Alias.get_index();

            // sample from a grid
            const int idx_b = grid_sampling(trial_number, idx_a);
            if (idx_b > -1) samples.push_back({idx_a, idx_b});
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_sampling = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " join sampling time: " << time_sampling << "[msec]\t";
        std::cout << " #trials: " << trial_number << "\t";
        std::cout << " success ratio: " << (double)sample_size / trial_number << "\t";
        std::cout << " time per sample: " << (time_sampling / trial_number) * 1000 << "[microsec]\n";
        time_total += time_sampling;
        std::cout << " total time: " << time_total << "[msec]\n\n";
    }

};