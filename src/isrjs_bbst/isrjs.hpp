#include "../utils/weighted_sampling.hpp"
#include <unordered_set>
#include <array>
#include "bbst.hpp"


class isrjs
{
    // grid (each cell has count, i.e., #points in the cell)
    std::unordered_map<std::string, std::vector<unsigned int>> grid;
    std::unordered_map<std::string, std::vector<point>> grid_x_sorted, grid_y_sorted;
    std::unordered_map<std::string, bucket_bst> grid_bbst;

    // arrays for maintaining (upper-bound of) range counts
    std::vector<unsigned int> upper_bound_array;

    std::vector<std::array<unsigned int, 9>> alias_local;

    // alias structure
    alias Alias;

    pcg32 mt;

    // cell key computation
    std::string get_key(const unsigned int x_key, const unsigned int y_key, const int cell_number)
    {
        std::string key = "";
        switch (cell_number)
        {
            case 1:
                key = std::to_string(x_key - 1) + "+" + std::to_string(y_key - 1);
                break;

            case 2:
                key = std::to_string(x_key - 1) + "+" + std::to_string(y_key);

                break;
            case 3:
                key = std::to_string(x_key - 1) + "+" + std::to_string(y_key + 1);

                break;
            case 4:
                key = std::to_string(x_key) + "+" + std::to_string(y_key - 1);

                break;
            case 5:
                key = std::to_string(x_key) + "+" + std::to_string(y_key);

                break;
            case 6:
                key = std::to_string(x_key) + "+" + std::to_string(y_key + 1);

                break;
            case 7:
                key = std::to_string(x_key + 1) + "+" + std::to_string(y_key - 1);

                break;
            case 8:
                key = std::to_string(x_key + 1) + "+" + std::to_string(y_key);

                break;
            case 9:
                key = std::to_string(x_key + 1) + "+" + std::to_string(y_key + 1);

                break;
            
            default:
                break;
        }

        return key;
    }

    // grid mapping
    void map_grid()
    {
        /******************/
        /* pre-processing */
        /******************/
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        const unsigned int size = dataset_b.size();
        point* point_set = new point[size];
        for (unsigned int i = 0; i < size; ++i)
        {
            point_set[i].idx = i;
            point_set[i].x = dataset_b[i].first;
            point_set[i].y = dataset_b[i].second;
        }

        // sort by x-coordinate
        std::sort(point_set, point_set + size, compare_x);

        // sort dataset_a (for duplication check)
        std::sort(dataset_a.begin(), dataset_a.end());

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_preprocess = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " Pre-processing time: " << time_preprocess << "[msec]\n";

        /****************/
        /* Grid mapping */
        /****************/
        mem = process_mem_usage();
        start = std::chrono::system_clock::now();

        const unsigned int size_b = dataset_b.size();
        for (unsigned int i = 0; i < size_b; ++i)
        {
            // compute key
            std::string key = std::to_string((unsigned int)(point_set[i].x / range)) + "+" + std::to_string((unsigned int)(point_set[i].y / range));
            grid[key].push_back(point_set[i].idx);
            grid_x_sorted[key].push_back(point_set[i]);
            grid_y_sorted[key].push_back(point_set[i]);
        }

        end = std::chrono::system_clock::now();
        time_mapping = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " grid mapping time: " << time_mapping << "[msec]\t";
        time_total += time_mapping;

        mem = process_mem_usage() - mem;
        std::cout << " Memory: " << mem << "[MB]\t grid size: " << grid.size() << "\n";

        // delete unnecessary space
        delete [] point_set;
    }

    // build bbforst
    void build_bbforst()
    {
        double _mem = process_mem_usage();
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // building bbforst
        bucket_bst temp;
        auto grid_iterator = grid.begin();
        while (grid_iterator != grid.end())
        {
            // sort grid_y_sorted
            std::sort(grid_y_sorted[grid_iterator->first].begin(), grid_y_sorted[grid_iterator->first].end(), compare_y);

            // init bbst of this grid cell
            grid_bbst[grid_iterator->first] = temp;

            // build bbsts for this grid cell
            grid_bbst[grid_iterator->first].build(grid_iterator->second);

            ++grid_iterator;
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        double time_build = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        time_total += time_build;
        std::cout << " bbforest building time: " << time_build << "[msec]\t";
        _mem = process_mem_usage() - _mem;
        mem += _mem;
        std::cout << " Memory: " << mem << "[MB]\n";
    }

    // add count in the cell
    unsigned int add_count(const std::string &key, const unsigned int case_number, const float x, const float y)
    {
        unsigned int count = 0;
        point p;

        switch (case_number)
        {
            case 1:  // left bottom
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    count = itr->second.approximate_range_counting(1, x - range, y - range);
                }

                break;
            }
            case 2:  // left
            {
                auto itr = grid_x_sorted.find(key);
                if (itr != grid_x_sorted.end())
                {
                    p.x = x - range;
                    auto it = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_x);
                    count = itr->second.end() - it;
                }

                break;
            }
            case 3:  // left top
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    count = itr->second.approximate_range_counting(3, x - range, y + range);
                }

                break;
            }
            case 4:  // bottom
            {
                auto itr = grid_y_sorted.find(key);
                if (itr != grid_y_sorted.end())
                {
                    p.y = y - range;
                    auto it = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y);
                    count = itr->second.end() - it;
                }

                break;
            }
            case 5:  // center
            {
                auto itr = grid.find(key);
                if (itr != grid.end()) count = itr->second.size();

                break;
            }
            case 6:  // top
            {
                auto itr = grid_y_sorted.find(key);
                if (itr != grid_y_sorted.end())
                {
                    p.y = y + range;
                    auto it = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_y);
                    count = it - itr->second.begin();
                }

                break;
            }
            case 7:  // right bottom
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    count = itr->second.approximate_range_counting(7, x + range, y - range);
                }

                break;
            }
            case 8:  // right
            {
                auto itr = grid_x_sorted.find(key);
                if (itr != grid_x_sorted.end())
                {
                    p.x = x + range;
                    auto it = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_x);
                    count = it - itr->second.begin();
                }

                break;
            }
            case 9:  // right top
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    count = itr->second.approximate_range_counting(9, x + range, y + range);
                }

                break;
            }
            default:
                break;
        }

        return count;
    }

    // upper-bounding
    void upper_bounding()
    {
        double _mem = process_mem_usage();
        unsigned int cnt = 0;

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        const unsigned int size_a = dataset_a.size();
        upper_bound_array.resize(size_a);
        alias_local.resize(size_a);

        for (unsigned int i = 0; i < size_a; ++i)
        {
            if (i > 0 && dataset_a[i - 1] == dataset_a[i])  // duplication check
            {
                upper_bound_array[i] = upper_bound_array[i - 1];
                alias_local[i] = alias_local[i - 1];
                ++cnt;
            }
            else
            {
                // compute key
                const unsigned int x_key = (unsigned int)(dataset_a[i].first / range);
                const unsigned int y_key = (unsigned int)(dataset_a[i].second / range);

                unsigned int cnt_1 = add_count(get_key(x_key, y_key, 1), 1, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_1;
                alias_local[i][0] = upper_bound_array[i];

                unsigned int cnt_2 = add_count(get_key(x_key, y_key, 2), 2, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_2;
                alias_local[i][1] = upper_bound_array[i];

                unsigned int cnt_3 = add_count(get_key(x_key, y_key, 3), 3, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_3;
                alias_local[i][2] = upper_bound_array[i];

                unsigned int cnt_4 = add_count(get_key(x_key, y_key, 4), 4, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_4;
                alias_local[i][3] = upper_bound_array[i];

                unsigned int cnt_5 = add_count(get_key(x_key, y_key, 5), 5, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_5;
                alias_local[i][4] = upper_bound_array[i];

                unsigned int cnt_6 = add_count(get_key(x_key, y_key, 6), 6, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_6;
                alias_local[i][5] = upper_bound_array[i];

                unsigned int cnt_7 = add_count(get_key(x_key, y_key, 7), 7, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_7;
                alias_local[i][6] = upper_bound_array[i];

                unsigned int cnt_8 = add_count(get_key(x_key, y_key, 8), 8, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_8;
                alias_local[i][7] = upper_bound_array[i];

                unsigned int cnt_9 = add_count(get_key(x_key, y_key, 9), 9, dataset_a[i].first, dataset_a[i].second);
                upper_bound_array[i] += cnt_9;
                alias_local[i][8] = upper_bound_array[i];
            }
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_upperbounding = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " upper-bounding time: " << time_upperbounding << "[msec]\t";
        time_total += time_upperbounding;

        unsigned long long all = 0;
        for (unsigned int i = 0; i < size_a; ++i) all += upper_bound_array[i];
        std::cout << " upper-bound sum: " << all << "\t" << "skip count: " << cnt << "\t";

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
        Alias.build(upper_bound_array);

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_alias = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " alias building time: " << time_alias << "[msec]\t";
        time_total += time_alias;

        _mem = process_mem_usage() - _mem;
        mem += _mem;
        std::cout << " Memory: " << mem << "[MB]\n";
    }

    // random sampling
    std::pair<bool, unsigned int> random_sampling_from_cell(const std::string &key, const unsigned int case_number, unsigned int range_count, const float x, const float y)
    {
        unsigned int count = 0;
        unsigned int idx_b = 0;
        point p;
        bool flag = 0;

        switch (case_number)
        {
            case 1: // left bottom
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    std::pair<bool, unsigned int> res = itr->second.random_sampling(1, range_count, x - range, y - range);
                    if (res.first)
                    {
                        flag = 1;
                        idx_b = res.second;
                    }
                }

                break;
            }
            case 2: // left
            {
                auto itr = grid_x_sorted.find(key);
                if (itr != grid_x_sorted.end())
                {
                    p.x = x - range;
                    auto it = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_x);
                    count = itr->second.end() - it;

                    if (count > 0)
                    {
                        std::uniform_int_distribution<> rnd_weight(itr->second.size() - count, itr->second.size() - 1);
                        unsigned int rnd = rnd_weight(mt);
                        idx_b = itr->second[rnd].idx;
                        flag = 1;
                    }
                }

                break;
            }
            case 3: // left top
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    std::pair<bool, unsigned int> res = itr->second.random_sampling(3, range_count, x - range, y + range);
                    if (res.first)
                    {
                        flag = 1;
                        idx_b = res.second;
                    }
                }

                break;
            }
            case 4: // bottom
            {
                auto itr = grid_y_sorted.find(key);
                if (itr != grid_y_sorted.end())
                {
                    p.y = y - range;
                    auto it = std::lower_bound(itr->second.begin(), itr->second.end(), p, compare_y);
                    count = itr->second.end() - it;

                    if (count > 0)
                    {
                        std::uniform_int_distribution<> rnd_weight(itr->second.size() - count, itr->second.size() - 1);
                        unsigned int rnd = rnd_weight(mt);
                        idx_b = itr->second[rnd].idx;
                        flag = 1;
                    }
                }

                break;
            }
            case 5: // center
            {
                auto itr = grid.find(key);
                if (itr != grid.end())
                {
                    std::uniform_int_distribution<> rnd_weight(0, itr->second.size() - 1);
                    unsigned int rnd = rnd_weight(mt);
                    idx_b = itr->second[rnd];
                    flag = 1;
                }

                break;
            }
            case 6: // top
            {
                auto itr = grid_y_sorted.find(key);
                if (itr != grid_y_sorted.end())
                {
                    p.y = y + range;
                    auto it = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_y);
                    count = it - itr->second.begin();

                    if (count > 0)
                    {
                        std::uniform_int_distribution<> rnd_weight(0, count - 1);
                        unsigned int rnd = rnd_weight(mt);
                        idx_b = itr->second[rnd].idx;
                        flag = 1;
                    }
                }

                break;
            }
            case 7: // right bottom
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    std::pair<bool, unsigned int> res = itr->second.random_sampling(7, range_count, x + range, y - range);
                    if (res.first)
                    {
                        flag = 1;
                        idx_b = res.second;
                    }
                }

                break;
            }
            case 8: // right
            {
                auto itr = grid_x_sorted.find(key);
                if (itr != grid_x_sorted.end())
                {
                    p.x = x + range;
                    auto it = std::upper_bound(itr->second.begin(), itr->second.end(), p, compare_x);
                    count = it - itr->second.begin();

                    if (count > 0)
                    {
                        std::uniform_int_distribution<> rnd_weight(0, count - 1);
                        unsigned int rnd = rnd_weight(mt);
                        idx_b = itr->second[rnd].idx;
                        flag = 1;
                    }
                }

                break;
            }
            case 9: // right top
            {
                auto itr = grid_bbst.find(key);
                if (itr != grid_bbst.end())
                {
                    std::pair<bool, unsigned int> res = itr->second.random_sampling(9, range_count, x + range, y + range);
                    if (res.first)
                    {
                        flag = 1;
                        idx_b = res.second;
                    }
                }

                break;
            }
            default:
                break;
        }

        return {flag, idx_b};
    }

    // range sampling
    void range_sampling()
    {
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // sampling
        std::vector<std::pair<unsigned int, unsigned int>> samples;
        while (samples.size() < sample_size)
        {
            ++iteration_count;

            // get a sample from dataset_a
            const unsigned int idx_a = Alias.get_index();

            // get a weighted random cell
            std::uniform_int_distribution<> rnd_weight(0, alias_local[idx_a][8]);
            int rnd = rnd_weight(mt);
            int case_number = 0;
            for (unsigned int i = 0; i < 9; ++i)
            {
                if (alias_local[idx_a][i] >= rnd)
                {
                    case_number = i;
                    break;
                }
            }
            ++case_number;

            // compute key
            const unsigned int x_key = (unsigned int)(dataset_a[idx_a].first / range);
            const unsigned int y_key = (unsigned int)(dataset_a[idx_a].second / range);
            const std::string key = get_key(x_key, y_key, case_number);

            // get range count
            unsigned int range_count = alias_local[idx_a][case_number - 1];
            if (case_number > 1) range_count -= alias_local[idx_a][case_number - 2];

            // random sampling
            std::pair sample = random_sampling_from_cell(key, case_number, range_count, dataset_a[idx_a].first, dataset_a[idx_a].second);
            if (sample.first) samples.push_back({idx_a, sample.second});
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_sampling = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000;
        std::cout << " join sampling time: " << time_sampling << "[msec]\t iteration count: " << iteration_count << "\t time per iteration: " << (time_sampling / iteration_count) * 1000 << "[microsec]\n";
        time_total += time_sampling;
    }

public:

    // constructor
    isrjs()
    {
        // bucket_bst bbst;
    }

    // join sampling
    void join_sampling()
    {
        // grid mapping
        map_grid();

        // bbforest building
        build_bbforst();

        // upper bounding
        upper_bounding();

        // build alias based on upper-bounds
        build_alias();

        // join sampling
        range_sampling();
        
        std::cout << " total time: " << time_total << "[msec]\n\n";

        // result output
        output_result();
    }

};