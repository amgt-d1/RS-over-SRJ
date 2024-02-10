#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include "pcg_random.hpp"


typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;

    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);

    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}



class alias
{
private:

    /* member variables */
    double avg = 0;
    unsigned int size = 0;
    std::vector<unsigned int> array;
    std::vector<unsigned int> indices;

    pcg32 mt;
    std::uniform_int_distribution<> rnd_idx;
    std::uniform_int_distribution<> rnd_ths;

public:

    // constructor
    alias(){}
    alias(std::vector<unsigned int> &arr) { build(arr); }

    // print alias
    void print_alias()
    {
        std::cout << " index:\n [";
        for (unsigned int i = 0; i < size; ++i)
        {
            std::cout << indices[i];
            if (i < size - 1)
            {
                std::cout << ", ";
            }
            else
            {
                std::cout << "]\n\n";
            }
        }

        std::cout << " threshold:\n [";
        for (unsigned int i = 0; i < size; ++i)
        {
            std::cout << array[i];
            if (i < size - 1)
            {
                std::cout << ", ";
            }
            else
            {
                std::cout << "]\n\n";
            }
        }
    }

    // build alias
    void build(std::vector<unsigned int> &arr)
    {
        // init array, indices & average
        array = arr;
        size = arr.size();
        for (unsigned int i = 0; i < size; ++i)
        {
            indices.push_back(i);
            avg += array[i];
        }
        avg /= size;

        // init small & large
        std::vector<unsigned int> small, large;
        for (unsigned int i = 0; i < size; ++i)
        {
            if (array[i] <= avg)
            {
                small.push_back(i);
            }
            else
            {
                large.push_back(i);
            }
        }

        while (small.size() && large.size())
        {
            const unsigned int i = small.back();
            small.pop_back();
            const unsigned int j = large.back();
            indices[i] = j;
            array[j] = array[j] - (avg - array[i]);
            if (array[j] <= avg)
            {
                small.push_back(j);
                large.pop_back();
            }
        }

        // random generator initialization
        std::uniform_int_distribution<>::param_type param_idx(0, size - 1);
        rnd_idx.param(param_idx);
        std::uniform_int_distribution<>::param_type param_ths(0, avg);
        rnd_ths.param(param_ths);
    }

    // return weighted random index
    unsigned int get_index()
    {
        // get random index & threshold
        const unsigned int idx = rnd_idx(mt);
        if (array[idx] > rnd_ths(mt)) return idx;
        return indices[idx];
    }

};