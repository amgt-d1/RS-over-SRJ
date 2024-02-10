#ifndef UTILS_H
#define UTILS_H
#define CSV_IO_NO_THREAD

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <cfloat>
#include "csv.h"


// dataset info
unsigned int dataset_id = 0;
std::vector<std::pair<float, float>> dataset, dataset_a, dataset_b;
float scalability = 0;

// sample size
unsigned int sample_size = 0;

// memory
double mem = 0;

// running time
double time_preprocess = 0;
double time_mapping = 0;
double time_upperbounding = 0;
double time_alias = 0;
double time_sampling = 0;
double time_total = 0;

unsigned int iteration_count = 0;

// query info.
float range = 0;


// compute memory usage
double process_mem_usage()
{
    double resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
            >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
            >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    resident_set = rss * page_size_kb;

	return resident_set / 1000;
}

// parameter input
void input_parameter()
{
    std::ifstream ifs_dataset_id("parameter/dataset_id.txt");
    std::ifstream ifs_sample_size("parameter/sample_size.txt");
    std::ifstream ifs_range("parameter/range.txt");
    std::ifstream ifs_scalability("parameter/scalability.txt");

    if (ifs_dataset_id.fail())
    {
        std::cout << " dataset_id.txt does not exist." << std::endl;
        std::exit(0);
    }
    if (ifs_sample_size.fail())
    {
        std::cout << " sample_size.txt does not exist." << std::endl;
        std::exit(0);
    }
    if (ifs_range.fail())
    {
        std::cout << " range.txt does not exist." << std::endl;
        std::exit(0);
    }
    if (ifs_scalability.fail())
    {
        std::cout << " scalability.txt does not exist." << std::endl;
        std::exit(0);
    }

    while (!ifs_dataset_id.eof()) { ifs_dataset_id >> dataset_id; }
    while (!ifs_sample_size.eof()) { ifs_sample_size >> sample_size; }
    while (!ifs_range.eof()) { ifs_range >> range; }
    while (!ifs_scalability.eof()) { ifs_scalability >> scalability; }
}

// data input
void input_data()
{
    std::mt19937 mt;
    std::uniform_real_distribution<> rnd_prob(0, 1.0);

    // point variable
    std::pair<float, float> point;
    double x, y;
    std::pair<float, float> domain_min = {FLT_MAX, FLT_MAX}, domain_max = {FLT_MIN, FLT_MIN};
    bool f = 1;

    // dataset input
    std::string f_name = "../dataset/";
    if (dataset_id == 1) f_name += "castreet_header.csv";
    if (dataset_id == 2) f_name += "imis_header.csv";
    if (dataset_id == 3) f_name += "nyc_taxi_header.csv";
    if (dataset_id == 4) f_name += "foursquare_header.csv";

    io::CSVReader<2> in(f_name);
    in.read_header(io::ignore_extra_column, "x", "y");

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    while(in.read_row(x, y))
    {
        point.first = x;
        point.second = y;

        if (f)
        {
            domain_min = point;
            domain_max = point;
            f = 0;
        }

        if (domain_min.first > point.first) domain_min.first = point.first;
        if (domain_min.second > point.second) domain_min.second = point.second;
        if (domain_max.first < point.first) domain_max.first = point.first;
        if (domain_max.second < point.second) domain_max.second = point.second;

        if (rnd_prob(mt) <= scalability) dataset.push_back(point);
	}

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    double time_io = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << " Data file I/O time: " << time_io << "[msec]\n";

    // normalization
    const unsigned int size = dataset.size();
    for (unsigned int i = 0; i < size; ++i)
    {
        if (domain_min.first >= 0)
        {
            dataset[i].first -= domain_min.first;
        }
        else
        {
            dataset[i].first += std::abs(domain_min.first);
        }

        if (domain_min.second >= 0)
        {
            dataset[i].second -= domain_min.second;
        }
        else
        {
            dataset[i].second += std::abs(domain_min.second);
        }

        dataset[i].first /= (domain_max.first - domain_min.first);
        dataset[i].second /= (domain_max.second - domain_min.second);

        dataset[i].first *= 10000;
        dataset[i].second *= 10000;
    }

    std::pair<float, float> domain_min_norm = dataset[0], domain_max_norm = dataset[0];
    for (unsigned int i = 1; i < size; ++i)
    {
        if (domain_min_norm.first > dataset[i].first) domain_min_norm.first = dataset[i].first;
        if (domain_min_norm.second > dataset[i].second) domain_min_norm.second = dataset[i].second;
        if (domain_max_norm.first < dataset[i].first) domain_max_norm.first = dataset[i].first;
        if (domain_max_norm.second < dataset[i].second) domain_max_norm.second = dataset[i].second;
    }

    // shuffle
    std::shuffle(dataset.begin(), dataset.end(), mt);

    // divide into two sets
    for (unsigned int i = 0; i < size; ++i)
    {
        if (rnd_prob(mt) <= 0.5)
        {
            dataset_a.push_back(dataset[i]);
        }
        else
        {
            dataset_b.push_back(dataset[i]);
        }
    }

    // change smaller dataset
    if (dataset_a.size() > dataset_b.size())
    {
        std::vector<std::pair<float, float>> temp = dataset_b;
        dataset_b = dataset_a;
        dataset_a = temp;
    }

    /* show input parameters */
    std::cout << " -------------------------------\n";
    std::cout << " dataset ID: " << dataset_id << "\n";
    std::cout << " scalability [%]: " << scalability * 100 << "\n";
    std::cout << " cardinality_a: " << dataset_a.size() << "\n";
    std::cout << " cardinality_b: " << dataset_b.size() << "\n";
    std::cout << " sample size: " << sample_size << "\n";
    std::cout << " range: " << range << "\n";
    std::cout << " domain [" << domain_min.first << ", " << domain_min.second << "] x [" << domain_max.first << ", " << domain_max.second << "]\n";
    std::cout << " normalized domain [" << domain_min_norm.first << ", " << domain_min_norm.second << "] x [" << domain_max_norm.first << ", " << domain_max_norm.second << "]\n";
    std::cout << " -------------------------------\n\n";
}

// result output
void output_result()
{
    std::string file_name = "result/";
    if (dataset_id == 1) file_name += "1_CaStreet/";
    if (dataset_id == 2) file_name += "2_IMIS/";
    if (dataset_id == 3) file_name += "3_NYC/";
    if (dataset_id == 4) file_name += "4_Foursquare/";
    file_name += "id(" + std::to_string(dataset_id) + ")_range(" + std::to_string(range) + ")_sample-size(" + std::to_string(sample_size) + ")_scalability(" + std::to_string(scalability) + ").csv";

    std::ofstream file;
    file.open(file_name.c_str(), std::ios::out | std::ios::app);

    if (file.fail())
    {
        std::cerr << " cannot open the output file." << std::endl;
        file.clear();
        return;
    }

    file << "Total time [msec]" << "," << "Memory [MB]" << "," << "Pre-processing time [msec]" << "," << "Mapping time [msec]" << "," << "Upper-bounding time [msec]" << "," << "Alias building time [msec]" << "," << "Sampling time [msec]" << "," << "#iterations" << "," << "Success ratio" << "\n";
    file << time_total << "," << mem << "," << time_preprocess << "," << time_mapping << "," << time_upperbounding << "," << time_alias << "," << time_sampling << "," << iteration_count << "," << (double)sample_size /iteration_count << "\n\n";
}

#endif
