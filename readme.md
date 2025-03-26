## Introduction
* This repository provides implementations of the algorithms introduced in `Random Sampling over Spatial Range Joins`.

[![MIT licensed](https://img.shields.io/badge/license-MIT-yellow.svg)](https://github.com/amgt-d1/IRS-interval/blob/main/license.txt)

## Requirement
* Linux OS (Ubuntu).
   * The others have not been tested.
* `g++ 11.3.0` (or higher version).

## How to use
* Parameter configuration can be done via txt files in `parameter` directory.
* Dataset should be stored in `dataset` directory.
	* We assign a unique dataset ID for each dataset. You can freely assign it.
	* In `input_data()` of `utils.hpp`, you can freely write codes for reading your dataset.
	* The file format is `x,y`.
    * For example, if you want to test on CaStreet, set `dataset_id = 0`  in `parameter` directory.
    * Foursquare, IMIS, and NYC cannot be uploaded (due to size limitation), whereas CaStreet.zip needs unzip.
* Compile: `g++ -o xxx.out -O3 main.cpp`.
* Run: `./xxx.out`

## Citation
If you use our implementation, please cite the following paper.
``` 
@inproceedings{amagata2024independent,  
    title={Random Sampling over Spatial Range Joins},  
    author={Amagata, Daichi},  
    booktitle={ICDE},  
    pages={xxx-xxx},  
    year={2025}  
}
``` 
