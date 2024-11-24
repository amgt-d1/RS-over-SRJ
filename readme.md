## Introduction
* This repository provides implementations of the algorithms introduced in `Random Sampling over Spatial Range Joins`.

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

