
# BRKGA-QL: An Adaptive and near Parameter-free BRKGA using Reinforcement Learning

This is an implementation of the Biased Random-key Genetic Algorithm with Q-Learning (BRKGA-QL) to solve combinatorial optmization problems.

The C++ code of this algorithm has been designed to be easy of reuse. Users can only implement specific methods (decoder and local search). 


## References

When using this algorithm in academic studies, please refer to the following works:

[1] Chaves, A.A. and Lorena, L.H.N. (2020)
An Adaptive and near Parameter-free BRKGAusing Reinforcement Learning, Arcxi. 
(Available in https://...  in technical report form).


## Scope

This code has been designed to solve the Traveling Salesman Problem (TSP). Users need to configure only Decoder.cpp and LocalSearch.cpp


## Running the algorithm

* Enter the Program directory: `cd Program`
* Run the make command: `make all`
* Run the BRKGA_QL: `./runTest`

* Or compile via terminal: `g++ -std=c++11 -o rutTest BRKGA_QL.cpp Decoder.cpp LocalSearch.cpp Read.cpp -O3 `


## Code structure

The code structure is documented in [1] and organized in the following manner:

* **SPECIFIC_CODE:**
    * **Decoder.cpp**: Contains the decoders of the BRKGA-QL.
    * **LocalSearch.cpp**: Includes the local search functions, including the RVND.
    * **Read.h**: Stores the instance data.

* **GENERAL_CODE:**
    * **BRKGA_QL.cpp**: Contains all of the BRKGA-QL algorithm's population mechanisms and the main function to start the algorithm.
    * **Data.h**: Represents the data structures of BRKGA, QL and problem specific.
    * **Define.h**: Stores the method variables.
    * **Output.h**: Stores the outputs functions, including the best solution found and statistical analysis.

## File arqProblems.csv is the input data problem and each line consists in:

- Instance Name
- Run mode (0 = debug, prints in the screen; 1 = run, prints in files)
- Local search module (0 = not available; 1 = available)
- Maximum rumtime (in seconds)
- Maximum number of runs
- Number of threads used in OpenMP

Users need to create a folder named "Instances/TSP", where the berlin52.tsp instance and others must be; Users also need to create a folder named "Results" where the results files are writed.
