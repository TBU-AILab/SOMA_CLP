# SOMA_CLP
Implementation of the SOMA-CLP in C++

Full paper: TBA

## Use in Linux
### Compilation
cmake CMakeLists.txt
make

### Run
./m_SOMA_CLP_CEC21 1 10 100 200000 30 0
#### CLI parameters:
* Test function ID <1, 10>
* Dimension size {10, 20}
* NP
* maxFES
* Number of repetitions
* Binary parameter for CEC21 test functions <0, 7>

### Output files
Based on the example above:

    info_SOMA_CLP_CEC21_d10_t1_c0.txt
    SOMA_CLP_CEC21_d10_t1_c0/

Info file contains basic statistics, and the directory contains the results of each run.
