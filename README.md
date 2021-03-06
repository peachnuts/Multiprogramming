# Multiprogramming
This repository contains the implementation of our multiprogramming mapping algorithms, including SRB, qubit mapping partition (GSP and QHSP), and mapping transition. HA algorithm is also included to evalutate the independent executions.

## Installation
The repository is a Python package. You can install by cloning with `git` and then using Python's package manager `pip`:

``` shell
git clone https://github.com/peachnuts/Multiprogramming.git
pip install multiprogramming==0.0.1
import Multiprogramming
```
## How to use?
The user can find an example in [`src/main.py`](https://github.com/peachnuts/Multiprogramming/blob/main/src/main.py), which includes a list of different circuits.

Another example of finding the ground state energy of Deuteron using VQE can be found in [`src/VQE/deuteron_test.py`](https://github.com/peachnuts/Multiprogramming/blob/main/src/VQE/deuteron_test.py).

Here are some notes:

1. The `crosstalk_properties` should be obtained before the multiprogramming algorithm, using Simultaneous Randomized Benchmarking (SRB). An example of performing SRB is shown in [`src/crosstalk/srb_hardware.py`](https://github.com/peachnuts/Multiprogramming/blob/main/src/crosstalk/srb_hardware.py). The optimization methods to parallelize SRB experiments of multiple CNOT pairs is shown in [`src/crosstalk/neighbour_pair.py`](https://github.com/peachnuts/Multiprogramming/blob/main/src/crosstalk/neighbour_pair.py).
The format of the `crosstalk_properties` should be like, for example: 
``` python
crosstalk_properties = {(2,3):{(5,8) : 0.0317} # E(g_i|g_j) = 0.0317 with g_i = CX(2,3), g_j = CX(5,8)
``` 

2. The user should create a list of circuits that are supposed to be executed simultaneously.

3. The function [`main.circuits_schedule`](https://github.com/peachnuts/Multiprogramming/blob/9c0069ffb1d69f9d648300dd6e1c2f180914a287/src/main.py#L81) is the entry to start the multiprogramming algorithm. The last parameter is the partition method selected by the user. It can be `partition_hardware_heuristic` (QHSP) or `partition_hardware` (GSP). 

4. The quantum register size of the benchmarks should be equal to the hardware qubit numbers.

## Contact information
If you have any questions, feel free to contact [siyuan.niu@lirmm.fr](mailto:siyuan.niu@lirmm.fr?subject=[GitHub]%20Multiprogramming).
