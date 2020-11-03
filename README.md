# tree-zig-zag
This repository contains zig-zag and Metropolis-Hastings samplers for the Kingman coalescent under finite or infinite sites mutation, along with a hybrid sampler combining both zig-zag and Metropolis-Hastings dynamics.
Specifications of all six methods are available in the preprint at <https://arxiv.org/abs/2004.08807>.

## Dependencies

tree-zig-zag has been tested on Kubuntu 18.04.
It requires the following:
- the g++ compiler (tested on version 7.5.0),
- the Gnu Scientific Library (tested on version 2.4),
- R (tested on version 4.0.3 using RStudio 1.3.1073).

The `post_processing.R` script requires two further non-core libraries:
- mcmcse (tested on version 1.4-1) available at <https://cran.r-project.org/web/packages/mcmcse/index.html>,
- tictoc (tested on version 1.0) available at <https://cran.r-project.org/web/packages/tictoc/index.html>.

## Tutorial
To replicate the analysis of <https://arxiv.org/abs/2004.08807>, execute the following steps in order.

### 1. Compile the simulators

Navigate to each of the following folders in turn:
- `Metropolis/Finite_sites/`
- `Metropolis/Infinite_sites/`
- `Zig-zag/Finite_sites/`
- `Zig-zag/Infinite_sites/`

Within each folder, call `make` to compile the corresponding simulator.
Compilation should be virtually instantaneous.

### 2. Run the simulations

Navigate to the root of the tree-zig-zag project and call 

    sh run.sh 2> acceptancerates-runtimes.txt

to execute the `run.sh` shell script.
This script runs each of the simulations specified in the `Configs` folder, and stores the resulting MCMC output into a new folder called `Results`.
The newly created `acceptancerates-runtimes.txt` file will contain the run times of each simulation, as well as the empirical acceptance rates of all Metropolis-Hastings moves.

**Warning**: The serial runtime of `run.sh` is around 100 hours on a mid-range laptop. However, all calls within the script are independent, and can be executed in parallel. The runtimes of the individual calls vary from 20 seconds to 30 hours.

### 3. Generate plots and effective sample sizes

Set the `working_dir` variable according to the commented instructions at the top of the `post_processing.R` script, and then execute the script. It creates PNG trace plots and txt files containing estimates of effective sample sizes into your R working directory. Filenames specify which plot corresponds to which simulation, and effective sample sizes are for mutation rates in the left column, tree heights in the right, and the three rows correspond to the zig-zag, hybrid, and Metropolis-Hastings methods in that order.

The script should execute in around 8 minutes on a mid-range laptop.

## Data sets and data simulation

All six data sets used in the simulation study are available in the `Data` folder.
Each row in each data file is an observed type.
The number in the last column specifies the number of times the type appeared in the sample, while the other columns specify the pattern of mutations carried by the type.

Four of the six data sets were simulated specifically for this study. The code for simulating finite sites data is available at `Metropolis/Finite_sites/`, and that for infinite sites data is at `Metropolis/Infinite_sites/`.
To compile a data simulator, navigate to the folder for the desired mutation model and call `make gen`. 
The compiled data simulator can then be called from the command line as

    ./gen_data <sample size> <number of sites> <mutation rate>

in the case of the finite sites model, and

    ./gen_data <sample size> <mutation rate>

in the case of the infinite sites model.
The output of either call is in the same format as the stored data files in the `Data` folder.
