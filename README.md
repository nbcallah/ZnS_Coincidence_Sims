# ZnS_Coincidence_Sims

## Synopsis
These simulations create synthetic ZnS detector pulses that can be analyzed to measure rate-dependent effects. Pulses are created from resampled data from the 2016-2017 run cycle. Estimates are made of the background and the efficiency as a function of rate and the shift in lifetime is estimated from those measurements.

## Inserting new Coincidence Routine
The code is mostly complete; new coincidence functions can be tested by inserting into the code. The function `sumCtsCorrected` should be replaced. Synthetic runs are provided as a vector of `struct event`, defined as follows:

```
struct event {
	double realtime;
	int ch;
};
```

The function should look like:

```
double sumCtsCorrected(std::vector<event> data, int parameter1, float parameter2, ...) {
    double cts;
    for ( ; ; ) {
        ...
        cts += 1;
    }
    ..
    // Add in DT correction
    cts += deadtime_correction;
    ...
    // Subtract Pileup correction
    cts -= pileup_correction;
    return cts;
}
```

# Usage
The program should be compileable by using `make` on the commandline. If mpic++/mpicc are available, a paralellized version will be made. If not, a single-core version will be made. After compilation, the simulation can be run via

```
mpirun ./randomGenerator-2016-2017_Efficiency --precision=0.0001 | ./calcCorrection.py
 -- or --
./randomGenerator-2016-2017_Efficiency --precision=0.0001 | ./calcCorrection.py
```
A precision of 0.0001 is sufficient to get a shift in lifetime precision of approximately 0.01 s. Status updates should be printed via stderr, and after termination of the simulation, an estimate of the shift in lifetime and approximate error of the shift is printed on the command line from the python program.