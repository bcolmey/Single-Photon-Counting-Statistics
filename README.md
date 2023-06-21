# Single-Photon-Counting-Statistics
A Python script simulating a single photon counting experiment using a Hanbury-Brown and Twiss (HBT) interferometer. The code simulates photon emission from a semiconductor quantum dot and its subsequent detection and processing, including calculations of the source's photon intensity autocorrelation function g²(τ). Realistic experimental parameters such as detector response functions, detector dead time, dark counts, and multiphoton emission are included. This program may serve as a digital companion to physical experimentation, helping to predict, interpret, and cross-verify experimental results. A full description of the code can be found here:
https://medium.com/@benjamincolmey/simulating-single-photon-emission-a-statistical-approach-in-python-aa348cc8119c

The entirety of the code can be found here and was originally undertaken in collaboration with James Godfrey at Queen’s University during the summer of 2020 while Covid collectively kept us at home and out of the lab.

## Running

To run this code simply clone this repository and run the photon_counting_experiment.py script with Python (the numpy and matplotlib modules are required):
 
```
$ git clone https://github.com/bcolmey/Single-Photon-Counting-Statistics
$ cd Single-Photon-Counting-Statistics
$ python run_me.py 
```

## Examples:

Here are some examples of plots produced by the program. 

### QD single photon emission 
![Alt text](https://github.com/bcolmey/Single-Photon-Counting-Statistics/blob/main/plots/emission_distribution.jpeg?raw=true "Title")

### QD multiphoton emission
![Alt text](https://github.com/bcolmey/Single-Photon-Counting-Statistics/blob/main/plots/multiphoton-emission.jpeg?raw=true "Title")

### Histogram of coincidence events for a pulsed QD source emitting single photons only
![Alt text](https://github.com/bcolmey/Single-Photon-Counting-Statistics/blob/main/plots/single-photon_histogram.jpeg?raw=true "Title")

### Gaussian jitter effect on detection times
![Alt text](https://github.com/bcolmey/Aharonov-Bohm-Space-Charge-Effects-in-Python/blob/main/Animations/Gaussian_jitter.jpeg?raw=true "Title")

### Histogram of coincidence events for a pulsed source of bi-photons pairs
![Alt text](https://github.com/bcolmey/Single-Photon-Counting-Statistics/blob/main/plots/bi-photon_histogram.jpeg?raw=true "Title")

