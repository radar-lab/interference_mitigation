## Project Title
Automotive Radar Interference Mitigation Using Adaptive Noise Canceller

Please cite this work as:

[F. Jin and S. Cao, "Automotive Radar Interference Mitigation Using Adaptive Noise Canceller," in IEEE Transactions on Vehicular Technology, vol. 68, no. 4, pp. 3747-3754, April 2019, doi: 10.1109/TVT.2019.2901493.](https://ieeexplore.ieee.org/abstract/document/8651538)

## Project Description
Interference among frequency modulated continues wave (FMCW) automotive radars can either increase the noise floor, which occurs in the most cases, or generate a ghost target in rare situations. To address the increment of noise floor due to interference, we proposed a low calculation cost method using adaptive noise canceller to increase the signal-to-interference ratio. In a quadrature receiver, the interference in the positive half of frequency spectrum is correlated to that in the negative half of frequency spectrum, whereas the beat frequencies from real targets are always present in the positive frequency. Thus, we estimated the power of the negative frequency as an indication of interference, and fed the positive frequency and negative frequency components into the primary and reference channel of an adaptive noise canceller, respectively. The least mean square algorithm was used to solve for the optimum filter solution. As a result, both the simulation and experiment showed a good interference mitigation performance.

## Project Source Code
To reproduce the figures presented in the paper, please refer to the main.m file in related folders.  
For the experiment data collection, we were using Texas Instruments AWR1243BOOST mmWave radar sensor and Texas Instruments TSW1400 as the data acquisition board.
