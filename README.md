# EITSourceSeparation
Repository of algorithms presented in [upcoming paper] for source-separation of respiratory and cardiac signals from raw thoracic electrical impedance tomography (EIT) data. 

## Source Separation Algorithms
In the paper, the following existing algorithms were implemented as reference:
1. Bandpass filtering (./algorithms/bandFilt.m)
2. Complete ensemble empirical mode decomposition (EMD) with adaptive noise (CEEMDAN) (./algorithms/ceemdanFilt.m)
3. Dynamic principal component analysis (PCA) (./algorithms/dPCA.m)

while the following novel approaches were suggested for more accurate, robust, and real-time separation of sources:
1. Harmonic regression (./algorithms/harmRegWin.m)
2. Optimal harmonic filtering (./algorithms/optFilt.m)
    1. Modelled output
    2. Filterd output
	
Here is a visual comparison of the algorithm performance shown in the paper.
<img src="https://github.com/user-attachments/assets/fcec7dee-e4c7-4ca4-81a0-8774b1bc499d" width=80% height=80%>

## Data Synthesizer
Example data is provided in ./data/. Real and artificial data from the lung and heart regions of interest is provided. Artificial generation is performed by generating a time-series of heart- and respiratory-rate signals which is then used to extrapolate pre-define waveforms for the cardiac and respiratory signals into full-length recordings. The synthesizer is found in ./generator/ where the set of differential equations describing the states of three coupled oscillators (heart, lungs, and Mayer waves) is defined and solved.

Here is an example of simulated waveforms for cardiac (red) and respiratory (blue) signals, as well as their total signal (black).
<img src="https://github.com/user-attachments/assets/728c7879-b987-4ba8-b8eb-41000c788775" width=50% height=50%>

## Showcasing Scripts
To test and compare the source-separation algorithms, check ./showcaseSeparation.m.

To test artifical generation of data, check ./generateDate.m.
