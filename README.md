## SHPB_Processing.py - An open-source Python algorithm for correcting stress wave dispersion in split-Hopkinson pressure bar experiments

#### DESCRIPTION:
When processing signals from split-Hopkinson pressure bar (SHPB) experiments, it is frequently presumed that longitudinal stress waves in the pressure bars travel one-dimensionally at a common velocity of c0. Hence, measurements recorded at the strain gauges are commonly simply translated to the end of the bar using a suitable time delay. In reality, stress waves travel at a certain phase velocity, cp, which varies with frequency, bar diameter, one-dimensional wave speed and Poisson’s ratio. As the frequency of a wave rises, the phase velocity drops, resulting in signal dispersion as it propagates down the bar. The dispersion of the stress pulse is followed by a frequency-dependent fluctuation in stress and strain throughout the bar cross-section. Therefore, a signal recorded on the surface of the bar at some distance from the specimen will not accurately represent the stresses the specimen was subjected to, and hence cannot be used to objectively measure the specimen response.


*SHPB_Processing.py* is the main algorithm required to process all the SHPB test data, from which the subroutines dispersion.py and dispersion_factors.py are called. The algorithm was inspired by a Matlab script created by Barr (2016). 

The functionalities, architecture and applied examples of the open-source algorithm *SHPB_Processing.py* are detailed in the paper by A.Van Lerberghe et al., currently going through the journal review process. 

The subroutine *dispersion.py*, in *process_SHPB.py*, uses an adaptation of Tyas and Pope’s (2005) dispersion-correction approach to verify that the inferred axial stress and strain data appropriately represent the specimen behaviour, as specified below:

1.	The fast Fourier transform (FFT) is used to convert the time-domain strain signal to the frequency domain.
2.	Using Bancroft’s equation, a correction is made to the phase angle of each frequency component to account for dispersion over the distance between the strain gauge and the bar end. 
3.	The amplitude of each frequency component is corrected using the factors m1, m2, v_ratio and norm_freq, which account for strain and Young’s modulus fluctuation over the bar cross section, respectively. Davies’ investigation of radial effects in a cylindrical pressure bar yielded these results.
4.	The inverse FFT is used to transform the signal back into the time domain.

To save time, dispersion.py utilises a subroutine *dispersion_factors.py*, which includes a precalculated, normalised look-up table of phase velocity, m1, m2, v_ratio and norm_freq. *SHPB_Processing.p*y, includes a lookup table for a Poisson’s ratio of 0.29. More tables can be generated using the calculation method outlined in Tyas and Pope (2005) and the Python script *phase_velocity.py*, (Van Lerberghe, A., Barr, A. D. (2023)), available on GitHub and ORDA, see links below.

#### FILES INCLUDED:
- *SHPB_Processing.py*: Includes the main python function, with the documentation on the use of the function included in the file as comments. It requires *dispersion.py*, *dispersion_factors.py*, and dispersion_factors folder to run.
- *dispersion.py*: A Python function, with the documentation on the use of the function included in the file as comments. It requires *dispersion_factors.py* to run.
- *dispersion_factors.py*: A Python function, with the documentation on the use of the function included in the file as comments. It requires the folder dispersion_factors.
- 'dispersion_factors' folder: A folder containing pre-calculated values of normalised frequency, normalised velocity and factors m1 and m2 for a material with a Poisson's ratio of 0.29 (i.e. stainless steel). This was calculated using the algorithm *phase_velocity.py* (Van Lerberghe, A., Barr, A. D. (2023)), available on GitHub and ORDA, see links below.

#### REFERENCES:
- Tyas, A., Pope, D.J., (2005). Full correction of first-mode Pochhammer–Chree dispersion effects in experimental pressure bar signals. Measurement science and technology, 16(3), p.642.

#### MATLAB SOFTWARE:
- Barr, A. D. (2016) *dispersion.m* - A Matlab script for phase angle and amplitude correction of pressure bar signals. University of Sheffield.\
Software ORDA link: [https://doi.org/10.15131/shef.data.3996876.v1]

#### PYTHON SOFTWARE:
- Van Lerberghe, A., Barr, A. D. (2023) *SHPB_Processing.py* - An open-source Python algorithm for correcting stress wave dispersion in split-Hopkinson pressure bar experiments. University of Sheffield.\
Software ORDA link: [https://doi.org/10.15131/shef.data.24570565]
- Van Lerberghe, A., Barr, A. D. (2023) *phase_velocity.py* - A Python algorithm for calculating frequency-dependent phase velocity and radial variation of elastic waves in cylindrical bars. University of Sheffield.\
Software ORDA link: [https://doi.org/10.15131/shef.data.22010999]\
Software GitHub link: [https://github.com/ArthurVL-maker/Phase_velocity.git]

#### AUTHORS:
Arthur Van Lerberghe <avanlerberghe1@sheffield.ac.uk>, Kin Shing Oswald Li <ksoli1@sheffield.ac.uk> & Andrew D. Barr <a.barr@sheffield.ac.uk>.
