# dispersion_factors.py
# ----------------------------------------------------------------------------------------------------------
# Calculate corrections to amplitude and phase angle to account for dispersion at a particular frequency.
# The dispersion factors are calculated using 'phase_velocity.py' for a specific Poisson's ratio.
# 'phase_velocity.py' is available on GitHub and ORDA (Van Lerberghe, A., Barr, A. D. (2023)), see links below.
# It was inspired by a Matlab script created by Barr (2023), see link below.

# REQUIRES:
# - dispersion_factors: 4 .pickle files containing pre-calculated vectors for a particular Poisson's ratio (0.29).
#                       - norm_freqs: normalised frequencies (f*a/c0).
#                       - v_ratios: normalised velocities (c/c0).
#                       - m1: amplitude factor m1.
#                       - m2: normalised amplitude factor m2 (m2/E)).

# INPUTS:
# - f: Frequency, Hz.
# - a: Bar radius, m.
# - c0: One-dimensional wave velocity of the bar, m/s.
# - z: Distance to apply correction over, m.

# OUTPUTS:
# - angle_mod: Phase angle correction, rad.
# - m1: Correction for variation in response across bar cross-section.
# - m2: Correction for variation in ratio of axial stress and axial strain (dynamic Young's Modulus).

# MATLAB SOFTWARE:
# - Barr, A. D. (2023) phasevelocity.m - A Matlab script to calculate the frequency-dependent phase velocity and
# radial variation of elastic waves in cylindrical bars. University of Sheffield.
# Software ORDA link: [https://doi.org/10.15131/shef.data.21982604.v1]

# PYTHON SOFTWARE:
# - Van Lerberghe, A., Barr, A. D. (2023) phase_velocity.py - A Python algorithm for calculating
# frequency-dependent phase velocity and radial variation of elastic waves in cylindrical bars. University of Sheffield.
# Software ORDA link: [https://doi.org/10.15131/shef.data.22010999]
# Software GitHub link: [https://github.com/ArthurVL-maker/Phase_velocity.git]

# AUTHORS:
# Arthur Van Lerberghe (<avanlerberghe1@sheffield.ac.uk>) & Andrew D. Barr (<a.barr@sheffield.ac.uk>).
# ----------------------------------------------------------------------------------------------------------
# Imported modules:
from pathlib import Path
import pandas as pd
import numpy as np

# Load dispersionFactors files containing velocity data:
m1_data = pd.read_pickle(Path('dispersion_factors/m1.pickle'))
m2_data = pd.read_pickle(Path('dispersion_factors/m2.pickle'))
norm_freqs_data = pd.read_pickle(Path('dispersion_factors/norm_freqs.pickle'))
v_ratios_data = pd.read_pickle(Path('dispersion_factors/v_ratios.pickle'))


def dispersion_factors(f, a, c0, z):
    # Calculate normalised frequency:
    norm_freq = f * a / c0

    # Find change in phase angle:
    phase_velocity = np.interp(norm_freq, norm_freqs_data.iloc[:, 0], v_ratios_data.iloc[:, 0]) * c0  # Interpolate phase velocity value.
    angle_mod = (2 * np.pi * f * z) / phase_velocity  # Change in phase angle at norm_freq.

    # Find amplitude factors m1 & m2:
    m1 = np.interp(norm_freq,  norm_freqs_data.iloc[:, 0], m1_data.iloc[:, 0])  # Interpolated value of m1 at norm_freq.
    m2 = np.interp(norm_freq,  norm_freqs_data.iloc[:, 0], m2_data.iloc[:, 0])  # Interpolated value of m2/E at norm_freq.

    return [angle_mod, m1, m2]
