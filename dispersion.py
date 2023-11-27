# dispersion.py
# ----------------------------------------------------------------------------------------------------------
# First-mode dispersion correction of a finite arbitrary signal in a cylindrical bar.

# OPERATION:
# - Finds FFT of the signal.
# - Corrects phase velocity and amplitude of each frequency using method described by Tyas & Pope (2005).
# - Reconstructs signal using IFFT.
# - Frequencies above fa/c0 = 0.2619 stripped (d/L = 0.6), due to limitations of m1 correction.

# INPUTS:
# - x: Zero-padded strain signal in time domain (1xN numeric).
# - fs: Sampling frequency, Hz.
# - a: Bar radius, m.
# - c0: One-dimensional wave velocity of the bar, m/s.
# - E: Young's modulus of the bar, GPa.
# - z: Distance to correction over, m (+ve in direction of propagation).

# OUTPUTS:
# - x_strain: Dispersion-corrected strain signal.
# - x_stress: Dispersion-corrected stress signal, MPa

# REFERENCES:
# - Tyas, A., Pope, D.J., (2005). Full correction of first-mode Pochhammer–Chree dispersion effects in experimental
# pressure bar signals. Measurement science and technology, 16(3), p.642.

# AUTHORS:
# Arthur Van Lerberghe (<avanlerberghe1@sheffield.ac.uk>) & Andrew D. Barr (<a.barr@sheffield.ac.uk>).
# ----------------------------------------------------------------------------------------------------------
# Imported modules:
import numpy as np

# Imported function:
from dispersion_factors import dispersion_factors


def dispersion(x, fs, a, c0, E, z):
    # Input signal:
    n = len(x)  # Number of elements in signal.
    f = np.arange(0, n-1) * (fs/n)  # FFT frequencies, Hz.
    fmax = 0.2619 * c0/a  # Max correctable frequency due to factor m1 limitations, Hz.

    # FFT the signal:
    X = np.fft.fft(x)
    XStrain = np.array(X)  # Create copy for strain correction.
    XStress = np.array(X)  # Create copy for stress correction.

    # Phase shift, adjust magnitude of frequency components:
    number_of_bins = len(X)

    if number_of_bins % 2 == 0:
        # n is even:
        positive_bins = np.arange(1, number_of_bins/2)  # Positive frequency bins.
        nyquist_bin = [number_of_bins/2]  # Nyquist frequency bin.
        bins_to_edit = np.concatenate((positive_bins, nyquist_bin))  # Total bins to edit individually.
        negative_bins = np.arange(number_of_bins/2+1, number_of_bins)  # Negative frequency bins.
    else:
        # n is odd:
        positive_bins = np.arange(1, (number_of_bins+1)/2)  # Positive frequency bins.
        bins_to_edit = np.array(positive_bins)  # Total bins to edit individually.
        negative_bins = np.arange((number_of_bins+1)/2+1, number_of_bins)  # Negative frequency bins.

    for b in bins_to_edit.astype(int):
        if f[b] <= fmax:
            # Find phase shift and factors m1 and m2 for current frequency:
            [angle_mod, m1, m2] = dispersion_factors(f[b], a, c0, z)
            # Apply shift and factors m1 to obtain corrected strain:
            XStrain[b] = m1 * np.abs(X[b]) * np.exp(1j * (np.angle(X[b])-angle_mod))
            # Apply phase shift and factors m1 & m2 to obtain corrected stress [/E]:
            XStress[b] = m1 * m2 * np.abs(X[b]) * np.exp(1j * (np.angle(X[b])-angle_mod))
        else:
            # Above fMax zero X data [Apply perfect low-pass filter].
            XStrain[b] = 0
            XStress[b] = 0

    # Correct negative bins by taking complex conjugate of positive bins:
    XStrain[negative_bins.astype(int)] = np.conj(XStrain[positive_bins.astype(int)[::-1]])
    XStress[negative_bins.astype(int)] = np.conj(XStress[positive_bins.astype(int)[::-1]])

    # Convert the corrected frequency components back into the time domain:
    x_strain = np.real(np.fft.ifft(XStrain))  # Corrected strain.
    x_stress = np.real(np.fft.ifft(XStress))*E*1000  # Corrected stress, MPa.

    return [x_strain, x_stress]
