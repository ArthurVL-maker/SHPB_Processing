# SHPB_Processing.py
# ----------------------------------------------------------------------------------------------------------------------------
# Signal processing function to correct and evaluate raw voltage signals optimised for high-strain rate SHPB testing.

# CODE DEVELOPED BY: Kin Shing Oswald Li, Arthur Van Lerberghe & Andrew D. Barr
# EMAIL: ksoli1@sheffield.ac.uk, avanlerberghe1@sheffield.ac.uk & a.barr@sheffield.ac.uk

# INPUTS: 
# csv_path = File path to CSV file containing oscilloscope columns in the form of
        # Time, CH1, CH2, CH3... (supports combined oscilloscopes).
# sample_data = List containing 3 elements that correspond to length (mm), mass (g), and dry mass (g) for the initial sample;
        # i.e. [initial_length, mass, dry_mass].
# confinement = Specify confinement of specimen.
        # 'None', 'Ring', or 'Reservoir' only.
        # confinement defaults to 'None' if not specified.
# signal_channels = Specify oscilloscope channels for input, output bars and any confinement mechanism.
        # i.e. [in_bar_gauge_channel, out_bar_gauge_channel, ring_gauge_channel OR reservoir_gauge_channel].
        # Input bar channel 7, output bar channel 8, and ring gauge channel 5 would correspond to list of [7, 8, 5].
        # signal_channels defaults to [7, 8, 5] if not specified.
# signal_amp = Specify oscilloscope signal amplification factors for input, output bars and any confinement mechanism.
        # i.e. [in_bar_gauge_amp, out_bar_gauge_amp, ring_gauge_amp].
        # Input bar amp 10, output bar amp 100, and ring gauge amp 10 would correspond to list of [10, 100, 10].
        # signal_amp defaults to [10, 100, 10] if not specified.
# disp_correction = Specify if dispersion correction or simple timeshift to process pulse.
        # True for dispersion correction, False for simple timeshift.
        # If disp_correction = True, ensure dispersion.py and dispersion_factors.py are attached.
        # disp_correction defaults to True if not specified.
# alignment = Specify the alignment mode for aligning stress waves at front and back of sample interface.
        # 'start' aligns the start of incident and transmitted pulse.
        # 'end' aligns the end of incident and transmitted pulse.
        # 'mid' aligns the median time of the pulse of the incident and transmitted pulse.
        # Integer/float values greater than 1 aligns the peaks of the incident and transmitted pulse to specific time (ms).
        # Float values greater than 0 and less than 1 aligns the incident and transmitted pulse based on a specific fraction of the maximum value of each pulse.
        # alignment defaults to 'start' if not specified.
# speedtrap = Specify if speed trap is employed and whether velocity of striker bar is determined.
        # True or False.
        # speedtrap is True if not specified (striker velocity is calculated).

# Other constants used as inputs for signal processing are defined in the 'INPUTS' section within the function, amend values if required.

# OUTPUTS:
# Generates processed data based on test data.
# if savedata=True (within 'INPUT' section):
        # Option to save processed data as CSV files within folder.
        # Log to show history for processing specific data set.

# ----------------------------------------------------------------------------------------------------------------------------
# Imported modules:
import pandas as pd
import numpy as np
import statistics
import warnings
import shutil
import time
import math
import os


def SHPB_Processing(csv_path, sample_data, confinement='None', signal_channels=[7, 8, 5], signal_amp=[10, 100, 10], disp_correction=True, alignment='start', speedtrap=True):
    time_start = time.time()  # Start timer
    warnings.filterwarnings("ignore")  # Ignore error warnings
    # ------------------------------------
    # INPUTS
    savedata = False  # Save data to CSV files after processing. Do not save to reduce computation time if not required.
    
    # Sample inputs:
    sample_diameter = 25  # mm
    
    # Speed trap inputs:
    if speedtrap is not False:
        speedtrap_distance = 50  # Distance between speed traps, mm
        speedtrap_trigger_voltage = 24  # Voltage when speed traps are activated, V
        speedtrap_front_channel = 1  # Speed trap front oscilloscope channel (i.e. 1 if CH1 is speed trap front).
        speedtrap_back_channel = 2  # Speed trap back oscilloscope channel (i.e. 2 if CH2 is speed trap back).
    
    # Incident bar inputs - Steel SS-25:
    in_bar_density = 7666  # Bar density, kg/m**3.
    in_bar_diameter = 25  # Bar diameter, mm.
    in_bar_wave_speed = 5376  # Bar wave speed, m/s.
    in_bar_gauge_factor = 123  # Input bar gauge factor.
    in_bar_gauge_voltage = 4  # Input bar signal voltage, V.
    in_bar_gauge_offset = 1000  # Distance from strain gauge to sample face, mm.

    # Transmitter bar inputs - Steel SS-25:
    out_bar_density = 7677  # Bar density, kg/m**3.
    out_bar_diameter = 25  # Bar diameter, mm.
    out_bar_wave_speed = 5305  # Bar wave speed, m/s.
    out_bar_gauge_factor = 127  # Output bar gauge factor. 
    out_bar_gauge_voltage = 4  # Output bar signal voltage, V.
    out_bar_gauge_offset = 500  # Distance from strain gauge to sample face, mm.

    # Confining ring inputs:
    confinement_type = str(confinement)  # Specimen confinement 'None'/'Ring'/'Reservoir'
    if confinement_type == 'Ring':
        ring_outside_diameter = 35  # Outside diameter, mm.
        ring_inside_diameter = 25  # Inside diameter, mm.
        ring_length = 5  # Length, mm.
        ring_gauge_factor = 124  # Gauge factor.
        ring_gauge_voltage = 8  # Signal voltage, V.
        ring_youngs_modulus = 206  # 200 - Young's modulus, GPa.
    
    # Steel reservoir and pressure transducer inputs:
    elif confinement_type == 'Reservoir':
        reservoir_fluid_wave_speed = 1482  # Wave speed of water 1482 m/s.
        reservoir_thickness = 7.5  # Thickness of fluid annulus at transducer, mm.
        reservoir_gauge_factor = 2.90  # Reservoir transducer calibration, mV/MPa.
        reservoir_gauge_voltage = 10  # Reservoir transducer voltage, V.
    else:
        confinement_type = 'None'

    # ------------------------------------
    # FILE AND LOG
    csv_name = os.path.basename(csv_path).split('.')[0]  # Name of raw data CSV file.
    
    with open(os.path.join('.', f'{csv_name}_log.txt'), 'w', encoding='utf-8') as file: 
        file.write("Please refer to SHPB_Processing.py for more info.\n"
                   'Code written by: Kin Shing Oswald Li, Arthur Van Lerberghe, Andrew D. Barr')

    def print_save(text):  # Function to print text and also write console log.
        print(text)
        with open(os.path.join('.', f'{csv_name}_log.txt'), 'a', encoding='utf-8') as file:
            print(text, file=file)

    # ------------------------------------
    # READING RAW DATA
    # Sample data:
    sample_initial_length = sample_data[0]  # mm
    sample_mass = sample_data[1]  # g
    sample_dry_mass = sample_data[2]  # g
    sample_initial_volume = sample_initial_length * math.pi * ((in_bar_diameter/2)**2) * 10**(-3)  # Sample initial volume, cm^3.
    
    if any(sd < 0 for sd in (sample_initial_length, sample_mass, sample_dry_mass)):
        raise ValueError('Invalid sample_data - negative value found.')
    if sample_dry_mass > sample_mass:
        raise ValueError('Invalid sample_data - dry_mass greater than mass.')
    
    in_bar_gauge_channel = signal_channels[0]  # Input bar oscilloscope channel.
    out_bar_gauge_channel = signal_channels[1]  # Output bar oscilloscope channel (i.e. 8 if CH8 is the input bar).
    if confinement_type == 'Ring':
        ring_gauge_channel = signal_channels[2]  # Confining ring oscilloscope channel (i.e. 5 if CH5 is the confining ring).
    elif confinement_type == 'Reservoir':
        reservoir_gauge_channel = signal_channels[2]  # Reservoir pressure transducer oscilloscope channel (i.e. 6 if CH6 is the reservoir pressure transducer).

    if any(ch < 0 for ch in signal_channels):
        raise ValueError('Invalid oscilloscope channel - negative value found.')
        
    in_bar_gauge_amp = signal_amp[0]  # Input bar signal amplification.
    out_bar_gauge_amp = signal_amp[1]  # Output bar signal amplification.
    if confinement_type == 'Ring':
        ring_gauge_amp = signal_amp[2]  # Signal amplification.
    
    if any(amp < 0 for amp in signal_amp):
        raise ValueError('Invalid oscilloscope amplification - negative value found.')
    
    # CSV file format (i.e. Relative time, Channel 1, Channel 2, Channel 3, Channel 4 ...):
    raw_data = pd.read_csv(csv_path, sep=';', skiprows=9, header=None, nrows=50000)  # Read csv file.
    time_base = raw_data.iloc[1:3, 0]  # First two time values, s.
    in_bar_gauge_signal = raw_data[in_bar_gauge_channel].iloc[1:50000]  # V.
    out_bar_gauge_signal = raw_data[out_bar_gauge_channel].iloc[1:50000]  # V.
   
    # ------------------------------------
    # Write initial log file:
    print_save('-' * 72 + '\n'
                + f'{time.strftime("%H:%M:%S %d-%m-%Y", time.localtime(time_start))}\n'
                + '-' * 72 + '\n'
                + f'PROCESSING DATA FROM:\n {csv_name}\n'
                + '-' * 72 + '\n'
                'FILE PATH:\n'
                + f' {csv_path}\n'
                'SAMPLE DATA USED:\n'
                + f' Length: {sample_data[0]}mm\n'
                + f' Wet mass: {sample_data[1]}g\n'
                + f' Dry mass: {sample_data[2]}g\n'
                'CONFINEMENT TYPE:\n'
                + f" {confinement_type if confinement_type=='Ring' or confinement_type=='Reservoir' else 'None'}\n"
                + '-' * 72 + '\n'
                + 'GAUGE INPUT SETTINGS:\n'
                + f' Incident channel: {in_bar_gauge_channel}\n'
                + f' Incident amplification: {in_bar_gauge_amp}\n'
                + f' Transmitted channel: {out_bar_gauge_channel}\n'
                + f' Transmitted amplification: {out_bar_gauge_amp}\n'
                + (f" Ring channel: {ring_gauge_channel}\n Ring amplification: {ring_gauge_amp}\n" if confinement_type == 'Ring' else '')
                + (f" Reservoir channel: {reservoir_gauge_channel}\n" if confinement_type == 'reservoir' else '')
                + '-' * 72)
    
    # ------------------------------------
    # FINDING STRIKER VELOCITY VIA SPEED TRAP
    # Finding index that activates speed trap
    time_step = time_base[2] - time_base[1]  # Oscilloscope time step, s.
    if speedtrap is not False:
        speedtrap_front_trigger = np.argmax(raw_data[speedtrap_front_channel].iloc[1:10000] > speedtrap_trigger_voltage)  # Index that activates front speed trap
        speedtrap_back_trigger = np.argmax(raw_data[speedtrap_back_channel].iloc[1:10000] > speedtrap_trigger_voltage)  # Index that activates back speed trap
    
    # Determining velocity based on time difference
        speedtrap_difference = abs(speedtrap_front_trigger - speedtrap_back_trigger) * time_step  # Difference in time between the two triggers, s.
        speedtrap_velocity = (speedtrap_distance/speedtrap_difference) / (10**3)  # Velocity of striker bar, m/s.
        
        print_save('STRIKER BAR VELOCITY:\n'
                    + ' ' + f'{round(speedtrap_velocity, 1)}m/s\n'
                    + '-' * 72)
    
    # ------------------------------------
    # PROCESSING BAR AXIAL SIGNALS
    # Bar strain gauge signals:
    in_bar_gauge_zero = statistics.mean(in_bar_gauge_signal.iloc[1000: 3000])  # Mean input bar "no signal" voltage, V.
    out_bar_gauge_zero = statistics.mean(out_bar_gauge_signal.iloc[1000: 3000])  # Mean output bar "no signal" voltage, V.

    # Strains:
    in_bar_strain = ((in_bar_gauge_signal - in_bar_gauge_zero) * 2) / (in_bar_gauge_factor * in_bar_gauge_voltage * in_bar_gauge_amp)  # Input bar strain assuming half wheatstone bridge.
    in_bar_youngs_modulus = (in_bar_wave_speed ** 2) * (in_bar_density / (10 ** 9))  # Input bar Young's modulus.
    in_bar_stress = in_bar_strain * in_bar_youngs_modulus * 1000  # Input bar stress
    
    out_bar_strain = ((out_bar_gauge_signal - out_bar_gauge_zero) * 2) / (out_bar_gauge_factor * out_bar_gauge_voltage * out_bar_gauge_amp)  # Output bar strain assuming half wheatstone bridge.
    out_bar_youngs_modulus = (out_bar_wave_speed ** 2) * (out_bar_density / (10 ** 9))  # Output bar Young's modulus.
    out_bar_stress = out_bar_strain * out_bar_youngs_modulus * 1000  # Output bar stress
    
    # ------------------------------------
    # DETECTING PULSES
    # Pulse triggers:
    incident_trigger_strain = 10*math.pow(10, math.ceil(math.log10(abs(max(in_bar_strain.iloc[1000: 3000])))))  # Find maximum absolute strain rounded up to nearest magnitude of 10 to use to indicate start of incident pulse.
    transmitted_trigger_strain = 10*math.pow(10, math.ceil(math.log10(abs(max(out_bar_strain.iloc[1000: 3000])))))  # Find maximum absolute strain rounded up to nearest magnitude of 10 to use to indicate start of transmitted pulse.
    if max(in_bar_strain.iloc[:20000]) < 2*incident_trigger_strain or max(out_bar_strain.iloc[:20000]) < 2*transmitted_trigger_strain:  # Raise error if no pulse is found that is at least 2x larger than the trigger strain.
        raise IndexError('Unable to detect pulse - check bar inputs or raw signal.')
    
    print_save('ABSOLUTE STRAIN TO TRIGGER PULSES:\n'
                f' Incident: {incident_trigger_strain}\n'
                f' Transmitted: {transmitted_trigger_strain}\n'
                + '-' * 72)
    
    # Finding incident pulse:
    incident_trigger = np.where(abs(in_bar_strain) > incident_trigger_strain)[0][0]  # Find when incident wave first has a value larger than incident_trigger_strain.
    if in_bar_strain[incident_trigger] < 0:
        in_bar_strain = -in_bar_strain  # If incident wave is negative, invert signal.
    incident_start = np.where(np.array(in_bar_strain.iloc[0:incident_trigger]) * np.array(in_bar_strain.iloc[1:incident_trigger+1]) < 0)[0][-1]  # Find last change of sign before trigger (start of incident pulse).
    incident_end = np.where((np.array(in_bar_strain.iloc[incident_start:-1]) * np.array(in_bar_strain.iloc[incident_start + 1:])) < 0)[0][1] + incident_start  # Find the next change of sign after trigger (end of incident pulse).
    incident_end = np.where((np.array(in_bar_strain.iloc[incident_start:-1]) * np.array(in_bar_strain.iloc[incident_start + 1:])) < 0)[0][1] + incident_start  # Find the next change of sign after trigger (end of incident pulse).
    incident_length = incident_end - incident_start  # Length of the incident pulse.

    # Finding reflected pulse:
    reflected_start = incident_start + round(2*in_bar_gauge_offset/(1000*in_bar_wave_speed*time_step))  # Find start of reflected wave based on wave speed.
    reflected_end = reflected_start + incident_length  # Length of the reflected pulse.
    
    # Finding transmitted pulse:
    transmitted_trigger = np.where(abs(out_bar_strain.iloc[incident_end:]) > transmitted_trigger_strain)[0][0] + incident_end  # Find when transmitted wave first has a value larger than transmitted_trigger_strain.
    if out_bar_strain[transmitted_trigger] < 0:
        out_bar_strain = -out_bar_strain  # If transmitted wave is negative invert signal.
    transmitted_start = np.where(np.array(out_bar_strain.iloc[0:transmitted_trigger]) * np.array(out_bar_strain.iloc[1:transmitted_trigger+1]) < 0)[0][-1]  # Find last change of sign before trigger (start of incident pulse).
    transmitted_end = np.where((np.array(out_bar_strain.iloc[transmitted_start:-1]) * np.array(out_bar_strain.iloc[transmitted_start + 1:])) < 0)[0][1] + transmitted_start  # Find the next change of sign after trigger (end of transmitted pulse).

    # ------------------------------------
    # TIME SHIFTING OR DISPERSION CORRECTION:
    # Creating signal cut off-length stress waves
    signal_cut_off = max(reflected_end, transmitted_end) + incident_length
    
    in_bar_incident = np.concatenate((np.zeros(incident_start), np.conj(np.array(in_bar_strain.iloc[incident_start:reflected_start+1])), np.zeros(signal_cut_off-reflected_start-1)))  # Concatenation of zero arrays with incident pulse.
    in_bar_reflected = np.concatenate((np.zeros(reflected_start), np.conj(np.array(in_bar_strain.iloc[reflected_start:reflected_end+1])), np.zeros(signal_cut_off-reflected_end-1)))  # Concatenation of zero array with reflected pulse.
    out_bar_transmitted = np.concatenate((np.zeros(transmitted_start), np.conj(np.array(out_bar_strain.iloc[transmitted_start:transmitted_end+1])), np.zeros(signal_cut_off-transmitted_end-1)))  # Concatenation of zero array with transmitted pulse.

    # Dispersion correction - see documentation for dispersion.py:
    if disp_correction is not False:
        print_save('PROCESSING WITH DISPERSION CORRECTION.' + '\n' + '-' * 72)
        fs = 1 / time_step  # Sampling frequency, Hz
        from dispersion import dispersion
        [in_bar_incident_strain, in_bar_incident_stress] = dispersion(in_bar_incident, fs, in_bar_diameter/2000, in_bar_wave_speed, in_bar_youngs_modulus, in_bar_gauge_offset/1000)
        [in_bar_reflected_strain, in_bar_reflected_stress] = dispersion(in_bar_reflected, fs, in_bar_diameter/2000, in_bar_wave_speed, in_bar_youngs_modulus, -in_bar_gauge_offset/1000)
        [out_bar_transmitted_strain, out_bar_transmitted_stress] = dispersion(out_bar_transmitted, fs, out_bar_diameter/2000, out_bar_wave_speed, out_bar_youngs_modulus, -out_bar_gauge_offset/1000)
    
    # Simple timeshift analysis without dispersion correction:
    else:
        print_save('PROCESSING WITH SIMPLE TIMESHIFT, NOT DISPERSION CORRECTION.' + '\n' + '-' * 72)
        in_bar_shift = round(((in_bar_gauge_offset/1000) / in_bar_wave_speed) / time_step)
        out_bar_shift = round(((out_bar_gauge_offset/1000) / out_bar_wave_speed) / time_step)

        in_bar_incident_strain = np.concatenate((np.array(in_bar_incident[-1-in_bar_shift:]), np.array(in_bar_incident[:-in_bar_shift])))
        in_bar_reflected_strain = np.concatenate((np.array(in_bar_reflected[in_bar_shift:]), np.array(in_bar_reflected[:in_bar_shift])))
        out_bar_transmitted_strain = np.concatenate((np.array(out_bar_transmitted[out_bar_shift:]), np.array(out_bar_transmitted[:out_bar_shift])))

        in_bar_incident_stress = in_bar_incident_strain * in_bar_youngs_modulus * 1000
        in_bar_reflected_stress = in_bar_reflected_strain * in_bar_youngs_modulus * 1000
        out_bar_transmitted_stress = out_bar_transmitted_strain * out_bar_youngs_modulus * 1000 

    # ------------------------------------
    # PROCESSING AXIAL STRESS AND STRAIN AT SPECIMEN INTERFACE
    # Find new start of incident pulse as start of sample stress/strain:
    incident_trigger_new = np.where(abs(in_bar_incident_strain) > incident_trigger_strain)[0][0]  # Find the new position of incident pulse.
    incident_start_new = np.where(in_bar_incident_strain[0:incident_trigger_new-1] * in_bar_incident_strain[1:incident_trigger_new] < 0)[0][-1]  # Use new start of incident pulse as start of sample interface analysis.

    # Bar displacement and sample strain:
    in_bar_displacement = np.zeros(signal_cut_off)  # Zero array placeholder for input bar displacement.
    out_bar_displacement = np.zeros(signal_cut_off)  # Zero array placeholder for output bar displacement.
    sample_strain = np.zeros(signal_cut_off)  # Zero array placeholder for sample strain.
    in_bar_displacement_alt = np.zeros(signal_cut_off)  # Zero array placeholder for input bar displacement.
    sample_strain_alt = np.zeros(signal_cut_off)  # Zero array placeholder for sample strain.

    for i in range(incident_start_new, signal_cut_off):
        in_bar_displacement[i] = in_bar_displacement[i-1] + ((in_bar_incident_strain[i] - in_bar_reflected_strain[i]) * 1000 * time_step * in_bar_wave_speed)  # Cumulative input bar displacement.
        out_bar_displacement[i] = out_bar_displacement[i-1] + (out_bar_transmitted_strain[i] * 1000 * time_step * out_bar_wave_speed)  # Cumulative output bar displacement.
        sample_strain[i] = (in_bar_displacement[i] - out_bar_displacement[i]) / sample_initial_length  # Sample axial strain.
        in_bar_displacement_alt[i] = in_bar_displacement[i-1] + ((in_bar_incident_strain[i]) * 1000 * time_step * in_bar_wave_speed)  # Cumulative input bar displacement, mm.
        sample_strain_alt[i] = 2 * (in_bar_displacement_alt[i] - out_bar_displacement[i]) / sample_initial_length  # Sample axial strain.
        if sample_strain[i] > 0.05 and abs(sample_strain[i] - sample_strain[i-1]) < 0.0001:  # Cut off when sample is fully strained i.e. when sample strain begins to flatten off (set to when strain difference between timesteps is 0.0001).
            in_bar_displacement = np.trim_zeros(in_bar_displacement, 'b')
            out_bar_displacement = np.trim_zeros(out_bar_displacement, 'b')
            sample_strain = np.trim_zeros(sample_strain, 'b')
            in_bar_displacement_alt = np.trim_zeros(in_bar_displacement_alt, 'b')
            sample_strain_alt = np.trim_zeros(sample_strain_alt, 'b')
            break
    
    # Final strain as end of sample stress/strain:
    sample_end = len(sample_strain) - 1  # Ending index based on when sample has been fully strained.

    # Define new starting index depending on alignment input:
    if alignment == 'end':  # Alignment mechanism for transmitted pulse based on detection of the ending of transmitted pulse.
        transmitted_end_new = np.where((out_bar_transmitted_strain[np.argmax(out_bar_transmitted_strain):-1]*out_bar_transmitted_strain[np.argmax(out_bar_transmitted_strain)+1:]) < 0)[0][0] + np.argmax(out_bar_transmitted_strain)
        transmitted_start_new = transmitted_end_new - (sample_end - incident_start_new)
    elif alignment == 'mid':  # Alignment mechanism for transmitted pulse on centering transmitted and incident/reflected pulses.
        transmitted_trigger_new = np.where(abs(out_bar_transmitted_strain[incident_start_new:]) > transmitted_trigger_strain)[0][0] + incident_start_new
        transmitted_start_new = np.where(out_bar_transmitted_strain[incident_start_new:transmitted_trigger_new-1] * out_bar_transmitted_strain[incident_start_new+1:transmitted_trigger_new] < 0)[0][-1] + incident_start_new
        transmitted_end_new = np.where((out_bar_transmitted_strain[np.argmax(out_bar_transmitted_strain):-1]*out_bar_transmitted_strain[np.argmax(out_bar_transmitted_strain)+1:]) < 0)[0][0] + np.argmax(out_bar_transmitted_strain)
        transmitted_start_new = round(transmitted_start_new + ((incident_start_new + sample_end)/2) - ((transmitted_start_new + transmitted_end_new)/2))
        transmitted_end_new = transmitted_start_new + (sample_end - incident_start_new)
    else: 
        transmitted_trigger_new = np.where(abs(out_bar_transmitted_strain[incident_start_new:]) > transmitted_trigger_strain)[0][0] + incident_start_new
        transmitted_start_new = np.where(out_bar_transmitted_strain[:transmitted_trigger_new-1] * out_bar_transmitted_strain[1:transmitted_trigger_new] < 0)[0][-1]
        transmitted_end_new = transmitted_start_new + (sample_end - incident_start_new)
        if isinstance(alignment, (int, float)) and 1 <= alignment < (sample_end - incident_start_new):  # Alignment mechanism based for transmitted pulse based on maximums of transmitted and incident/reflected pulses set to a defined value.
            incident_start_new = incident_start_new + np.argmax(in_bar_incident_strain[incident_start_new:sample_end+1]) - alignment  # Set peak of incident pulse to align with stress_peak input.
            transmitted_start_new = transmitted_start_new + np.argmax(out_bar_transmitted_strain[transmitted_start_new:transmitted_end_new+1]) - alignment  # Set peak of transmitted pulse to align with stress_peak input.
            transmitted_end_new = transmitted_start_new + (sample_end - incident_start_new)
        elif isinstance(alignment, float) and 0 < alignment < 1:  # Alignment mechanism for transmitted pulse based on a proportion of the maximums of transmitted and incident/reflected pulses (alignment value should be a decimal between 0-1)
            transmitted_trigger_new = np.where(abs(out_bar_transmitted_strain[transmitted_start_new:transmitted_end_new]) > (alignment * max(abs(out_bar_transmitted_strain[transmitted_start_new:transmitted_end_new]))))[0][0]
            incident_trigger_new = np.where(abs(in_bar_incident_strain[incident_start_new:sample_end]) > (alignment * max(abs(in_bar_incident_strain[incident_start_new:sample_end]))))[0][0]
            transmitted_start_new = transmitted_start_new - (incident_trigger_new - transmitted_trigger_new)
            transmitted_end_new = transmitted_start_new + (sample_end - incident_start_new) 
        else:
            alignment = 'start'
    
    # Sample front stress:
    in_bar_incident_strain = in_bar_incident_strain[incident_start_new:sample_end+1]  # Redefining incident strain to be within sample boundaries.
    in_bar_reflected_strain = in_bar_reflected_strain[incident_start_new:sample_end+1]  # Redefining reflected strain to be within sample boundaries.
    in_bar_incident_stress = in_bar_incident_stress[incident_start_new:sample_end+1]  # Redefining incident stress to be within sample boundaries.
    in_bar_reflected_stress = in_bar_reflected_stress[incident_start_new:sample_end+1]  # Redefining reflected stress to be within sample boundaries.
    
    stress_factor = ((in_bar_diameter/2)**2) / ((sample_diameter/2)**2)  # Stress factor to adjust for difference in diameter between sample and pressure bars.
    sample_front_stress = stress_factor * (in_bar_incident_stress + in_bar_reflected_stress)  # Stress at incident bar specimen face, MPa.
    
    # Sample back stress:
    out_bar_transmitted_strain = out_bar_transmitted_strain[transmitted_start_new+10:transmitted_end_new+1+10]  # Redefining transmitted strain to be within sample boundaries.
    out_bar_transmitted_stress = out_bar_transmitted_stress[transmitted_start_new+10:transmitted_end_new+1+10]  # Redefining transmitted stress to be within sample boundaries.
    
    sample_back_stress = stress_factor * out_bar_transmitted_stress  # Stress at transmitter bar specimen face, MPa.
    
    # Sample stress and strain:
    sample_strain = sample_strain[incident_start_new:sample_end+1]  # New sample strain set within bounds of sample stress/strain analysis.
    sample_length = (1 - sample_strain) * sample_initial_length  # Sample length.

    sample_mid_stress = (sample_front_stress + sample_back_stress)/2  # Mean axial specimen stress, MPa.
    
    # Sample axial strain rate:
    rel_time = time_step * np.arange(0, sample_end-incident_start_new+1)  # Relative time, s.
    sample_strain_rate = np.zeros((2, len(sample_strain)))
    for i in range(0, len(sample_strain)-1):
        sample_strain_rate[0, i] = (rel_time[i] + rel_time[i+1])/2
        sample_strain_rate[1, i] = ((sample_length[i] - sample_length[i+1]) / sample_length[i]) / time_step
    
    sample_strain_rate_1 = sample_strain_rate[0]
    sample_strain_rate_2 = sample_strain_rate[1]  # Sample strain rate
    
    in_bar_displacement = in_bar_displacement[incident_start_new:sample_end+1]  # Redefining incident bar displacement to be within sample boundaries.
    in_bar_displacement_alt = in_bar_displacement_alt[incident_start_new:sample_end+1]  # Redefining incident bar displacement to be within sample boundaries.
    out_bar_displacement = out_bar_displacement[incident_start_new:sample_end+1]  # Redefining transmitted bar displacement to be within sample boundaries.
    
    print_save("SAMPLE AXIAL RESULTS:\n"
                # f" Maximum sample stress: {round(max(sample_mid_stress),2)}MPa\n"
                f" Maximum sample strain: {round(max(sample_strain), 3)*100}%\n"
                f" Maximum strain rate: {round(max(sample_strain_rate_2),0)}s⁻¹\n"
                + '-' * 72)
    
    strain_difference = in_bar_incident_strain - (out_bar_transmitted_strain - in_bar_reflected_strain)  # Difference in strain to check condition for stress equilibrium (strainI = strainT - strainR).
    
    # ------------------------------------
    # PROCESSING RADIAL STRESSES
    # Confining ring strain input:
    if confinement_type == 'Ring':
        # Processing raw data from confining ring input:
        ring_gauge_signal = raw_data[ring_gauge_channel].iloc[1:50000]  # Confining ring signal, V.
        ring_gauge_zero = statistics.mean(ring_gauge_signal.iloc[:1000])  # Mean confining ring "no signal" voltage, V.
        ring_radial_strain = np.array((ring_gauge_signal - ring_gauge_zero) * 4 / (ring_gauge_factor * ring_gauge_voltage * ring_gauge_amp))  # Confining ring strain assuming quarter wheatstone bridge.
       
        # Finding radial pulse:
        ring_trigger_strain = math.pow(10, math.ceil(math.log10(abs(max(ring_radial_strain[1000: 3000])))))  # Find maximum absolute strain rounded up to nearest magnitude of 10 to use to indicate start of radial pulse.
        ring_pulse_trigger = np.where(abs(ring_radial_strain[incident_start_new:]) > ring_trigger_strain)[0][0] + incident_start_new  # Find when transmitted wave first has a value larger than transmitted_trigger_strain.
        if ring_radial_strain[np.argmax(abs(ring_radial_strain[incident_start_new:transmitted_end_new]))] < 0:
            ring_radial_strain = -ring_radial_strain
        if max(ring_radial_strain) < 1.1*ring_trigger_strain:  # Raise error if no pulse is found that is at least 10% larger than the trigger strain.
            raise IndexError('Unable to detect radial pulse - check ring inputs or raw signal.')
        radial_start = np.where(ring_radial_strain[incident_start_new:ring_pulse_trigger-1] * ring_radial_strain[incident_start_new+1:ring_pulse_trigger] < 0)[0][-1] + incident_start_new  # Find last change of sign before trigger (start of radial pulse).
        radial_end = radial_start + (sample_end - incident_start_new)  # Find end of radial pulse based on sample pulse length.
        
        # Sample radial stress and strain:
        ring_thick_walled_pipe_factor = (((ring_outside_diameter/2)**2) - ((ring_inside_diameter/2)**2)) / (2*(ring_inside_diameter/2)**2)  # Ratio of internal radial stress on the specimen to circumferential stress in the ring.
        sample_radial_strain = ring_radial_strain[radial_start:radial_end+1]  # Sample radial strain.
        sample_radial_stress = (ring_thick_walled_pipe_factor * (ring_youngs_modulus * 1000) * sample_radial_strain) * (ring_length / sample_length)  # Radial stress from the ring, MPa.
        sample_volume = sample_initial_volume * (1-sample_strain)  # Soil sample volume, cm^3.
        sample_density = sample_mass / sample_volume  # Sample density, Mg/m^3.
        sample_dry_density = sample_dry_mass / sample_volume  # Sample dry density, Mg/m^3.
        sample_mean_stress = (sample_mid_stress + 2 * sample_radial_stress) / 3
        
        print_save('SAMPLE RADIAL RESULTS VIA CONFINING RING STRAIN:\n'
                    f' Absolute radial strain trigger: {ring_trigger_strain}\n'
                    f' Maximum radial stress: {round(max(sample_radial_stress),2)}MPa\n'
                    f' Maximum radial strain: {round(max(sample_radial_strain)*100,5)}%\n'
                    f' Maximum mean stress: {round(max(sample_mean_stress), 2)}MPa\n'
                    f' Change in volume: {round(sample_initial_volume - sample_volume[-1], 3)}cm³\n'
                    f' Change in density: {round(sample_density[-1] - sample_density[0], 3)}Mg/m³\n'
                    f' Change in dry density: {round(sample_dry_density[-1] - sample_dry_density[0], 3)}Mg/m³\n'
                    + '-' * 72)
    
    # Reservoir with pressure transducer input:
    elif confinement_type == 'Reservoir':
        # Processing raw data from pressure transducer input:
        reservoir_gauge_signal = raw_data[reservoir_gauge_channel].iloc[1:50000]  # Pressure transducer signal, V.
        reservoir_gauge_zero = statistics.mean(reservoir_gauge_signal[:1000])  # Mean pressure transducer "no signal" voltage, V.
        reservoir_stress = ((reservoir_gauge_signal - reservoir_gauge_zero) * 1000) / reservoir_gauge_factor  # Pressure transducer stress, MPa.
        
        # Timeshifting radial stress based on reservoir travel time:
        reservoir_transit = (reservoir_thickness / 1000) / reservoir_fluid_wave_speed  # Time for pulse to travel through reservoir fluid, s.
        reservoir_time_steps = round(reservoir_transit / time_step)  # Timeshift in oscilloscope timesteps.
        reservoir_radial_stress = np.array(reservoir_stress[reservoir_time_steps:-1], reservoir_stress[:reservoir_time_steps - 1])
        reservoir_radial_stress = np.concatenate((reservoir_radial_stress, np.full(reservoir_time_steps+1, reservoir_radial_stress[-1])))
        
        # Finding radial pulse:
        # reservoir_trigger_stress = math.pow(10, math.ceil(math.log10(abs(max(reservoir_radial_stress[:4000]))))) # Find maximum absolute stress rounded up to nearest magnitude of 10 to use to indicate start of radial pulse. ***(which one?)
        reservoir_trigger_stress = abs(max(reservoir_radial_stress[incident_start_new-1000: incident_start_new]))  # Find maximum absolute stress to indicate start of radial pulse. ***(which one?)
        reservoir_pulse_trigger = np.where(abs(reservoir_radial_stress[incident_start_new:]) > reservoir_trigger_stress)[0][0] + incident_start_new  # Find when transmitted wave first has a value larger than transmitted_trigger_strain.
        if max(reservoir_radial_stress) < 0:
            reservoir_radial_stress = -reservoir_radial_stress
        if max(reservoir_radial_stress) < 2*reservoir_trigger_stress:  # Raise error if no pulse is found that is at least 2x larger than the trigger stress.
            raise IndexError('Unable to detect radial pulse - check reservoir inputs or raw signal.')
        
        radial_start = reservoir_pulse_trigger  # Start of radial pulse.
        radial_end = radial_start + (sample_end - incident_start_new)  # Find end of radial pulse based on sample pulse length.
        
        # Sample radial stress:
        sample_radial_stress = reservoir_radial_stress[radial_start:radial_end+1]  # Sample radial stress.
        
        print_save('SAMPLE RADIAL RESULTS VIA WATER RESERVOIR PRESSURE:\n'
                    f' Absolute radial stress trigger: {reservoir_trigger_stress}\n'
                    f' Maximum radial stress: {round(max(sample_radial_stress),2)}MPa\n'
                    + '-' * 72)

        radial_start = radial_start + (np.argmax(sample_radial_stress) - np.argmax(sample_mid_stress))  # Aligning maximum of radial pulse with axial mid stress.
        radial_end = radial_start + sample_end - incident_start_new  # New end of radial pulse.
        sample_radial_stress = reservoir_radial_stress[radial_start:radial_end+1]  # New sample radial stress aligned with pulse length of axial mid stress.

    else:
        print_save('NO RADIAL MEASUREMENT SELECTED:\n'
                    ' No radial stress and strain can be obtained.\n'
                    ' No volume and density data can be obtained.\n'
                    + '-' * 72)
    
    # ------------------------------------
    # EXPORTING AND SAVING PROCESSED DATA    
    if savedata is True:
        if speedtrap == False:
            folder_path = f"./Processed Data/{confinement_type if confinement_type=='Ring' or confinement_type=='Reservoir' else 'None'}/{csv_name}"  # File path to find processed data.
        else:
            folder_path = f"./Processed Data/{confinement_type if confinement_type=='Ring' or confinement_type=='Reservoir' else 'None'}/{round(speedtrap_velocity)}ms/{csv_name}"
        os.makedirs(folder_path, exist_ok=True)
        
        # Saving sample results:
        sample_results = {
            'Relative time (s)': rel_time,
            'Front stress (MPa)': sample_front_stress,
            'Back Stress (MPa)': sample_back_stress,
            'Mid Stress (MPa)': sample_mid_stress,
            'Strain': sample_strain,
            'Strain rate (s⁻¹)': sample_strain_rate_2,
            'Length (mm)': sample_length, }
        
        # Saving radial sample results if confinement is selected:
        if confinement_type == 'Ring':
            sample_results.update({
                'Radial stress (MPa)': sample_radial_stress,
                'Radial strain': sample_radial_strain,
                'Mean stress (MPa)': sample_mean_stress,
                'Volume (cm³)': sample_volume,
                'Density (Mg/m³)': sample_density,
                'Dry density (Mg/m³)': sample_dry_density})
        elif confinement_type == 'Reservoir':
            sample_results.update({
                'Radial stress (MPa)': sample_radial_stress})
        sample_results = pd.DataFrame(sample_results)
        sample_results.to_csv(os.path.join(folder_path, f'{csv_name}_sample_results.csv'), index=False)
    
        # Saving constants used as inputs:
        constants = {
            'Sample initial length': sample_initial_length,
            'Sample mass': sample_mass,
            'Sample dry mass': sample_dry_mass,
            'Sample diameter': sample_diameter,
            'Sample initial volume': sample_initial_volume,
            'In-bar density': in_bar_density,
            'In-bar diameter': in_bar_diameter,
            'In-bar wavespeed': in_bar_wave_speed,
            'In-bar gauge channel': in_bar_gauge_channel,
            'In-bar gauge factor': in_bar_gauge_factor,
            'In-bar gauge amp': in_bar_gauge_amp,
            'In-bar gauge voltage': in_bar_gauge_voltage,
            'In-bar gauge offset': in_bar_gauge_offset,
            'In-bar Youngs modulus': in_bar_youngs_modulus,
            'Out-bar density': out_bar_density,
            'Out-bar diameter': out_bar_diameter,
            'Out-bar wavespeed': out_bar_wave_speed,
            'Out-bar gauge channel': out_bar_gauge_channel,
            'Out-bar gauge factor': out_bar_gauge_factor,
            'Out-bar gauge amp': out_bar_gauge_amp,
            'Out-bar gauge voltage': out_bar_gauge_voltage,
            'Out-bar gauge offset': out_bar_gauge_offset,
            'Out-bar Youngs modulus': out_bar_youngs_modulus}
        
        # Saving radial constants if confinement is selected:
        if confinement_type == 'Ring':
            constants.update({
                'Ring outside diameter': ring_outside_diameter,
                'Ring inside diameter': ring_inside_diameter,
                'Ring length': ring_length,
                'Ring gauge channel': ring_gauge_channel,
                'Ring gauge factor': ring_gauge_factor,
                'Ring gauge amp': ring_gauge_amp,
                'Ring gauge voltage': ring_gauge_voltage,
                'Ring Youngs modulus': ring_youngs_modulus})
        elif confinement_type == 'Reservoir':
            constants.update({
                'Reservoir fluid wave speed': reservoir_fluid_wave_speed,
                'Reservoir thickness': reservoir_thickness,
                'Reservoir gauge channel': reservoir_gauge_channel,
                'Reservoir gauge factor': reservoir_gauge_factor,
                'Reservoir gauge voltage': reservoir_gauge_voltage})
        constants = pd.DataFrame(constants.items(), columns=['Constant', 'Value'])
        constants.to_csv(os.path.join(folder_path, f'{csv_name}_constants.csv'), index=False)
        
        # Saving incident and transmitter bar signal data:
        signal_results = {
            'Gauge signal (In Bar)': in_bar_gauge_signal,
            'Strain (In Bar)': in_bar_strain,
            'Stress (In Bar)': in_bar_stress,
            'Gauge signal (Out Bar)': out_bar_gauge_signal,
            'Strain (Out Bar)': out_bar_strain,
            'Stress (Out Bar)': out_bar_stress}
        
        # Saving radial signal data if confinement is selected:
        if confinement_type == 'Ring' or confinement_type == 'Reservoir':
            if confinement_type == 'Ring':
                signal_results.update({
                    'Gauge signal (Ring)': ring_gauge_signal,
                    'Strain (Ring)': ring_radial_strain})
            else:
                signal_results.update({
                    'Gauge signal (Reservoir)': reservoir_gauge_signal,
                    'Stress (Reservoir)': reservoir_stress,
                    'Adjusted stress (reservoir) at interface': reservoir_radial_stress})
        signal_results = pd.DataFrame(signal_results)
        signal_results.to_csv(os.path.join(folder_path, f'{csv_name}_signal_results.csv'), index=False)
        
        # Saving incident, reflected, and transmitted pulse data:
        pulse_results = {
            'Incident pulse strain': in_bar_incident_strain,
            'Incident pulse stress': in_bar_incident_stress,
            'Reflected pulse strain': in_bar_reflected_strain,
            'Reflected pulse stress': in_bar_reflected_stress,
            'Displacement (In bar)': in_bar_displacement,
            'Displacement alt (In bar)': in_bar_displacement_alt,
            'Transmitted pulse strain': out_bar_transmitted_strain,
            'Transmitted pulse stress': out_bar_transmitted_stress,
            'Displacement (Out Bar)': out_bar_displacement}
        pulse_results = pd.DataFrame(pulse_results)
        pulse_results.to_csv(os.path.join(folder_path, f'{csv_name}_pulse_results.csv'), index=False)
            
        print_save(f'PROCESSED DATA SAVED SUCCESSFULLY\n'
                   f" CSV file names: sample_results, constants, signal_results, pulse_results\n"
                   f' CSV files saved to folder: {os.path.abspath(folder_path)}\n'
                   f' Console log saved to log.txt\n'
                   + '-' * 72)
    else:
        print_save('PROCESSED DATA HAS NOT BEEN SAVED\n'
                   'savedata=False\n'
                   + '-' * 72)
    # ------------------------------------
    print_save(f'ALL DATA PROCESSED SUCCESSFULLY\n'
                f' To rerun, enter:\n'
                f" SHPB_Processing(r'{csv_path}', {sample_data}, '{confinement_type}', signal_channels=[{in_bar_gauge_channel}, {out_bar_gauge_channel}{', '+ str(ring_gauge_channel) if confinement_type=='Ring' else ''}{', '+ str(reservoir_gauge_channel) if confinement_type=='Reservoir' else ''}], signal_amp=[{in_bar_gauge_amp}, {out_bar_gauge_amp}{', '+ str(ring_gauge_channel) if confinement_type=='Ring' else ''}], disp_correction={False if disp_correction==False else True}, alignment={alignment if isinstance(alignment, (int, float)) else repr(alignment)}, speedtrap={False if speedtrap==False else True})\n"
                + '-' * 72 + '\n'
                + f'TIME REQUIRED: {round(time.time() - time_start, 3)}s.\n'
                + '-' * 72)
    
    if savedata == True:
        with open(os.path.join(folder_path, f'{csv_name}_log.txt'), 'a', encoding='utf-8') as file:
            file.write('-' * 8 + "Hope you're doing okay. Stay positive and keep vibing!!" + '-' * 9)
        
        shutil.move(f'./{csv_name}_log.txt', os.path.join(folder_path, f'{csv_name}_log.txt'))
    else:
        os.remove(f'{csv_name}_log.txt')
    
    # ------------------------------------
    # RETURN STATEMENTS
    # Change, add or remove return statements as required:
      
    if confinement_type == 'Ring' or confinement_type == 'Reservoir':
        return [sample_mid_stress, sample_front_stress, sample_back_stress, in_bar_stress, out_bar_stress, sample_radial_stress, sample_strain, sample_strain_rate_2]
        # return
    else:
        return [sample_mid_stress, sample_front_stress, sample_back_stress, in_bar_stress, out_bar_stress, sample_strain, sample_strain_rate_2]
        # return
