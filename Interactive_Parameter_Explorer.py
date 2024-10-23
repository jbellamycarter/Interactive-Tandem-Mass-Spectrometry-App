## INTERACTIVE PARAMETER EXPLORER ##
"""
This script is used for the visualisation of Mass Spectrometry (MS) parameters. 
It imports the necessary libraries, defines functions for peak detection, centroid calculation, 
averages spectra, interpolates collision induced dissociation energy, retrieves fragments and includes 
functionalities to load and process mzML data. 
"""
import streamlit as st
import numpy as np
import requests
import io
from bokeh.models import ColumnDataSource, LabelSet, HoverTool, Range1d
from bokeh.plotting import figure
from pyteomics import mzml, mass, parser
from scipy import signal


## FUNCTIONS ##

def peak_detection(spectrum, threshold=5, prominence=0.8, width=2, distance=4, centroid=False):
    # Calculates relative intensity threshold for peak detection
    relative_threshold = spectrum['intensity array'].max() * (threshold / 100)
    if centroid:
        # For the centroid mode: Indices where the intensity exceeds the relative threshold are returned 
        peaks = np.where(spectrum['intensity array'] > relative_threshold)[0]
        return peaks
    else:
        # For non-centroid mode: scipy's find_peak for detection of peaks with specified propperties is used
        peaks, properties = signal.find_peaks(
            spectrum['intensity array'], 
            width=width,
            height=relative_threshold, 
            distance=distance,
            prominence=prominence
            )
        return peaks, properties

def get_centroid(spectrum, peaks, properties):
    # Initialises array to hold centroid values
    _centroids = np.zeros_like(peaks, dtype='float32')
    # Iterates over each peak detected 
    for i, peak in enumerate(peaks):
        # Retrieves right and left indices of peak's range from the dictionary: properties
        left_ip = int(properties['left_ips'][i])
        right_ip = int(properties['right_ips'][i])
        peak_range = range(left_ip, right_ip + 1)
        
        mz_value_range = spectrum['m/z array'][peak_range]
        intensity_range = spectrum['intensity array'][peak_range]

        # Calculates centroid as the weighted average of m/z values based on intensity 
        _centroids[i] = np.sum(mz_value_range * intensity_range) / np.sum(intensity_range)

    return _centroids


def average_spectra(spectra, filter_string=None, bin_width=None):
    # Initialises an array to hold reference m/z values from the first spectrum 
    reference_scan = np.unique(spectra[0]['m/z array'])
    # Bin width set for interpolation; if bin width not provided, default set to the minimum difference between m/z values
    if bin_width is None:
        bin_width = np.min(np.diff(reference_scan))

    # Determine m/z range from the first spectrum's scan window 
    scan_min = spectra[0]['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
    scan_max = spectra[0]['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']

    # Creates a new m/z array for the averaged spectrum with specified bin width
    reference_mz_value = np.arange(scan_min, scan_max, bin_width)
    # Initialises an array for the accumulation of the summed intensities across all spectra 
    merge_intensity = np.zeros_like(reference_mz_value)

    # Iterate over each spectrum to accumulate interpolated intensities
    for scan in spectra:
        temp_mz = scan['m/z array'] # Temporary mz 
        temp_intensity = scan['intensity array'] # Temporary intensity
        merge_intensity += np.interp(reference_mz_value, temp_mz, temp_intensity, left=0, right=0)

    # Gives mean average of accumulated intensities by dividing by the number of spectra 
    mean_intensity = merge_intensity / len(spectra)

    avg_spectrum = spectra[0].copy()
    avg_spectrum['m/z array'] = reference_mz_value
    avg_spectrum['intensity array'] = mean_intensity

    avg_spectrum['scanList']['scan'][0]['filter string'] = "AV: {:.2f}-{:.2f}; {}".format(
        spectra[0]['scanList']['scan'][0]['scan start time'], 
        spectra[-1]['scanList']['scan'][0]['scan start time'], filter_string)

    return avg_spectrum

def interpolate_spectra(spectra, target_energies, energies=[0, 5, 10, 15, 20]):
    # Extracts intensity arrays from each spectrum in the input list 
    intensity_arrays = [spectrum['intensity array'] for spectrum in spectra]
    # Initialises a dictionary to hold interpolated intensity values 
    interpolated_spectra = {target_energy: [] for target_energy in target_energies}

    # Iterate over each index in the intensity arrays and extract intensity values for the current index from each spectrum 
    for i in range(len(intensity_arrays[0])):
        _intensities = [intensity_arrays[index][i] for index in range(len(energies))]
        
        # Interpolate intensity values for each target array
        for target_energy in target_energies:
            #if target_energy < energies[0] or target_energy > energies[-1]:
            #    raise ValueError(f"Target energy {target_energy} not in the available range of energies in the spectra.")
            
            # Performs interpolation to estimate the intensity at target energy 
            interpolated_intensity = np.interp(target_energy, energies, _intensities)
            interpolated_spectra[target_energy].append(interpolated_intensity)
    
    return {key: np.array(value) for key, value in interpolated_spectra.items()}

# Standard amino acid masses with additional modifications for specific post-translational modifications 
aa_mass = mass.std_aa_mass
aa_mass['p'] = 79.966331  # phosphorylation (STY)
aa_mass['ox'] = 15.994915  # oxidation (MW)
aa_mass['d'] = 0.984016  # deamidation (NQ)
aa_mass['am'] = -0.984016  # amidation (C-term)

def get_fragments(sequence, selected_charge_state, peaks_data, ion_types=('b', 'y')):
    fragments = []
    _sequence = parser.parse(sequence) 

    pep_length = len(_sequence)
    print(f"Length of Peptide: {pep_length}")

    # Defines neutral losses which occur during fragmentation
    neutral_losses = {
        '': 0, 
        '-H2O': -18.01528, 
        '-NH3': -17.02655
    }

    def is_in_peaks_data(mass):
        tolerance = 0.2 
        return any(abs(mass - peak) <= tolerance for peak in peaks_data)

    # Iterate over each position in peptide sequence for fragment generation 
    for pos in range(1, pep_length):
        for ion_type in ion_types:
       
            # Determine sequence for fragmentation based on ion type (N-terminal or C-terminal)
            if ion_type[0] in ('a', 'b', 'c'):
                seq = ''.join(_sequence[:pos])
            elif ion_type[0] in ('x', 'y', 'z'):
                seq = ''.join(_sequence[-pos:])
            
            # Fragments generated for each charge state
            for charge in range(1, selected_charge_state + 1):
                for loss, mass_diff in neutral_losses.items():
                    # Calculate fragment mass with potential neutral loss
                    _mass = mass.fast_mass2(seq, ion_type=ion_type, charge=charge, aa_mass=aa_mass) + (mass_diff / charge)  

                    # Determine ion label based on ion type and neutral loss and charge state 
                    ion_label = ion_type + str(pos) + loss + "+"*charge  

                    # Check if calculated mass is close to any observed peaks
                    if is_in_peaks_data(_mass):
                            fragments.append({'sequence': seq, 'ion': ion_label, 'm/z': _mass, 'type': ion_type})
                            print(f"Annotated fragment: {ion_label}, m/z: {_mass}")

    # Handle precursor ion annotation
    for charge in range(1, selected_charge_state + 1):
        for loss, mass_diff in neutral_losses.items():
            # Calculate fragment mass with potential neutral loss
            seq = ''.join(_sequence)
            _mass = mass.fast_mass2(seq, ion_type="M", charge=charge, aa_mass=aa_mass) +  (mass_diff / charge)

            # Determine ion label based on ion type and neutral loss
            ion_label = "M" + loss + "+"*charge

            if is_in_peaks_data(_mass):
                    fragments.append({'sequence': seq, 'ion': ion_label, 'm/z': _mass, 'type': "M"})
                    print(f"Annotated fragment: {ion_label}, m/z: {_mass}")
        
    return fragments

def load_predefined_data(peptide, charge_state, resolution, energy_ramp, isolation=None):
    file_map = {
        ('MRFA', '1+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/05Mar2024_MJ_MRFA_1%2B_collision_energy_ramp_enhanced_01.mzML',
        ('MRFA', '1+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/05Mar2024_MJ_MRFA_1%2B_collision_energy_ramp_turbo_01.mzML',
        ('MRFA', '1+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/05Mar2024_MJ_MRFA_1%2B_collision_energy_ramp_zoom_01.mzML',
        ('MRFA', '2+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/05Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_zoom_01.mzML',
        ('MRFA', '2+', 'Zoom', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_zoom_iso_2.mzML',
        ('MRFA', '2+', 'Zoom', 'Iso 3', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_zoom_iso_3.mzML',
        ('MRFA', '2+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/05Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_turbo_01.mzML',
        ('MRFA', '2+', 'Turbo', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_turbo_iso_2.mzML',
        ('MRFA', '2+', 'Turbo', 'Iso 3', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_turbo_iso_3.mzML',
        ('MRFA', '2+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/05Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_normal_01.mzML',
        ('MRFA', '2+', 'Normal', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_normal_iso_2.mzML',
        ('MRFA', '2+', 'Normal', 'Iso 3', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_normal_iso_3.mzML',
        ('MRFA', '2+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/05Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_enhanced_01.mzML',
        ('MRFA', '2+', 'Enhanced', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_enhanced_iso_2.mzML',
        ('MRFA', '2+', 'Enhanced', 'Iso 3', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_2%2B_collision_energy_ramp_enhanced_iso_3.mzML',

        ('Bradykinin', '2+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_enhanced.mzML',
        ('Bradykinin', '2+', 'Enhanced', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_enhanced_defined.mzML',
        ('Bradykinin', '2+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_turbo.mzML',
        ('Bradykinin', '2+', 'Turbo', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_turbo_defined.mzML',
        ('Bradykinin', '2+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_zoom.mzML',
        ('Bradykinin', '2+', 'Zoom', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_zoom_defined.mzML',
        ('Bradykinin', '2+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_normal_01.mzML',
        ('Bradykinin', '2+', 'Normal', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_2%2B_collision_energy_ramp_normal_defined.mzML',
        ('Bradykinin', '3+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_enhanced.mzML',
        ('Bradykinin', '3+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_normal.mzML',
        ('Bradykinin', '3+', 'Normal', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_normal_defined.mzML',
        ('Bradykinin', '3+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_turbo.mzML',
        ('Bradykinin', '3+', 'Turbo', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_turbo_defined.mzML',
        ('Bradykinin', '3+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_zoom.mzML',
        ('Bradykinin', '3+', 'Zoom', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_zoom_defined.mzML',
        ('Bradykinin', '3+', 'Enhanced', 'Default', 'Defined'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/19Mar2024_MJ_bradykinin_3%2B_collision_energy_ramp_enhanced_defined.mzML', 

        ('Substance_P', '2+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_2%2B_collision_energy_ramp_enhanced.mzML',
        ('Substance_P', '2+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_2%2B_collision_energy_ramp_normal.mzML',
        ('Substance_P', '2+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_2%2B_collision_energy_ramp_turbo.mzML',
        ('Substance_P', '2+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_2%2B_collision_energy_ramp_zoom.mzML',
        ('Substance_P', '3+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_3%2B_collision_energy_ramp_enhanced.mzML',
        ('Substance_P', '3+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_3%2B_collision_energy_ramp_normal.mzML',
        ('Substance_P', '3+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_3%2B_collision_energy_ramp_turbo.mzML',
        ('Substance_P', '3+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Substance_P/20Mar2024_MJ_subp_3%2B_collision_energy_ramp_zoom.mzML',
        
        ('GRGDS', '1+', 'Enhanced', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_1%2B_collision_energy_ramp_enhanced_02.mzML',
        ('GRGDS', '1+', 'Normal', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_1%2B_collision_energy_ramp_normal_02.mzML',
        ('GRGDS', '1+', 'Turbo', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_1%2B_collision_energy_ramp_turbo_02.mzML',
        ('GRGDS', '1+', 'Zoom', 'Iso 2', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_1%2B_collision_energy_ramp_zoom_02.mzML',
        ('GRGDS', '2+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_2%2B_collision_energy_ramp_enhanced.mzML',
        ('GRGDS', '2+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_2%2B_collision_energy_ramp_normal.mzML',
        ('GRGDS', '2+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_2%2B_collision_energy_ramp_turbo.mzML',
        ('GRGDS', '2+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_2%2B_collision_energy_ramp_zoom.mzML',

        ('SDGRG', '1+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_1%2B_collision_energy_ramp_enhanced.mzML',
        ('SDGRG', '1+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_1%2B_collision_energy_ramp_normal.mzML',
        ('SDGRG', '1+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_1%2B_collision_energy_ramp_turbo.mzML',
        ('SDGRG', '1+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_1%2B_collision_energy_ramp_zoom.mzML',
        ('SDGRG', '2+', 'Enhanced', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_2%2B_collision_energy_ramp_enhanced.mzML',
        ('SDGRG', '2+', 'Normal', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_2%2B_collision_energy_ramp_normal.mzML',
        ('SDGRG', '2+', 'Turbo', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_conc_2%2B_collision_energy_ramp_turbo.mzML',
        ('SDGRG', '2+', 'Zoom', 'Default', 'Centre'): 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_conc_2%2B_collision_energy_ramp_zoom.mzML',

    }

    selected_file_url = file_map.get((peptide, charge_state, resolution, energy_ramp, "Centre"))

    if selected_file_url:
        file_response = requests.get(selected_file_url)
        # Convert response content to a BytesIO object
        raw_data = io.BytesIO(file_response.content)
        reader = mzml.read(raw_data, use_index=True)

        scan_energy_list = {}
    # Process each scan in the mzML file 
    reader.reset()
    for scan in reader:
        x = scan['index']
        if scan['ms level'] == 1: 
            continue  # Skip MS1 scans

        # Checks scan for precursor information and checks collision energy is available in precursor info
        if 'precursorList' in scan and 'precursor' in scan['precursorList'] and len(scan['precursorList']['precursor']) > 0:
            if 'activation' in scan['precursorList']['precursor'][0] and 'collision energy' in scan['precursorList']['precursor'][0]['activation']:
                _energy = scan['precursorList']['precursor'][0]['activation']['collision energy']
                if _energy not in scan_energy_list:
                    scan_energy_list[_energy] = []
                # Appends the scan index to the list for energy if not already present 
                scan_energy_list[_energy].append(x)
    
    return reader, scan_energy_list
    
# Allows user to use own data 
@st.cache_data
def load_data(raw_file):
    reader = mzml.read(raw_file, use_index=True)
    scan_energy_list = {}
    reader.reset()
    for scan in reader:
        x = scan['index']
        if scan['ms level'] == 1:
            continue  # Skip MS1 scans
        if 'precursorList' in scan and 'precursor' in scan['precursorList'] and len(scan['precursorList']['precursor']) > 0:
            if 'activation' in scan['precursorList']['precursor'][0] and 'collision energy' in scan['precursorList']['precursor'][0]['activation']:
                _energy = scan['precursorList']['precursor'][0]['activation']['collision energy']
                if _energy not in scan_energy_list:
                    scan_energy_list[_energy] = []
                scan_energy_list[_energy].append(x)
    
    return reader, scan_energy_list

## APP LAYOUT ##

# Set up the configuration for the Streamlit app 
st.set_page_config(page_title="MS/MS Parameter Explorer", 
                   layout="wide", 
                   menu_items={'about': "This application is a parameter explorer for mzML mass spectrometry data. Written by Kiah Slater and Dr Jedd Bellamy-Carter; data generated by Mitch Jones. Loughborough University."})

st.sidebar.title("MS/MS Parameter Explorer")
st.sidebar.markdown("This is an interactive parameter explorer for tandem mass spectrometry data stored in .mzML data format")

# Defines options for different peptides 
peptide_options = {
    "MRFA": {
        "sequence": "MRFA",
        "charge_states": ["1+", "2+"],
        "resolutions": {
            "1+": ["Turbo", "Enhanced", "Zoom"],
            "2+": ["Turbo", "Normal", "Enhanced", "Zoom"]
        },
        "isolation": {
            "1+": ["Default"],
            "2+": ["Default", "Iso 2", "Iso 3"]
        }
    },
    "GRGDS": {
        "sequence": "GRGDS",
        "charge_states": ["1+", "2+"],
        "resolutions": ["Turbo", "Normal", "Enhanced", "Zoom"],
        "isolation": {
            "1+": ["Default"],
            "2+": ["Default"]
        }
    },
    "SDGRG": {
        "sequence": "SDGRG",
        "charge_states": ["1+", "2+"],
        "resolutions": ["Turbo", "Normal", "Enhanced", "Zoom"],
        "isolation": ["Default"]
    },
    "Bradykinin": {
        "sequence": "RPPGFSPFR",
        "charge_states": ["2+", "3+"],
        "resolutions": ["Turbo", "Normal", "Enhanced", "Zoom"],
        "isolation": ["Default"]
    },
    "Substance_P": {
        "sequence": "RPKPQQFFGLamM",
        "charge_states": ["2+", "3+"],
        "resolutions": ["Turbo", "Normal", "Enhanced", "Zoom"],
        "isolation": ["Default"]
    }
}

        ## SIDEBAR LAYOUT ##

# Streamlit sidebar widgets
use_predefined_data = st.sidebar.checkbox(
    "Use Predefined Data", 
    value=True, 
    help="Toggle to use predefined data")

# Selectbox for choosing peptide
selected_peptide = st.sidebar.selectbox("Select Peptide", 
                                        # List of peptides obtained from the keys of the 'peptide_options' dictionary
                                        list(peptide_options.keys()),
                                        help = "Small selection of real peptides") 

def get_options(selected_peptide, option_type):
    return peptide_options[selected_peptide][option_type]

# Selectbox for choosing charge state 
selected_charge_state = st.sidebar.selectbox(
    "Select Charge State", 
    get_options(selected_peptide, "charge_states"),
    help = "Charge states of peptide with reasonable intensity after ESI.")

def get_resolutions(selected_peptide, selected_charge_state):
    if selected_peptide == "MRFA":
        # Exception for handling 'MRFA'
        if selected_charge_state in peptide_options[selected_peptide]["resolutions"]:
            return peptide_options[selected_peptide]["resolutions"][selected_charge_state]
        else:
            return []  # Return an empty list if no valid options are found
    else:
        return peptide_options[selected_peptide]["resolutions"]

# Selectbox for choosing resolution 
selected_resolution = st.sidebar.selectbox(
    "Select Resolution", 
    get_resolutions(selected_peptide, selected_charge_state),
    help="Available mass analyser resolutions. Turbo = 125000 Da/sec (FWHM ⩽3 Th), Normal = 16666 Da/sec (FWHM ⩽0.7 Th), Enhanced = 5000 Da/sec (FWHM ⩽0.45 Th), Zoom = 1111 Da/sec (FWHM ⩽0.3 Th)")

def get_isolation_options(selected_peptide, selected_charge_state):
    if selected_peptide in peptide_options and selected_charge_state in peptide_options[selected_peptide]["isolation"]:
        return peptide_options[selected_peptide]["isolation"][selected_charge_state]
    else:
        return ["Default"]  # Default option if no valid options are found

# Selectbox for choosing energy collision ramp
selected_isolation_options = get_isolation_options(selected_peptide, selected_charge_state)
if selected_isolation_options:
    selected_isolation = st.sidebar.selectbox("Select Isolation Parameters", selected_isolation_options)
else:
    selected_isolation = None

## TAB LAYOUT ##

# Creates tabs in the Streamlit app for displaying the spectrum and instructions 
spectrum_tab, instructions_tab = st.tabs(["Spectrum", "Instructions"])

# Instructions tab content
with instructions_tab:  
    # URLs for instructional images 
    image_predefined_data = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Data%20Selection.png'
    image_user_data = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Drag%20and%20drop%20data.png'
    image_parameter_selection = 'https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/raw/main/Instruction%20images/Parameter%20selection%20.png'
    image_setting_selection = 'https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/raw/main/Instruction%20images/Settings%20selection%20.png'
    image_plot_expansion = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/View%20fullscreen.png'
    image_zoom_function =  'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Zoom%20function.png'
    image_hover_function = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Hover%20function.png'

    st.header("Instructions")
    st.markdown("Instructions for the use of the Interactive Parameter Explorer")
    st.subheader("Data Selection")
    st.write("""
- To begin the user can either use the predefined data by selecting the toggle box or 
             use their own data in the mzML format via the drag and drop box. 
- Once the data has been selected, the spectra will be plotted automatically.
             """)
    
    # Creates two columns to display images side by side 
    col1, col2 = st.columns(2)
    with col1: 
        # Displays image and caption for using the predefined data 
        st.image(image_predefined_data, caption='Data Selection: Using data already available', width=300)
    with col2:
        st.image(image_user_data, caption='Data Selection: User utilising their own data', width=300)

    st.subheader("Changing the parameters")
    st.write("""
    - Both the settings and parameters can be changed. 
    - The parameters can be changed via the sidebar, which allows the peptide, charge state, resolution, energy ramp and isolation to be changed for the predefined data. """)
    st.image(image_parameter_selection, caption='Parameter Selection', width=300)
    
    ("""
    - The changing of parameters via the sidebar is not applicable if the user chooses to upload their own data. 
    - Once the parameters have been selected or the users data inserted, the settings of the plot can be altered. 
    - The Settings:
        - Allows for a collision induced dissociation (CID) energy to be selected, which displays the spectrum with the selected collision energy. 
        - Allows for m/z labels to be selected to label the spectrum. The threshold of these labels can be selected too.
        - Checkbox can be selected to annotate the fragments with the ions within the spectrum. 
             """)
    st.image(image_setting_selection, caption='The settings available for selection', width=300)
    
    st.subheader("Plot interactivity")
    st.write(""" 
    - Various features allow for the plot to be explored. These features include:

        - The expansion of the plot to a full screen.""")
    st.image(image_plot_expansion, caption='Button for spectrum plot expansion', width=800)
    
    ("""
        - The ability to drag the cursor of a section of the plot to zoom in for exploration.""")
    st.image(image_zoom_function, caption='The cursor drag zoom function', width=600)
    
    ("""
        - A hover tool allows the cursor, when hovering over a peak, to display the m/z, intensity and centroid data regarding each peak.
             """)
    st.image(image_hover_function, caption='The hover function in use', width=600)
    
# Spectrum tab contents 
with spectrum_tab:
    st.markdown("Explore the parameters influencing the spectra, over a series of scans. Select instructions tab for help.")


# Load predefined data based on selected parameters
if use_predefined_data:
    reader, scan_filter_list = load_predefined_data(
        selected_peptide, 
        selected_charge_state, 
        selected_resolution, 
        selected_isolation, 
        isolation=None)
else:
    # Allows user to upload their own mzML file
    raw_file = st.sidebar.file_uploader(
        "Choose a file", 
        type=['mzml'], 
        key="rawfile", 
        help="Choose an mzML file for exploration.")
    if raw_file is not None:
        reader, scan_filter_list = load_data(raw_file)
    else: 
        reader = None

# Initialise the plot labels variables to True
labels_on = True 
label_ions = True

# Streamlit layout for displaying the spectrum tab 
with spectrum_tab: 
    scol1, scol2 = st.columns([0.2, 0.5])
    with scol1:
        if reader is not None:
            st.markdown("### Settings")

            st.write("You have selected the ", selected_charge_state, "of", selected_peptide, ", this has a monoisotopic precursor m/z of ", "{:.2f}".format(mass.fast_mass2(peptide_options[selected_peptide]['sequence'], ion_type="M", charge=int(selected_charge_state[0]), aa_mass=aa_mass)))

            # Defines available collision energies 
            _available_energies = [0, 5, 10, 15, 20]
            available_energies = [e for e in _available_energies if e in scan_filter_list]

            scan_filter = st.number_input(
                "Set Collision Energy", 
                min_value=0.0,
                max_value=40.0, 
                value=0.0, step=0.5,
                format="%0.1f",
                help="Filter scans by collision energy, from 0 to 40 V in 0.5 V steps")

            if scan_filter in scan_filter_list:
                # Get range of scans for the selected collision energy 
                scan_range = (scan_filter_list[scan_filter][0], scan_filter_list[scan_filter][-1])
                spectra = [reader[i] for i in range(scan_range[0], scan_range[1] + 1)]
                selected_scan = average_spectra(spectra, filter_string=scan_filter)
            else:
                # If energy is not in the list, then interpolate the spectra for available energies
                spectra = [average_spectra(reader[scan_filter_list[energy]]) for energy in available_energies]
                interpolated_spectra = interpolate_spectra(spectra, [scan_filter], energies=available_energies)
                selected_scan = {
                    'm/z array': spectra[0]['m/z array'],
                    'intensity array': interpolated_spectra[scan_filter],
                    'scanList': {'scan': [{'scanWindowList': {'scanWindow': [{'scan window lower limit': min(spectra[0]['m/z array']),
                                                                             'scan window upper limit': max(spectra[0]['m/z array'])}]}}]}
                }

            label_threshold = st.number_input(
                "Label Threshold (%)", 
                min_value=0, 
                value=2, 
                help="Label peaks with intensity above threshold% of maximum.")
            labels_on = st.checkbox("Show m/z labels", help="Display all peak labels on plot.", value=True)
            label_ions = st.checkbox("Annotate Fragments", help="Display all fragment labels on plot", value=False)
            
            # Plot spectrum function
            def plot_spectrum(selected_scan, labels_on, label_ions, selected_peptide):
                # Create a Bokeh figure for the plot 
                _spectrum_plot = figure(
                    y_axis_label='intensity',
                    x_axis_label='m/z',
                    tools='pan,box_zoom,xbox_zoom,reset,save',
                    active_drag='xbox_zoom'
    )
                # Format the y-axis and x-axis on the plot 
                _spectrum_plot.left[0].formatter.use_scientific = True
                _spectrum_plot.left[0].formatter.power_limit_high = 0
                _spectrum_plot.left[0].formatter.precision = 1
                _spectrum_plot.y_range.start = 0

                # Set x-axis range based on scan window 
                if 'scanList' in selected_scan and 'scan' in selected_scan['scanList'] and len(selected_scan['scanList']['scan']) > 0:
                    if 'scanWindowList' in selected_scan['scanList']['scan'][0]:
                        min_mz = selected_scan['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
                        max_mz = selected_scan['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']
                        _spectrum_plot.x_range = Range1d(min_mz, max_mz, bounds="auto")

                _spectrum_plot.line(selected_scan['m/z array'], selected_scan['intensity array'], line_width=1.5, color='black')

                # Peak detection and centroid calculation
                _peaks, _properties = peak_detection(selected_scan, threshold=5, centroid=False)
                _peak_centroids = get_centroid(selected_scan, _peaks, _properties)

                # Create ColumnDataSource for peaks
                peaks_data = {
                    'x': selected_scan['m/z array'][_peaks],
                    'y': selected_scan['intensity array'][_peaks],
                    'cent': ["%.2f" % x for x in _peak_centroids]
    }
                peaks_source = ColumnDataSource(data=peaks_data if _peaks.size > 0 else {'x': [], 'y': [], 'cent': []})
                r = _spectrum_plot.circle('x', 'y', size=6, source=peaks_source, color='skyblue', alpha=0.6)

                # Configures hover tool 
                hover_tool = HoverTool(tooltips=[
                    ("m/z", "@x{0.00}"),
                    ("intensity", "@y{0.00}"),
                    ("centroid", "@cent{0.00}")
                ], renderers=[r])
                _spectrum_plot.add_tools(hover_tool)

                # Conditionally add labels if labels_on is True
                if labels_on:
                    labels = LabelSet(x='x', 
                                      y='y',  
                                      source=peaks_source, 
                                      text='cent',
                                      text_font_size='8pt', 
                                      text_color='black')
                    _spectrum_plot.add_layout(labels)
                
                # Conditionally annotate spectrum with theoretical fragments if label_ions is True
                if label_ions:
                        charge_state_cleaned = int(selected_charge_state.rstrip('+'))  # Remove '+' and convert charge state to integer

                        fragments = get_fragments(peptide_options[selected_peptide]['sequence'], charge_state_cleaned, _peak_centroids)
                           
                        # Annotate the spectrum with theoretical fragments
                        ions_data = {
                            'x': [frag['m/z'] for frag in fragments],
                            'y': [selected_scan['intensity array'][np.argmin(np.abs(selected_scan['m/z array'] - frag['m/z']))] * 1.05 for frag in fragments],
                            'ion_type': [frag['ion'] for frag in fragments]
                        }
                        ions_source = ColumnDataSource(data=ions_data)

                        ion_labels = LabelSet(x='x', 
                                              y='y', 
                                              source=ions_source, 
                                              text='ion_type', 
                                              text_font_size='8pt', 
                                              text_color='blue', 
                                              y_offset=10,
                                              x_offset=1.5)
                        _spectrum_plot.add_layout(ion_labels)

                        print(_spectrum_plot)        

                return _spectrum_plot

    with scol2:
        if reader is not None:
            _spectrum_plot = plot_spectrum(selected_scan, labels_on, label_ions, selected_peptide)
            st.bokeh_chart(_spectrum_plot, use_container_width=True)
