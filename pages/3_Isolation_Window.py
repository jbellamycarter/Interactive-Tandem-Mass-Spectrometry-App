## ISOLATION WINDOW ##
"""
This script is used to visualise the m/z range that should be isolated for further analysis, this is the isolation window. 
Again imported are the necessary libraries. The functions defined in this script are average spectra, download files, plot the spectrum and
importantly, the tophat filter. One specific MRFA file is used in this script. 
"""

import streamlit as st
import numpy as np 
import requests
import os
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure 
from bokeh.palettes import Spectral11
from pyteomics import mzml 


## FUNCTIONS ##

# Function for downloading of MRFA file
def load_MRFA_file(url, MRFA_filename):
    # Retrieve specific file content 
    file_response = requests.get(url)
    # Checks if the file request was successful by checking status code is below 400
    if file_response.ok:
        # Opens file in binary write mode 
        with open(MRFA_filename, 'wb') as file:
            file.write(file_response.content)
    else:
        st.error("Unable to open file")

# MRFA file pathway and URL 
MRFA_file_path = "12Mar2024_MJ_ultramark1621.mzML"
MRFA_url = "https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_ultramark1621.mzML"

# Checks if file exists 
if not os.path.isfile(MRFA_file_path):
    load_MRFA_file(MRFA_url, MRFA_file_path)

# Average spectra function 
def average_spectra(spectra, filter_string=None, bin_width=None):
    reference_scan = np.unique(spectra[0]['m/z array'])
    if bin_width is None:
        bin_width = np.min(np.diff(reference_scan))

    scan_min = spectra[0]['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
    scan_max = spectra[0]['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']

    reference_mz_value = np.arange(scan_min, scan_max, bin_width) 
    merge_intensity = np.zeros_like(reference_mz_value)

    for scan in spectra:
        temp_mz = scan['m/z array'] # Temporary mz 
        temp_intensity = scan['intensity array'] # Temporary intensity
        merge_intensity += np.interp(reference_mz_value, temp_mz, temp_intensity, left=0, right=0)

    mean_intensity = merge_intensity / len(spectra)

    avg_spectrum = spectra[0].copy()
    avg_spectrum['m/z array'] = reference_mz_value
    avg_spectrum['intensity array'] = mean_intensity

    avg_spectrum['scanList']['scan'][0]['filter string'] = "AV: {:.2f}-{:.2f}; {}".format(
        spectra[0]['scanList']['scan'][0]['scan start time'], 
        spectra[-1]['scanList']['scan'][0]['scan start time'], filter_string)

    return avg_spectrum

# Plots spectrum 
def plot_spectrum(MRFA_file, isolation_width, isolation_centre):
    reader = mzml.read(MRFA_file, use_index=True)
    spectra = list(reader)

    avg_spectrum = average_spectra(spectra)

    mz_values = avg_spectrum["m/z array"]
    _intensities = avg_spectrum["intensity array"]

    tophat_filter = tophat(isolation_width, isolation_centre, mz_values)
    intensities_filtered = _intensities * tophat_filter

    source_original = ColumnDataSource(
        data=dict(
            _mzs=mz_values, 
            intensities=_intensities))
    source_filtered = ColumnDataSource(
        data=dict(
            _mzs=mz_values, 
            filtered_intensities=intensities_filtered))

    _plot = figure(title="Full MRFA MS1 Scan with Isolation Window",
               x_axis_label='m/z', 
               y_axis_label='Intensity',
               x_range=(1000,1300),
               tools='pan,box_zoom,xbox_zoom,reset,save',
               active_drag='xbox_zoom')
    
    _plot.line('_mzs', 'intensities', 
               source=source_original, legend_label='Original Intensity', color=Spectral11[0], alpha=0.3, line_width=2)
    _plot.line('_mzs', 'filtered_intensities', 
               source=source_filtered, legend_label='Filtered Intensity', color=Spectral11[1], line_dash='dashed', line_width=2)

    _plot.legend.location = "top_left"
    _plot.legend.background_fill_color = 'pink'

    _plot.left[0].formatter.use_scientific = True
    _plot.left[0].formatter.power_limit_high = 0
    _plot.left[0].formatter.precision = 1
    _plot.y_range.start = 0

    return _plot

# Tophat filter 
def tophat(centre, width, x_values):
    # Stores the filter output, initialised to zero
    y_values = np.zeros_like(x_values)
    # Sets the lower boundary of tophat filter 
    lowerl = centre - (width / 2)
    # Sets the upper boundary of tophat filter
    upperl = centre + (width / 2)  
    # Applies the filter 
    y_values[(x_values > lowerl) & (x_values < upperl)] = 1
    
    return y_values


st.title("Isolation Window Visualisation")

# Isolation and Instruction tab 
Isolation_tab, Instruction_tab = st.tabs(["Isolation Window", "Instructions"])

with Instruction_tab:
    st.header("Instructions")
    st.markdown("*The isolation window is range of m/z values selected for further MS analysis.")
    st.markdown("Instructions for the use of the Isolation Window:")
    st.write("Upon selecting this Isolation Window page, a spectrum plot is generated of the full MS1 scan of MRFA with a default centrepoint of 1215.0 and width of 15.0.")
    st.write("To use this page, the parameters that can be adapted are the centrepoint and the width of the isolation window.")
    st.write("Simply use the addition and subtraction buttons on the input box to select a number for each of these parameters to identify how the adjustment of these parameters influences the isolation window.")

with Isolation_tab:
    st.sidebar.title("Isolation Window Visualisation")
    st.sidebar.markdown("Explore the impact of the isolation centre and width on the isolation window.")

    # Configures the sidebar settings 
    isolation_centre = st.sidebar.number_input(
        "Select the centrepoint value (m/z) of the isolation window",
        min_value=1.0,
        max_value=2000.0,
        value=1215.0,
        step=0.1,
        format="%0.1f")

    isolation_width = st.sidebar.number_input(
        "Select the width of the isolation window",
        min_value=1.0,
        max_value=100.0,
        value=15.0,
        step=0.1,
        format="%0.1f")

    bokeh_figure = plot_spectrum(MRFA_file_path, isolation_centre, isolation_width)
    st.bokeh_chart(bokeh_figure)

