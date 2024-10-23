## PEPTIDE EXPLORER ##
"""
This script is used to explore a full MS scan of each peptide using the data provided. 
It importst the necessary libraries, defines the functions for peak detection, centroid calculation and
retrieves fragments. 
"""
import streamlit as st
import numpy as np
import pandas as pd
import requests
import io
from bokeh.models import HoverTool, ColumnDataSource, LabelSet
from bokeh.plotting import figure
from pyteomics import mzml, mass, parser
from scipy import signal

## FUNCTIONS ##

# Define the peak detection function
def peak_detection(spectrum, threshold=5, prominence=0.8, width=2, distance=4, centroid=False):
    relative_threshold = spectrum['intensity array'].max() * (threshold / 100)
    if centroid:
        peaks = np.where(spectrum['intensity array'] > relative_threshold)[0]
        return peaks
    else:
        peaks, properties = signal.find_peaks(
            spectrum['intensity array'], 
            width=width,
            height=relative_threshold, 
            distance=distance,
            prominence=prominence
            )
        return peaks, properties
    
# Define the function to get centroid values
def get_centroid(spectrum, peaks, properties):
    _centroids = np.zeros_like(peaks, dtype='float32')
    mz_array = spectrum['m/z array']
    intensity_array = spectrum['intensity array']

    for i, peak in enumerate(peaks):
        left_ip = int(properties['left_ips'][i])
        right_ip = int(properties['right_ips'][i])
        peak_range = range(left_ip, right_ip + 1)

        mz_value_range = mz_array[peak_range]
        intensity_range = intensity_array[peak_range]

        _centroids[i] = np.sum(mz_value_range * intensity_range) / np.sum(intensity_range)
    return _centroids


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

    neutral_losses = {
        '': 0, 
        '-H2O': -18.01528, 
        '-NH3': -17.02655
    }

    def is_in_peaks_data(mass):
        tolerance = 0.2 
        return any(abs(mass - peak) <= tolerance for peak in peaks_data)
   
    for pos in range(1, pep_length):
        for ion_type in ion_types:
       
            if ion_type[0] in ('a', 'b', 'c'):
                seq = ''.join(_sequence[:pos])
            elif ion_type[0] in ('x', 'y', 'z'):
                seq = ''.join(_sequence[-pos:])
            
         
            for charge in range(1, selected_charge_state + 1):
                for loss, mass_diff in neutral_losses.items():
                    # Calculate fragment mass with potential neutral loss, need to account for charge state with losses 
                    _mass = mass.fast_mass2(seq, ion_type=ion_type, charge=charge, aa_mass=aa_mass) + (mass_diff / charge)  

                    # Determine ion label based on ion type and neutral loss
                    ion_label = ion_type + str(pos) + loss + "+"*charge # adds charge to ion labels

                    # Check if calculated mass is close to any observed peaks
                    if is_in_peaks_data(_mass):
                            fragments.append({'sequence': seq, 'ion': ion_label, 'm/z': _mass, 'type': ion_type})
                            print(f"Annotated fragment: {ion_label}, m/z: {_mass}")

    # Precursor ion annotation
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


# Defines the function to load mzML data
def load_mzml_data(peptide):
    file_path = {
        'MRFA': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_full_scan_enhanced.mzML',
        'Bradykinin': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin.mzML',
        'GRGDS': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_full_ms.mzML',
        'SDGRG': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_full_ms.mzML'
    }
    file_url = file_path.get(peptide)
    if file_url:
        file_response = requests.get(file_url)
        mzml_data = io.BytesIO(file_response.content)
        spectra = mzml.read(mzml_data)
        return spectra
    else:
        st.error("No mzML data available for the chosen peptide.")
        return None
    
# Plot spectrum function 
def plot_spectrum(spectrum, show_labels, label_ions):
    mz_values = spectrum['m/z array']
    intensity_values = spectrum['intensity array']

    _peaks, _properties = peak_detection(spectrum, threshold=5, centroid=False)
    peak_centroids = get_centroid(spectrum, _peaks, _properties)

    source = ColumnDataSource(data=dict(
        x=mz_values[_peaks],
        y=intensity_values[_peaks],
        cent=["%.2f" % x for x in peak_centroids]
    ))

    TOOLTIPS = [
        ("m/z", "@x{0.00}"),
        ("intensity", "@y{0.00}"),
        ("centroid", "@cent{0.00}")
    ]

    _plot = figure(title="", x_axis_label="m/z", y_axis_label="Intensity", tooltips=None,
               tools='pan,box_zoom,xbox_zoom,reset,save', active_drag='xbox_zoom')

    _plot.line(mz_values, intensity_values, line_width=1.5, color='black')
    r = _plot.circle('x', 'y', size=5, source=source, color='red')

    if show_labels:
        labels = LabelSet(x='x', 
                          y='y', 
                          source=source, 
                          text='cent', 
                          text_font_size='8pt', 
                          text_color='black')
        _plot.add_layout(labels)

    if label_ions:
        fragments = get_fragments(peptide_options[selected_peptide]['sequence'], 3, peak_centroids)

        ions_data = {
            'x': [frag['m/z'] for frag in fragments],
            'y': [spectrum['intensity array'][np.argmin(np.abs(spectrum['m/z array'] - frag['m/z']))] * 1.05 for frag in fragments],
            'ion_type': [frag['ion'] for frag in fragments]
        }
        ions_source = ColumnDataSource(data=ions_data)

        ion_labels = LabelSet(x='x', 
                              y='y', 
                              source=ions_source, 
                              text='ion_type', 
                              text_font_size='8pt', 
                              text_color='blue', 
                              y_offset=8,
                              x_offset=1.5)
        _plot.add_layout(ion_labels)

    # Defines range for x and y axes of the plot 
    # Sets x-axis range from the minimum to the maximum m/z values 
    _plot.x_range.start = min(mz_values)
    _plot.x_range.end = max(mz_values)
    # Setting the y-axis range from 0 to 10% above the maximum intensity value 
    _plot.y_range.start = 0
    _plot.y_range.end = max(intensity_values) * 1.1

    hover_tool = HoverTool(tooltips=TOOLTIPS, renderers=[r])
    _plot.add_tools(hover_tool)
    # Displays the plot in the Streamlit app with width that adapts to the container 
    st.bokeh_chart(_plot, use_container_width=True)

# Streamlit layout, creation of tabs 
Spectrum_tab, Instruction_tab = st.tabs(["Spectrum", "Instructions"])

# Instruction tab content
with Instruction_tab:  

    image_peptide_selection = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/peptide%20selection.png'
    image_full_peptide_plot = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/plotting%20of%20peptide.png'
    image_plot_expansion = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/View%20fullscreen.png'
    image_zoom_function =  'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Zoom%20function.png'
    image_hover_function = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Hover%20function.png'

    st.header("Instructions")
    st.markdown("Instructions for the use of the Peptide Explorer")
    st.subheader("Peptide Selection")
    st.write("""
- Within this interactive peptide explorer, the entire spectra for a single peptide can be explored. 
- The user can select which peptide they would like to investigate. 
- Once selected the spectra is automatically plotted.
             """)

    col1, col2 = st.columns(2)
    with col1: 
        st.image(image_peptide_selection, caption='Data Selection: Using data already available', width=300)
    with col2:
        st.image(image_full_peptide_plot, caption='Automatic plotting of the full spectrum specific to each peptide', width=900)    

    st.subheader("Plot interactivity")
    st.write("""
    - Once plotted, the only settings that can be changed is whether the user wants to 'Show m/z Labels', 
             which can be selected via the checkbox. 
             A checkbox can be selected to annotate the fragments with the ions within the spectrum. 
    - The plot has various interactive features that allows the user to undertake a thorough exploration. These features include:

    - The expansion of the plot to a full screen.""")
    st.image(image_plot_expansion, caption='Button for spectrum plot expansion', width=800)
    
    ("""
        - The ability to drag the cursor of a section of the plot to zoom in for exploration.""")
    st.image(image_zoom_function, caption='The cursor drag zoom function', width=600)
    
    ("""
        - A hover tool allows the cursor, when hovering over a peak, to display the m/z, intensity and centroid data regarding each peak.
             """)
    st.image(image_hover_function, caption='The hover function in use', width=600)


## PAGE LAYOUT ## 

# Sidebar configuration for the Streamlit app
st.sidebar.title("Interactive Peptide Explorer")
st.sidebar.markdown("Select a peptide to explore its mass spectrometry data.")

# Options for different peptides with associated sequence 
peptide_options = {
    'MRFA': {'sequence': 'MRFA'},
    'Bradykinin': {'sequence': 'RPPGFSPFR'},
    'GRGDS': {'sequence': 'GRGDS'},
    'SDGRG': {'sequence': 'SDGRG'}
}

selected_peptide = st.sidebar.selectbox("Select Peptide", 
                                        list(peptide_options.keys()))
show_labels = st.sidebar.checkbox("Show m/z Labels", value=False)
label_ions = st.sidebar.checkbox("Annotate Fragments", value=False)

# Main content for Spectrum tab 
with Spectrum_tab:
    st.markdown("Explore a full MS scan for different peptides. Select instructions for help.")
    spectra = load_mzml_data(selected_peptide)
    if spectra is not None:
        # Get the first spectrum from the loaded data or None if no spectra 
        first_spectrum = next(spectra, None)  
        # Checks whether spectrum contains necessary data arrays 
        if first_spectrum and 'm/z array' in first_spectrum and 'intensity array' in first_spectrum:
            plot_spectrum(first_spectrum, show_labels, label_ions)
        else:
            st.error("The mzML data for the chosen peptide is not in the correct format.")
        