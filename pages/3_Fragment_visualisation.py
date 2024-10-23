## FRAGMENT VISUALISATION 2D ##
"""
This script allows for the 2D visualisation of the fragments found in each peptide. Include in the visualisation is 
a scatter plot and a table containing m/z value, the ion type and the specfic ion.
"""
import streamlit as st
import numpy as np
import pandas as pd
import requests
import io
from bokeh.models import ColumnDataSource, HoverTool, LegendItem
from bokeh.plotting import figure
from bokeh.palettes import Category10
from pyteomics import mzml, parser, mass
from scipy import signal

## FUNCTIONS ##

# Peak detection function
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

# Centroid calculation function
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

# Fragment annotation function with isolation window
def get_fragments(sequence, selected_charge_state, peaks_data, isolation_window, ion_types=('b', 'y'), calc_neutral_losses=False):
    fragments = []
    _sequence = parser.parse(sequence)
    pep_length = len(_sequence)
    print(f"Peptide length: {pep_length}")

    neutral_losses = { '' : 0}

    if calc_neutral_losses:
        neutral_losses['-H20'] = -18.01528
        neutral_losses['-NH3'] = -17.02655

    def is_in_peaks_data(mass):
        tolerance = 0.2
        lower_bound = isolation_window[0] - tolerance
        upper_bound = isolation_window[1] + tolerance
        return any(lower_bound <= peak <= upper_bound for peak in peaks_data)

    for pos in range(1, pep_length):
        for ion_type in ion_types:
            if ion_type[0] in ('a', 'b', 'c'):
                seq = ''.join(_sequence[:pos])
            elif ion_type[0] in ('x', 'y', 'z'):
                seq = ''.join(_sequence[-pos:])
            
            for charge in range(1, selected_charge_state + 1):
                for loss, mass_diff in neutral_losses.items():
                    _mass = mass.fast_mass2(seq, ion_type=ion_type, charge=charge, aa_mass=aa_mass) + (mass_diff / charge)
                    ion_label = '[' + ion_type + str(pos) + loss + ']' + '+'*charge

                    if is_in_peaks_data(_mass):
                        fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': ion_type, 'start': _mass - 0.5, 'end': _mass + 0.5})
                        print(f"Annotated fragment: {ion_label}, m/z: {_mass}")


    # Precursor ion annotation
    for charge in range(1, selected_charge_state + 1):
        for loss, mass_diff in neutral_losses.items():
            seq = ''.join(_sequence)
            _mass = mass.fast_mass2(seq, ion_type="M", charge=charge, aa_mass=aa_mass) + (mass_diff / charge)
            ion_label = "M" + loss + "+"*charge

            if is_in_peaks_data(_mass):
                fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': "M", 'start': _mass - 0.5, 'end': _mass + 0.5})
                print(f"Annotated fragment: {ion_label}, m/z: {_mass}")
    
    return fragments

# Function to load mzML data
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
            if file_response.ok:
                mzml_data = io.BytesIO(file_response.content)
                spectra = list(mzml.read(mzml_data))
                return spectra
            else:
                print(f"Unable to retrieve file.")
                return None
    else:
        print("mzML data unavailable for the selected peptide.")
        return None


# Plot fragments function
def plot_fragments(fragments, sequence):
    if not fragments:
        st.write("Unable to plot: No fragments found.")
        return
    
    # Extracts m/z values, ion labels, ion types and the peptide sequence from the fragments
    mz_values = [fragment['m/z'] for fragment in fragments]
    ion_labels = [fragment['ion'] for fragment in fragments]
    ion_types = [fragment['type'] for fragment in fragments]
    pep_sequence = [len(fragment['seq']) for fragment in fragments]

    ion_start = [fragment.get('start', fragment['m/z'] - 0.5) for fragment in fragments]
    ion_end = [fragment.get('end', fragment['m/z'] + 0.5) for fragment in fragments]

    # Define colors for different ion types
    _colours = {'b': 'blue', 'y': 'red', 'M': 'orange', '': 'pink'}
    colour_values = [_colours.get(ion_type, 'black') for ion_type in ion_types]

    # ColumnDataSource is created from the fragment data
    fragment_data = ColumnDataSource(data=dict(
        mz_values=mz_values,
        ion_labels=ion_labels,
        ion_start=ion_start,
        ion_end=ion_end,
        ion_types=ion_types,
        pep_sequence=pep_sequence,
        color_values=colour_values
    ))
    
    # Creates the Bokeh figure to be plotted 
    _plot = figure(
        title="Fragment Ion Visualisation",
        x_axis_label='m/z',
        y_axis_label='Ion',
        y_range=fragment_data.data['ion_labels'],
        tools='pan,box_zoom,xbox_zoom,reset,save',
        active_drag='xbox_zoom'
    )

    # Add horizontal bars for each fragment ion
    h_bars = _plot.hbar(
        y='ion_labels',
        left='ion_start',
        right='ion_end',
         source=fragment_data,
        height=0.7,
        color='color_values'
    )

    # Add HoverTool for displaying tooltips
    hover_tool = HoverTool(
        tooltips=[
            ("m/z", "@mz_values"),
            ("Ion", "@ion_labels"),
            ("Type", "@ion_types"),
            ("Sequence", "@pep_sequence")
        ], 
        renderers=[h_bars]
    )
    _plot.add_tools(hover_tool)
 
    offset = 4

    # Add the sequence text labels to the plot
    for i, (mz, ion_label, pep_seq) in enumerate(zip(mz_values, ion_labels, pep_sequence)):
        _plot.text(
            x=[mz + offset], 
            y=[ion_label], 
            text=[sequence[:pep_seq]], 
            text_font_size="8pt", 
            text_align="left", 
            text_baseline="middle"
        )

    # Set x-axis ticks and labels
    if mz_values:
        max_mz_value = int(np.ceil(max(mz_values)))  # Round up to the nearest integer
        x_ticks = list(range(0, max_mz_value + 1, 10))  # Create ticks at every 10 units
        _plot.xaxis.ticker = x_ticks
        _plot.xaxis.major_label_overrides = {tick: str(tick) for tick in x_ticks}
    else:
        _plot.xaxis.ticker = [0]

    _plot.xaxis.major_label_orientation = 1

    # Create legend items for each ion type
    ion_legend_items = []
    unique_ions = set(ion_types)  # To avoid duplicate legend items: get unique ion types 
    for ion_type in unique_ions:
        colour = _colours.get(ion_type, 'pink')
        
        # Create a temporary plot for the legend item
        legend_source = ColumnDataSource(data=dict(x=[0], y=[0]))
        
        # Create a glyph for the legend item
        ion_legend_glyph = _plot.scatter(
            x='x',
            y='y',
            source=legend_source,
            color=colour,
            size=8,
            legend_label=ion_type
        )
        
        ion_legend_items.append(LegendItem(label=f"{ion_type}", renderers=[ion_legend_glyph]))

    # Configure the aesthetic of the legend
    _plot.legend.title = 'Ion Type'
    _plot.legend.location = 'bottom_right'
 
    st.bokeh_chart(_plot, use_container_width=True)

fragment_mass_tab, instructions_tab = st.tabs(["Fragment Masses", "Instructions"])

with fragment_mass_tab:
    peptide_options = {
        'MRFA': 'MRFA',
        'Bradykinin': 'RPPGFSPFR',
        'GRGDS': 'GRGDS',
        'SDGRG': 'SDGRG'
    }

    selected_peptide_name = st.sidebar.selectbox(
        "Select a peptide sequence", 
        list(peptide_options.keys())
    )

    peptide_sequence = peptide_options[selected_peptide_name]
    st.sidebar.write(f"Sequence of Selected Peptide: {peptide_sequence}")

    # User-selectable slider to select charge state
    selected_charge_state = st.sidebar.slider(
        "Select Precursor Charge State", 
        min_value=1, 
        max_value=3, 
        value=1,
        help="Will calculate fragments with up to this charge state.")

    calc_neutral_losses = st.sidebar.checkbox("Include neutral losses", value=False)

    isolation_window = (0.0, 2000.0)
   
    # Loads the mzML data
    mzml_data = load_mzml_data(selected_peptide_name)
    if mzml_data is not None:
        spectrum = mzml_data[0]
        _peaks, _properties = peak_detection(spectrum)

        # Retrieves centroid m/z values for detected peaks 
        peaks_data = get_centroid(spectrum, _peaks, _properties)

        # Annotates fragments based on peaks data, peptide sequence, selected charge state and isolation window
        fragments = get_fragments(peptide_sequence, selected_charge_state, peaks_data, isolation_window, calc_neutral_losses=calc_neutral_losses)
        
        # Pass both fragments and peptide_sequence to plot_fragments
        plot_fragments(fragments, peptide_sequence)
        
        # Convert fragments data to a DataFrame for display 
        df_fragments = pd.DataFrame(fragments)
        # Display the fragments in a table, with columns as specified 
        st.table(df_fragments[['m/z', 'type', 'ion']])

with instructions_tab:
    st.header("Instructions")
    st.write("""
             
             Explore the fragments through visualisation. 
            
             View the fragmentation pattern in a 2D Scatter plot.
                
                ### How to Use:
                1. Select the type of visualisation you would like to see from the sidebar, i.e. 2D or 3D. 
                2. Select the peptide you would like to view.
                3. Interact with the plots to explore the fragments, peptides and final protein structure.
                4. Use the available controls to customise visualisation. 
                
             Important: Neutral loss ions 
              Loss of neutral fragments during fragmentation in MS. This neutral fragment loss is reflected via a change in the remaining ion's mass. 
              These neutral fragments include Water and Ammonia:
              * Loss of Water (H$_2$O): Loss of H$_2$O, often observed during peptide fragmentation, indicates the dehydration of backbone structures or side chains,
              or the cleavage of peptide bonds. 
              * Loss of Ammonia (NH$_3$): Loss of NH$_3$, again observed during peptide fragmentation, indicates the presence of the basic residues 
             like histidine, lysine or arginine. 
                """)