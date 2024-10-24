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
from bokeh.models import ColumnDataSource, LabelSet, HoverTool, LegendItem
from bokeh.plotting import figure
from bokeh.palettes import Category10
from pyteomics import mzml, parser, mass
from scipy import signal

## FUNCTIONS ##

aa_mass = mass.std_aa_mass.copy()
aa_mass['p'] = 79.966331  # phosphorylation (STY)
aa_mass['ox'] = 15.994915  # oxidation (MW)
aa_mass['d'] = 0.984016  # deamidation (NQ)
aa_mass['am'] = -0.984016  # amidation
#aa_mass['-NH2'] = -0.984016 # amidation (C-term)

# Fragment annotation function with isolation window
def get_fragments(sequence, selected_charge_state, isolation_window, ion_types=('b', 'y'), calc_neutral_losses=False):
    fragments = []
    _sequence = parser.parse(sequence)

    if parser.is_term_mod(_sequence[0]):
        n_term = _sequence.pop(0)
        _sequence[0] = n_term + _sequence[0]
    
    if parser.is_term_mod(_sequence[-1]):
        c_term = _sequence.pop(-1)
        _sequence[-1] = _sequence[-1] + c_term

    print(_sequence)
    pep_length = len(_sequence)
    print(f"Peptide length: {pep_length}")

    neutral_losses = { '' : 0}

    if calc_neutral_losses:
        neutral_losses['-H20'] = -18.01528
        neutral_losses['-NH3'] = -17.02655

    for pos in range(1, pep_length):
        for ion_type in ion_types:
            if ion_type[0] in ('a', 'b', 'c'):
                _start = 0
                _end = pos
                
            elif ion_type[0] in ('x', 'y', 'z'):
                _start = pep_length - pos
                _end = pep_length
            
            seq = ''.join(_sequence[_start:_end])
            print(seq)
            
            for charge in range(1, selected_charge_state + 1):
                for loss, mass_diff in neutral_losses.items():
                    _mass = mass.fast_mass2(seq, ion_type=ion_type, charge=charge, aa_mass=aa_mass) + (mass_diff / charge)
                    ion_label = '[' + ion_type + str(pos) + loss + ']' + '+'*charge

                    fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': ion_type, 'start': _start, 'end': _end})
                    print(f"Annotated fragment: {ion_label}, m/z: {_mass}")


    # Precursor ion annotation
    for charge in range(1, selected_charge_state + 1):
        for loss, mass_diff in neutral_losses.items():
            seq = ''.join(_sequence)
            _mass = mass.fast_mass2(seq, ion_type="M", charge=charge, aa_mass=aa_mass) + (mass_diff / charge)
            ion_label = "M" + loss + "+"*charge

            fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': "M", 'start': 0, 'end': pep_length})
            print(f"Annotated fragment: {ion_label}, m/z: {_mass}")
    
    return fragments

# Plot fragments function
def plot_fragments(fragments, sequence):
    if not fragments:
        st.write("Unable to plot: No fragments found.")
        return
    
    # Extracts m/z values, ion labels, ion types and the peptide sequence from the fragments
    mz_values = [fragment['m/z'] for fragment in fragments]
    ion_labels = [fragment['ion'] for fragment in fragments]
    ion_types = [fragment['type'] for fragment in fragments]
    pep_sequence = [fragment['seq'] for fragment in fragments]

    ion_start = [fragment['m/z'] - 0.5 for fragment in fragments]
    ion_end = [fragment['m/z'] + 0.5 for fragment in fragments]

    # Define colors for different ion types
    _colours = {'b': 'steelblue', 'y': 'firebrick', 'M': 'orange', '': 'pink'}
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

    frag_scatter = _plot.scatter(
        x='mz_values',
        y='ion_labels',
        color='color_values',  # Reference the color values from the data source
        size=8,
        source=fragment_data,
        legend_field='ion_types'
    )    

    # Add HoverTool for displaying tooltips
    hover_tool = HoverTool(
        tooltips=[
            ("m/z", "@mz_values"),
            ("Ion", "@ion_labels"),
            ("Type", "@ion_types"),
            ("Sequence", "@pep_sequence")
        ],
        renderers = [frag_scatter])
    _plot.add_tools(hover_tool)

    # Add the sequence text labels to the plot
    for i, (mz, ion_label, pep_seq) in enumerate(zip(mz_values, ion_labels, pep_sequence)):
        _plot.text(
            x=[mz + 4], 
            y=[ion_label], 
            text=[ion_label], 
            text_font_size="8pt", 
            text_align="left", 
            text_baseline="middle"
        )

    # Configure the aesthetic of the legend
    _plot.legend.title = 'Ion Type'
    _plot.legend.location = 'bottom_right'
    _plot.legend.click_policy = 'hide'

    return _plot

def plot_fragment_coverage(fragments, sequence):
    mz_values = [fragment['m/z'] for fragment in fragments]
    ion_seq = [fragment['seq'] for fragment in fragments]
    ion_labels = [fragment['ion'] for fragment in fragments]
    ion_types = [fragment['type'] for fragment in fragments]
    ion_start = [fragment['start'] - 0.5 for fragment in fragments]
    ion_end = [fragment['end'] - 0.5 for fragment in fragments]

    # Define colors for different ion types
    colors = {'b': 'steelblue', 'y': 'firebrick', 'M': 'orange', '': 'black'}
    color_values = [colors.get(ion_type, 'black') for ion_type in ion_types]

    # Create a ColumnDataSource from the fragment data
    fragment_data = ColumnDataSource(data=dict(
        mz_values=mz_values,
        ion_seq=ion_seq,
        ion_start=ion_start,
        ion_end=ion_end,
        ion_labels=ion_labels,
        ion_types=ion_types,
        color_values=color_values  # Add color values to the data source
    ))

    p = figure(x_axis_label='Residue',
               y_axis_label='Ion',
               y_range=fragment_data.data['ion_labels'])
    p.hbar(y='ion_labels', left='ion_start', right='ion_end', height=0.6, color='color_values', source=fragment_data)

    p.xaxis.ticker = list(range(0, len(sequence)))
    p.xaxis.major_label_overrides = dict(zip(range(len(sequence)), list(sequence)))

    hover = HoverTool(tooltips=[
        ("m/z", "@mz_values"),
        ("Ion", "@ion_labels"),
        ("Type", "@ion_types"),
        ("Sequence", "@ion_seq")
    ])
    p.add_tools(hover)

    st.bokeh_chart(p, use_container_width=True)

fragment_mass_tab, fragment_coverage_tab, instructions_tab = st.tabs(["Fragment Masses", "Fragment Coverage", "Instructions"])

with fragment_mass_tab:
    st.header("Fragment m/z visualisation")
    peptide_options = {
        'MRFA': 'MRFA',
        'Bradykinin': 'RPPGFSPFR',
        'GRGDS': 'GRGDS',
        'SDGRG': 'SDGRG',
        'Substance P': 'RPKPQQFFGLM-NH2',
        'Custom': ''
    }


    selected_peptide_name = st.sidebar.selectbox(
        "Select a peptide sequence", 
        list(peptide_options.keys())
    )

    if selected_peptide_name is 'Custom':
        peptide_sequence = st.sidebar.text_input("Enter 1-letter peptide sequence", value='', help="One-letter amino acid sequence. ModX support available: _ox = oxidation_")
        peptide_sequence = "".join(peptide_sequence.split())
    else:
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

    # Annotates fragments based on peaks data, peptide sequence, selected charge state and isolation window
    fragments = get_fragments(peptide_sequence, selected_charge_state, isolation_window, calc_neutral_losses=calc_neutral_losses)
    
    # Pass both fragments and peptide_sequence to plot_fragments
    _frag_plot = plot_fragments(fragments, peptide_sequence)
    st.bokeh_chart(_frag_plot, use_container_width=True)
    
    # Convert fragments data to a DataFrame for display 
    df_fragments = pd.DataFrame(fragments)
    # Display the fragments in a table, with columns as specified 
    st.dataframe(df_fragments[['m/z', 'type', 'ion']])

with fragment_coverage_tab:
    st.header("Fragment coverage visualisation")

    plot_fragment_coverage(fragments, peptide_sequence)

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