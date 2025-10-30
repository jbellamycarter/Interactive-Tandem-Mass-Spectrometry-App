## PEPTIDE AND PROTEIN VISUALISATION ##
"""
This script creates a plot for visualisation of both the peptide Bradykinin and the protein Ubiquitin. 
Only Streamlit and Streamlit components are imported as libraries. With the function defined to generate a 
JSmol viewer using JavaScript. Alongside the show function to allow the importation of this script when called upon 
by the navigation menu. 
"""

import streamlit as st
from streamlit_molstar import st_molstar_rcsb

st.title("3D Protein Structure Visualization")
st.markdown("Explore the 3D protein structure of Bradykinin as a peptide and the full protein structure of Ubiquitin.") 
st.write("With structures directly taken from RCSB Protein Data Bank to provide a visual image of both a peptide structure and a protein structure.")

peptide_options = {'Bradykinin':'6F3V', 'Ubiquitin':'1UBQ', 'Custom': ''}
selected_peptide = st.sidebar.selectbox("Select Peptide", peptide_options.keys())

rcsb_id = st.sidebar.text_input("RCSB PDB ID", value=peptide_options[selected_peptide], max_chars=4, help="Enter 4-character PDB code.")
st_molstar_rcsb(rcsb_id)


