## Protein Digester ##
"""
This script generates tryptic peptides from intact protein sequences.
"""
import streamlit as st
import numpy as np
import pandas as pd
from pyteomics import mzml, parser, mass

aa_mass = mass.std_aa_mass.copy()
aa_mass['p'] = 79.966331  # phosphorylation (STY)
aa_mass['ox'] = 15.994915  # oxidation (MW)
aa_mass['d'] = 0.984016  # deamidation (NQ)
aa_mass['am'] = -0.984016  # amidation


## Interface ##
st.title("Protein digester")
st.markdown("Generate *tryptic* peptide sequences from intact protein sequence.") 
st.write("Enter the proten sequence and adjust parameters in the sidebar to suit your experiment.")

st.markdown("Once you've generated your peptides below, copy the relevant mass values and go to [MASCOT Peptide Mass Fingerprint](https://www.matrixscience.com/cgi/search_form.pl?FORMVER=2&SEARCH=PMF).")

protein_sequence = st.text_input("Enter 1-letter protein sequence", value='', help="One-letter amino acid sequence. ModX support available: _ox = oxidation_")
protein_sequence = "".join(protein_sequence.split())

num_missed_cleavages = st.sidebar.slider("Select number of missed cleavages",
                                        min_value=0,
                                        max_value=2,
                                        value=1,
                                        help="Accounts for inefficient cleavage by protease.")

min_peptide_length = st.sidebar.slider("Minimum peptide length",
                                        min_value=1,
                                        max_value=10,
                                        value=5,
                                        help="Shortest peptide expected.")



peptides = pd.DataFrame(parser.xcleave(protein_sequence, 'trypsin', num_missed_cleavages),
                        columns=["Start", "Sequence"])
peptides["Start"] += 1
peptides["Length"] = [parser.length(x) for x in peptides["Sequence"]]

peptides["Mass (Da)"] = ["{0:.2f}".format(mass.fast_mass2(seq, ion_type="M", charge=0, aa_mass=aa_mass)) for seq in peptides["Sequence"]]
peptides["m/z (z=1)"] = ["{0:.2f}".format(mass.fast_mass2(seq, ion_type="M", charge=1, aa_mass=aa_mass)) for seq in peptides["Sequence"]]
peptides["m/z (z=2)"] = ["{0:.2f}".format(mass.fast_mass2(seq, ion_type="M", charge=2, aa_mass=aa_mass)) for seq in peptides["Sequence"]]
peptides["m/z (z=3)"] = ["{0:.2f}".format(mass.fast_mass2(seq, ion_type="M", charge=3, aa_mass=aa_mass)) for seq in peptides["Sequence"]]

peptides.style.bar(subset=["Mass (Da)", "m/z (z=1)"], color='#d65f5f')

st.dataframe(peptides[(peptides["Length"] >= min_peptide_length)],
             hide_index=True)
