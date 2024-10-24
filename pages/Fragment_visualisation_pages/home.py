## FRAGMENT VISUALISATION HOME ##
"""
This script displays the necessary information and instructions to help the user utilise the fragment visualisation 
menu and navigate its features. 
"""
import streamlit as st 

## INSTRUCTIONS ##
# Displays introductory text and instructions for the fragment visualisation features within the app
def show(): 
    st.title("Fragment Visualisation")
    st.write("""
             
             Explore the fragments through visualisation. 
            
              Use the sidebar to navigate between visualisations: 
                - **2D Fragment Visualisation**: View the fragmentation pattern in a 2D Scatter plot.
                - **Peptide and Protein Visualisation**: View Bradykinin in peptide structure and the full protein structure of Ubiquitin.
                
                ### How to Use:
                1. Select the type of visualisation you would like to see from the sidebar, i.e. 2D or 3D. 
                2. Select the peptide you would like to view.
                3. Interact with the plots to explore the fragments, peptides and final protein structure.
                4. Use the available controls to customise visualisation. 
                
              For the best view of each plot, expand to full screen. 
             
              * Bradykinin and Ubiquitin plots taken directly from the Research Collaboratory for Structural Bioinformatics Protein Data Bank (RCSB PDB). 
              
             
             Important: Neutral loss ions 
              Loss of neutral fragments during fragmentation in MS. This neutral fragment loss is reflected via a change in the remaining ion's mass. 
              These neutral fragments include Water and Ammonia:
              * Loss of Water (H$_2$O): Loss of H$_2$O, often observed during peptide fragmentation, indicates the dehydration of backbone structures or side chains,
              or the cleavage of peptide bonds. 
              * Loss of Ammonia (NH$_3$): Loss of NH$_3$, again observed during peptide fragmentation, indicates the presence of the basic residues 
             like histidine, lysine or arginine. 
                """)