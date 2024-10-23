## PEPTIDE AND PROTEIN VISUALISATION ##
"""
This script creates a plot for visualisation of both the peptide Bradykinin and the protein Ubiquitin. 
Only Streamlit and Streamlit components are imported as libraries. With the function defined to generate a 
JSmol viewer using JavaScript. Alongside the show function to allow the importation of this script when called upon 
by the navigation menu. 
"""

import streamlit as st
import streamlit.components.v1 as components 
from streamlit_molstar import st_molstar_rcsb

# # Function for the generation of the HTML for the JSmol viewer
# def generate_jsmol_html(pdb_id, display_mode="cartoon"):
#     return f"""
#     <!-- Incorporates JSmol JavaScript library -->
#     <script type="text/javascript" src="https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js"></script>
    
#     <!-- Initialises JSmol viewer with options for configuration --> 
#     <script type="text/javascript">

#         // Generates new JSmol applet 
#         Jmol.getApplet("JSmolApplet1", {{
#             height: 700,
#             width: 600,
#             use: "HTML5",                   // HTML5 used for rendering 
#             j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",  // File pathway for the JSmol JavaScript
#             script: "load https://files.rcsb.org/view/{pdb_id}.pdb; {display_mode} only; color structure;"    // Loads and applies visualisation settings to the PDB file
#         }});
#     </script>

#     <!-- Container defined to hold JSmol viewer --> 
#     <div id="JSmolApplet1" style="display: flex;"></div>
#     """

# # Streamlit feature to visualize peptide/protein structure in 3D using JSmol viewer
def show():
    st.title("3D Protein Structure Visualization")
    st.markdown("Explore the 3D protein structure of Bradykinin as a peptide and the full protein structure of Ubiquitin.") 
    st.write("With structures directly taken from RCSB Protein Data Bank to provide a visual image of both a peptide structure and a protein structure.")

#     peptide_options = {"Bradykinin", "Ubiquitin"}
#     selected_peptide = st.sidebar.selectbox("Select Peptide", peptide_options)

#     # Sidebar to alllow user to select style of visualisation 
#     display_mode = st.sidebar.selectbox("Select the style of visualisation:",["Cartoon", "Spacefill"]
#     )

#     # PDB ID for specific peptide/protein
#     peptide_pdb = {
#         "Ubiquitin": '1UBQ',
#         "Bradykinin": '6F3V', 
        
#     }
#     # Gets PDB ID when user selects peptide from dictionary above 
#     pdb_id = peptide_pdb[selected_peptide]

#     # Generates HTML JSmol viewer of selected peptide with chosen style 
#     JSmol_viewer = generate_jsmol_html(pdb_id, display_mode)
#     # Embedding of the JSmol viewer into Streamlit app 
#     st.components.v1.html(JSmol_viewer, height=650)
    st_molstar_rcsb('1UBQ')

if __name__ == "__main__":
    show()


