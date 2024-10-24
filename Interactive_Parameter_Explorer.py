import streamlit as st

st.set_page_config(page_title="Home", 
                   layout="wide", 
                   menu_items={'about': "This application is a parameter explorer for mzML mass spectrometry data. Written by Kiah Slater and [Dr Jedd Bellamy-Carter](https://www.lboro.ac.uk/departments/chemistry/staff/jedd-bellamycarter); data generated by Mitch Jones. Loughborough University."})

st.title("Interactive Peptide MS Explorer")

"## Loughborough Proteomics"
"This is an interactive application developed for teaching proteomics at Loughborough University by [Dr Jedd Bellamy-Carter](https://www.lboro.ac.uk/departments/chemistry/staff/jedd-bellamycarter)."

"""
The sidebar on the left has several different interactive pages to explore different aspects of peptide mass spectrometry. These contain real peptide MS data and simulate a simplified instrument interface. All data given were collected on a Thermo LTQ mass spectrometer.
"""

"""
## Contributors
* This project was conceptualised and lead by Dr Jedd Bellamy-Carter (Chemistry, Loughborough University).
* App development and deployment was done by Kiah Slater _(MSc Data Science student)_ and Dr Jedd Bellamy-Carter.
* Data collection was performed by Mitch Jones _(Natural Sciences Part C student)_ and Dr Jedd Bellamy-Carter.
"""

"""
### Code
This app is written in Python using the Streamlit library and deployed to the Streamlit Community Cloud. 
The source code for this app can be found on GitHub at: [https://github.com/jbellamycarter/Interactive-Tandem-Mass-Spectrometry-App](https://github.com/jbellamycarter/Interactive-Tandem-Mass-Spectrometry-App)
"""