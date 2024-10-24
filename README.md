# Interactive Mass Spectrometry App for Proteomics

## Overview 

This project is an interactive web application developed and designed for exploring mass spectrometry data. The app provides tools for visualising spectrum plots, fragments, peptide/protein structures and the isolation window using a variety of interactive plots. The application is built using Streamlit and Bokeh and offers a user-friendly interface to interact and analyse proteomic data. 

## Table of Contents 

1. [Main Page: Interactive Parameter Explorer](#interactive-parameter-explorer)
2. [Peptide Explorer](#peptide-explorer)
3. [Isolation Window](#isolation-window)
4. [Fragment Visualisation](#fragment-visualisation)
    - [2D Fragment Visualisation](#2d-fragment-visualisation)
    - [Protein and Peptide Visualisation](#protein-and-peptide-visualisation)
5. [Requirements](#requirements)
6. [How to Run](#how-to-run)
7. [Usage](#usage)
8. [Data](#data)
9. [Contributing](#contributing)
10. [License](#license)
11. [Acknowledgements](#acknowledgements)

## Interactive Parameter Explorer

The **Interactive Parameter Explorer** is the main page of the application. Users are able to interactively adjust parameters of provided data and observe real-time updates on Spectra. As well as the ability to upload their own data for exploration. Users can explore various settings and their effects on the data. 

For the Interactive Parameter Explorer code, see the [Interactive_Parameter_Explorer.py](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/blob/main/Interactive_Parameter_Explorer.py) file. 

## Peptide Explorer 

The **Peptide Explorer** page provides the user with tools to visualise and analyse one full MS2 scan for each peptide. Users can select peptides and observe the fragmentation of the peptide via a spectrum plot. 

For the code of the Peptide Explorer, see the [Peptide_Explorer.py](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/blob/main/pages/Peptide_Explorer.py) file. 

## Isolation Window 

The **Isolation Window** page allows users to adjust the isolation window settings (centre and width of the isolation window) to observe how these adjustments affect the visualisation of the MS data. 

For the code of the Isolation Window, see the [Isolation_Window.py](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/blob/main/pages/Isolation_Window.py) file. 

## Fragment Visualisation

The **Fragment Visualisation** page contains two sub-sections: 

### 2D Fragment Visualisation 
This sections allows users to visualise peptide fragments in 2D. The key features include:
- Interactive plots for fragment ions linking the ions to the peptide sequence. 
- Hover tooltips for detailed fragment information. 
- The customisation of the charge state parameter. 

### Protein and Peptide Visualisation
This section provides visualisations specifically for the Bradykinin peptide and the full protein structure of Ubiquitin. 

For the code of the Fragment Visualisation, see the [Fragment_visualisation.py](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/blob/main/pages/Fragment_visualisation.py) file. 
For the code of the 2D Fragment Visualisation, see the [fragment_visualisation_2d.py](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/blob/main/pages/Fragment_visualisation_pages/fragment_visualisation_2d.py) file. 
For the code of the Peptide and Protein Visualisation, see the[Peptide_and_protein_visualisation.py](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/blob/main/pages/Fragment_visualisation_pages/Peptide_and_protein_visualisation.py) file. 

## Requirements 

- Python 3.x 
- Streamlit 
- Bokeh
- NumPy
- Pandas
- SciPy
- Pynumpress
- Pyteomics
- Requests
- io 
- os 

Install the required packages using pip:

```bash 
pip install -r requirements.txt 
```

## How to Run 
1. Clone the Repository:

```bash
git clone https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-.git
cd Interactive-Tandem-Mass-Spectrometry-App-
```

2. Run the Streamlit App:
```bash
streamlit run Interactive_Parameter_Explorer.py
```

3. Open in Browser: 
```bash
Access the application via the URL provided by Streamlit. 
``` 

## Usage 
- Interactive Parameter Explorer: Adjust parameters and observe spectra changes in real-time. 
- Peptide Explorer: Select peptide and analyse fragmentation pattern. 
- Isolation Window: Modify and view the isolation window to identify the m/z range for further analysis. 
- Fragment Visualisation: Use the 2D and peptide and protein visualisation tools to visualise ion fragmentation and protein structure. 

## Data
The application uses mzML files for MS data. The data files used by the application are all included in the repository. You can find and download them from the following directory:

- [Data Directory](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/tree/main/Data)

## Contributions 
Contributions are welcome! Please fork the repository and submit pull requests with your changes. Follow the steps below: 
1. Fork the Repository.
2. Clone Your Fork:
```bash
git clone <fork-url>
```
3. Create a Branch:
```bash
git checkout -b feature-branch
```
4. Make the Changes.
5. Push to Your Fork:
```bash
git add
git commit -m "Description of changes"
git push origin feature-branch 
```
6. Create a Pull Request. 

## License
This project is licensed under the MIT License- see the [LICENSE](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/blob/main/__init__.py) file for details. 

## Acknowledgments 
This is a fork of the [initial app created by Kiah Slater](https://github.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-).
Bokeh for interactive plotting. 
Streamlit for building the web app interface. 
Pyteomics for handling mass spectrometry data. 


