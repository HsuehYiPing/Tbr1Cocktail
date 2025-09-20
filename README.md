# Dietary supplement cocktails improve social behaviors and neural ensembles of multiple autism models

Tzyy-Nan Huang<sup>1</sup>, Ming-Hui Lin<sup>1,2</sup>, Tsan-Ting Hsu<sup>1</sup>, Chen-Hsin Yu<sup>1</sup>,
and Yi-Ping Hsueh<sup>1,2</sup>

Author affiliations:  
<sup>1</sup> Institute of Molecular Biology, Academia Sinica, Taipei, 11529, Taiwan, ROC  
<sup>2</sup> Molecular and Cell Biology, Taiwan International Graduate Program, Institute of Molecular Biology, Academia Sinica and Graduate Institute of Life Sciences, National Defense Medical Center, Taipei 11529, Taiwan, ROC 

Autism spectrum disorders (ASD) are heterogeneous neurodevelopmental disorders caused by both genetic variations and environmental factors. Since synaptic function is particularly susceptible to ASD-linked conditions, we hypothesized that a cocktail of various nutrients that benefit synaptic function, such as zinc, BCAA and serine, would ameliorate various ASD-related phenotypes. We explored the effects of three different concentrations of supplement cocktails on three previously established mouse models (i.e. _Tbr1_<sup>+/–</sup>, _Nf1_<sup>+/–</sup> and _Cttnbp2_ M120I mice), revealing that both one-week and two-month treatments elicited beneficial effects on social behaviors without noticeable side effects. The supplement cocktails affected the synaptic proteomes of mutant mouse brains. In vivo calcium imaging further revealed that one week of supplementation with a cocktail corrected hyperactivity and hyperconnectivity of basolateral amygdalar neurons in _Tbr1_<sup>+/–</sup> mice. Thus, nutrient cocktails containing zinc, BCAA and serine improve synaptic function, modulate neuronal activity and ensembles, and control social behaviors in multiple ASD contexts.

# Installation
From the top-level directory of the project, install the required Python packages using:

```bash
pip install -r requirements.txt
```

> This will install all the dependencies listed in `requirements.txt`. -e .


# Data preparation
Download the dataset from the release page and arrange the files according to the data structure shown below.

## Data structure
Some neuronal data may be missing due to temporary hardware connection failures during recording. Corresponding `_blank_frame.csv` files are provided in such cases.
The files `*_acc.csv` and `Multi-session_registration_*.csv` were generated using [CaImAn](https://github.com/flatironinstitute/CaImAn) (A Python toolbox for large-scale Calcium Imaging Analysis; Giovannucci et al., 2019).

```
./Calcium_imaging_data/
├── Water/
│   ├── NO.1/
│   │   ├── OE/
│   │   │   ├── *_acc.csv
│   │   │   ├── *_behavior_period.csv
│   │   │   └── *_blank_frame.csv     # included only when some neuronal data are missing
│   │   └── RSI/
│   │       ├── *_acc.csv
│   │       ├── *_behavior_period.csv
│   │       └── *_blank_frame.csv
│   ├── NO.2/
│   │   └── ...
│   └── NO.n/
│       └── ...
├── Cocktail/
│   ├── NO.1/
│   │   ├── OE/
│   │   │   ├── *_acc.csv
│   │   │   ├── *_behavior_period.csv
│   │   │   └── *_blank_frame.csv
│   │   └── RSI/
│   │       ├── *_acc.csv
│   │       ├── *_behavior_period.csv
│   │       └── *_blank_frame.csv
│   ├── NO.2/
│   │   └── ...
│   └── NO.n/
│       └── ...
└── Multi-session_registration_*.csv     # e.g., Multi-session_registration_NO.1.csv
```

## CSV File Descriptions

| CSV Filename Pattern        | Description |
| :-------------------------- | :---------- |
| *_acc.csv                    | Records neuronal activity data (generated using CaImAn)|
| *_behavior_period.csv        | Contains behavioral period annotations |
| *_blank_frame.csv            | Lists blank frames or missing data (included only when some neuronal data are missing)|
| `Multi-session_registration_*.csv` | Contains multi-session registration data (generated using CaImAn; one file per mouse) |
        
# Usage
The table below summarizes the purpose of each notebook and indicates the corresponding figures in [this paper](https://www.biorxiv.org/content/10.1101/2025.05.29.656761v1.abstract).

| Use case | Notebook | Notes |
| :----- | :----- | :----- |
| Population-level analysis of the correlation between behavior and neuronal activity | [`Neural_ensemble_analysis.py`](Cell_type_identification.py)  | Related to Fig. 2 and Supplementary Fig. 5, 6 |
| Identification of cell types at single-cell level | [`Cell_type_identification.py`](Cell_type_identification.py)  | Related to Fig. 3A, 4B |
| Parsing cell types across sessions  | [`parsing_cell_type_to_links.py`](parsing_cell_type_to_links.py) | - |
| Functional neural network analysis | [`NetworkX_analysis.ipynb`](Cell_type_identification.py)  | Related to Fig. 4 and Supplementary Fig. 7, 8 |

# Acknowledgements
  We thank Academia Sinica Common Mass Spectrometry Facilities for Proteomics and Protein Modification Analysis located at the Institute of Biological Chemistry, Academia Sinica (supported by AS-CFII-108-107), the Transgenic Core Facility located at the Institute of Molecular Biology, Academia Sinica (supported by AS-CFII-108-104), the Mass Spectrometry Facility of the Genomics Core at the Institute of Molecular Biology, Dr. John O’Brien for English editing, and members of Y.-P.H.’s laboratory who relabeled samples for blind experiments. 
