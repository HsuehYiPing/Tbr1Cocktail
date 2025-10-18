# Low-dose mixtures of dietary nutrients ameliorate behavioral deficits in multiple mouse models of autism

Tzyy-Nan Huang<sup>1</sup>, Ming-Hui Lin<sup>1,2</sup>, Tsan-Ting Hsu<sup>1</sup>, Chen-Hsin Yu<sup>1</sup>,
and Yi-Ping Hsueh<sup>1,2</sup>

Author affiliations:  
<sup>1</sup> Institute of Molecular Biology, Academia Sinica, Taipei, 11529, Taiwan, Republic of China 
<sup>2</sup> Molecular and Cell Biology, Taiwan International Graduate Program, Institute of Molecular Biology, Academia Sinica and Graduate Institute of Life Sciences, National Defense Medical Center, Taipei 11529, Taiwan, Republic of China

Autism spectrum disorders (ASD) are a group of heterogeneous, behaviorally defined neurodevelopmental conditions influenced by both genetic and environmental factors. Here, we show that supplementation of low-dose multiple nutrients—an important environmental factor—can modulate synaptic proteomes, reconfigure neural ensembles, and improve social behaviors in mice. We first used _Tbr1_<sup>+/–</sup> mice, a well-established model of ASD, to investigate the effect of nutrient cocktails containing zinc, branched-chain amino acids (BCAA), and serine, all of which are known to regulate synapse formation and activity. Supplementation of nutrient cocktails for seven days altered total proteomes by increasing synapse-related proteins. Our results further revealed that _Tbr1_ haploinsufficiency promotes hyperactivation and hyperconnectivity of basolateral amygdala (BLA) neurons, enhancing the activity correlation between individual neurons and their corresponding ensembles. Nutrient supplementation normalized the activity and connectivity of the BLA neurons in _Tbr1_<sup>+/–</sup> mice during social interactions. We further show that although a low dose of individual nutrients did not alter social behaviors, treatment with supplement mixtures containing low-dose individual nutrients improved social behaviors and associative memory of _Tbr1_<sup>+/–</sup> mice, implying a synergistic effect of combining low-dose zinc, BCAA, and serine. Moreover, the supplement cocktails also improved social behaviors in _Nf1_<sup>+/–</sup> and _Cttnbp2_<sup>+/M120I</sup> mice, two additional ASD mouse models. Thus, our findings reveal aberrant neural connectivity in the BLA of _Tbr1_<sup>+/–</sup> mice and indicate that dietary supplementation with zinc, BCAA, and/or serine offers a safe and accessible approach to mitigate neural connectivity and social behaviors across multiple ASD models.

# Installation
From the top-level directory of the project, create the Conda environment with Python 3.9:

```bash
conda create -n calcium_env python=3.9
```

After creating the environment, activate it:

```bash
conda activate calcium_env

```

Then install all required Python packages using:

```bash

pip install -r requirements.txt
```


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
| Population-level analysis of the correlation between behavior and neuronal activity | [`Neural_ensemble_analysis.ipynb`](Neural_ensemble_analysis.ipynb)  | Related to Fig. 4 and Supplementary Fig. 6, 7 |
| Identification of cell types at single-cell level | [`Cell_type_identification.py`](Cell_type_identification.py)  | Related to Fig. 5A, 5B |
| Parsing cell types across sessions  | [`parsing_cell_type_to_links.py`](parsing_cell_type_to_links.py) | - |
| Functional neural network analysis | [`NetworkX_analysis.ipynb`](NetworkX_analysis.ipynb)  | Related to Fig. 6 and Supplementary Fig. 8, 9 |

# Acknowledgements
  We thank Academia Sinica Common Mass Spectrometry Facilities for Proteomics and Protein Modification Analysis located at the Institute of Biological Chemistry, Academia Sinica (supported by AS-CFII-108-107), the Mass Spectrometry Facility of the Genomics Core at the Institute of Molecular Biology, Dr. John O’Brien for English editing, and members of Y.-P.H.’s laboratory who relabeled samples for blind experiments.
