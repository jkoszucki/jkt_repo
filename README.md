## Overview
This repository contains code used to process the data and generate all figures used in the PhD thesis titled "XXX" by Janusz Koszucki. Code is organised to reflect scripts that are used in each of the corresponding chapters.

## Repository Structure

```
.
├── config
│   └── config.yml      # Master configuration file
├── README.md
└── scripts
    ├── analysis        # Data processing scripts
    ├── figures         # Scripts generating figures
    └── helpers         # Helper scripts
```

## Config

Each script loads class `config.py` which load the `config.yml` master configuration file and return the dictonary with all necessary parameters and paths. Each script across the repository loads dictionary with all the paths from this class.

Master configuration file:<br>
`input_dir` - input data<br>
`analysis_dir` - output folder for `scripts/analysis` organised per chapter.<br>
`figures_dir` - output folder for `scripts/figures` organised per chapter.<br>


## Scripts

#### analysis

Contains the code used to process the data used in the PhD thesis titled “XXX” by Janusz Koszucki. Each folder `chapterX` corresponds to analysis executed for this chapter of the PhD Thesis. 


#### figures

Contains the code used to generate all figures in the PhD thesis titled “XXX” by Janusz Koszucki. The analysis integrated analysis of Klebsiella nuclear magnetic resonance capsule polysaccharide (CPS) structures, their loci (K-loci) and prophage receptor binding proteins.



## Details

`/Users/januszkoszucki/Projects/thesis-claude/other/docs` - Documents submitted to doctoral school for extension of the PhD Thesis (eg, the schedule and plan analyses), PhD theses of other people along with the documents.
`/Users/januszkoszucki/Projects/thesis-claude/other/literature` - Reference publications which will be used to insert references in the PhD thesis.

