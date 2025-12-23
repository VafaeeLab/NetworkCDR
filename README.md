# Network-Based Combinatorial Drug Repositioning with Protein Context

Drug repurposing holds vast potential for speeding up clinical outcomes by utilising drugs that are known to have low toxicity. Drugs can also be combined to increase efficacy and reduce toxicity. This project aims to find drug combinations through network-based repurposing utilising human protein-protein interactions, improving on previous attempts by also including contextual protein information such as the expression level of proteins in different cell-types and cell-lines. 

This is a project for the UNSW Summer Vacation Research Scholarship program.


Personal Repo: https://github.com/Enough-Said/svrs  

VafaeeLab Repo: https://github.com/VafaeeLab/NetworkCDR

### TODO
- Many of the diseases in verifiedCombinations do not have exact matches in the graph so need fuzzy string matching, some sort of disease ID, or manually rename diseases in verifiedCombinations to match disease names in the graph
- Need a numerical value representing the predictive ability before and after filtration (ideally with a p-value or other statistical result)

### Repository Structure
```bash
├── CDCDB/                          # Full dataset from CDCDB (Not included)
├── README.md                       # Readme file containing information pertaining to the project
├── analysisScripts/                # Directory containing scripts for utilising the clean data for graph creation, filtration, usage, and analysis
│   ├── buildGraph.r                # Contains functions for creating the graph from a dataset
│   ├── examplePipeline.r           # An example of how to use the functions created in this project
│   ├── explore.r                   # Some simple exploration into the graph and its properties
│   ├── filterEdges.r               # Initial functions for filtering edges (Outdated)
│   ├── filterNodes.r               # Initial functions for filtering nodes (Outdated)
│   ├── filters.r                   # New filtration function that is better optimised
│   ├── findDistance.r              # Functions for determining the distance between drugs and diseases in the graph
│   ├── findDrugs.r                 # Functions for finding drug combinations for given diseases based on the graph
│   ├── methodTesting.r             # Script to get numerical and statistical results on the efficacy of this model (Incomplete)
│   └── visualisation.r             # Functions to visualise the graph or parts of the graph
├── clean/                          # Cleaned datasets that are ready for usage in this project
│   ├── baseGraph.rds               # The base graph with no filtration
│   ├── cellLine.tsv                # Cell line RNA expression data from protein atlas
│   ├── cellType.json               # Cell type RNA expression data from protein atlas
│   ├── disease.json                # Diseases and associated genes from DisGeNET (from enrichR site)
│   ├── drugs.json                  # Drugs and target genes from Open Targets
│   ├── ppi.json                    # Protein-protein interactions from stringdb (Only known and confident PPI)
│   ├── splitData.txt               # Which data files are too large and need splitting to upload into github
│   ├── subcell.json                # Subcellular RNA expression data from protein atlas
│   ├── tissue.json                 # RNA expression by tissue from protein atlas
│   └── verifiedCombinations.RData  # Known drug combinations from CDCDB
├── dataScripts/                    # Directory of scripts used to acquire, manipulate, and clean data
│   ├── buildBaseGraph.r            # Script to build the base graph
│   ├── cleanCombs.r                # Cleaning the drug combination data from CDCDB
│   ├── combineData.sh              # Used to recombine data split by splitData.sh for use in the project
│   ├── dataAcquisition.sh          # Some shell scripts used during data acquisition
│   ├── makeFormatConsistent.r      # Cleaning data and ensuring all the data uses the same consistent IDs
│   └── splitData.sh                # Used to split larger data files into smaller ones to fit into the github size limit
└── rawdata/                        # Raw datasets from various sources (Not included)
```

### Main References
- https://doi.org/10.1016/j.patter.2021.100325
- https://platform.opentargets.org/
- https://maayanlab.cloud/Enrichr/#libraries
- https://string-db.org/
- https://www.proteinatlas.org/
