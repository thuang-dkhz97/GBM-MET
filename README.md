# GBM-MET

---

## Repository structure

```
BlastRNAPredict-main/
├── src/
│   ├── process-2.6/                # R scripts for processing the data.
│   │   ├── combine-data-1.R                 # Initial preprocessing, such as removal of high missingness features
│   │   ├── combine-data-2.R               # Application of missforest
│   │   ├── combine-data-3.R                    # Generation of data structures for downstream analysis.
│   ├── pipeline-2.7/        # Contains code for the ML pipeline.
│   │   ├── initial-selection.R          # External validation (test set)
│   │   │   ├── gbm-v-met.R  # Run the gbm-v-met feature selection pipeline
│   │   │   ├── gbm-v-hc.R   # Run the gbm-v-hc feature selection pipeline
│   │   │   └── met-v-hc.R   # Run the met-v-hc feature selection pipeline.
│   │   ├── filtration-after-initial-selection/       # Contains code to extract the top 27-30 features from each FSM
│   │   │   └── filtration.R
│   |   ├── second-selection/     # Contains code to rerun the predictive modelling after filtering initial feature selections
|   |   |   ├── gbm-v-met.R  # Run the gbm-v-met predictive model without feature selection
│   │   │   ├── gbm-v-hc.R   # Run the gbm-v-hc predictive model without feature selection
│   │   │   └── met-v-hc.R   # Run the met-v-hc predictive model without feature selection
│   |   ├── final-filtration/     # Contains code to extract the features of the best method
|   |   |   └── final-selection.R
│   ├── DE_pipeline-2.8/    # Run the DE pipeline to get a final set of features
│   |   ├── compute-feature-ranks/  # Obtain the ranking of features based on DE and ML
|   |   |   ├── cm_rf.R  # Helper module that uses an RF to help compute feature rankings
│   │   │   └── compute-feature-rank.R   # Primary code to conduct feature ranking
│   |   ├── RFE/  # Conduct recursive feature selection based on feature rankings
|   |   |   └── RFE.R 
│   |   ├── post-RFE-select-by-AUC/  # Scan through different thresholds for RFE and select best performing features
|   |   |   └── post-rfe.R 
│   ├── figures/    # Code to generate specific figures
│   ├── tables/    # Code to generate specific tables
│   ├── nested/    # Run nested validation pipeline
│   |   ├── preprocess/  # Preprocess proteomics
|   |   |   ├── nested-fold-creation.R  # Creates initial fold separations
│   │   │   └── nested-selection-initial.R   # Initialisation code - preprocessing, etc.
│   |   ├── preprocess-RNA/  # Preprocess transcriptomics
|   |   |   ├── nested-fold-creation.R  # Creates initial fold separations
│   │   │   └── nested-selection-initial.R   # Initialisation code - preprocessing, etc.
│   |   ├── new-implementation/  # Nested validation for proteomics
|   |   |   ├── de.R  # Helper function for differential expression analysis
|   |   |   ├── models.R  # Helper functions for running predictive models
|   |   |   ├── selections.R  # Helper functions for running feature selection methods
|   |   |   ├── nested-pipeline.R  # Main code to be called to run nested validation
|   |   |   ├── nested.sh  # Shell scripts to run nested-pipeline.R with command line arguments.
│   │   │   └── helper.R  
│   |   ├── new-implementation-RNA/  # Nested validation for transcriptomics
|   |   |   ├── de.R  # Helper function for differential expression analysis
|   |   |   ├── models.R  # Helper functions for running predictive models
|   |   |   ├── selections.R  # Helper functions for running feature selection methods
|   |   |   ├── nested-pipeline.R  # Main code to be called to run nested validation
|   |   |   ├── nested.sh  # Shell scripts to run nested-pipeline.R with command line arguments.
│   │   │   └── helper.R  
│   |   └── markdown/  # Code to generate markdowns for visual rundown of results
|   |       ├── *.Rmd  # Code to help generate markdowns
│   │       └── *.html # The markdown documents.  
├── .gitignore
└── README.md
```

---
