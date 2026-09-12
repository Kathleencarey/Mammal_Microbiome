# Data and code to reproduce the analyses from: Carey K.A., Cove M.V., Hui-Ling L., Baruzzi C. Similarity between Skin and Gut Microbiomes in Pine Woodland Small Mammals. In Review. 
This repository provides the data and R code to reproduce the analyses from: Carey K.A., Cove M.V., Hui-Ling L., Baruzzi C. Similarity between Skin and Gut Microbiomes in Pine Woodland Small Mammals. In Review.

## Data Folder Description 
1. The dataset `Bacteria_Reads.csv` contains data on the number of bacterial taxa reads for each small mammal microbiome sample, with the full taxonomic ranking available. Data was obtained from Zymo Research. This dataset was used for the bacterial perMANOVA analysis, bacterial differential abundance analysis, and the bacterial observed richness and diveristy comparison between the bacterial skin and gut microbiome.

2. The dataset `Fungi_Reads.csv` contains data on the number of fungal taxa reads for each small mammal microbiome sample, with the full taxonomic ranking available. Data was obtained from Zymo Research. This dataset was used for the fungal perMANOVA analysis, fungal differential abundance analysis, and the fungal observed richness and diveristy comparison between the bacterial skin and gut microbiome.

3. The dataset `Mammal_Capture_Data.csv` contains information on each skin and gut microbiome sample collected and was used in all analyses. 

## Code Folder Description
1. The file  `Mammal_Bacteria_Community_Composition.R` runs the code for the bacterial perMANOVA analysis and the bacterial differential abundance analysis between the skin and gut microbiome. It also runs the associated NMDS figure.
   
2. The file `Mammal_Bacteria_Richness_Diversity.R` runs the code for the observed bacterial richness and diversity between the skin and gut microbiome. It also plots these observations in ggplot2 by individual small mammal sample pair.
   
3. The file `Mammal_Fungi_Community_Composition.R` runs the code for the fungal perMANOVA analysis and the fungal differential abundance analysis between the skin and gut microbiome. It also runs the associated NMDS figure.
   
4. The file `Mammal_Fungi_Richness_Diversity.R` runs the code for the observed fungal richness and diversity between the skin and gut microbiome. It also plots these observations in ggplot2 by individual small mammal sample pair.
