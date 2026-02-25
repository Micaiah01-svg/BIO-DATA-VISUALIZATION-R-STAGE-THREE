# BIO-DATA-VISUALIZATION-R-STAGE-THREE
**Cell Type Snapshot Explorer (scRNA‑seq)
AUTHOR: MICAIAH ADEOLUWA ADEDEJI
A lightweight Shiny app for exploring UMAP embeddings, marker gene expression, and gene specificity scores across cell types.**

**How to Run the App**
1. Install required packages:install.packages(c("shiny", "ggplot2", "dplyr", "viridis", "readr"))
2. Run directly from the project folder: shiny::runApp()

Or, launch the deployed app via the shinyapps.io link: https://adeoluwa.shinyapps.io/Stage_three/

**Gene Specificity Score**
Each gene is ranked by a simple specificity score: diff = mean_in − mean_out
Where:
mean_in: average expression of the gene within the selected cell type
mean_out: average expression outside that cell type
A larger diff means the gene is more specific (highly expressed inside the group and low outside).

**Marker Gene Selection**
The app automatically identifies a marker gene for the selected cell type by:
1.Computing diff for every gene
2.Ranking genes by descending diff
3.Selecting the top gene (highest specificity score) as the marker
This marker gene is then used to color the UMAP plot and summarize expression patterns.
