# Single cell visualizations
Package of single cell vizualization methods 

## Install 
To download in R run the following commands

```R
library('devtools')
install_github('galelab/scViz')
library(scViz)
```
## How to execute pipeline for heatmap

```R
folder <- "./test_DE_data/" # test data a part of repository

data <- format_data(folder=folder)

clusters <- sc_heatmap(data$data)
```

### Output files for heatmap
----------------
* heatmap.png - heatmap of DE genes in each comparison 
<img src="./output_example/heatmap.png?raw=true" width="300"/>


