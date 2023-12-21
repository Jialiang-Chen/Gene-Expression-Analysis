# Gene-Expression-Analysis
An analysis on the differentially expressed genes between ERBB2+ and other breast cancer tumours

## Code Explanation
### 0. preparation 
This part loads equired libraries, sets the work directroy, and set the thresholds/cut-offs.  
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L3-L45

### 1. read files and match patient ids
The data files are loaded and processed for metadata construction. 
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L47-L87


### 2. create metadata
The metadata is created based on CNA level, and then DESeq object is built. Quality checks and data processing are performed in this part to ensure the reliability of the following analysis.
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L88-L146


### 3. differential expression analysis
The RNASeq data are normalized, and tested for differential expression. Further analysis are obtained with graphs. 
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L147-L291

### 4. pathway enrichment analysis
Pathway enrichment analysis is conducted on up-regulated and down-regulated genes separately, and are explored with various plotting functions. A PCA plot is obtained using vst values.
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L292-L338

### 5. principal component analysis
Principal componet analysis with vst values.
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L339-L360

### 6. gene expression cluster
Gene expression analysis results are further used for clustering.
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L361-L439

### 7. survival analysis
Overall survival analysis comparing subtypes and gene expression levels. 
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/73008123c7ea65d69ff0e659578e80dde44023e0/Gene%20Expression%20Analysis.R#L440-L531


