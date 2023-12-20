# Gene-Expression-Analysis
An analysis on the differentially expressed genes between ERBB2+ and other breast cancer tumours

## Code Explanation
### 0. preparation 
This part loads equired libraries, sets the work directroy, and set the thresholds/cut-offs.  
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/03bc61cbed9008d7d72ca3d61c60c2d0bdb81406/Gene%20Expression%20Analysis.R#L3-L38

### 1. read files and match patient ids
The data files are loaded and processed for metadata construction. 
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/03bc61cbed9008d7d72ca3d61c60c2d0bdb81406/Gene%20Expression%20Analysis.R#L41-L80

### 2. create metadata
The metadata is created based on CNA level, and then DESeq object is built. Quality checks and data processing are performed in this part to ensure the reliability of the following analysis.
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/03bc61cbed9008d7d72ca3d61c60c2d0bdb81406/Gene%20Expression%20Analysis.R#L82-L139

### 3. differential expression analysis
The RNASeq data are normalized, and tested for differential expression. Further analysis are obtained with graphs. 
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/03bc61cbed9008d7d72ca3d61c60c2d0bdb81406/Gene%20Expression%20Analysis.R#L141-L316

### 4. pathway enrichment analysis
Pathway enrichment analysis is conducted on up-regulated and down-regulated genes separately, and are explored with various plotting functions. A PCA plot is obtained using vst values.
https://github.com/Jialiang-Chen/Gene-Expression-Analysis/blob/03bc61cbed9008d7d72ca3d61c60c2d0bdb81406/Gene%20Expression%20Analysis.R#L317-L383
