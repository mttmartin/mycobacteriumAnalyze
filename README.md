# mycobacteriumAnalyze

Analyzes protein and associated gene data for specific Mycobacterium species.

This is a small R package that encapsulates functionality from UniProt.ws and clusterProfiler for these specific organisms and data format.

## Installation
You can install this package directly from within R using the devtools package.
```
> library(devtools)
> install_github("mttmartin/mycobacteriumAnalyze") 
```

You can also clone this repository with:
```
$ git clone https://github.com/mttmartin/mycobacteriumAnalyze.git
```

### Dependencies

The package depends on clusterProfiler and org.Mabscessus.eg.db.

## Usage

The most straight forward method for analysis is to use the "" function provided by this package. This function can take a CSV file with protein and gene information and output tables for KEGG and GO enrichment.

The get_KEGG_enrichment and get_GO_enrichment functions can also be used to retrieve clusterProfiler result objects. These can then be used with clusterProfiler to conduct further analysis such as creating various types of figures(see clusterProfiler documentation for more details on its capabilities). 

## Further Documentation

There is a PDF reference manual located within the doc directory. This file documents all of the functionality exposed to users. 

