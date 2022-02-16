# ReporTree

ReporTree was developed at the Bioinformatics Unit of the National Institute of Health Dr. Ricardo Jorge (INSA) with the aim of:
- obtaining **genetic clusters at any partition level(s)** of a minimum spanning tree derived from cg/wgMLST allele matrix or of a newick tree (e.g. SNP-scaled tree)
- obtaining **summary reports with the statistics/trends** (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for the derived genetic clusters or for any other provided grouping variable (e.g. clade, lineage, ST, vaccination status, etc.)
- identifying regions of **cluster stability** (i.e. cg/wgMLST partition threshold ranges in which cluster composition is similar)


_Note: ReporTree relies on the usage of programs/modules of other developers. DO NOT FORGET TO ALSO CITE THEM!_

## Implementation
ReporTree is implemented in python 3.6 and comprises four modules available in standalone mode that are orchestrated by _reportree.py_:
- _partitioning_grapetree.py_ (only run when an allele matrix is provided)   
This script reconstructs the minimum spanning tree of a user-provided cg/wgMLST allele-matrix using a modified version of [GrapeTree](https://github.com/insapathogenomics/GrapeTree), and cuts this tree at all (or any user-specified) thresholds. The resulting genetic clusters are reported in a single matrix file (*_partitions.tsv).


- _partitioning_treecluster.py_ (only run when a tree is provided)    
This script takes advantage of [TreeCluster](https://github.com/niemasd/TreeCluster) to cut a provided tree for all the user-specified clustering methods and respective thresholds. It then reports all the resulting genetic clusters in a partitions matrix file (*_partitions.tsv).


- _comparing_partitions_v2.py_ (only run when the user requests “stability regions”)    
This is a modified version of comparing_partitions.py that analyzes the cluster congruence between subsequent partitions of a given clustering method (using the Adjusted Wallace coefficient ([Carriço et al. 2006](https://journals.asm.org/doi/10.1128/JCM.02536-05))), and identifies the partition threshold ranges for which cluster composition is similar (i.e. the “stability regions”). ReporTree then includes the minimum threshold of each of the “stability regions” in the summary report.


- _metadata_report.py_      
This script takes a metadata table as input and provides a separated summary report (*_summary.tsv) for each variable specified by the user. If a partitions file is also provided, this script includes the genetic clusters in a new metadata file (*_metadata_w_partitions.tsv) and creates a summary report of the metadata associated with each of the genetic clusters.


## Installation and dependencies

Dependencies:
- [TreeCluster](https://github.com/niemasd/TreeCluster) v1.0.3
- [A modified version of GrapeTree](https://github.com/insapathogenomics/GrapeTree)
- [Pandas](https://pandas.pydata.org)
- [Ete3](http://etetoolkit.org)

Installation:
```bash
conda create -n reportree -c etetoolkit -c anaconda -c bioconda ete3 scikit-learn pandas grapetree=2.1 treecluster=1.0.3
git clone ReporTree
cd ReporTree/scripts/
git clone https://github.com/insapathogenomics/GrapeTree
git clone https://github.com/insapathogenomics/ComparingPartitions
```
