# ReporTree

**Genomics-informed pathogen surveillance** strengthens public health decision-making, thus playing an important role in infectious diseases’ prevention and control. A pivotal outcome of genomics surveillance is the **identification of pathogen genetic clusters/lineages and their characterization in terms of geotemporal spread or linkage to clinical and demographic data**. This task usually relies on the visual exploration of (large) phylogenetic trees (eg. Minimum Spanning Trees (MST) for bacteria or rooted SNP-based trees for viruses). As this may be a non-trivial, non-reproducible and time consuming task, we aimed to develop a flexible pipeline that facilitates the detection of genetic clusters and their linkage to epidemiological data. 


_ReporTree can help you to:_      
- obtain **genetic clusters at any partition level(s)** of a newick tree (e.g. SNP-scaled tree) or an MST derived from cg/wgMLST allele matrix
- obtain **summary reports with the statistics/trends** (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for the derived genetic clusters or for any other provided grouping variable (e.g. clade, lineage, ST, vaccination status, etc.)
- identify regions of **cluster stability** (i.e. cg/wgMLST partition/threshold ranges in which cluster composition is similar)


_Note: ReporTree relies on the usage of programs/modules of other developers. DO NOT FORGET TO ALSO CITE THEM!_


## Implementation

ReporTree is implemented in python 3.6 and comprises four modules available in standalone mode that are orchestrated by _reportree.py_:
- _partitioning_grapetree.py_ (only run when an allele matrix is provided)   
This script reconstructs the minimum spanning tree of a user-provided cg/wgMLST allele-matrix using a modified version of [GrapeTree](https://github.com/insapathogenomics/GrapeTree), and cuts this tree at all (or any user-specified) thresholds. The resulting genetic clusters are reported in a single matrix file (partitions.tsv).


- _partitioning_treecluster.py_ (only run when a newick tree is provided)    
This script takes advantage of [TreeCluster](https://github.com/niemasd/TreeCluster) to cut a newick tree for all the user-specified clustering methods and respective thresholds. It then reports all the resulting genetic clusters in a partitions matrix file (partitions.tsv).


- _comparing_partitions_v2.py_ (only run when the user requests “stability regions”)    
This is a modified and automated version of [comparing_partitions.py](https://github.com/jacarrico/ComparingPartitions) that analyzes the cluster congruence between subsequent partitions of a given clustering method (using the Adjusted Wallace coefficient ([Carriço et al. 2006](https://journals.asm.org/doi/10.1128/JCM.02536-05))) and identifies the stability regions. ReporTree then includes the minimum threshold of each “stability region” in the summary report.


- _metadata_report.py_      
This script takes a metadata table as input (.tsv) and provides a separated summary report for each variable specified by the user. If a partitions matrix file is also provided, this script includes the genetic clusters in a new metadata file (metadata_w_partitions.tsv) and creates a summary report of the metadata associated with each of the genetic clusters.

![reporTree](https://user-images.githubusercontent.com/87025402/154258478-f64fae75-31c6-49a6-a6ed-54948f4acb8d.png)


## Input files
Metadata table in .tsv format (column 'sequence' is mandatory and columns should not have blank spaces)           
**AND (optionally)**        

Newick tree which will be used to obtain genetic clusters       
**OR**      
Allele matrix which will be used to obtain genetic clusters from a MST      
**OR**     
Partitions table (i.e. matrix with genetic clusters) in .tsv format (column 'sequence' is mandatory and columns should not have blank spaces)      


## Main output files
1. partitions.tsv - genetic clusters obtained for each user-selected partition threshold (only generated when a newick file or an allele matrix is provided)
2. metrics.tsv - metrics resulting from the cluster congruence analysis, with indication of the Adjusted Wallace and the Ajusted Rand coefficients for each comparison of subsequent partitions, and the Simpson's Index of Diversity for each partition.
3. stableRegions.tsv - partition ranges for which Adjusted Wallace coefficient is higher than the cut-off defined by the user
4. metadata_w_partitions.tsv - initial metadata information with additional columns comprising information on the genetic clusters at different partitions
5. partitions_summary.tsv - summary report with the statistics/trends (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for the derived genetic clusters present in partitions.tsv
6. variable_summary.tsv - summary report with the statistics/trends (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for any (and as many) grouping variable present in metadata_w_partitions.tsv (such as, clade, lineage, ST, vaccination status, etc.)


## Installation and dependencies

Dependencies:
- [TreeCluster](https://github.com/niemasd/TreeCluster) v1.0.3
- [A modified version of GrapeTree](https://github.com/insapathogenomics/GrapeTree)
- [Pandas](https://pandas.pydata.org)
- [Ete3](http://etetoolkit.org)

Installation:
```bash
conda create -n reportree -c etetoolkit -c anaconda -c bioconda ete3 scikit-learn pandas grapetree=2.1 treecluster=1.0.3 python=3.6
git clone https://github.com/insapathogenomics/ReporTree
cd ReporTree/scripts/
git clone https://github.com/insapathogenomics/GrapeTree
git clone https://github.com/insapathogenomics/ComparingPartitions
```

## Usage

```bash
ReporTree input/output file specifications

  -m METADATA, --metadata METADATA
                        [MANDATORY] Metadata file in .tsv format (column 'sequence' is mandatory)
  -t TREE, --tree TREE  [OPTIONAL] Input tree
  -a ALLELE_PROFILE, --allele-profile ALLELE_PROFILE
                        [OPTIONAL] Input allele profile matrix
  -p PARTITIONS, --partitions PARTITIONS
                        [OPTIONAL] Partitions file in .tsv format (column 'sequence' is mandatory) - 'partition' represents any variable that is not in the metadata
  -out OUTPUT, --output OUTPUT
                        [OPTIONAL] Tag for output file name (default = ReporTree)
  --list                [OPTIONAL] If after your command line you specify this option, ReporTree will list all the possible columns that you can use as input in '--columns_summary_report'. NOTE!! The
                        objective of this argument is to help you with the input of '--columns_summary_report'. So, it will not run reportree.py main functions!!

Partition minimum unit:
  Minimum unit/distance between partition thresholds

  -d DIST, --dist DIST  Distance unit by which partition thresholds will be multiplied (example: if -d 10 and -thr 5,8,10-30, the minimum spanning tree will be cut at 50,80,100,110,120,...,300. If -d
                        10 and --method-threshold avg_clade-2, the avg_clade threshold will be set at 20). This argument is particularly useful for non- SNP-scaled trees. Currently, the default is 1,
                        which is equivalent to 1 allele distance or 1 SNP distance. [1.0]

Partitioning with TreeCluster:
  Specifications to cut the tree with TreeCluster [only if a tree file is provided]

  --method-threshold METHOD_THRESHOLD
                        List of TreeCluster methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write the method name
                        (e.g. root_dist). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. root_dist-10). To get clustering at a specific range, indicate the range
                        with a hyphen (e.g. root_dist-2-10). Default: root_dist,avg_clade-1 (List of possible methods: avg_clade, leaf_dist_max, leaf_dist_min, length, length_clade, max, max_clade,
                        root_dist, single_linkage, single_linkage_cut, single_linkage_union) Warning!! So far, ReporTree was only tested with avg_clade and root_dist!
  --support SUPPORT     [OPTIONAL: see TreeCluster github for details] Branch support threshold
  --root-dist-by-node   [OPTIONAL] Set only if you WANT to cut the tree with root_dist method at each tree node distance to the root (similar to root_dist at all levels but just for informative
                        distances)

Partitioning with GrapeTree:
  Specifications to get and cut minimum spanning trees derived from cg/wgMLST allele data [only if an allele profile file is provided]

  --method GRAPETREE_METHOD
                        "MSTreeV2" [DEFAULT] Alternative:"MSTree"
  --missing HANDLER     ONLY FOR MSTree. 0: [DEFAULT] ignore missing data in pairwise comparison. 1: remove column with missing data. 2: treat missing data as an allele. 3: use absolute number of
                        allelic differences.
  --wgMLST              [EXPERIMENTAL: see GrapeTree github for details] a better support of wgMLST schemes
  --n_proc NUMBER_OF_PROCESSES
                        Number of CPU processes in parallel use. [5]
  -thr THRESHOLD, --threshold THRESHOLD
                        Partition threshold for clustering definition. Different thresholds can be comma-separated (e.g. 5,8,16). Ranges can be specified with a hyphen (e.g. 5,8,10-20). If this option
                        is not set, the script will perform clustering for all the values in the range 1 to max
  --subset              Reconstruct the minimum spanning tree using only the samples that correspond to the filters specified at the '--filter' argument

ReporTree metadata report:
  Specific parameters to report clustering/grouping information associated to metadata

  --columns_summary_report COLUMNS_SUMMARY_REPORT
                        Columns (i.e. variables of metadata) to get statistics for the derived genetic clusters or for other grouping variables defined in --metadata2report (comma-separated). If the
                        name of the column is provided, the different observations and the respective percentage are reported. If 'n_column' is specified, the number of the different observations is
                        reported. For example, if 'n_country' and 'country' are specified, the summary will report the number of countries and their distribution (percentage) per cluster/group.
                        Exception: if a 'date' column is in the metadata, it can also report first_seq_date, last_seq_date, timespan_days. Default =
                        n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days [the order of the list will be the order of the columns in the report]
  --partitions2report PARTITIONS2REPORT
                        Columns of the partitions table to include in a joint report (comma-separated). Other alternatives: 'all' == all partitions; 'stability_regions' == first partition of each
                        stability region as determined by comparing_partitions_v2.py. Warning!! 'stability_regions' can only be inferred when partitioning TreeCluster or GrapeTree is run for all
                        possible thresholds or when a similar partitions table is provided (i.e. sequential partitions obtained with the same clustering method) [all]
  --metadata2report METADATA2REPORT
                        Columns of the metadata table for which a separated summary report must be provided (comma-separated)
  -f FILTER_COLUMN, --filter FILTER_COLUMN
                        [OPTIONAL] Filter for metadata columns to select the samples to analyze. This must be specified within quotation marks in the following format 'column< >operation< >condition'
                        (e.g. 'country == Portugal'). When more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When
                        filters include more than one column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are
                        important in this argument, so, do not leave spaces before and after commas/semicolons.
  --frequency-matrix FREQUENCY_MATRIX
                        [OPTIONAL] Metadata column names for which a frequency matrix will be generated. This must be specified within quotation marks in the following format 'variable1,variable2'.
                        Variable1 is the variable for which frequencies will be calculated (e.g. for 'lineage,iso_week' the matrix reports the percentage of samples that correspond to each lineage per
                        iso_week). If you want more than one matrix you can separate the different requests with semicolon (e.g. 'lineage,iso_week;country,lineage'). If you want a higher detail in
                        your variable2 and decompose it into two columns you use a colon (e.g. lineage,country:iso_week will report the percentage of samples that correspond to each lineage per
                        iso_week in each country)
  --count-matrix COUNT_MATRIX
                        [OPTIONAL] Same as '--frequency-matrix' but outputs counts and not frequencies
  --mx-transpose        [OPTIONAL] Set ONLY if you want that the variable1 specified in '--frequency-matrix' or in '--count-matrix' corresponds to the matrix first column.
  --metadata-4-grapetree
                        Output an additional metadata file with the header ready for GrapeTree visualization. Set only if you WANT the file

Stability regions:
  Congruence analysis of cluster composition at all possible partitions to determine regions of cluster stability (automatically run if you set --partitions2report 'stability_regions'). WARNING! This option is planned to handle sequential partitions obtained with the same clustering method, such as a partitions table derived from cg/wgMLST data (from 1 to max allele threshold). Use it at your own risk, if you provide your own partitions table.

  -AdjW ADJUSTEDWALLACE, --AdjustedWallace ADJUSTEDWALLACE
                        Threshold of Adjusted Wallace score to consider an observation for method stability analysis [0.99]
  -n N_OBS, --n_obs N_OBS
                        Minimum number of sequencial observations that pass the Adjusted Wallace score to be considered a 'stability region' (i.e. a threshold range in which cluster composition is
                        similar) [5]
  -o ORDER, --order ORDER
                        [Set only if you provide your own partitions table] Partitions order in the partitions table (0: min -> max; 1: max -> min) [0]
```


#### Note on the '--method-threshold' argument

As above-mentioned, reporTree uses [TreeCluster](https://github.com/niemasd/TreeCluster) to obtain clustering information from a SNP-distance tree. To provide some flexibility, it uses the argument '--method-threshold' to run this program for all the combinations of method-threshold that the user needs:
- Clustering at all possible thresholds of a single method (from 1 SNP to the maximum distance of the tree) -> set only the method (e.g. root_dist)
- Clustering at a specific threshold for a single method -> set the method and use a hyphen to indicate the threshold (e.g. root_dist-10)
- Clustering at all possible thresholds of a range for a single method -> set the method and use a hyphen to indicate the threshold range (e.g. root_dist-10-30)
- Clustering using two or more methods and/or different thresholds -> set the different requirements separated by a comma (e.g. root_dist,avg_dist-10,avg_dist-20-30)


#### Note on the '--partitions2report' argument

This argument is used to select the columns of the partitions table that will be incorporated into the metadata table. If you use ReporTree to obtain the partitions table, the column names specification follows the same rules as the '--method-threshold' or the '--threshold' argument, depending on whether you provided a newick tree or an allele matrix, respectively.


#### Note on the columns for summary reports

This argument is used to select the columns that will be provided in the summary report.
-	To take the most profit of ReporTree, we recommend that you include the column 'date' in your metadata. This column must follow the format YYYY-MM-DD. If you only provide YYYY, it will assume YYYY-01-01!
-	If a 'date' column is provided in the metadata, ReporTree will determine and provide in the new metadata table the columns:
    - iso_year
    - iso_week_nr 
    - iso_week
-	While for nominal or categorical variables ReporTree can provide in the summary report the number of observations (‘n_column’) or the frequency of each observation, for the 'date' column this script can provide:
    - first_seq_date
    - last_seq_date
    - timespan_days
-	The columns of the summary reports are defined by the ‘--columns_summary_report’ argument. To know the columns that you can include based on your metadata table and the outputs of ReporTree, write the full command line and add the argument ‘--list’ in the end. This will give you a list of the possible columns that your summary reports can include, i.e. that you can request in ‘--columns_summary_report’.


## Examples

## Citation

If you run ReporTree, please do not forget to cite this page.

Also, ReporTree relies on the work of other developers. So, depending on the functionalities you use, there are other tools that you must cite:     
- Ete3: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4868116/pdf/msw046.pdf (all tools)     
- Grapetree: http://www.genome.org/cgi/doi/10.1101/gr.232397.117 (if you provided an allele matrix)      
- TreeCluster: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6705769/pdf/pone.0221068.pdf (if you provided a newick tree)      
- ComparingPartitions: https://journals.asm.org/doi/10.1128/jcm.02536-05?permanently=true (if you requested "stability_regions")      
