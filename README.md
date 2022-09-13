# ReporTree

<p align="center">
  <img width="300" height="180" src=https://user-images.githubusercontent.com/19263468/189874644-00ff1b8b-1d1d-4d69-80df-0bd5893a9df2.png>
</p>


**Genomics-informed pathogen surveillance** strengthens public health decision-making, thus playing an important role in infectious diseases’ prevention and control. A pivotal outcome of genomics surveillance is the **identification of pathogen genetic clusters/lineages and their characterization in terms of geotemporal spread or linkage to clinical and demographic data**. This task usually relies on the visual exploration of (large) phylogenetic trees (e.g. Minimum Spanning Trees (MST) for bacteria or rooted SNP-based trees for viruses). As this may be a non-trivial, non-reproducible and time consuming task, we developed **ReporTree, a flexible pipeline that facilitates the detection of genetic clusters and their linkage to epidemiological data**. 


_ReporTree can help you to:_      
- obtain **genetic clusters at any threshold level(s)** of a tree, SNP or cg/wgMLST allele matrix, VCF files, sequence alignment, or distance matrix
- obtain **summary reports with the statistics/trends** (e.g., timespan, location, cluster/group composition, age distribution etc.) for the derived genetic clusters or for any other provided grouping variable (e.g., clade, lineage, ST, vaccination status, etc.)
- obtain **count/frequency matrices** for the derived genetic clusters or for any other provided grouping variable
- identify **regions of cluster stability** (i.e., threshold ranges in which cluster composition is similar), a key step for nomenclature design


In summary, ReporTree facilitates and accelerates the production of surveillance-oriented reports, thus contributing to a sustainable and efficient public health genomics-informed pathogen surveillance.


_2022.09.13 - NEWS!! ReporTree can take a list 
_2022.06.25 - ReporTree can take VCF files as input!_


_Note: this tool relies on the usage of programs/modules of other developers. DO NOT FORGET TO ALSO CITE THEM!_


## Implementation

ReporTree is implemented in python 3.8 and comprises six main modules available in standalone mode that are orchestrated by _reportree.py_ (see details in [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/1.-Implementation#reportree-modular-organization)).


## Input files

Metadata table in .tsv format (column should not have blank spaces)           
**AND (optionally)**        

Newick tree which will be used to obtain genetic clusters       
**OR**      
Allele/SNP profile matrix which will be used to obtain genetic clusters from a MST      
**OR**  
Sequence alignment which will be converted into a profile and used to obtain genetic clusters  
**OR**  
VCF which will be converted into a profile and used to obtain genetic clusters    
**OR**     
Distance matrix which will be used to obtain genetic clusters       
**OR**  
Partitions table (i.e. matrix with genetic clusters) in .tsv format (columns should not have blank spaces)       


_In the following table we summarize the different options that ReporTree provides to determine genetic clusters, as well as the different types of file that each of them can take as input:_

<p align="center">
  <img width="800" alt="Captura de ecrã 2022-06-25, às 11 22 47" src="https://user-images.githubusercontent.com/19263468/175769475-5fc6efb3-266b-4302-b85d-0400d5671df0.png">
</p>


## Main output files

- metadata_w_partitions.tsv - initial metadata information with additional columns comprising information on the genetic clusters at different partitions

_TIP: Users can interactively visualize and explore the ReporTree derived clusters by uploading this metadata_w_partitions.tsv table together with either the original newick tree (e.g. rooted SNP-scaled tree) or the dendrogram resulting from hierarchical clustering at [auspice.us](https://auspice.us) or the MST resulting from GrapeTree at [GrapeTree](https://github.com/achtman-lab/GrapeTree). With these tools your dataset is visualised from the client-side in the browser._

- partitions_summary.tsv - summary report with the statistics/trends (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for the derived genetic clusters present in partitions.tsv (note: singletons are not reported in this file but indicated in metadata_w_partitions.tsv)
- variable_summary.tsv - summary report with the statistics/trends (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for any (and as many) grouping variable present in metadata_w_partitions.tsv (such as, clade, lineage, ST, vaccination status, etc.)
- partitions.tsv - genetic clusters obtained for each user-selected partition threshold
- freq_matrix.tsv - frequencies of grouping variable present in metadata_w_partitions.tsv (e.g. lineage, ST, etc.) across another grouping variable (e.g. iso_week, country, etc.)
- count_matrix.tsv - counts of a grouping variable present in metadata_w_partitions.tsv (e.g. lineage, ST, etc.) across another grouping variable (e.g. iso_week, country, etc.)
- metrics.tsv - metrics resulting from the cluster congruence analysis, with indication of the Adjusted Wallace and the Ajusted Rand coefficients for each comparison of subsequent partitions, and the Simpson's Index of Diversity for each partition.
- stableRegions.tsv - partition ranges for which Adjusted Wallace coefficient is higher than the cut-off defined by the user (useful to study cluster stability and infer possible nomenclature) 
- Newick file with the dendrogram resulting of the hierarchical clustering analysis or with the minimum spanning tree of GrapeTree


## Installation and dependencies

Dependencies:
- [TreeCluster](https://github.com/niemasd/TreeCluster) v1.0.3
- [A modified version of GrapeTree](https://github.com/insapathogenomics/GrapeTree)
- [Biopython](https://biopython.org)
- [Pandas](https://pandas.pydata.org)
- [Ete3](http://etetoolkit.org)
- [cgmlst-dists](https://github.com/tseemann/cgmlst-dists)
- [vcf2mst](https://github.com/vmixao/vcf2mst)

Installation:
```bash
conda create -n reportree -c anaconda -c bioconda -c etetoolkit python=3.8 biopython=1.77 pandas=1.1.3 numpy=1.19.2 grapetree=2.1 treecluster=1.0.3 ete3 scikit-learn cgmlst-dists
git clone https://github.com/insapathogenomics/ReporTree
cd ReporTree/scripts/
git clone https://github.com/insapathogenomics/GrapeTree
git clone https://github.com/insapathogenomics/ComparingPartitions
git clone https://github.com/vmixao/vcf2mst.git
```


## Usage

```bash
  -h, --help            show this help message and exit

ReporTree:
  ReporTree input/output file specifications

  -a ALLELE_PROFILE, --allele-profile ALLELE_PROFILE
                        [OPTIONAL] Input allele/SNP profile matrix (tsv format)
  -align ALIGNMENT, --alignment ALIGNMENT
                        [OPTIONAL] Input multiple sequence alignment (fasta format)
  -vcf VCF, --vcf VCF   [OPTIONAL] Single-column list of VCF files (txt format). This file must comprise the full PATH to each vcf file.
  -var VARIANTS, --variants VARIANTS
                        [OPTIONAL] Input table (tsv format) with sample name in the first column and a comma-separated list of variants in the second column with the following regular expression:
                        '\w(\d+)\w'
  -d_mx DISTANCE_MATRIX, --distance_matrix DISTANCE_MATRIX
                        [OPTIONAL] Input pairwise distance matrix (tsv format)
  -t TREE, --tree TREE  [OPTIONAL] Input tree (newick format)
  -p PARTITIONS, --partitions PARTITIONS
                        [OPTIONAL] Partitions file (tsv format) - 'partition' represents the threshold at which clustering information was obtained
  -m METADATA, --metadata METADATA
                        [MANDATORY] Metadata file (tsv format). To take the most profit of ReporTree functionalities, you must provide this file.
  -out OUTPUT, --output OUTPUT
                        [OPTIONAL] Tag for output file name (default = ReporTree)
  --list                [OPTIONAL] If after your command line you specify this option, ReporTree will list all the possible columns that you can use as input in '--columns_summary_report'. NOTE!! The
                        objective of this argument is to help you with the input of '--columns_summary_report'. So, it will not run reportree.py main functions!!

Analysis details:
  Analysis details

  --analysis ANALYSIS   Type of clustering analysis (options: grapetree, HC, treecluster). If you provide a tree, genetic clusters will always be obtained with treecluster. If you provide a distance
                        matrix, genetic clusters will always be obtained with HC. If you provide any other input, it is MANDATORY to specify this argument.
  --subset              Obtain genetic clusters using only the samples that correspond to the filters specified in the '--filter' argument (only valid for analysis == grapetree or HC)
  -d DIST, --dist DIST  Distance unit by which partition thresholds will be multiplied (example: if -d 10 and -thr 5,8,10-30, the minimum spanning tree will be cut at 50,80,100,110,120,...,300. If -d 10
                        and --method-threshold avg_clade-2, the avg_clade threshold will be set at 20). This argument is particularly useful for non-SNP-scaled trees. Currently, the default is 1, which
                        is equivalent to 1 allele distance or 1 SNP distance. [1.0]
  --site-inclusion N_CONTENT
                        [OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum proportion of samples per site without missing data (e.g. '--site-inclusion 1.0' will only
                        keep loci/positions without missing data, i.e. a core alignment/profile; '--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment
                        positions/loci (i.e. columns)! [default: 0.0 - content of missing data is not considered during matrix/alignment cleaning].

Processing profiles:
  Processing profiles

  --loci-called LOCI_CALLED
                        [OPTIONAL] Minimum proportion of loci/positions called for SNP/allele matrices (e.g. '--loci-called 0.95' will only keep in the profile matrix samples with > 95% of
                        alleles/positions, i.e. <= 5% missing data). Applied after '--site-inclusion' argument! Code for missing data: 0.

Alignment processing:
  Alignment processing

  --sample-ATCG-content ATCG_CONTENT
                        [OPTIONAL] Minimum proportion (0 to 1) of ATCG in informative sites of the alignment per sample (e.g. '--sample-ATCG-content 1.0' will only keep samples without N's or any non-
                        ATCG code in informative sites)
  --remove-reference    Set only if you want to remove the reference sequence of the alignment (reference name must be provided with the argument '--reference').
  --use-reference-coords
                        Set only if you want that column names in the final alignment matrix represent the reference coordinates (reference name must be provided with the argument '--reference') Note:
                        Depending on the alignment size, this argument can make alignment processing very slow!
  -r REFERENCE, --reference REFERENCE
                        [OPTIONAL] Name of reference sequence. Required if '--remove-reference' and/or '--use-reference-coords' specified.

Partitioning with GrapeTree:
  Specifications to get and cut minimum spanning trees

  --method GRAPETREE_METHOD
                        "MSTreeV2" [DEFAULT] Alternative:"MSTree (goeBURST)"
  --missing HANDLER     ONLY FOR MSTree. 0: [DEFAULT] ignore missing data in pairwise comparison. 1: remove column with missing data. 2: treat missing data as an allele. 3: use absolute number of
                        allelic differences.
  --n_proc NUMBER_OF_PROCESSES
                        Number of CPU processes in parallel use. [5]
  -thr THRESHOLD, --threshold THRESHOLD
                        Partition threshold for clustering definition. Different thresholds can be comma-separated (e.g. 5,8,16). Ranges can be specified with a hyphen (e.g. 5,8,10-20). If this option
                        is not set, the script will perform clustering for all the values in the range 0 to max. Note: Threshold values are inclusive, i.e. '-thr 7' will consider samples with <= 7
                        differences as belonging to the same cluster!
  --matrix-4-grapetree  Output an additional allele profile matrix with the header ready for GrapeTree visualization. Set only if you WANT the file!
  --wgMLST              [EXPERIMENTAL] a better support of wgMLST schemes (check GrapeTree github for details).

Partitioning with HC:
  Specifications to genetic clusters with hierarchical clustering

  --HC-threshold HCMETHOD_THRESHOLD
                        List of HC methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write the method name (e.g.
                        single). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. single-10). To get clustering at a specific range, indicate the range with a hyphen
                        (e.g. single-2-10). Note: Threshold values are inclusive, i.e. '--HC-threshold single-7' will consider samples with <= 7 differences as belonging to the same cluster! Default:
                        single (Possible methods: single, complete, average, weighted, centroid, median, ward)

Partitioning with TreeCluster:
  Specifications to cut the tree with TreeCluster

  --method-threshold METHOD_THRESHOLD
                        List of TreeCluster methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write the method name
                        (e.g. root_dist). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. root_dist-10). To get clustering at a specific range, indicate the range
                        with a hyphen (e.g. root_dist-2-10). Default: root_dist,avg_clade-1 (List of possible methods: avg_clade, leaf_dist_max, leaf_dist_min, length, length_clade, max, max_clade,
                        root_dist, single_linkage, single_linkage_cut, single_linkage_union) Warning!! So far, ReporTree was only tested with avg_clade and root_dist!
  --support SUPPORT     [OPTIONAL: see TreeCluster github for details] Branch support threshold
  --root-dist-by-node   [OPTIONAL] Set only if you WANT to cut the tree with root_dist method at each tree node distance to the root (similar to root_dist at all levels but just for informative
                        distances)

ReporTree metadata report:
  Specific parameters to report clustering/grouping information associated to metadata

  --columns_summary_report COLUMNS_SUMMARY_REPORT
                        Columns (i.e. variables of metadata) to get statistics for the derived genetic clusters or for other grouping variables defined in --metadata2report (comma-separated). If the
                        name of the column is provided, the different observations and the respective percentage are reported. If 'n_column' is specified, the number of the different observations is
                        reported. For example, if 'n_country' and 'country' are specified, the summary will report the number of countries and their distribution (percentage) per cluster/group.
                        Exception: if a 'date' column is in the metadata, it can also report first_seq_date, last_seq_date, timespan_days. Check '--list' argument for some help. Default =
                        n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days [the order of the list will be the order of the columns in the report]
  --partitions2report PARTITIONS2REPORT
                        Columns of the partitions table to include in a joint report (comma-separated). Other alternatives: 'all' == all partitions; 'stability_regions' == first partition of each
                        stability region as determined by comparing_partitions_v2.py. Note: 'stability_regions' can only be inferred when partitioning TreeCluster or GrapeTree is run for all possible
                        thresholds or when a similar partitions table is provided (i.e. sequential partitions obtained with the same clustering method) [all]. Check '--list' argument for some help
  --metadata2report METADATA2REPORT
                        Columns of the metadata table for which a separated summary report must be provided (comma-separated)
  -f FILTER_COLUMN, --filter FILTER_COLUMN
                        [OPTIONAL] Filter for metadata columns to select the samples to analyze. This must be specified within quotation marks in the following format 'column< >operation< >condition'
                        (e.g. 'country == Portugal'). When more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When
                        filters include more than one column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are
                        important in this argument, so, do not leave spaces before and after commas/semicolons.
  --sample_of_interest SAMPLE_OF_INTEREST
                        Comma-separated list of samples of interest for which summary reports will be created. If none provided, only the summary reports comprising all samples will be generated.
  --frequency-matrix FREQUENCY_MATRIX
                        [OPTIONAL] Metadata column names for which a frequency matrix will be generated. This must be specified within quotation marks in the following format 'variable1,variable2'.
                        Variable1 is the variable for which frequencies will be calculated (e.g. for 'lineage,iso_week' the matrix reports the percentage of samples that correspond to each lineage per
                        iso_week). If you want more than one matrix you can separate the different requests with semicolon (e.g. 'lineage,iso_week;country,lineage'). If you want a higher detail in your
                        variable2 and decompose it into two columns you use a colon (e.g. lineage,country:iso_week will report the percentage of samples that correspond to each lineage per iso_week in
                        each country)
  --count-matrix COUNT_MATRIX
                        [OPTIONAL] Same as '--frequency-matrix' but outputs counts and not frequencies
  --mx-transpose        [OPTIONAL] Set ONLY if you want that the variable1 specified in '--frequency-matrix' or in '--count-matrix' corresponds to the matrix first column.

Stability regions:
  Congruence analysis of cluster composition at all possible partitions to determine regions of cluster stability (automatically run if you set --partitions2report 'stability_regions'). WARNING! This option is planned to handle sequential partitions obtained with the same clustering method, such as a partitions table derived from cg/wgMLST data (from 1 to max allele threshold). Use it at your own risk, if you provide your own partitions table.

  -AdjW ADJUSTEDWALLACE, --AdjustedWallace ADJUSTEDWALLACE
                        Threshold of Adjusted Wallace score to consider an observation for method stability analysis [0.99]
  -n N_OBS, --n_obs N_OBS
                        Minimum number of sequencial observations that pass the Adjusted Wallace score to be considered a 'stability region' (i.e. a threshold range in which cluster composition is
                        similar) [5]
  -o ORDER, --order ORDER
                        [Set only if you provide your own partitions table] Partitions order in the partitions table (0: min -> max; 1: max -> min) [0]
  --keep-redundants     Set ONLY if you want to keep all samples of each cluster of the most discriminatory partition (by default redundant samples are removed to avoid the influence of cluster size)
```


#### Note on the '--method-threshold' argument

ReporTree uses [TreeCluster](https://github.com/niemasd/TreeCluster) to obtain clustering information from a SNP-distance tree. To provide some flexibility, it uses the argument '--method-threshold' to run this program for all the combinations of method-threshold that the user needs:
- Clustering at all possible thresholds of a single method (from 1 SNP to the maximum distance of the tree) -> set only the method (e.g. root_dist)
- Clustering at a specific threshold for a single method -> set the method and use a hyphen to indicate the threshold (e.g. root_dist-10)
- Clustering at all possible thresholds of a range for a single method -> set the method and use a hyphen to indicate the threshold range (e.g. root_dist-10-30)
- Clustering using two or more methods and/or different thresholds -> set the different requirements separated by a comma (e.g. root_dist,avg_dist-10,avg_dist-20-30)


#### Note on the '--HC-threshold' argument
'--HC-threshold' is the equivalent to '--method-threshold' for the hierarchical clustering analysis. As such, the input of this argument follows the same structure as mentioned in the previous paragraph. Possible methods: single, complete, average, weighted, centroid, median, ward.


#### Note on the '--partitions2report' argument

This argument is used to select the columns of the partitions table that will be incorporated into the metadata table. If you use ReporTree to obtain the partitions table, the column names specification follows the same rules as the '--method-threshold' or the '--threshold' argument, depending on whether you provided a newick tree or an allele matrix, respectively.


#### Note on the columns for summary reports

This argument is used to select the columns that will be provided in the summary report.
-	To take the most profit of ReporTree, we recommend that you include the column 'date' in your metadata. This column must follow the format YYYY-MM-DD. If you only provide YYYY or YYYY-MM, it will assume YYYY-01-01!
-	If a 'date' column is provided in the metadata, ReporTree will determine and provide in the new metadata table the columns:
    - iso_year
    - iso_week_nr 
    - iso_week
-	While for nominal or categorical variables ReporTree can provide in the summary report the number of observations (‘n_column’) or the frequency of each observation, for the 'date' column this script can provide:
    - first_seq_date
    - last_seq_date
    - timespan_days
-	The columns of the summary reports are defined by the ‘--columns_summary_report’ argument. To know the columns that you can include based on your metadata table and the outputs of ReporTree, write the full command line and add the argument ‘--list’ in the end. This will give you a list of the possible columns that your summary reports can include, i.e. that you can request in ‘--columns_summary_report’.


#### Note on the '--frequency-matrix' and '--count-matrix' arguments

These arguments take as input two variables separated by a comma (variable1,variable2). Frequencies will be calculated for variable1 and by default the different observations of this variable (e.g. different lineages if variable 1 == 'lineage') will correspond to different columns, while the observations of variable 2 will correspond to different rows. To transpose this you can use '--mx-transpose'.

_TIP: If you want, you can split variable 2 in up to two variables. To this end you can indicate them separated by colon (e.g. lineage,country:iso_week)_


### How to run ReporTree with a newick tree as input?

```bash
reportree.py -m metadata.tsv -t tree.nw -out output -d dist --method-threshold method1,method2-threshold --root-dist-by-node --columns_summary_report columns,summary,report --partitions2report all --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2
```


### How to run ReporTree with an allele matrix as input and obtain the genetic clusters from a MST?

```bash
reportree.py -m metadata.tsv -a allele_matrix.tsv -out output -d dist --method MSTreeV2 -thr 5,8,15 --matrix-4-grapetree --columns_summary_report columns,summary,report --partitions2report all --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2 -AdjW adjusted_wallace -n n --analysis grapetree
```

### How to run ReporTree with an alignment as input and obtain the genetic clusters with hierarchical clustering?

```bash
reportree.py -m metadata.tsv -align alignment.fasta -out output -d dist --HC-threshold single --columns_summary_report columns,summary,report --partitions2report all --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2 -AdjW adjusted_wallace -n n --analysis HC
```


## Examples

### Routine surveillance - viral pathogen (e.g. SARS-CoV-2) - [click here](https://github.com/insapathogenomics/ReporTree/wiki/4.-Examples#routine-surveillance---viral-pathogen-eg-sars-cov-2)

ReporTree is currently applied to generate weekly reports about SARS-CoV-2 variant circulation in Portugal. In [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/4.-Examples#routine-surveillance---viral-pathogen-eg-sars-cov-2), we give some examples on how to rapidly generate key surveillance metrics taking as input metadata tables (tsv format) and rooted divergence (SNP) trees (newick format) provided for download in regular Nextstrain (auspice) builds, such as those maintained by the National Institute of Health Dr. Ricardo Jorge, Portugal (INSA) at https://insaflu.insa.pt/covid19/.


### Outbreak detection - bacterial foodborne pathogen (e.g. _Listeria monocytogenes_) - [click here](https://github.com/insapathogenomics/ReporTree/wiki/4.-Examples#outbreak-detection---bacterial-foodborne-pathogen-eg-listeria-monocytogenes)

ReporTree can facilitate the routine surveillance and outbreak investigation of bacterial pathogens, such as foodborne pathogens. In [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/4.-Examples#outbreak-detection---bacterial-foodborne-pathogen-eg-listeria-monocytogenes), we provide a simple example of the usage of ReporTree to rapidly identify and characterize potential Listeriosis outbreaks. With a single command, ReporTree builds a MST from cgMLST data and **automatically extracts genetic clusters at three high resolution levels (<5, <8, <15 allelic differences)**, and provides comprehensive reports about the sample collection (e.g. ST sequence count/frequency per year, etc).  


### Large-scale genetic clustering and linkage to antibiotic resistance data (e.g. _Neisseria gonorrhoeae_) - [click here](https://github.com/insapathogenomics/ReporTree/wiki/4.-Examples#large-scale-genetic-clustering-and-linkage-to-antibiotic-resistance-data-eg-neisseria-gonorrhoeae)

ReporTree can enhance genomics surveillance and quickly identify/characterize genetic clusters from large datasets. In [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/4.-Examples#large-scale-genetic-clustering-and-linkage-to-antibiotic-resistance-data-eg-neisseria-gonorrhoeae), with a single command line, we reproduce part of the extensive genomics analysis performed by [Pinto et al., 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208699/pdf/mgen-7-481.pdf) over 3,791 _N. gonorrhoeae_ genomes from isolates collected across Europe.


## Citation

If you run ReporTree, please cite our preprint publication:    
[Mixão V, Pinto M, Gomes JP, Borges V (2022) ReporTree: a surveillance-oriented tool to strengthen the linkage between pathogen genetic clusters and epidemiological data. _Research Square_. doi: 10.21203/rs.3.rs-1404655/v1](https://www.researchsquare.com/article/rs-1404655/v1)

Also, ReporTree relies on the work of other developers. So, depending on the functionalities you use, there are other tools that you must cite:
1. Grapetree: http://www.genome.org/cgi/doi/10.1101/gr.232397.117 (if you requested a grapetree analysis)
2. TreeCluster: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6705769/pdf/pone.0221068.pdf (if you provided a newick tree)
3. vcf2mst: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-08112-0 (if you provided a vcf or a list of variants)
4. ComparingPartitions: https://journals.asm.org/doi/10.1128/jcm.02536-05?permanently=true (if you requested "stability_regions")
5. Adjusted Wallace and cluster stability: https://www.biorxiv.org/content/10.1101/299347v1 (if you requested "stability_regions")
6. Ete3: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4868116/pdf/msw046.pdf (if you provided a newick tree)     
7. cgmlst-dists: https://github.com/tseemann/cgmlst-dists (if you requested a HC analysis and did not provide a distance matrix)


## Funding

This work was supported by funding from the European Union’s Horizon 2020 Research and Innovation programme under grant agreement No 773830: One Health European Joint Programme.
