# ReporTree

![Test Pass:](https://github.com/insapathogenomics/ReporTree/actions/workflows/python-app.yml/badge.svg)

<p align="center">
  <img width="300" height="180" src=https://user-images.githubusercontent.com/19263468/189874644-00ff1b8b-1d1d-4d69-80df-0bd5893a9df2.png>
</p>


**Genomics-informed pathogen surveillance** strengthens public health decision-making, thus playing an important role in infectious diseases’ prevention and control. A pivotal outcome of genomics surveillance is the **identification of pathogen genetic clusters/lineages and their characterization in terms of geotemporal spread or linkage to clinical and demographic data**. This task usually relies on the visual exploration of (large) phylogenetic trees (e.g. Minimum Spanning Trees (MST) for bacteria or rooted SNP-based trees for viruses). As this may be a non-trivial, non-reproducible and time consuming task, we developed **ReporTree, a flexible pipeline that facilitates the detection of genetic clusters and their linkage to epidemiological data**. 


_ReporTree can help you to:_      
- obtain **genetic clusters at any threshold level(s)** of a tree, SNP or cg/wgMLST allele matrix, VCF files, sequence alignment, or distance matrix
- obtain **summary reports with the statistics/trends** (e.g., timespan, location, cluster/group composition, age distribution etc.) for the derived genetic clusters or for any other provided grouping variable (e.g., clade, lineage, ST, vaccination status, etc.)
- obtain **count/frequency matrices** for the derived genetic clusters or for any other provided grouping variable
- identify the phylogenetic context of **samples of interest** through a zoom-in on their clusters and/or through an in-depth analysis with the closest related samples
- maintain **cluster nomenclature** between runs and generate **hierarchical codes** at your levels of interest  
- identify **regions of cluster stability** (i.e., threshold ranges in which cluster composition is similar), a key step for pathogen-specific nomenclature design


In summary, ReporTree facilitates and accelerates the production of surveillance-oriented reports, thus contributing to a sustainable and efficient public health genomics-informed pathogen surveillance.

_Note: this tool relies on the usage of programs/modules of other developers. DO NOT FORGET TO ALSO [CITE](https://github.com/insapathogenomics/ReporTree/edit/main/README.md#citation) THEM!_

## News!
#### 2024.04.02 - ReporTree v2.4.1
We release a new version of ReporTree that is compatible with [SPREAD](https://github.com/genpat-it/spread), an extended version of GrapeTree. With this new version of ReporTree:
1. A new column named "category" is added to the _metadata_w_patitions.tsv_ in order to tag the samples of interest
2. A new optional argument ('--unzip') can be povided in order to output the results of '--zoom-cluster-of-interest' and '--subtree-of-interest' in unzipped format
3. Whenever '--zoom-cluster-of-interest' or '--subtree-of-interest' arguments are used, a file *zooms.txt* is provided for automated visualization of the zoom-in trees/subtrees with SPREAD


## Implementation

ReporTree is implemented in python 3.8 and comprises six main scripts available in standalone mode that are orchestrated by _reportree.py_ (see details in [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/1.-Implementation#reportree-modular-organization)).


## Input files

Metadata table in .tsv format (please avoid headers with blank spaces)           
**AND (optionally)**        

Newick tree which will be used to obtain genetic clusters       
**OR**      
Allele/SNP profile matrix which will be used to obtain genetic clusters from a MST      
**OR**  
Sequence alignment which will be converted into a profile and used to obtain genetic clusters  
**OR**  
VCF which will be converted into a profile and used to obtain genetic clusters    
**OR**   
TSV with a list of mutations per samples which will be converted into a profile and used to obtain genetic clusters    
**OR**     
Distance matrix which will be used to obtain genetic clusters       
**OR**  
Partitions table (i.e. matrix with genetic clusters) in .tsv format (columns should not have blank spaces)       

_**Nomenclature only:** If you want to maintain cluster names between ReporTree runs, you should also provide the partitions table of the previous run._

In the following table we summarize the different options that ReporTree provides to determine genetic clusters, as well as the different types of file that each of them can take as input:

<p align="center">
<img width="1038" alt="Captura de ecrã 2023-03-28, às 18 26 09" src="https://user-images.githubusercontent.com/19263468/228320038-424ca291-d1be-4b32-9222-40f8758547d4.png">
</p>


## Main output files

- metadata_w_partitions.tsv - initial metadata information with additional columns comprising information on the genetic clusters at different partitions and nomenclature code

_TIP: Users can interactively visualize and explore the ReporTree derived clusters by uploading this metadata_w_partitions.tsv table together with either the original newick tree (e.g. rooted SNP-scaled tree) or the dendrogram resulting from hierarchical clustering at [auspice.us](https://auspice.us) or the MST resulting from GrapeTree at [GrapeTree](https://github.com/achtman-lab/GrapeTree) or [SPREAD](https://github.com/genpat-it/spread). With these tools your dataset is visualised from the client-side in the browser._

- partitions_summary.tsv - summary report with the statistics/trends (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for the derived genetic clusters present in partitions.tsv (note: singletons are not reported in this file but indicated in metadata_w_partitions.tsv)
- SAMPLES_OF_INTEREST_partitions_summary.tsv - similar to partitions_summary.tsv but exclusively for the samples of interest     
- variable_summary.tsv - summary report with the statistics/trends (e.g. timespan, location range, cluster/group size and composition, age distribution etc.) for any (and as many) grouping variable present in metadata_w_partitions.tsv (such as, clade, lineage, ST, vaccination status, etc.)
- partitions.tsv - genetic clusters obtained for each user-selected partition threshold
- nomenclature_changes.tsv - record of all cluster alterations in terms of nomenclature and composition
- freq_matrix.tsv - frequencies of grouping variable present in metadata_w_partitions.tsv (e.g. lineage, ST, etc.) across another grouping variable (e.g. iso_week, country, etc.)
- count_matrix.tsv - counts of a grouping variable present in metadata_w_partitions.tsv (e.g. lineage, ST, etc.) across another grouping variable (e.g. iso_week, country, etc.)
- metrics.tsv - metrics resulting from the cluster congruence analysis, with indication of the Adjusted Wallace and the Ajusted Rand coefficients for each comparison of subsequent partitions, and the Simpson's Index of Diversity for each partition.
- stableRegions.tsv - partition ranges for which Adjusted Wallace coefficient is higher than the cut-off defined by the user (useful to study cluster stability and infer possible nomenclature) 
- Newick file with the dendrogram resulting of the hierarchical clustering analysis or with the minimum spanning tree of GrapeTree
- .zip - compressed folders with the output files of a high resolution analysis of the clusters with samples of interest and of the N closest samples to the samples of interest
- zooms.txt - list of output folders containing the outputs of '--zoom-cluster-of-interest' and/or '--subtree-of-interest'. This file is required for visualization of the zoom-in trees/subtrees with SPREAD.

The following figure presents the different input/output workflow possibilities provided by ReporTree:

<p align="center">
<img width="700" alt="Figure1_v3" src="https://github.com/insapathogenomics/ReporTree/assets/19263468/9047702b-5ad9-41d5-9a9e-a00e1d3add79">
</p>

_Note: This is the [Figure 1](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1/figures/1) of [ReporTree publication](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1). For details please check the respective caption._


## Nomenclature

ReporTree cluster names follow a regular expression:
- cluster_N - for each genetic cluster (i.e. >= 2 samples) found at each threshold of the current analysis (e.g. cluster_1)
- singleton_N - for each sample that does not belong to a cluster at a given threshold of the current analysis (e.g. singleton_1)

To facilitate routine surveillance and cluster monitoring over time, the user can provide the 'partitions.tsv' from the previous run in the '--nomenclature-file' argument, and ReporTree will use this information to (re)name the clusters in the current run. Below, we show a summary of the behavior of the “Cluster Nomenclature System” in some of the most common situations in a routine surveillance scenario:

![Nomenclature_2](https://github.com/insapathogenomics/ReporTree/assets/19263468/eea36fe7-113e-47cb-9eba-912e0d36f57a)

_Details about these and other situations are presented in the following table:_

<img width="1111" alt="Captura de ecrã 2023-04-06, às 14 42 41" src="https://user-images.githubusercontent.com/19263468/230396669-94077006-3593-4a78-bc6d-dd1fe1598092.png">


Of note, to increase the flexibility of the nomenclature system, ReporTree also allows the users to change the regular expression for cluster nomenclature (i.e., starting with “cluster_” or “singleton_”) by other nomenclature of interest (e.g., other official codes for outbreaks, genogroups, etc.), which will be kept afterwards. If the cluster/singleton names do not follow ReporTree's regular expression, the rules mentioned in the above table will be applied, except for the particular situation in which the name a former singleton does not follow ReporTree's regular expression. In this case, if a sigleton named as 'mycode' integrates a cluster only with new samples, this new cluster will be named 'mycode'.


## Installation and dependencies

### Dependencies:
- [TreeCluster](https://github.com/niemasd/TreeCluster) v1.0.3
- [A modified version of GrapeTree](https://github.com/insapathogenomics/GrapeTree)
- [Biopython](https://biopython.org)
- [Pandas](https://pandas.pydata.org)
- [Ete3](http://etetoolkit.org)
- [cgmlst-dists](https://github.com/genpat-it/cgmlst-dists)
- [vcf2mst](https://github.com/genpat-it/vcf2mst)
- [snp-sites](https://github.com/sanger-pathogens/snp-sites)

### Installation with conda
```bash
git clone https://github.com/insapathogenomics/ReporTree
cd ReporTree/scripts
git clone https://github.com/insapathogenomics/GrapeTree.git
git clone https://github.com/insapathogenomics/ComparingPartitions.git
git clone https://github.com/insapathogenomics/vcf2mst.git
cd ..
conda env create --name reportree --file=reportree_env.yml
```

_Note: If you are using Mac OS X, please use reportree_env_osx.yml to create the conda environment!_

  
Activate conda environment
```bash
conda activate reportree
```

Run pytest to check that your installation was well succeeded
```bash
pytest
```

Run ReporTree:
```bash
python reportree.py -h
```

### Installation with Docker
```bash
docker pull insapathogenomics/reportree:v2.4.1
```

Run ReporTree:
```bash
docker run insapathogenomics/reportree:v2.4.1 reportree.py -h
```

## Usage

```bash
optional arguments:
  -h, --help            show this help message and exit

Version:
  ReporTree version

  -v, --version         Print version and exit

ReporTree:
  ReporTree input/output specifications

  -a ALLELE_PROFILE, --allele-profile ALLELE_PROFILE
                        Input allele/SNP profile matrix (tsv format)
  -align ALIGNMENT, --alignment ALIGNMENT
                        Input multiple sequence alignment (fasta format)
  -d_mx DISTANCE_MATRIX, --distance_matrix DISTANCE_MATRIX
                        Input pairwise distance matrix (tsv format)
  -t TREE, --tree TREE  Input tree (newick format)
  -p PARTITIONS, --partitions PARTITIONS
                        Partitions file (tsv format) - 'partition' represents the threshold at which clustering information was obtained
  -m METADATA, --metadata METADATA
                        [MANDATORY] Metadata file (tsv format). To take the most profit of ReporTree functionalities, you must provide this file.
  -vcf VCF, --vcf VCF   Single-column list of VCF files (txt format). This file must comprise the full PATH to each vcf file.
  -var VARIANTS, --variants VARIANTS
                        Input table (tsv format) with sample name in the first column and a comma-separated list of variants in the second column with the following regular expression: '\w(\d+)\w'
  -out OUTPUT, --output OUTPUT
                        [OPTIONAL] Tag for output file name (default = ReporTree)
  --list                [OPTIONAL] If after your command line you specify this option, ReporTree will list all the possible columns that you can use as input in '--columns_summary_report'. To obtain
                        information about the partition name for other arguments('--frequency-matrix' and/or '--count-matrix'), please also indicate the type of analysis. NOTE!! The objective of this
                        argument is to help you with the input of some other arguments. So, it will not run reportree.py main functions!!

Analysis details:
  Analysis details

  --analysis ANALYSIS   Type of clustering analysis (options: grapetree, HC, treecluster). If you provide a tree, genetic clusters will always be obtained with treecluster. If you provide a distance
                        matrix, genetic clusters will always be obtained with HC. If you provide any other input, it is MANDATORY to specify this argument.
  --subset              [OPTIONAL] Obtain genetic clusters using only the samples that correspond to the filters specified in the '--filter' argument (only valid for analysis == grapetree or HC)
  -d DIST, --dist DIST  [OPTIONAL] Distance unit by which partition thresholds will be multiplied (example: if -d 10 and -thr 5,8,10-30, the minimum spanning tree will be cut at
                        50,80,100,110,120,...,300. If -d 10 and --method-threshold avg_clade-2, the avg_clade threshold will be set at 20). This argument is particularly useful for non-SNP-scaled trees.
                        Currently, the default is 1, which is equivalent to 1 allele distance or 1 SNP distance. [1.0]

Cleaning missing data:
  Remove loci/positions and samples based on missing content

  --missing-code MISSING_CODE
                        [OPTIONAL] Code representing missing data. If different from '0' or 'N', please try to avoid a IUPAC character (even in lower-case). [default: 'N' when '-align' provided and '0'
                        for other inputs]
  --site-inclusion N_CONTENT
                        [OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum proportion of samples per site without missing data (e.g. '--site-inclusion 1.0' will only
                        keep loci/positions without missing data, i.e. a core alignment/profile; '--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment
                        loci/positions (i.e. columns)! [default: 0.0].
  --loci-called LOCI_CALLED
                        [OPTIONAL - only works for matrices] Minimum proportion of loci/positions called for SNP/allele matrices (e.g. '--loci-called 0.95' will only keep in the profile matrix samples
                        with > 95% of alleles/positions, i.e. <= 5% missing data). Applied after '--site-inclusion' argument! [default: 0.0]
  --sample-ATCG-content ATCG_CONTENT
                        [OPTIONAL - only works for alignment] Minimum proportion (0 to 1) of ATCG in informative sites of the alignment per sample (e.g. '--sample-ATCG-content 1.0' will only keep
                        samples without Ns or any non-ATCG code in informative sites) [default: 0 - keep all samples]

Alignment processing:
  Alignment processing

  --remove-reference    Set only if you want to remove the reference sequence of the alignment (reference name must be provided with the argument '--reference').
  --use-reference-coords
                        Set only if you want that column names in the final alignment matrix represent the reference coordinates (reference name must be provided with the argument '--reference') Note:
                        Depending on the alignment size, this argument can make alignment processing very slow!
  -r REFERENCE, --reference REFERENCE
                        [OPTIONAL] Name of reference sequence. Required if '--remove-reference' and/or '--use-reference-coords' specified.
  --get-position-correspondence POS_CORR
                        [OPTIONAL] Request a .tsv with position correspondence between any sequences of your alignment. These should be indicated separated by a comma (e.g. seqA,seqB). To get the
                        position coordinates of all sequences just write 'all'.
  --position-list POS_LIST
                        [OPTIONAL] .tsv file with the positions of interest to be reported when '--get-position-correspondence' is requested. Each column should correspond to the positions of a sequence
                        and the sequence name should be indicated in the header. If this file is not provided, all positions of the alignment will be reported.

Partitioning with GrapeTree:
  Specifications to get and cut minimum spanning trees

  --method GRAPETREE_METHOD
                        "MSTreeV2" [DEFAULT] Alternative:"MSTree (goeBURST)"
  --missing HANDLER     ONLY FOR MSTree. 0: [DEFAULT] ignore missing data in pairwise comparison. 1: remove column with missing data. 2: treat missing data as an allele. 3: use absolute number of
                        allelic differences.
  --n_proc NUMBER_OF_PROCESSES
                        Number of CPU processes in parallel use. [5]
  -thr THRESHOLD, --threshold THRESHOLD
                        [OPTIONAL] Partition threshold for clustering definition (integer). Different thresholds can be comma-separated (e.g. 5,8,16). Ranges can be specified with a hyphen separating
                        minimum and maximum (e.g. 5,8,10-20). If this option is not set, the script will perform clustering for all the values in the range 0 to max. If you prefer to exclusively use the
                        '-pct_thr' argument, please set '-thr none'. Note: Threshold values are inclusive, i.e. '-thr 7' will consider samples with <= 7 differences as belonging to the same cluster!
  -pct_thr PCT_THRESHOLD, --pct_threshold PCT_THRESHOLD
                        [OPTIONAL] Similar to 'thr' but values are indicated as the proportion of differences to the final allelic schema size or number of informative positions, e.g. '-pct_thr 0.005'
                        corresponds to a threshold of 5 allelic/SNP differences in a matrix with 1000 loci/sites under analysis). Different values can be comma-separated (e.g. 0.005,0.01,0.1). Ranges
                        CANNOT be specified. This option is particularly useful for dynamic wgMLST analysis for which the size of the schema under analysis is contigent on dataset diversity. Note: This
                        argument can be specified even if you used the '-thr' argument.
  --matrix-4-grapetree  Output an additional allele profile matrix with the header ready for GrapeTree visualization. Set only if you WANT the file!
  --wgMLST              [EXPERIMENTAL] a better support of wgMLST schemes (check GrapeTree github for details).

Partitioning with HC:
  Specifications to genetic clusters with hierarchical clustering

  --HC-threshold HCMETHOD_THRESHOLD
                        [OPTIONAL] List of HC methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write the method name
                        (e.g. single). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. single-10). To get clustering at a specific range, indicate the range with a
                        hyphen separating minimum and maximum (e.g. single-2-10). If you prefer to exclusively use the '--pct-HC-threshold' argument, please set '--HC-threshold none'. Note: Threshold
                        values are inclusive, i.e. '--HC-threshold single-7' will consider samples with <= 7 differences as belonging to the same cluster! Default: single (Possible methods: single,
                        complete, average, weighted, centroid, median, ward)
  --pct-HC-threshold PCT_HCMETHOD_THRESHOLD
                        [OPTIONAL] Similar to '--HC-threshold' but the partition threshold for cluster definition is set as the proportion of differences to the final allelic schema size or number of
                        informative positions, e.g. '--pct-HC-threshold single-0.005' corresponds to a threshold of 5 allelic/SNP differences in a matrix with 1000 loci/sites under analysis. Ranges
                        CANNOT be specified.

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
  -root ROOT, --root ROOT
                        Set root of the input tree. Specify the leaf name to use as output. Alternatively, write 'midpoint', if you want to apply midpoint rooting method.

ReporTree cluster nomenclature:
  Cluster nomenclature instructions

  --nomenclature-file NOMENCLATURE
                        [OPTIONAL] Intended usage: provide a .tsv file with the nomenclature information to be used (normaly the '*_partitions.tsv' file of a previous ReporTree run) with the following
                        structure: First column should have the sample names; Subsequent columns have the partitions at any/all levels. Columns matching the column names of the partitions requested in
                        the current run (e.g. MST-7x1.0) will be used to name the current clusters at the respective levels. For more details on how to use this nomenclature system please visit our
                        github! Alternative usage (at your own risk): you can provide a .tsv file with just a single column with a grouping variable that does not match any partition requested in the
                        current run (e.g. ST). In this case, all clusters will be named according to this column (e.g. ST1.1,ST1.2,ST2.1...).
  --nomenclature-code-levels CODE_LEVELS
                        [OPTIONAL and only available for --analysis grapetree or HC] This argument allows getting a nomenclature code combining cluster information at different hierarchical levels (e.g.
                        if '150,30,7' is provided when the analysis is GrapeTree, a code combining the cluster names at these levels will be generated: C3-C2-C1). The order of levels indicated in the
                        command-line will be kept in the code (in the example, C3 indicates that the sample belongs to cluster_3 at MST-150). Of note, if a '--nomenclature-file' is provided in
                        subsequent ReporTree runs using the same method and thresholds, the nomenclature code will be kept. You can also add one metadata variable (e.g. Country) to get an extra layer to
                        the code (e.g. C3-C2-C1-Portugal). Partition thresholds can be indicated in this argument following the same rules as the arguments '-thr' and '-pct_thr' for GrapeTree or '--HC-
                        threshold' and '--pct-HC-threshold' for HC.

ReporTree metadata report:
  Specific parameters to report clustering/grouping information associated to metadata

  --columns_summary_report COLUMNS_SUMMARY_REPORT
                        Columns (i.e. variables of metadata) to get statistics for the derived genetic clusters or for other grouping variables defined in --metadata2report (comma-separated). If the
                        name of the column is provided, the different observations and the respective percentage are reported. If 'n_column' is specified, the number of the different observations is
                        reported. For example, if 'n_country' and 'country' are specified, the summary will report the number of countries and their distribution (percentage) per cluster/group.
                        Exception: if a 'date' column is in the metadata, it can also report first_seq_date, last_seq_date, timespan_days. Check '--list' argument for some help. Default =
                        n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days [the order of the list will be the order of the columns in the report]
  --partitions2report PARTITIONS2REPORT
                        Thresholds for which clustering information will be included in a joint report (comma-separated). Other alternatives: 'all' == all partitions; 'stability_regions' == first
                        partition of each stability region as determined by comparing_partitions_v2.py. Note: 'stability_regions' can only be inferred when partitioning TreeCluster or GrapeTree or HC is
                        run for all possible thresholds or when a similar partitions table is provided (i.e. sequential partitions obtained with the same clustering method) [all]. Partition thresholds
                        can be indicated in this argument following the same rules as the arguments '-thr' and '-pct_thr' for GrapeTree or '--HC-threshold' and '--pct-HC-threshold' for HC or '--method-
                        threshold' for TreeCluster.
  --metadata2report METADATA2REPORT
                        Columns of the metadata table for which a separated summary report must be provided (comma-separated)
  -f FILTER_COLUMN, --filter FILTER_COLUMN
                        [OPTIONAL] Filter for metadata columns to select the samples to analyze. This must be specified within quotation marks in the following format 'column< >operation< >condition'
                        (e.g. 'country == Portugal'). When more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When
                        filters include more than one column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are
                        important in this argument, so, do not leave spaces before and after commas/semicolons.
  --sample_of_interest SAMPLE_OF_INTEREST
                        [OPTIONAL] List of samples of interest for which summary reports will be created. This list can be a comma-separated list in the command line, or a comma-separated list in a
                        file, or a list in the first column of a tsv file. No headers should be provided in the input files. If nothing is provided, only the summary reports comprising all samples will
                        be generated.
  --zoom-cluster-of-interest ZOOM
                        [OPTIONAL and only available for --analysis grapetree or HC] Repeat the analysis using only the samples that belong to each cluster of the samples of interest at a given distance
                        threshold. This argument takes as input a comma-separated list of partitions for which you want the zoom-in. Partition thresholds can be indicated in this argument following the
                        same rules as the arguments '-thr' and '-pct_thr' for GrapeTree or '--HC-threshold' and '--pct-HC-threshold' for HC. This argument requires that a metadata table was provided
                        with '-m'. Default: no zoom-in.
  --subtree-of-interest SUBTREE
                        [OPTIONAL and only available for --analysis grapetree or HC] Repeat the analysis using the n closest samples of each sample of interest. This argument takes as input a comma-
                        separated list of n's, corresponding to the number of closest samples you want to include for the samples of interest. This argument requires that a metadata table was provided
                        with '-m'. Default: no subtree.
  --unzip               [OPTIONAL and only available for --analysis grapetree or HC] Provide the outputs of '--zoom-cluster-of-interest' and '--subtree-of-interest' in unzipped format.
  --frequency-matrix FREQUENCY_MATRIX
                        [OPTIONAL] Metadata column names for which a frequency matrix will be generated. This must be specified within quotation marks in the following format 'variable1,variable2'.
                        Variable1 is the variable for which frequencies will be calculated (e.g. for 'lineage,iso_week' the matrix reports the percentage of samples that correspond to each lineage per
                        iso_week). If you want more than one matrix you can separate the different requests with semicolon (e.g. 'lineage,iso_week;country,lineage'). If you want a higher detail in your
                        variable2 and decompose it into two columns you use a colon (e.g. lineage,country:iso_week will report the percentage of samples that correspond to each lineage per iso_week in
                        each country)
  --count-matrix COUNT_MATRIX
                        [OPTIONAL] Same as '--frequency-matrix' but outputs counts and not frequencies
  --mx-transpose        [OPTIONAL] Set ONLY if you want that the variable1 specified in '--frequency-matrix' or in '--count-matrix' corresponds to the matrix first column.
  --pivot               [OPTIONAL] Set ONLY if you want an additional table for each count/frequency matrix in pivot format.

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


#### Note on the '--metadata' argument

Please avoid white spaces in the metadata columns. If you still want to have blank spaces, ReporTree can deal with them but please remember to replace them by "_" everytime you want to write them in the command line.


#### Note on the '--missing-code' argument

By default, ReporTree assumes 0 as the code for missing data in SNP/allele matrices and N as the code for missing data in multi-sequence alignments. You can use this argument to mention which code you used. If you left an empty space in your matrix to indicate missing data, please write 'empty' in this argument.


#### Note on the '--method-threshold' argument

ReporTree uses [TreeCluster](https://github.com/niemasd/TreeCluster) to obtain clustering information from a SNP-distance tree. To provide some flexibility, it uses the argument '--method-threshold' to run this program for all the combinations of method-threshold that the user needs:
- Clustering at all possible thresholds of a single method (from 1 SNP to the maximum distance of the tree) -> set only the method (e.g. root_dist)
- Clustering at a specific threshold for a single method -> set the method and use a hyphen to indicate the threshold (e.g. root_dist-10)
- Clustering at all possible thresholds of a range for a single method -> set the method and use a hyphen to indicate the threshold range (e.g. root_dist-10-30)
- Clustering using two or more methods and/or different thresholds -> set the different requirements separated by a comma (e.g. root_dist,avg_dist-10,avg_dist-20-30)


#### Note on the '--HC-threshold' argument

'--HC-threshold' is the equivalent to '--method-threshold' for the hierarchical clustering analysis. As such, the input of this argument follows the same structure as mentioned in the previous paragraph. Possible methods: single, complete, average, weighted, centroid, median, ward.


#### Note on the '--pct_threshold' and the '--pct-HC-threshold' arguments

These arguments take as input the distance thresholds as a proportion of the number of alleles/sites used in the analysis. This is particularly useful for wgMLST and assembly-based approaches in which not all loci/positions of the schema/alignment are used, being difficult to the user to predict the size of the final core. For convenience, the proportion will be translated in absolute allele/SNP distances.


#### Note on the '--nomenclature-file' argument

The intended usage to this argument is to provide a partitions table of a previous run. ReporTree will then use this information to (re)name the clusters of the current run. However, we are quite flexible with this input. For instance, you can modify the name of particular clusters (e.g., Outbreak1) and this new name will be kept onwards. Alternative usage: you can provide a .tsv file with just a single column with a grouping variable that does not match any partition requested in the current run (e.g. ST). In this case, the generated clusters at all levels will take this variable as the prefix for cluster naming (e.g., ST1.1, ST1.2, ST2.1, ...).


#### Note on the '--partitions2report' argument

This argument is used to select the columns of the partitions table that will be incorporated into the metadata table. This arguments follows the same rules as the threshold argument you used. Please note that there are other expressions that you can use, namely, 'all' for all partitions (default) and 'stability_regions' to run comparing_partitions_v2.py and report the first partition of stable regions. Run '--list' for more help.


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


#### Simple ReporTree command line with an allele matrix as input to obtain the genetic clusters from a MST:

```bash
reportree.py -m metadata.tsv -a allele_matrix.tsv -out output --method MSTreeV2 -thr all --matrix-4-grapetree --columns_summary_report columns,summary,report --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2 --analysis grapetree
```


#### Simple ReporTree command line with an alignment as input to obtain the genetic clusters with hierarchical clustering (single-linkage):

```bash
reportree.py -m metadata.tsv -align alignment.fasta -out output -d dist --HC-threshold single --columns_summary_report columns,summary,report --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2 --analysis HC
```


#### Simple ReporTree command line with a newick tree to obtain the genetic clusters:

```bash
reportree.py -m metadata.tsv -t tree.nw -out output -d dist --method-threshold method1,method2-threshold --columns_summary_report columns,summary,report --partitions2report all --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2
```


#### More advanced ReporTree command line with an cgMLST allele matrix (e.g., _Listeria monocytogenes_) as input to obtain the genetic clusters from a MST and keep nomenclature from previous run (and generate hierarchical codes):

```bash
reportree.py -m metadata.tsv -a allele_matrix.tsv -out output --loci-called 0.95 --method MSTreeV2 -thr all --matrix-4-grapetree --columns_summary_report columns,summary,report --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2 --analysis grapetree --nomenclature-file previous_run_partitions.tsv --nomenclature-code 150,30,7
```


#### More advanced ReporTree command line with an wgMLST allele matrix (e.g., _Salmonella enterica_) as input to obtain the genetic clusters from a MST with levels defined as proportion, keep nomenclature from previous run (and generate hierarchical codes) and provide in-depth dynamic analysis (extended core) for samples of interest:

```bash
reportree.py -m metadata.tsv -a allele_matrix.tsv -out output --site-inclusion 0.99 --loci-called 0.95 --method MSTreeV2 -pct_thr 0.3,0.1,0.04 --matrix-4-grapetree --columns_summary_report columns,summary,report --metadata2report column1,column2 -f ‘column1 == observation;date > year-mm-dd’ --frequency-matrix variable1,variable2 --count-matrix variable1,variable2 --analysis grapetree --nomenclature-file previous_run_partitions.tsv --nomenclature-code 0.3,0.1,0.04 --sample_of_interest A,B,C --zoom-cluster-of-interest 0.04 --subtree-of-interest 20
```


## Examples

### Outbreak detection - bacterial foodborne pathogen (e.g. _Listeria monocytogenes_) - [click here](https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples#outbreak-detection---bacterial-foodborne-pathogen-eg-listeria-monocytogenes)

ReporTree can facilitate the routine surveillance and outbreak investigation of bacterial pathogens, such as foodborne pathogens. In [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples#outbreak-detection---bacterial-foodborne-pathogen-eg-listeria-monocytogenes), we provide a simple example of the usage of ReporTree to rapidly identify and characterize potential Listeriosis outbreaks. With a single command, ReporTree builds a MST from cgMLST data and **automatically extracts genetic clusters at three high resolution levels (<=4, <=7, <=14 allelic differences)**, and provides comprehensive reports about the sample collection (e.g. ST sequence count/frequency per year, etc).


### Large-scale genetic clustering and linkage to antibiotic resistance data (e.g. _Neisseria gonorrhoeae_) - [click here](https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples#large-scale-genetic-clustering-and-linkage-to-antibiotic-resistance-data-eg-neisseria-gonorrhoeae)

ReporTree can enhance genomics surveillance and quickly identify/characterize genetic clusters from large datasets. In [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples#large-scale-genetic-clustering-and-linkage-to-antibiotic-resistance-data-eg-neisseria-gonorrhoeae), with a single command line, we reproduce part of the extensive genomics analysis performed by [Pinto et al., 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208699/pdf/mgen-7-481.pdf) over 3,791 _N. gonorrhoeae_ genomes from isolates collected across Europe.


 ### Routine surveillance - viral pathogen (e.g. SARS-CoV-2) - [click here](https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples#routine-surveillance---viral-pathogen-eg-sars-cov-2)

ReporTree is currently applied to generate periodic reports about SARS-CoV-2 variant circulation in Portugal. In [ReporTree wiki](https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples#routine-surveillance---viral-pathogen-eg-sars-cov-2), we give some examples on how to rapidly generate key surveillance metrics taking as input metadata tables (tsv format) and rooted divergence (SNP) trees (newick format) provided for download in regular Nextstrain (auspice) builds, such as those maintained by the National Institute of Health Dr. Ricardo Jorge, Portugal (INSA) at https://insaflu.insa.pt/covid19/.


## Citation

If you run ReporTree, please cite the publication:    
[Mixão V, Pinto M, Sobral D, Di Pasquale A, Gomes JP, Borges V (2023) ReporTree: a surveillance-oriented tool to strengthen the linkage between pathogen genetic clusters and epidemiological data. _Genome Medicine_. doi: 10.1186/s13073-023-01196-1](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1)

Also, ReporTree relies on the work of other developers. So, depending on the functionalities you use, there are other tools that you must cite:
1. Grapetree: http://www.genome.org/cgi/doi/10.1101/gr.232397.117 (if you requested a grapetree analysis)
2. TreeCluster: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6705769/pdf/pone.0221068.pdf (if you provided a newick tree)
3. vcf2mst: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-08112-0 (if you provided a vcf or a list of variants)
4. ComparingPartitions: https://journals.asm.org/doi/10.1128/jcm.02536-05?permanently=true (if you requested "stability_regions")
5. Adjusted Wallace and cluster stability: https://www.biorxiv.org/content/10.1101/299347v1 (if you requested "stability_regions")
6. Ete3: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4868116/pdf/msw046.pdf (if you provided a newick tree)     
7. cgmlst-dists: https://github.com/tseemann/cgmlst-dists for original code and https://github.com/genpat-it/cgmlst-dists for improvements regarding memory efficiency 
8. snp-sites: https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056 (if you provided a multi-sequence alignment)


## Funding

This work was supported by funding from the European Union’s Horizon 2020 Research and Innovation programme under grant agreement No 773830: One Health European Joint Programme.
