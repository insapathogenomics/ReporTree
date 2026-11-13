
FROM conda/miniconda3

WORKDIR /app

# Create an isolated environment for ReporTree
RUN conda create -y -p /opt/reportree-env \
    -c etetoolkit \
    -c anaconda \
    -c bioconda \
    -c conda-forge \
    -c defaults \
    python=3.8 \
    biopython=1.77 \
    cgmlst-dists=0.4.0 \
    ete3 \
    grapetree=2.1 \
    numba=0.55.1 \
    numpy=1.19.5 \
    pandas=1.1.3 \
    scipy \
    scikit-learn \
    libgomp \
    snp-sites=2.5.1 \
    treecluster=1.0.3 \
    zip \
    pytest \
    git \
    && conda clean -afy

# Make the new environment available
ENV PATH="/opt/reportree-env/bin:/app/ReporTree:/app/ReporTree/scripts:/app/ReporTree/scripts/ComparingPartitions:/app/ReporTree/scripts/GrapeTree:${PATH}"

# Make OpenMP library available
ENV LD_LIBRARY_PATH="/opt/reportree-env/lib"

# Download ReporTree and dependencies
RUN git clone -b v2.7.1 https://github.com/insapathogenomics/ReporTree \
    && chmod 755 ReporTree/reportree.py \
    && chmod 755 ReporTree/scripts/partitioning_grapetree.py \
    && chmod 755 ReporTree/scripts/partitioning_HC.py \
    && chmod 755 ReporTree/scripts/partitioning_treecluster.py \
    && chmod 755 ReporTree/scripts/alignment_processing.py \
    && chmod 755 ReporTree/scripts/metadata_report.py \
    && cd ReporTree/scripts/ \
    && git clone https://github.com/insapathogenomics/ComparingPartitions.git \
    && git clone https://github.com/insapathogenomics/GrapeTree.git \
    && git clone https://github.com/insapathogenomics/vcf2mst.git \
    && cd ComparingPartitions \
    && chmod 755 comparing_partitions_v2.py
