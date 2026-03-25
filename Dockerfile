FROM continuumio/miniconda3:latest

# Prevent Python from writing pyc files and enable unbuffered logging
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# Install Python packages and bioinformatics tools via conda
RUN conda install -y -c conda-forge -c bioconda \
    python=3.11 \
    lxml \
    pandas \
    numpy \
    biopython \
    requests \
    mafft \
    iqtree \
    && conda clean -afy

# Set the Python path so local modules like 'utils' can be imported easily
ENV PYTHONPATH="/app/src"

CMD ["bash"]