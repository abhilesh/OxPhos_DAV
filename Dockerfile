FROM continuumio/miniconda3:latest

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# 1. Create the environment
RUN conda create -y -n oxphos_dav --solver=libmamba -c conda-forge -c bioconda -c defaults \
    python=3.11 \
    cxx-compiler \
    pytest \
    lxml \
    pandas \
    pyarrow \
    numpy \
    biopython \
    scipy \
    requests \
    macse \
    iqtree \
    r-base \
    r-ape \
    r-phytools \
    && conda clean -afy

# 2. Force Docker to use this environment for all subsequent RUN commands during the build
SHELL ["conda", "run", "-n", "oxphos_dav", "/bin/bash", "-c"]

COPY setup.py /app/
COPY src /app/src/

# 3. Install pip dependencies, including the local package and aaindex
RUN pip install --no-cache-dir -e . aaindex

# 4a. Install DCA / tree-simulation packages
#     pydca:    plmDCA (pseudolikelihood) for phylogenetically-corrected coevolution
#     pyvolve:  sequence evolution simulation on gene trees (WAG model)
#     evcouplings: independent DCA implementation for top-50 cross-validation
# Install pyvolve and evcouplings normally
RUN pip install --no-cache-dir pyvolve evcouplings

# Force install pydca without its ancient dependencies (uses the modern Conda SciPy instead)
RUN pip install --no-cache-dir --no-deps pydca

# 4b. CRITICAL FIX: Append the activation command to .bashrc to override Miniconda's default behavior
RUN echo "conda activate oxphos_dav" >> ~/.bashrc

CMD ["/bin/bash"]
