FROM continuumio/miniconda

COPY containers/ersa2/ersa ersa

RUN conda env create -f ersa/ersa2.yaml && \
    conda clean --all -y

RUN echo "source activate ersa2" > ~/.bashrc

ENV PATH "$PATH:/ersa"

ENV PATH /opt/conda/envs/ersa2/bin:${PATH}

ENTRYPOINT ["python"]