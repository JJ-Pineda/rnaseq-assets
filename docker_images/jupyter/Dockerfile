FROM --platform=linux/amd64 jupyter/datascience-notebook:x86_64-ubuntu-22.04

USER root

WORKDIR /home/jovyan

RUN wget -nv https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip && \
    mv awscli-exe-linux-x86_64.zip awscliv2.zip && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws

RUN chown jovyan -R /home/jovyan

USER jovyan

#RUN conda create -n my-new-python-env python=3.9 -y && \
#    source activate my-new-python-env && \
#    conda install ipykernel && \
#    python -m ipykernel install --user --name my-new-python-env --display-name "Python (My-New-Python-Env)" && \
#    source deactivate

RUN mkdir /home/jovyan/javier

WORKDIR /home/jovyan/javier
