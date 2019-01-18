FROM jupyter/datascience-notebook

MAINTAINER Guilherme Garcia <pectinidae@gmail.com>

USER root

RUN apt-get update

RUN apt-get -y install libglu-dev libfreetype6-dev

COPY requirements.R .

RUN R CMD BATCH requirements.R

RUN rm requirements.R

RUN usermod -u 1003 jovyan

USER jovyan

WORKDIR /home/jovyan/app

COPY requirements.txt .

RUN conda install --file=requirements.txt

# RUN rm requirements.txt

# RUN jupyter notebook --generate-config

RUN echo \
    "c.NotebookApp.token = '4fcbc72db3dadfc4585c374a5c83b566b9a86944b1489e79'" \
    >> /home/jovyan/.jupyter/jupyter_notebook_config.py

RUN echo \
    "c.NotebookApp.password = u'sha1:0bff055a03ce:c09cde8bab532481515d631285948252f72e8cc5'" \
    >> /home/jovyan/.jupyter/jupyter_notebook_config.py
