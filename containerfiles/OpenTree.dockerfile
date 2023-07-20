# Docker image with the software required for obtaining phylogenetic trees from Open Tree of Life
# https://github.com/McTavishLab/GBIF-Biodiverse-OpenTree

FROM python:3.10.12-bullseye

MAINTAINER vladimir.mikryukov@ut.ee

ENV LANG C.UTF-8
ENV PYTHONUNBUFFERED=true
SHELL [ "/bin/bash", "-c" ]


RUN git clone --depth 1 https://github.com/McTavishLab/GBIF-Biodiverse-OpenTree
# RUN pip install -r GBIF-Biodiverse-OpenTree/requirements.txt

## NB! For `python-opentree`, the `itol_annot` branch is used!
RUN pip install DendroPy \
    && pip install git+https://github.com/OpenTreeOfLife/python-opentree.git@itol_annot#egg=opentree \
    && chmod +x /GBIF-Biodiverse-OpenTree/scripts/*.py

    ## && pip install git+https://github.com/OpenTreeOfLife/python-opentree.git \
    ## && pip install git+https://github.com/OpenTreeOfLife/python-opentree.git@main#egg=opentree \

# ENV PATH=/GBIF-Biodiverse-OpenTree/scripts/:$PATH

# CMD [ "python", "GBIF-Biodiverse-OpenTree/scripts/induced_synth_subtree_from_csv.py" ]

CMD [ "python", "-u" ]
