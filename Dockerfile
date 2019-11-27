FROM gapsystem/gap-docker

MAINTAINER 

COPY --chown=1000:1000 . $HOME/1-knots

RUN sudo pip3 install ipywidgets RISE

RUN jupyter-nbextension install rise --user --py

RUN jupyter-nbextension enable rise --user --py

USER gap

WORKDIR $HOME/1-knots
