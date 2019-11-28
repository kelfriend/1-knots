FROM gapsystem/gap-docker

COPY --chown=1000:1000 . $HOME/1-knots
  
WORKDIR $HOME/1-knots
