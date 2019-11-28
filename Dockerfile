FROM gapsystem/gap-docker-master:20180820

COPY --chown=1000:1000 . $HOME/1-knots
