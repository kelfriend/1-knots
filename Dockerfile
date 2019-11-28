FROM gapsystem/gap-docker

COPY --chown=1000:1000 . $HOME/1-knots

RUN cd /gap-4.10.2/pkg \
  && wget -q http://hamilton.nuigalway.ie/Hap/hap1.23.tar.gz \
  && tar xzf hap1.23.tar.gz \
  && rm hap1.23.tar.gz \
  && cd Hap1.23 \
  && ./configure \ 
  && make \
  
WORKDIR $HOME/1-knots
