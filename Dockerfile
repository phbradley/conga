FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
  apt-get install -y wget bzip2 ca-certificates curl git && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y python3 && apt-get install -y python3-pip
RUN apt-get install -y inkscape --fix-missing
RUN pip3 install jedi==0.17.2
RUN pip3 install ipython==7.18.1
RUN pip3 install scanpy[leiden]
RUN mkdir /gitrepos && cd /gitrepos && git clone https://github.com/phbradley/conga
RUN cd /gitrepos/conga/tcrdist_cpp/ && make

CMD [ "/bin/bash" ]
