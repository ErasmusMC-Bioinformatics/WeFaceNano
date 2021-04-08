FROM ubuntu:latest

RUN apt-get update && apt-get install -y make curl wget git python3 unrar

ADD . /wefacenano
WORKDIR /wefacenano

RUN make install-conda
RUN make create-env
RUN make install
RUN make download-testdata

CMD make run
EXPOSE 8008
