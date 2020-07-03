FROM ubuntu:20.04
COPY ./ /genetIC

RUN apt-get update && apt-get install -y g++-9 libgsl-dev libfftw3-dev python3-numpy python3-pip
RUN pip3 install pynbody
RUN cd /genetIC && make clean && make

ENTRYPOINT ["/genetIC/genetIC"]