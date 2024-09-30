FROM ubuntu:24.04
COPY ./ /genetIC

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y g++-12 libgsl-dev libfftw3-dev python3-pip pkg-config libhdf5-serial-dev
RUN pip3 install git+https://github.com/pynbody/pynbody.git --break-system-packages # should use a venv but this is simpler for now
RUN cd /genetIC && make clean && make

ENTRYPOINT ["/genetIC/genetIC"]
