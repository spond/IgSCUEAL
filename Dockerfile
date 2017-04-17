# Use phusion/baseimage as base image
FROM phusion/baseimage:0.9.19

# Set environment variables the phusion way
RUN echo en_US.UTF-8 > /etc/container_environment/LANGUAGE
RUN echo en_US.UTF-8 > /etc/container_environment/LANG
RUN echo en_US.UTF-8 > /etc/container_environment/LC_ALL
RUN echo UTF-8 > /etc/container_environment/PYTHONIOENCODING

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

MAINTAINER Simon Frost <sdwfrost@gmail.com>

## Set a default user. Available via runtime flag `--user docker`
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& mkdir /home/docker/programs \
	&& addgroup docker staff

RUN apt-get update -qq && \
  DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
 	build-essential \
 	python3 \
	git \
	cmake \
	openmpi-bin \
	libopenmpi-dev

## Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## HyPhy
RUN cd /home/docker/programs && \
  git clone https://github.com/veg/hyphy.git && \
	cd hyphy && \
  git checkout v2.3-dev && \
  cmake . && \
  make MPI && \
  make install && \
  cd /home/docker/programs && \
  rm -rf hyphy


RUN cd /home/docker/programs && \
  git clone https://github.com/antibodyome/IgSCUEAL.git && \
  cd IgSCUEAL && \
  git checkout v2

VOLUME ["/home/docker"]
