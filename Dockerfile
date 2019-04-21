# To use this Dockerfile, run the following:
#
# docker build --tag ryanjoneil/tsppd-hybrid \
#              --build-arg UID=$(id -u) \
#              --build-arg GID=$(id -g) .

FROM ubuntu:19.04
MAINTAINER roneil1@gmu.edu

RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y curl

# Gecode
RUN apt-get install -y libmpfr-dev
RUN apt-get install -y qt5-default qtbase5-dev

RUN curl -L https://github.com/Gecode/gecode/archive/release-6.2.0.tar.gz | \
    tar xvz && \
    cd gecode-release-6.2.0/ && \
    ./configure --with-mpfr-lib --disable-examples --disable-flatzinc && \
    make install -j$(nproc --all) && \
    cd .. && \
    rm -rf gecode-release-6.2.0/

# Gurobi
RUN cd /opt
RUN curl -L https://packages.gurobi.com/8.1/gurobi8.1.1_linux64.tar.gz | \
    tar xvz && \
    cd gurobi811/linux64/src/build/ && \
    make

# User mapping
ARG UNAME=tsppd
ARG UID
ARG GID

RUN groupadd -g $GID -o $UNAME
RUN useradd -m -u $UID -g $GID -o -s /bin/bash $UNAME

USER $UNAME
ENV PATH="/opt/gurobi811/linux64/bin/"
ENV LD_LIBRARY_PATH="/opt/gurobi811/linux64/src/build/"

WORKDIR /home/$UNAME
