FROM ubuntu:18.04
MAINTAINER Ricard Campos <ricard.campos@coronis.es>

# =================================================
# Install all dependencies installable via apt-get
# =================================================
RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get -y install build-essential gcc-4.8 g++-4.8 libarpack2-dev libarpack++2-dev libeigen3-dev cmake wget
# Set the default compiler to be version 4.8
ENV CC=gcc-4.8
ENV CXX=g++-4.8

# ===============================================================
# Install the dependencies that need to be compiled from sources
# ===============================================================

# --- Levmar ---
RUN apt-get -y install libf2c2-dev
WORKDIR /opt
RUN wget users.ics.forth.gr/~lourakis/levmar/levmar-2.6.tgz
RUN tar -zxvf levmar-2.6.tgz
WORKDIR /opt/levmar-2.6
RUN make

# --- Boost 1.52 ---
WORKDIR /opt
RUN wget http://sourceforge.net/projects/boost/files/boost/1.52.0/boost_1_52_0.tar.gz
RUN tar -zxvf boost_1_52_0.tar.gz
WORKDIR /opt/boost_1_52_0
# Patch! see patches/Readme.md
COPY patches/boost_1_52_0/tools/build/v2/user-config.jam /opt/boost_1_52_0/tools/build/v2/user-config.jam
RUN ./bootstrap.sh
# Note the "exit 0", not all the libraries of boost will be compiled, but all the ones that we need do!
RUN ./b2 install; exit 0 

# --- CGAL 4.0.2 ---
WORKDIR /opt
RUN wget https://github.com/CGAL/cgal/archive/refs/tags/releases/CGAL-4.0.2.tar.gz
RUN tar -zxvf CGAL-4.0.2.tar.gz
WORKDIR /opt/cgal-releases-CGAL-4.0.2
# Patch! see patches/Readme.md
COPY patches/CGAL-4.0.2/Installation/src/CMakeLists.txt /opt/cgal-releases-CGAL-4.0.2/Installation/src/CMakeLists.txt
RUN mkdir build
WORKDIR /opt/cgal-releases-CGAL-4.0.2/build
RUN cmake -DWITH_CGAL_Core=ON ..
RUN make install

# =========================================================
# Compile the splats_surface_reconstruction project itself
# =========================================================

WORKDIR /opt/splats_surface_reconstruction
COPY ./3rd_party /opt/splats_surface_reconstruction/3rd_party
COPY ./cmake /opt/splats_surface_reconstruction/cmake
COPY ./include /opt/splats_surface_reconstruction/include
COPY ./src /opt/splats_surface_reconstruction/src
COPY ./demo /opt/splats_surface_reconstruction/demo
COPY CMakeLists.txt /opt/splats_surface_reconstruction
COPY Readme.md /opt/splats_surface_reconstruction
WORKDIR /opt/splats_surface_reconstruction/build
ENV LEVMAR_INC_DIR=/opt/levmar-2.6
ENV LEVMAR_LIB_DIR=/opt/levmar-2.6
RUN cmake -DSPLATS_BUILD_TESTS=False ..
RUN make install
ENV LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
WORKDIR /opt/splats_surface_reconstruction

# ========
# Cleanup
# ========

# Cleanup (remove unnecessary files from apt-get, makes the image a bit lighter...)
RUN apt-get -y autoremove &&\
    apt-get clean autoclean &&\
    rm -fr /var/lib/apt/lists/{apt,dpkg,cache,log} /tmp/* /var/tmp/*

# Remove all sources from the libraries
RUN rm /opt/levmar-2.6.tgz
RUN rm -rf /opt/levmar-2.6
RUN rm /opt/boost_1_52_0.tar.gz
RUN rm -rf /opt/boost_1_52_0
RUN rm /opt/CGAL-4.0.2.tar.gz
RUN rm -rf /opt/cgal-releases-CGAL-4.0.2


