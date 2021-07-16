FROM centos:7.9.2009

WORKDIR /src

RUN yum update -y && \
    yum install -y mercurial.x86_64 git && \
    yum install -y "@Development Tools" qt5-qtbase-devel.x86_64  qt5-qttools-devel.x86_64 && \
    yum install -y libjpeg-turbo-devel libtiff-devel java-1.8.0-openjdk-devel tcsh.x86_64

RUN hg clone http://bio3d.colorado.edu/imod/nightlyBuilds/IMOD /src

RUN cd /src && hg update -r IMOD_4-11-7

RUN cp sourcecode/raza.c /src/mrc/ && cp sourcecode/Makefile.mrc.IMOD_4-11-7 /src/mrc/

RUN mkdir -p /opt/local/imod-raza/4.11.7

ENV PATH="/usr/lib64/qt5/bin:${PATH}"
ENV QTDIR="/usr/lib64/qt5"
ENV QMAKESPEC="linux-g++-64" 

RUN cd /src && ./setup -i /opt/local/imod-raza/4.11.7 && make -j 2 && make install

ENV IMOD="/opt/local/imod-raza/4.11.7"
ENV PATH="${IMOD}/bin:${PATH}" 
ENV LD_LIBRARY_PATH="${IMOD}/lib:${LD_LIBRARY_PATH}" 
ENV MANPATH="${IMOD}/man:${MANPATH}" 