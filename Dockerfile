FROM centos:7.9.2009

WORKDIR /src

RUN yum update -y && \
    yum install -y mercurial.x86_64 git && \
    yum install -y "@Development Tools" qt5-qtbase-devel.x86_64  qt5-qttools-devel.x86_64 && \
    yum install -y libjpeg-turbo-devel libtiff-devel java-1.8.0-openjdk-devel tcsh.x86_64

# This step is to build from its repo, but it is SO SLOW
#RUN hg clone http://bio3d.colorado.edu/imod/nightlyBuilds/IMOD /src/IMOD
#RUN cd /src/IMOD && hg update -r IMOD_4-11-7

RUN git clone https://github.com/uqrali1/raza.git /src/raza
RUN tar -xvf /src/raza/IMOD_4-11-7.tar.gz -C /src
RUN mv /src/IMOD_4-11-7 /src/IMOD

RUN cp /src/raza/sourcecode/raza.c /src/IMOD/mrc/ && cp /src/raza/sourcecode/Makefile.mrc.IMOD_4-11-7 /src/IMOD/mrc/

RUN mkdir -p /opt/local/imod-raza/4.11.7

ENV PATH="/usr/lib64/qt5/bin:${PATH}"
ENV QTDIR="/usr/lib64/qt5"
ENV QMAKESPEC="linux-g++-64" 

RUN cd /src/IMOD && ./setup -i /opt/local/imod-raza/4.11.7 && make -j 2 && make install

ENV IMOD="/opt/local/imod-raza/4.11.7"
ENV PATH="${IMOD}/bin:${PATH}" 
ENV LD_LIBRARY_PATH="${IMOD}/lib:${LD_LIBRARY_PATH}" 
ENV MANPATH="${IMOD}/man:${MANPATH}" 
