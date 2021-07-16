*Note*: The only reason IMOD source code is included because cloing the repo is SO SLOW


This instruction is for:
* Centos 7.9.2009
* IMOD 4.11-7 (might work with higher version but did not test)


Steps to compile and install (copy from Dockerfile =)) )


    yum update -y
    
    yum install -y mercurial.x86_64 git
    
    yum install -y "@Development Tools" qt5-qtbase-devel.x86_64  qt5-qttools-devel.x86_64
    
    yum install -y mesa-libGL-devel mesa-libGLU-devel.x86_64 libjpeg-turbo-devel libtiff-devel java-1.8.0-openjdk-devel tcsh.x86_64

    mkdir /src

    tar -xvf IMOD_4-11-7.tar.gz -C /src

    mv /src/IMOD_4-11-7 /src/IMOD

    cp sourcecode/raza.c /src/IMOD/mrc/

    cp sourcecode/Makefile.mrc.IMOD_4-11-7 /src/IMOD/mrc/Makefile

    mkdir -p /opt/local/imod-raza/4.11.7

    export PATH="/usr/lib64/qt5/bin:${PATH}"
    export QTDIR="/usr/lib64/qt5"
    export QMAKESPEC="linux-g++-64" 

    cd /src/IMOD && ./setup -i /opt/local/imod-raza/4.11.7 && make -j 2 && make install
    export IMOD="/opt/local/imod-raza/4.11.7"
    export PATH="${IMOD}/bin:${PATH}" 
    export LD_LIBRARY_PATH="${IMOD}/lib:${LD_LIBRARY_PATH}" 
    export MANPATH="${IMOD}/man:${MANPATH}" 


Follow README.md to see how to run it