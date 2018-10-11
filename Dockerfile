# Generated by Neurodocker version 0.4.0
# Timestamp: 2018-08-01 20:21:01 UTC
# 
# Thank you for using Neurodocker. If you discover any issues
# or ways to improve this software, please submit an issue or
# pull request on our GitHub repository:
# 
#     https://github.com/kaczmarj/neurodocker

FROM centos:7

ARG DEBIAN_FRONTEND="noninteractive"

ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ND_ENTRYPOINT="/neurodocker/startup.sh"
RUN export ND_ENTRYPOINT="/neurodocker/startup.sh" \
    && yum install -y -q \
           bzip2 \
           ca-certificates \
           curl \
           localedef \
           unzip \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/* \
    && localedef -i en_US -f UTF-8 en_US.UTF-8 \
    && chmod 777 /opt && chmod a+s /opt \
    && mkdir -p /neurodocker \
    && if [ ! -f "$ND_ENTRYPOINT" ]; then \
         echo '#!/usr/bin/env bash' >> "$ND_ENTRYPOINT" \
    &&   echo 'set -e' >> "$ND_ENTRYPOINT" \
    &&   echo 'if [ -n "$1" ]; then "$@"; else /usr/bin/env bash; fi' >> "$ND_ENTRYPOINT"; \
    fi \
    && chmod -R 777 /neurodocker && chmod a+s /neurodocker

ENTRYPOINT ["/neurodocker/startup.sh"]

#ants
ENV ANTSPATH="/opt/ants-2.2.0" \
    PATH="/opt/ants-2.2.0:$PATH"
RUN echo "Downloading ANTs ..." \
    && mkdir -p /opt/ants-2.2.0 \
    && curl -fsSL --retry 5 https://dl.dropbox.com/s/2f4sui1z6lcgyek/ANTs-Linux-centos5_x86_64-v2.2.0-0740f91.tar.gz \
    | tar -xz -C /opt/ants-2.2.0 --strip-components 1

#fsl
ENV FSLDIR="/opt/fsl-5.0.11" \
    PATH="/opt/fsl-5.0.11/bin:$PATH" 
ENV FSLCLUSTER_MAILOPTS=n
ENV FSLLOCKDIR=
ENV FSLMACHINELIST=
ENV FSLMULTIFILEQUIT=TRUE
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLREMOTECALL=
ENV FSLTCLSH=/usr/bin/tclsh
ENV FSLWISH=/usr/bin/wish
ENV POSSUMDIR="/opt/fsl-5.0.11"
RUN yum install -y -q \
           bc \
           file \
           libGL \
           libGLU \
           libICE \
           libSM \
           libX11 \
           libXcursor \
           libXext \
           libXft \
           libXinerama \
           libXrandr \
           libXt \
           libgomp \
           libjpeg \
           libmng \
           libpng12 \
           wget \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/* \
    && echo "Downloading FSL ..." \
    && mkdir -p /opt/fsl-5.0.11 \
    && curl -fsSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-5.0.11-centos6_64.tar.gz \
    | tar -xz -C /opt/fsl-5.0.11 --strip-components 1 \
    && sed -i '$iecho Some packages in this Docker container are non-free' $ND_ENTRYPOINT \
    && sed -i '$iecho If you are considering commercial use of this container, please consult the relevant license:' $ND_ENTRYPOINT \
    && sed -i '$iecho https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence' $ND_ENTRYPOINT \
    && sed -i '$isource $FSLDIR/etc/fslconf/fsl.sh' $ND_ENTRYPOINT \
    && echo "Installing FSL conda environment ..." \
    && bash /opt/fsl-5.0.11/etc/fslconf/fslpython_install.sh -f /opt/fsl-5.0.11

#mrtrix3
#ENV PATH="/opt/mrtrix3-3.0/bin:$PATH"
#RUN echo "Downloading MRtrix3 ..." \
#    && mkdir -p /opt/mrtrix3-3.0 \
#    && curl -fsSL --retry 5 https://dl.dropbox.com/s/2g008aaaeht3m45/mrtrix3-Linux-centos6.tar.gz \
#    | tar -xz -C /opt/mrtrix3-3.0 --strip-components 1
ENV PATH=/opt/rh/devtoolset-4/root/usr/bin:$PATH:/usr/lib64/qt4/bin
RUN yum update -y \
    && yum install -y -q git gcc-c++ zlib-devel gsl-devel qt-devel \
    && yum install -y -q https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm \
    && yum install -y -q eigen3-devel \
    && yum install -y -q centos-release-scl centos-release-scl-rh \
    && yum install -y -q devtoolset-4-gcc-c++ \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/* \
    && cd /opt && git clone https://github.com/MRtrix3/mrtrix3.git \
    && cd /opt/mrtrix3 && ./configure &&  ./build
ENV PATH=/opt/mrtrix3/bin:$PATH

#c3d
ENV C3DPATH="/opt/convert3d-1.0.0" \
    PATH="/opt/convert3d-1.0.0/bin:$PATH"
RUN echo "Downloading Convert3D ..." \
    && mkdir -p /opt/convert3d-1.0.0 \
    && curl -fsSL --retry 5 https://sourceforge.net/projects/c3d/files/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz/download \
    | tar -xz -C /opt/convert3d-1.0.0 --strip-components 1

#HCP
RUN yum install -y -q \
           python-pip \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/* \
    && pip install numpy \
    && cd /opt && git clone https://github.com/Washington-University/HCPpipelines
