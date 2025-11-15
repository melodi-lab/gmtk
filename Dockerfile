FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV LANGUAGE=C.UTF-8

ARG uid
ARG username

# docker build -t gmtk:24.04 --build-arg uid=$UID --build-arg username="John Smith <jsmith@isp.com>" .
# docker run --rm -it -v $HOME/gmtk:/home/gmtk/gmtk gmtk:24.04

RUN apt-get -y update && \
    apt-get -y install apt-utils && \
    apt-get -y install linux-headers-generic && \
    apt-get -y install kmod pkg-config && \   
    apt-get -y install \
        sudo software-properties-common \
	git \
	less vim lynx wget gnupg \
	g++ clang lcov valgrind lldb gdb automake autoconf2.69 make cmake \
	bison flex help2man \
        texlive-latex-extra
RUN for f in /usr/bin/*2.69; do x=${f##*/}; ln -s $f /usr/local/bin/${x%%2.69}; done

RUN userdel -r ubuntu
RUN adduser --quiet --disabled-password -gecos GMTK,,, --uid ${uid} gmtk
RUN echo "gmtk ALL=(ALL:ALL) NOPASSWD:ALL" >>/etc/sudoers

USER gmtk
WORKDIR /home/gmtk
CMD ["bash"]
