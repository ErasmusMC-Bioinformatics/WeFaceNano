FROM continuumio/miniconda3:latest

RUN mkdir /root/WeFaceNano
WORKDIR /root/WeFaceNano
ENV QT_SELECT=5
ENV PERL5LIB="/root/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}"
ENV PERL_LOCAL_LIB_ROOT="/root/perl5${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"
ENV PERL_MB_OPT="--install_base \"/root/perl5\""
ENV PERL_MM_OPT="INSTALL_BASE=/root/perl5"
ENV PATH "/root/perl5/bin${PATH:+:${PATH}}"
ENV PATH "$PATH:/root/WeFaceNano/static"
ENV PATH "$PATH:/root/perl5/bin"
ENV PATH "$PATH:/root/.local/bin"
ENV PATH "$PATH:/root/Flye/bin"
ENV PATH "$PATH:/root/kmergenie-1.7051"
RUN export PATH=$PATH
RUN apt-get update --fix-missing
RUN apt-get -y install r-base perlbrew cpanminus ncbi-blast+ blast2 libgd-dev libbz2-dev cmake build-essential g++ qtbase5-dev libqt5svg5-dev libxml-libxml-perl
RUN pip install django==2.2.8 biopython NanoPlot tabulate cgecore --user
WORKDIR /root
RUN perlbrew init
RUN perlbrew install perl-5.26.1
RUN perlbrew use
RUN wget https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip
RUN unzip Bandage_Ubuntu_static_v0_8_1.zip
WORKDIR /root
RUN git clone https://github.com/rrwick/Porechop.git
WORKDIR /root/Porechop
RUN python3 setup.py install
WORKDIR /root
RUN git clone https://github.com/lh3/minimap2
WORKDIR /root/minimap2
RUN make
WORKDIR /root
RUN git clone https://github.com/lh3/miniasm
WORKDIR /root/miniasm
RUN make
WORKDIR /root
RUN git clone https://github.com/fenderglass/Flye
WORKDIR /root/Flye
RUN python setup.py build
WORKDIR /root
RUN git clone --recursive https://github.com/isovic/racon.git racon
WORKDIR /root/racon
RUN mkdir build
WORKDIR /root/racon/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make
WORKDIR /root
RUN wget http://kmergenie.bx.psu.edu/kmergenie-1.7051.tar.gz && tar -xvzf kmergenie-1.7051.tar.gz
WORKDIR /root/kmergenie-1.7051
RUN make
WORKDIR /root
RUN git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
RUN git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
WORKDIR /root/plasmidfinder
WORKDIR /root
RUN git clone https://github.com/Kzra/Simple-Circularise

# RUN cpan -u
# RUN cpan force install Getopt::Long Bio::SeqIO Bio::SearchIO Try::Tiny::Retry GD

COPY . /root/WeFaceNano

RUN chmod 777 /root/WeFaceNano/wefacenano_start.sh

EXPOSE 8008

WORKDIR /root/WeFaceNano
CMD ["/root/WeFaceNano/wefacenano_start.sh"]