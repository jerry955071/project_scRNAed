FROM python:3.9

# install gawk
RUN wget https://ftp.gnu.org/gnu/gawk/gawk-5.2.1.tar.gz && \
    tar -xzf gawk-5.2.1.tar.gz && \
    cd gawk-5.2.1 && \
    sh ./configure && \
    make && \
    make install

# install vawk (version 0.0.2)
RUN git clone https://github.com/cc2qe/vawk.git
ENV PATH=$PATH:/vawk

# set default command
CMD sh