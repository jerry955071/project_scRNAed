FROM python:2.7

# TODO: https://github.com/BioinfoUNIBA/REDItools/blob/master/README_1.md#download
RUN pip install pysam
RUN git clone https://github.com/BioinfoUNIBA/REDItools \
    && cd REDItools \
    && python setup.py install
