FROM python:3.7.3-stretch
WORKDIR /Magphi
COPY . .
# TODO - Change the python up top
# Install the python package (and executable)
RUN pip3 install .

# Override some of the dependencies with the hard-coded versions
RUN pip3 install -r requirements-dev.txt
RUN conda install -c bioconda bedtools=2.29.9
RUN conda install -c bioconda samtools=1.9
#bedtools==2.29.2
 # TODO - add docker container with bedtools - biocontainers?
#samtools==1.9
