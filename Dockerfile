FROM python:3.9.7-stretch
WORKDIR /Magphi
COPY . .
# TODO - Change the python up top
# Install the python package (and executable)
RUN pip3 install .

# Override some of the dependencies with the hard-coded versions
RUN pip3 install -r requirements-dev.txt

RUN apt-get update && apt-get install -y \
  bedtools \
  samtools \
  && rm -rf /var/lib/apt/lists/*
#bedtools==2.29.2
 # TODO - add docker container with bedtools - biocontainers?
#samtools==1.9
