FROM python:3.9.7-buster
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
  ncbi-blast+ \
  && rm -rf /var/lib/apt/lists/*
  # TODO - install BLAST+!


