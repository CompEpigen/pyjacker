FROM ubuntu:22.04
RUN apt update && apt install -y --no-install-recommends bzip2 python3-pip gcc g++ make zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev
RUN pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir pyjacker

