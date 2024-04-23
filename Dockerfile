FROM ubuntu

RUN apt-get update && apt-get install -y \
    g++ \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip && \
    pip3 install numpy

COPY testing /root/testing

WORKDIR /root/testing

SHELL ["bash", "-c"]
CMD ["./tester.sh"]

