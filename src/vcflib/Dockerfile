FROM ubuntu:22.04

ENV TZ=Asia/Taipei
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libvcflib-tools libvcflib-dev