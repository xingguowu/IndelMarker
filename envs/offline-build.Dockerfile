FROM rockylinux/rockylinux:8.10@sha256:e8a49c5403b687db05d4d67333fa45808fbe74f36e683cec7abb1f7d0f2338c6
ARG MINIFORGE_VERSION=26.7.2-0
ARG MINIFORGE_SHA256=281b0ac7d550802efc81af633225a5e6116d29ae72f3ab4eae7168c3931a4c05
RUN curl -fL "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-Linux-x86_64.sh" -o /tmp/Miniforge3.sh \
    && printf '%s  %s\n' "$MINIFORGE_SHA256" /tmp/Miniforge3.sh | sha256sum -c - \
    && bash /tmp/Miniforge3.sh -b -p /opt/build-conda \
    && rm /tmp/Miniforge3.sh
ENV PATH=/opt/build-conda/bin:$PATH
RUN python -m pip install --no-deps conda-pack==0.8.1 setuptools==80.9.0
WORKDIR /work
