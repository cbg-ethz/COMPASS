FROM ubuntu:latest

LABEL author="Etienne Sollier" \
      description="COMPASS"

RUN apt-get update && apt-get install -y make g++ graphviz

# Add the COMPASS source files to the container
WORKDIR /app
COPY . .

# Install COMPASS
RUN make
RUN cp COMPASS /usr/local/bin/

WORKDIR /data/