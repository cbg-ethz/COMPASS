# Teeny-tiny matplotlib image based on alpine
FROM ubuntu:latest

LABEL author="Etienne Sollier" \
      description="COMPASS"

RUN apk add --no-cache bash

# Add the COMPASS source files to the container
ADD . /usr/src/COMPASS
WORKDIR /usr/src/COMPASS

# Install COMPASS
RUN make

# Set up entrypoint and cmd for easy docker usage
ENTRYPOINT [ "COMPASS" ]
CMD [ "." ]