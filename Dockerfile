# Teeny-tiny matplotlib image based on alpine
FROM ubuntu:latest

LABEL author="Etienne Sollier" \
      description="COMPASS"

# Add the COMPASS source files to the container
ADD . /usr/src/COMPASS
WORKDIR /usr/src/COMPASS

# Install COMPASS
RUN make

# Set up entrypoint and cmd for easy docker usage
ENTRYPOINT [ "COMPASS" ]
CMD [ "." ]