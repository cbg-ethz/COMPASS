# Teeny-tiny matplotlib image based on alpine
FROM ubuntu:latest

LABEL author="Etienne Sollier" \
      description="COMPASS"

RUN apt-get update && apt-get install -y make g++ graphviz

# Add the COMPASS source files to the container
ADD . /app
WORKDIR /app
COPY . .

# Install COMPASS
RUN make

# Set up entrypoint and cmd for easy docker usage
ENTRYPOINT [ "/app/COMPASS" ]
CMD [ "." ]