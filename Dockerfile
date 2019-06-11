FROM rust:1.31

RUN useradd -ms /bin/bash rustacean
WORKDIR /home/rustacean
COPY . .
RUN chown -R rustacean:rustacean /home/rustacean

USER rustacean
