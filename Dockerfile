FROM rust:1.31

WORKDIR /usr/src/app
COPY . .

RUN cargo install --path .

CMD ["cargo", "run"]
