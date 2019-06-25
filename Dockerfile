FROM rust:1.31

WORKDIR /usr/src/bcl2fastr

COPY . .

RUN cargo install --path .

CMD ["cargo", "run"]
