FROM rust:1.46

ARG dev=0

WORKDIR /usr/src/bcl2fastr

COPY . .

RUN cargo install --path .
RUN bash -c 'if [ ${dev} == 1 ]; then cargo test; fi'

CMD ["cargo", "run"]
