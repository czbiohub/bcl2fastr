FROM rust:1.37

ARG dev=0

WORKDIR /usr/src/bcl2fastr

COPY . .

RUN cargo install --path .
RUN bash -c 'if [ ${dev} == 1 ]; then cargo install assert_cli; fi'

CMD ["cargo", "run"]
