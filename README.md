### Run Instructions

[![Build Status](https://travis-ci.org/czbiohub/bcl2fastr.svg?branch=master)](https://travis-ci.org/czbiohub/bcl2fastr)
[![Code Coverage](https://codecov.io/gh/czbiohub/bcl2fastr/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/bcl2fastr)

 - Build the Docker image:

   `docker build --tag=bcl2fastr-dev . # in bcl2fastr top directory`

 - Build the docs:

  The first command builds the documentation for the dependencies and only needs to be run occasionally. The second builds richer documentation for `bcl2fastr`.

   ```
   docker run -v `pwd`:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev cargo doc
   docker run -v `pwd`:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev cargo doc --no-deps --document-private-items
   ```

 - Develop:

   `docker run -it -v <your local bcl2fastr path>:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev`

 - Run tests:

   `docker run -it -v <your local bcl2fastr path>:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev cargo test`

 - Production (???):

   `docker run -it --rm --name bcl2fastr-dev bcl2fastr-dev`
