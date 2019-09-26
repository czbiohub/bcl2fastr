### Run Instructions

[![Build Status](https://travis-ci.org/czbiohub/bcl2fastr.svg?branch=master)](https://travis-ci.org/czbiohub/bcl2fastr)
[![Code Coverage](https://codecov.io/gh/czbiohub/bcl2fastr/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/bcl2fastr)

 - Build the Docker image:

   `docker build --tag=bcl2fastr-dev . # in bcl2fastr top directory`

 - Build the docs:

   ```
   docker run -v `pwd`:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev cargo doc --document-private-items
   ```

 - Develop:

   `docker run -it -v <your local bcl2fastr path>:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev`

 - Run tests:

   `docker run -it -v <your local bcl2fastr path>:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev cargo test`

 - Production (???):

   `docker run -it --rm --name bcl2fastr-dev bcl2fastr-dev`
