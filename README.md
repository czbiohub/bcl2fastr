### Run Instructions

	`docker build --tag=bcl2fastr-dev . # in bcl2fastr top directory`
	
	for developement:
	`docker run -it -v <your local bcl2fastr path>:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev`
	
	for testing:
	`docker run -it -v <your local bcl2fastr path>:/usr/src/bcl2fastr --rm --name bcl2fastr-dev bcl2fastr-dev cargo test`
	
	for production:
	`docker run -it --rm --name bcl2fastr-dev bcl2fastr-dev`
