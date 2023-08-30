all: v0.10.7-beta

v0.10.7-beta:
	docker build \
		-t entap/entap:v0.10.7-beta \
		--build-arg ENTAP_VERSION=v0.10.7-beta \
		.

local:
	docker build \
		-f ./docker/Dockerfile \
		-t entap/entap:local \
		--build-arg ENTAP_VERSION=local \
		.

push:
	docker push entap/entap:v0.10.7-beta

pull:
	docker pull entap/entap:v0.10.7-beta

clean:
	docker image rm -f entap/entap:v0.10.7-beta
	docker image rm -f entap/entap:local
