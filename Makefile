DEP_EXECUTABLE := ${GOPATH}/bin/dep
GO_EXECUTABLE := go
VERSION := $(shell git describe --abbrev=10 --dirty --always --tags)
DIST_DIRS := find * -type d -exec
VERSION_PACKAGE := github.com/matsengrp/goliketree/cmd.Version
NAME := goliketree
PACKAGE:=github.com/matsengrp/goliketree

all: dep build test 

dep:
	${DEP_EXECUTABLE} ensure

build:
	${GO_EXECUTABLE} build -o ${NAME} -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}

install:
	rm -f ${GOPATH}/bin/${NAME}
	${GO_EXECUTABLE} install -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}

test:
	${GO_EXECUTABLE} test ${PACKAGE}/...
