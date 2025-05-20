.PHONEY: build test

build:
	go build -o bin/camus ./cmd/camus/

test:
	go test -v ./...
