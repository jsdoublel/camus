BINARY_NAME := camus
MAIN_GO_FILE := camus.go
VERSION := $(shell git describe --tags --always --dirty || echo "dev")
LDFLAGS := -ldflags="-X 'main.Version=$(VERSION)'"

TARGETS := \
	linux/amd64 \
	linux/arm64 \
	windows/amd64 \
	windows/arm64 \
	darwin/amd64 \
	darwin/arm64

build: | bin
	@echo "Building for $$(go env GOOS)/$$(go env GOARCH)..."
	go build $(LDFLAGS) -o bin/$(BINARY_NAME) $(MAIN_GO_FILE)

all: $(TARGETS)

$(TARGETS): | bin
	$(eval GOOS := $(word 1,$(subst /, ,$@)))
	$(eval GOARCH := $(word 2,$(subst /, ,$@)))
	$(eval EXT := $(if $(findstring windows,$(GOOS)),.exe,))
	@echo "Building for $(GOOS)/$(GOARCH)..."
	GOOS=$(GOOS) GOARCH=$(GOARCH) go build $(LDFLAGS) -o bin/$(BINARY_NAME)-$(GOOS)-$(GOARCH)$(EXT) $(MAIN_GO_FILE)

bin:
	@mkdir -p bin

clean:
	rm bin/*

.PHONY: all build clean $(TARGETS)
