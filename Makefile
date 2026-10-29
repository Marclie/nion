# =============================================================================
# Makefile - Convenience Build & CI Test Runner for Nion
# =============================================================================
BUILD_DIR        ?= build
CMAKE_BUILD_TYPE ?= Release
TEST_LARGE_DEGREE ?= OFF
NPROCS           ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

.PHONY: all configure build test ci-test bench clean

all: build

configure:
	cmake -B $(BUILD_DIR) \
		-DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) \
		-DNION_BUILD_TESTS=ON \
		-DNION_BUILD_EXAMPLES=ON \
		-DTEST_LARGE_DEGREE=$(TEST_LARGE_DEGREE)

build: configure
	cmake --build $(BUILD_DIR) -j$(NPROCS)

test: build
	ctest --test-dir $(BUILD_DIR) --output-on-failure -V

ci-test: build
	./$(BUILD_DIR)/Tests/nion_ci_test

bench: build
	./$(BUILD_DIR)/Tests/nion

clean:
	rm -rf $(BUILD_DIR)
