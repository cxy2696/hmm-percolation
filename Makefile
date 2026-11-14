# Build: make          -> static libs + CLIs
#        make check    -> C core tests + C++ app tests
#        make reproduce-> full deterministic experiment set (scripts/reproduce.sh)
#        make clean

CC      ?= cc
CXX     ?= c++
PYTHON  ?= python3
JULIA   ?= julia
OPT     ?= -O2 -g
STD_C   ?= -std=c11
STD_CXX ?= -std=c++17
WARN    ?= -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude -Itests
CFLAGS   += $(STD_C) $(WARN) -fPIC -MMD -MP
CXXFLAGS += $(STD_CXX) $(WARN) -MMD -MP
LDLIBS   += -lm

BUILD := build
LIB_SRCS  := $(wildcard src/core/*.c) $(wildcard src/hmm/*.c) $(wildcard src/perc/*.c) \
             $(wildcard src/market/*.c) $(wildcard src/cross/*.c)
APP_SRCS  := $(wildcard src/app/*.cpp)
LIB_OBJS  := $(patsubst %.c,$(BUILD)/obj/%.o,$(LIB_SRCS))
APP_OBJS  := $(patsubst %.cpp,$(BUILD)/obj/%.o,$(APP_SRCS))
LIB       := $(BUILD)/libhp.a
APPLIB    := $(BUILD)/libhpapp.a

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHARED_LIB := $(BUILD)/libhp.dylib
SHARED_LDFLAGS := -dynamiclib
else
SHARED_LIB := $(BUILD)/libhp.so
SHARED_LDFLAGS := -shared
endif

APPS := $(BUILD)/hp_synth $(BUILD)/hp_regime $(BUILD)/hp_connectedness $(BUILD)/hp_perc_lab $(BUILD)/hp_strategies
EXAMPLES := $(BUILD)/ex_discrete_hmm $(BUILD)/ex_gaussian_hmm $(BUILD)/ex_percolation

TEST_C_OBJS := $(patsubst %.c,$(BUILD)/obj/%.o,$(wildcard tests/test_*.c))
TEST_APP_OBJS := $(patsubst %.cpp,$(BUILD)/obj/%.o,$(wildcard tests/test_app.cpp))

.PHONY: all check clean reproduce shared examples python-test python-example julia-test julia-example figures bindings-check
all: $(APPS) shared examples

$(BUILD)/obj/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(OPT) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/obj/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(OPT) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB): $(LIB_OBJS)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $^

$(SHARED_LIB): $(LIB_OBJS)
	@mkdir -p $(dir $@)
	$(CC) $(SHARED_LDFLAGS) -o $@ $^ $(LDLIBS)

shared: $(SHARED_LIB)

$(APPLIB): $(APP_OBJS) $(LIB)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $(APP_OBJS)

$(BUILD)/hp_%: apps/hp_%.cpp $(APPLIB) $(LIB)
	$(CXX) $(OPT) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(APPLIB) $(LIB) $(LDLIBS)

$(BUILD)/ex_%: examples/%.c $(LIB)
	$(CC) $(OPT) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(LIB) $(LDLIBS)

examples: $(EXAMPLES)

$(BUILD)/hp_tests: $(TEST_C_OBJS) $(LIB)
	$(CC) $(OPT) -o $@ $(TEST_C_OBJS) $(LIB) $(LDLIBS)

$(BUILD)/hp_tests_app: $(TEST_APP_OBJS) $(APPLIB) $(LIB)
	$(CXX) $(OPT) -o $@ $(TEST_APP_OBJS) $(APPLIB) $(LIB) $(LDLIBS)

check: $(BUILD)/hp_tests $(BUILD)/hp_tests_app
	./$(BUILD)/hp_tests
	./$(BUILD)/hp_tests_app

python-test: shared
	PYTHONPATH=bindings/python $(PYTHON) -m unittest discover -s bindings/python/tests -v

python-example: shared
	PYTHONPATH=bindings/python $(PYTHON) bindings/python/examples/regime_demo.py

julia-test: shared
	$(JULIA) --startup-file=no bindings/julia/test/runtests.jl

julia-example: shared
	$(JULIA) --startup-file=no bindings/julia/examples/regime_demo.jl

bindings-check: python-test

figures:
	MPLCONFIGDIR=build/mplcache build/venv/bin/python scripts/plot_results.py

reproduce: all
	./scripts/reproduce.sh

clean:
	rm -rf $(BUILD)

-include $(LIB_OBJS:.o=.d) $(APP_OBJS:.o=.d) $(TEST_C_OBJS:.o=.d) $(TEST_APP_OBJS:.o=.d)
