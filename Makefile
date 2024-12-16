CC := g++ # This is the main compiler
# CC := clang++ --analyze # and comment out the linker last line for sanity
SRCDIR := src
BUILDDIR := build
TARGET1 := bin/runner1
TARGET2 := bin/runner2
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -O2 -std=c++17 # -Wall
LIB := -fopenmp #-pthread -lmongoclient -L lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt
INC := -I include

# Separate main source files for each target
MAIN1 := src/main1.cpp
MAIN2 := src/main2.cpp

# Object files for each target
OBJS1 := $(filter-out $(BUILDDIR)/main2.o,$(OBJECTS))
OBJS2 := $(filter-out $(BUILDDIR)/main1.o,$(OBJECTS))

# Rules for each target
$(TARGET1): $(OBJS1) $(BUILDDIR)/main1.o
	@echo " Linking $@..."
	@echo " $(CC) $^ -o $@ $(LIB)"; $(CC) $^ -o $@ $(LIB)

$(TARGET2): $(OBJS2) $(BUILDDIR)/main2.o
	@echo " Linking $@..."
	@echo " $(CC) $^ -o $@ $(LIB)"; $(CC) $^ -o $@ $(LIB)

# Compile object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) bin/*"; $(RM) -r $(BUILDDIR) bin/*

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

# Spikes
ticket:
	$(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean
