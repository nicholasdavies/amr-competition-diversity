# Standard project makefile

# defines EXTDIR
include ./Makedefs

# OPTIONS

# Basic options: executable name, external modules, source directory, build directory
EXE := ./elegy
EXT := Randomizer/randomizer.cpp LuaFunc/luafunc.cpp Config/config.cpp
SRCDIR := .
OBJDIR := ./Build

# Compiler, language flags, preprocessor flags, compiler flags, linker flags
CPP := g++-9
LANG := -std=c++17
IFLAGS := -I . -I $(EXTDIR) -I /usr/local/include $(EXTRA_IFLAGS)
CFLAGS := -c -O3 -W -Wall -Werror -Wno-misleading-indentation -fopenmp
LFLAGS := -L/usr/local/lib -llua -lgsl -lgslcblas -lnlopt -fopenmp $(EXTRA_LFLAGS)

# Directories and files
SRC := $(wildcard $(SRCDIR)/*.cpp)
OBJ := $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o)))
DEPEND := $(OBJDIR)/depend.d
DEPSRC := $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/*.h) $(wildcard $(SRCDIR)/*.hpp)
SRCX := $(addprefix $(EXTDIR)/, $(EXT))
OBJX := $(addprefix $(OBJDIR)/,$(notdir $(SRCX:.cpp=.o)))
XRULES := $(OBJDIR)/xrules.make

# RULES

# Debugging, profiling, optimization, and cleaning
ifeq ($(DEBUG),1)
 CFLAGS := $(filter-out -s -fomit-frame-pointer -O3, $(CFLAGS)) -g -DDEBUG
endif

ifeq ($(PROFILE),1)
 CFLAGS := $(filter-out -s -fomit-frame-pointer, $(CFLAGS)) -pg
 LFLAGS := $(LFLAGS) -pg
endif

ifeq ($(NO_OPT),1)
 CFLAGS := $(filter-out -O3, $(CFLAGS))
endif

# Phony and default targets
.PHONY: clean veryclean so
default: $(DEPEND) $(XRULES) $(OBJX) $(OBJ) $(EXE)

# Dependencies and external module rules
$(DEPEND): $(DEPSRC)
	@$(CPP) -MM $(IFLAGS) $(LANG) $(SRC) $(SRCX) > $(DEPEND).x
	@sed -e "s/\([a-zA-Z0-9_]*\.o:\)/$(subst /,\/,$(OBJDIR)/)\1/" $(DEPEND).x > $(DEPEND)
	@$(RM) $(DEPEND).x

$(XRULES):
	@echo $$'$(foreach X,$(SRCX),$(OBJDIR)/$(notdir $(basename $(X))).o: $(X)\n\t$$(CPP) -o $$@ $$(IFLAGS) $$(CFLAGS) $$(LANG) $$<\n)' > $(XRULES)

# Include dependencies and external module compilation rules, unless we're cleaning
ifeq ($(MAKECMDGOALS),clean)
 CLEANING := 1
endif

ifeq ($(MAKECMDGOALS),veryclean)
 CLEANING := 1
endif

ifndef $(CLEANING)
 -include $(DEPEND)
 -include $(XRULES)
endif

# Compilation and linking
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CPP) -o $@ $(IFLAGS) $(CFLAGS) $(LANG) $<

$(EXE): $(OBJX) $(OBJ)
	$(CPP) -o $(EXE) $(OBJX) $(OBJ) $(LFLAGS)

# Cleaning
clean:
	$(RM) $(OBJX) $(OBJ) $(DEPEND) $(XRULES)

veryclean:
	$(RM) $(OBJX) $(OBJ) $(DEPEND) $(XRULES) $(EXE)

# Dynamic library
so: $(OBJX) $(OBJ)
	$(CPP) -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o $(EXE).so $(OBJ) $(OBJX) $(LFLAGS) -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
