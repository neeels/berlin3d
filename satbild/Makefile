####################################################
# dexample1
####################################################

convert-all:
	./ecw2jpg.sh -c src-ecw/*.ecw

erdas = erdas-sdk/Desktop_Read-Only

#include $(erdas)/examples/ecwjp2make.opt
CONFIG_TARGET=release
#CONFIG_TARGET=debug

UNICODE_DEFINES=-DUNICODE -D_UNICODE

EXTRA_CCFLAGS=$(UNICODE_DEFINES)

ifeq ($(CONFIG_TARGET),debug)
SDK_SUFFIX=d
CONFIG_DEFINES=-DDEBUG
TBB_SUFFIX=_debug
COMPILER_OPTIMISATION=
DEBUG_LEVEL=-g
else
SDK_SUFFIX=
CONFIG_DEFINES=-DNDEBUG
TBB_SUFFIX=
COMPILER_OPTIMISATION=-O3
DEBUG_LEVEL=
endif



#PLATFORM_CONF=x86
PLATFORM_CONF=x64

ifeq ($(PLATFORM_CONF),x86)
PLATFORM_DEFINES=-DX86 -D_X86 -D_X86=1 -m32
else
PLATFORM_DEFINES=-Dx64linux -m64
endif

# Source Dirs
SRC_FILES 	= ecw2raw.c
INCLUDE_DIR	:= $(erdas)/include

OBJ 		= ecw2raw.o
OUT_DIR		:= $(erdas)/bin/$(PLATFORM_CONF)/$(CONFIG_TARGET)
LIB_DIR		:= $(erdas)/lib/$(PLATFORM_CONF)/$(CONFIG_TARGET)

# Set Includes
CPP_FLAGS	= -I . \
		  -I$(INCLUDE_DIR)

# Set Libraries
LD_FLAGS	= -L$(LIB_DIR) 

# Compiler options
CXX 		= g++
CC 		= $(CXX)

CXXFLAGS	= $(DEBUG_LEVEL) $(COMPILER_OPTIMISATION) $(PLATFORM_DEFINES) $(EXTRA_CCFLAGS) $(CONFIG_DEFINES) $(CPP_FLAGS) $(LD_FLAGS)
CFLAGS		= $(CXXFLAGS)

# Set Dependencies
LDLIBS		= -lNCSEcw$(SDK_SUFFIX)


# Set the object outputs
O_FILES		= $(SRC_FILES:%.c=%.o)
O_FILES_OUT	=$(addprefix $(OUT_DIR)/,dexample1)

####################################################
# Additional make rules
####################################################


def: ecw2raw
#LD_LIBRARY_PATH=$(LIB_DIR) ./ecw2raw src-ecw/dop20_404_5810.ecw

erdas-sdk:
	@echo
	@echo "  * You need the ERDAS SDK from hexagon."
	@echo "  * http://go.hexagongeospatial.com/ERDASECW/JP2SDKv5.2.1%28Linux%29"
	@echo "  *"
	@echo "  * Place the ERDAS-yada-yada dir here and rename/link as erdas-sdk,"
	@echo "  * so that $(erdas) exists."
	@echo
	@false

depend:
	makedepend -- $(CXXFLAGS) -- -Y $(SRC_FILES)

clean:
	$(RM) $(O_FILES) $(O_FILES_OUT) core *.rpo

clean-all: clean
	$(RM) $(EXEC)

ecw2raw: ecw2raw.c erdas-sdk
	g++ -g -I. -I$(erdas)/include -o ecw2raw ecw2raw.c -L$(LIB_DIR) -lNCSEcw



# DO NOT DELETE
