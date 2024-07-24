SHELL = /usr/bin/env bash
default_target: install



# ----- Detect FFTW path ----- #

FFTW_PATH_Lnx := /usr/local
FFTW_PATH_Win := C:/MinGW/msys/1.0/local
FFTW_PATH_Vbx := ~/cp2k-master/tools/toolchain/install/fftw-3.3.8
FFTW_PATH_Ssh := ~/bin/fftw

ifneq "$(wildcard ${FFTW_PATH_Win} )" ""
  # If it exists, use the Windows path
  MSG = " Windows FFTW version"
  FFTW_PATH ?= $(FFTW_PATH_Win)
else ifneq "$(wildcard ${FFTW_PATH_Vbx} )" ""
  # If it exists, use the CP2K toolchain path
  MSG = " CP2K FFTW version"
  FFTW_PATH ?= $(FFTW_PATH_Vbx)
else ifneq "$(wildcard ${FFTW_PATH_Ssh} )" ""
  # If it exists, use the small cluster path
  MSG = " Cluster FFTW version"
  FFTW_PATH ?= $(FFTW_PATH_Ssh)
else ifneq "$(wildcard ${FFTW_PATH_Lnx} )" ""
  # If none of those exist, use the default path
  MSG = " Linux/Mac default FFTW version"
  FFTW_PATH ?= $(FFTW_PATH_Lnx)
else
    $(error FFTW not found. Please use: make ... FFTW_PATH=/custom/fftw/path)
endif

FFTW_LIB ?= $(FFTW_PATH)/lib
FFTW_INC ?= $(FFTW_PATH)/include
LFFTW := -lfftw3 -L $(FFTW_LIB) -I $(FFTW_INC)



# ----- Detect platform ----- #

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    PLATFORM = UNIX
    RM = rm -f
    MKDIR = mkdir -p
    CP = cp
    INSTALL_DIR = /usr/local/bin
else ifeq ($(UNAME_S),Darwin)
    PLATFORM = UNIX
    RM = rm -f
    MKDIR = mkdir -p
    CP = cp
    INSTALL_DIR = /usr/local/bin
else ifeq ($(OS),Windows_NT)
    PLATFORM = WINDOWS
    RM = del
    MKDIR = mkdir
    CP = copy
    INSTALL_DIR = C:\Scripts
else
    $(error Unsupported platform)
endif



# ----- Detect C++ compiler ----- #

# If available, use g++, otherwise, use default.
GCC_AVAILABLE := $(shell command -v g++ >/dev/null 2>&1 && echo yes || echo no)
ifeq ($(GCC_AVAILABLE),yes)
    COMPILER ?= g++
else
    COMPILER ?= $(CXX)
    COMPILER ?= $(CPP)
endif
CPPFLAGS := -O2 -std=c++11
OBJECTS  := auxiliary.o gfmdSheet.o atomicSheet.o interSheet.o contMech.o

# Warning flags can be changed by user
WFLAGS   ?= -Wno-unused-command-line-argument

# Current date (for backups)
NOW      := $(shell date +%F)
FOLDER   := backup-$(NOW)



# ----- Compilation Rules ----- #

gapconvert.exe: header.h gapconvert.cpp 
	@echo " Building gapconvert.exe..."
	@$(COMPILER) $(CPPFLAGS) $(WFLAGS) gapconvert.cpp -o gapconvert.exe

surf.exe: header.h statsSurf.cpp
	@echo " Building surf.exe..."
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@$(COMPILER) $(CPPFLAGS) $(WFLAGS) statsSurf.cpp $(LFFTW) -o surf.exe

params.exe: header.h statsParams.cpp
	@echo " Building params.exe..."
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@$(COMPILER) $(CPPFLAGS) $(WFLAGS) statsParams.cpp $(LFFTW) -o params.exe

range.exe: header.h extractRange.cpp
	@echo " Building range.exe..."
	@$(COMPILER) $(CPPFLAGS) $(WFLAGS) extractRange.cpp -o range.exe

convert.exe: header.h convertImg.cpp
	@echo " Building convert.exe..."
	@$(COMPILER) $(CPPFLAGS) $(WFLAGS) convertImg.cpp -o convert.exe

green: header.h GreensFunc2D.cpp GreensFunc3D.cpp
	@echo " Building GF2D.exe..."
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@$(COMPILER) $(CPPFLAGS) $(WFLAGS) GreensFunc2D.cpp $(LFFTW) -o GF2D.exe
	@echo " Building GF3D.exe..."
	@$(COMPILER) $(CPPFLAGS) $(WFLAGS) GreensFunc3D.cpp $(LFFTW) -o GF3D.exe

all: surf.exe params.exe range.exe convert.exe gapconvert.exe

# ----- Install Scripts for Running Executables from Anywhere ----- #

install_scripts_UNIX: scripts
	@$(MKDIR) $(INSTALL_DIR)
	@$(CP) cm_convert.sh $(INSTALL_DIR)/cm_convert
	@$(CP) cm_range.sh $(INSTALL_DIR)/cm_range
	@$(CP) cm_params.sh $(INSTALL_DIR)/cm_params
	@$(CP) cm_surf.sh $(INSTALL_DIR)/cm_surf
	@$(CP) cm_gapconvert.sh $(INSTALL_DIR)/cm_gapconvert
	@$(CP) cmjobs.sh $(INSTALL_DIR)/cmjobs
	@chmod +x $(INSTALL_DIR)/cm_convert
	@chmod +x $(INSTALL_DIR)/cm_range
	@chmod +x $(INSTALL_DIR)/cm_params
	@chmod +x $(INSTALL_DIR)/cm_surf
	@chmod +x $(INSTALL_DIR)/cm_gapconvert
	@chmod +x $(INSTALL_DIR)/cmjobs
	@echo " Installed scripts to $(INSTALL_DIR)."

install_scripts_WINDOWS: scripts
	@if not exist "$(INSTALL_DIR)" mkdir "$(INSTALL_DIR)"
	@$(CP) cm_convert.bat "$(INSTALL_DIR)\cm_convert.bat"
	@$(CP) cm_range.bat "$(INSTALL_DIR)\cm_range.bat"
	@$(CP) cm_params.bat "$(INSTALL_DIR)\cm_params.bat"
	@$(CP) cm_surf.bat "$(INSTALL_DIR)\cm_surf.bat"
	@$(CP) cm_gapconvert.bat "$(INSTALL_DIR)\cm_gapconvert.bat"
	@$(CP) cmjobs.bat "$(INSTALL_DIR)\cmjobs.bat"
	@echo " Installed bat files to $(INSTALL_DIR)."

scripts: convert.exe range.exe params.exe surf.exe gapconvert.exe
ifeq ($(PLATFORM),WINDOWS)
	@echo " Creating bat files..."
	@echo @echo off > cm_convert.bat
	@echo @echo off > cm_range.bat
	@echo @echo off > cm_params.bat
	@echo @echo off > cm_surf.bat
	@echo @echo off > cm_gapconvert.bat
	@echo @echo off > cmjobs.bat
	@echo $(CURDIR)\\convert.exe %* >> cm_convert.bat
	@echo $(CURDIR)\\range.exe %* >> cm_range.bat
	@echo $(CURDIR)\\params.exe %* >> cm_params.bat
	@echo $(CURDIR)\\surf.exe %* >> cm_surf.bat
	@echo $(CURDIR)\\gapconvert.exe %* >> cm_gapconvert.bat
	@echo python3 $(CURDIR)\\cmjobs.py %* >> cmjobs.bat
else
	@echo " Creating shell scripts..."
	@echo "#!/bin/bash" > cm_convert.sh
	@echo "#!/bin/bash" > cm_range.sh
	@echo "#!/bin/bash" > cm_params.sh
	@echo "#!/bin/bash" > cm_surf.sh
	@echo "#!/bin/bash" > cm_gapconvert.sh
	@echo "#!/bin/bash" > cmjobs.sh
	@echo $(CURDIR)/convert.exe \"\$$@\" >> cm_convert.sh
	@echo $(CURDIR)/range.exe \"\$$@\" >> cm_range.sh
	@echo $(CURDIR)/params.exe \"\$$@\" >> cm_params.sh
	@echo $(CURDIR)/surf.exe \"\$$@\" >> cm_surf.sh
	@echo $(CURDIR)/gapconvert.exe \"\$$@\" >> cm_gapconvert.sh
	@echo python3 $(CURDIR)/cmjobs.py \"\$$@\" >> cmjobs.sh
endif

install: install_scripts_$(PLATFORM)



# ----- Cleanup and backup rules ----- #

clean:
	@$(RM) *.o
	@$(RM) *.exe
	@$(RM) *.x
	@$(RM) contMech.e

bak:
	@$(MKDIR) $(FOLDER)
	@$(CP) *.cpp $(FOLDER)
	@$(CP) *.h $(FOLDER)
	@#if [ -f "params.in" ]; then $(CP) params.in $(FOLDER); fi
	@#if [ -f "equilPos0.in" ]; then $(CP) equilPos0.in $(FOLDER); fi
	@#if [ -f "konfig0E.real" ]; then $(CP) konfig0E.real $(FOLDER); fi
	@$(CP) Makefile $(FOLDER)
	@echo " Source code backed up in $(FOLDER)."
