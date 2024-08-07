SHELL = /usr/bin/env bash
default_target: all

COMPILER := g++
CPPFLAGS := -O2 -std=c++11 
#WFLAGS   := -Wno-unused-command-line-argument -Wno-shift-count-overflow -Wno-unused-result
OBJECTS  := auxiliary.o gfmdSheet.o atomicSheet.o interSheet.o contMech.o

# Current date (for backups)
NOW      := $(shell date +%F)
FOLDER   := backup-$(NOW)



# ------ MinGW Installation Instructions ------ #
# (from http://www.fftw.org/install/windows.html)

# Download latest FFTW3 tar file from http://www.fftw.org/download.html
# Go to the Download folder, unpack, open a MinGW Terminal and execute:
# ./configure
# make
# make install



# ----- Determine correct FFTW paths ----- #

FFTW_PATH_Lnx := /usr/local
FFTW_PATH_Win := C:/MinGW/msys/1.0/local
FFTW_PATH_Vbx := ~/cp2k-master/tools/toolchain/install/fftw-3.3.8
FFTW_PATH_Ssh := ~/bin/fftw

ifneq "$(wildcard ${FFTW_PATH_Win} )" ""
  # If it exists, use the Windows path
  MSG = " Windows FFTW version"
  FFTW_PATH = $(FFTW_PATH_Win)
else ifneq "$(wildcard ${FFTW_PATH_Vbx} )" ""
  # If it exists, use the CP2K toolchain path
  MSG = " CP2K FFTW version"
  FFTW_PATH = $(FFTW_PATH_Vbx)
else ifneq "$(wildcard ${FFTW_PATH_Ssh} )" ""
  # If it exists, use the small cluster path
  MSG = " Cluster FFTW version"
  FFTW_PATH = $(FFTW_PATH_Ssh)
else
  # If none of those exist, use the default path
  MSG = " Linux/Mac default FFTW version"
  FFTW_PATH = $(FFTW_PATH_Lnx)
endif

FFTW_LIB = $(FFTW_PATH)/lib
FFTW_INC = $(FFTW_PATH)/include
LFFTW := -lfftw3 -L $(FFTW_LIB) -I $(FFTW_INC)



# ----- Compilation Rules ----- #

gap:
	@echo " Building readGap.exe..."
	@$(COMPILER) $(CPPFLAGS) readGap.cpp -o readGap.exe


surf:
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@echo " Building surf.exe..."
	@$(COMPILER) $(CPPFLAGS) statsSurf.cpp $(LFFTW) -o surf.exe


params:
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@echo " Building params.exe..."
	@$(COMPILER) $(CPPFLAGS) statsParams.cpp $(LFFTW) -o params.exe


range:
	@echo " Building range.exe..."
	@$(COMPILER) $(CPPFLAGS) extractRange.cpp -o range.exe


convert:
	@echo " Building convert.exe..."
	@$(COMPILER) $(CPPFLAGS) convertImg.cpp -o convert.exe


green:
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@echo " Building GF2D.exe..."
	@$(COMPILER) $(CPPFLAGS) GreensFunc2D.cpp $(LFFTW) -o GF2D.exe
	@echo " Building GF3D.exe..."
	@$(COMPILER) $(CPPFLAGS) GreensFunc3D.cpp $(LFFTW) -o GF3D.exe


all: surf params range convert

# ----- Cleanup and backup rules ----- #

clean:
	@rm -f *.o
	@rm -f *.exe
	@rm -f *.x
	@rm -f contMech.e


bak:
	@mkdir -p $(FOLDER)
	@cp *.cpp $(FOLDER)
	@cp *.h $(FOLDER)
	@#if [ -f "params.in" ]; then cp -f params.in $(FOLDER); fi
	@#if [ -f "equilPos0.in" ]; then cp -f equilPos0.in $(FOLDER); fi
	@#if [ -f "konfig0E.real" ]; then cp -f konfig0E.real $(FOLDER); fi
	@cp Makefile $(FOLDER)
	@echo " Source code backed up in $(FOLDER)."
