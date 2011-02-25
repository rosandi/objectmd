#MPI_HOME = $(HOME)/mpi-$(shell uname -m)

# the environment variable $OMD_HOME (see $HOME/.omdconf*) has highest priority
# to change remove this variable from the shell
INSTALL_DIR = $(HOME)/omd-$(shell uname -m)

M_OMD_CXX = $(shell if [ -z "$(OMD_CXX)" ]; then echo g++; else echo $(OMD_CXX); fi)
