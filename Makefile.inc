# PYTHIA configuration file.
# Generated on Fri Apr 20 16:03:34 EDT 2018 with the user supplied options:

# Install directory prefixes.
PREFIX_BIN=/usr/local/bin
PREFIX_INCLUDE=/usr/local/include
PREFIX_LIB=/usr/local/lib
PREFIX_SHARE=/usr/local/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=
CXX_SHARED=
CXX_SONAME=-Wl,-install_name,@rpath/
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=./
EVTGEN_LIB=./

# FASTJET3 configuration.
FASTJET3_USE=true
FASTJET3_BIN=/usr/local/bin/
FASTJET3_INCLUDE=/usr/local/include/fastjet
FASTJET3_LIB=/usr/local/lib

# HEPMC2 configuration.
HEPMC2_USE=false
HEPMC2_BIN=
HEPMC2_INCLUDE=./
HEPMC2_LIB=./

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=./
LHAPDF5_LIB=./

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=./
LHAPDF6_LIB=./

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE=true
ROOT_BIN=${ROOTSYS}/bin/
ROOT_INCLUDE=${ROOTSYS}/include/
ROOT_LIB=${ROOTSYS}/lib/

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=./
BOOST_LIB=./

# PYTHON configuration.
PYTHON_USE=true
PYTHON_BIN=/usr/bin/
PYTHON_INCLUDE=/usr/include/python3.8
PYTHON_LIB=/usr/lib



# HOTVR configuration.
HOTVR_USE=true
HOTVR_BIN=./HOTVRContrib
HOTVR_INCLUDE=./HOTVRContrib
HOTVR_LIB=./HOTVRContrib
