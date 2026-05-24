#!/bin/sh

# Requirements to install on Mac:
# Install from homebrew: htslib boost
# Make sure that your environment has Homebrew's include and library paths are on your C compiler's search path
# For example, in your .zshrc
# export LIBRARY_PATH="$LIBRARY_PATH:$(brew --prefix)/lib"
# export CPATH="$(brew --prefix)/include"

# We are using dynamically linked libhts, which makes building much easier. We don't need libcrypto, libssl, or even pthreads.
make -j8 HTSLIB_LIB="-lhts" BOOST_LIB_IO="-lboost_iostreams" BOOST_LIB_PO="-lboost_program_options" DYN_LIBS=""
