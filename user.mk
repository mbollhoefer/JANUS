# 1. step. 
# which platform do we use?
# ILUPACK comes along with some standard settings for various compilers and
# platforms. Often enough it suffices to choose one from the list below.
# In case you do not find an acceptable configuration, you can change the
# associated "makefile.include*" configuration in "makefiles/"
# The ILUPACK libraries and some dependencies will be provided in the
# associated sub directory of "libs/"

# GNU-compiler-based options

# 64 BIT gcc/gfortran linux system
# PLATFORM=GNU64

# 64 BIT gcc/gfortran linux system with 64 bit integer
PLATFORM=GNU64_long

# Intel-compiler-based options

# 64 BIT icc/ifort linux system
# PLATFORM=INTEL64

# 64 BIT icc/ifort linux system with 64 bit integer
# PLATFORM=INTEL64_long


# APPLE MAC OSX 64 bit Intel-based gcc/gfortran linux system
#PLATFORM=MACOSX64

# APPLE MAC OSX 64 bit Intel-based wiht 64-bit integer gcc/gfortran linux system
#PLATFORM=MACOSX64_long



