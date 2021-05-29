MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY EXPRESSED 
OR IMPLIED. ANY USE IS AT YOUR OWN RISK. 

Using JANUS is free for non commercial applications 
(For commercial use, please contact the Matthias Bollhoefer
 m.bollhoefer@tu-bs.de).  

Please acknowledge, using the reference 

Matthias Bollhoefer, Olaf Schenk, Fabio Verbosio.
A High Performance Level-Block Approximate LU Factorization Preconditioner Algorithm.
Applied Numerical Mathematics 162:265-282, 2021.
DOI:10.1016/j.apnum.2020.12.023

For the software package you may use as footnote http://bilu.tu-bs.de

the contribution of this package in any scientific publication
dependent upon the use of the package. 
You shall use reasonable endeavors to notify the authors of 
the package of this publication. 



----------------------------------------------------------------------



You will find the documentation in JANUS/doc/ (a user guide giving
additional information to the ANM paper referenced above).

The user guide refers to several examples. The related sources are located
in JANUS/samples/



building the library
--------------------

To build the JANUS library, edit the "user.mk" file and adapt it to your
specific needs. If you are planning to use the MATLAB interface, 64 bit
long integer are required. Also, the GNU compilers are recommended for
interaction with MATLAB. Otherwise make your preferred choice.
After that, type "make".
JANUS will compute the library "libjanus.a" which will be located at
JANUS/lib/ArchComp/, "ArchComp" refers to your most favoured setting from
"user.mk".

Before you compile any sample code, you need to get the AMD library which is
part of SuiteSparse library of Tim Davis (Texas A&M UNiversity). Downloading and
compiling the SuiteSparse library will also provide libraries "libamd.a" and
"libsuitesparseconfig.a" (These are typically located in some subdirectory such
as SuiteSparse/AMD/Lib/ and SuiteSparse/SuiteSparse_config/). You could copy
these libraries (for simplicity) to JANUS/lib/ArchComp/ or refer differently
to them.

Furthermore you need your own BLAS and LAPACK library. Simple precompiled F77
sources are could be used, but these are likely to be slow. You better use some
vendor-specific libraries which highly optimized (if in doubt, ask your system
administrator for details).

Last but not least you need to get MC64 and MC21 from the HSL Mathematical
Software Library (https://www.hsl.rl.ac.uk/). For copyrights, terms of use, etc.,
we kindly refer to read carefully the their conditions and make sure that they
apply to you. You could compile these sources in the same way as the library has
been built including options for long integer, position independent code, memory
model, "-c" etc. Then copy the object files MC*.o to  JANUS/lib/ArchComp/
For example, if you are using the GNU compiler and you want to use 64 bit
long integer, then use a command such as
> gfortran -O -fPIC -m64 -fdefault-integer-8 -mcmodel=medium -c MC64D.f
This refers to the "GNU64_long" option

Other compilers and options are treated accordingly. You may take a look at
JANUS/src/makefiles for some suggestions.



compiling the main program
--------------------------

To compile any program using JANUS, for instance the examples in
JANUS/samples/ , go to the directory where the program is 
located and select from the "makefile" your most favoured sample
program. Please refer to the JANUS/samples/README for the various
sample programs and make your choice. After that type "make".

Of course, some of these sample progroams these are very simple examples whereas
other sample programs allow you to read in a matrix in Matrix-Market format and
try. 


Compiling the CMEX interfaces for MATLAB
----------------------------------------
When you have typed "make" in the JANUS/ main directory, CMEX executables will
be compiled as well. These executables along with the MATLAB entry points
will be found in JANUS/matlab. The main JANUS functions (beside several others)
are "janus.m" for the approximate block factorization and "janussolver.m" for
the embedded Krylov subspace solver. A README file as well as a MATLAB online
"help" function report about the details.


________________________________________________________________




However note that we can only distribute our own source codes. External
software is not included.




Here is a list of the following external software codes.


1. MC64
C COPYRIGHT (c) 1999 Council for the Central Laboratory
*                    of the Research Councils
CCCCC PACKAGE MC64A/AD
CCCCC AUTHORS Iain Duff (i.duff@rl.ac.uk) and Jacko Koster (jak@ii.uib.no)
CCCCC LAST UPDATE 20/09/99
CCCCC
C *** Conditions on external use ***
C
C The user shall acknowledge the contribution of this
C package in any publication of material dependent upon the use of
C the package. The user shall use reasonable endeavours to notify
C the authors of the package of this publication.
C
C The user can modify this code but, at no time
C shall the right or title to all or any part of this package pass
C to the user. The user shall make available free of charge
C to the authors for any purpose all information relating to any
C alteration or addition made to this package for the purposes of
C extending the capabilities or enhancing the performance of this
C package.
C
C The user shall not pass this code directly to a third party without
C the express prior consent of the authors.  Users wanting to licence
C their own copy of these routines should send email to hsl@aeat.co.uk
C
C None of the comments from the Copyright notice up to and including
C this one shall be removed or altered in any way.



2. MT-METIS
   MIT License (MIT)

   Copyright (c) 2013-2015, Regents of the University of Minnesota 

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

        -----------------------------------------------------------------

   This license only applies to the mt-Metis code, the software packages included
   with this distribution are distributed under their respective licenses:

   metis/LICENSE.txt
   bowstring/LICENSE
   domlib/LICENSE



3. AMD
   AMD is a set of routines for ordering a sparse matrix prior to Cholesky 
   factorization (or for LU factorization with diagonal pivoting). 

   Copyright (c) 2004-2006 by Timothy A. Davis, Patrick R. Amestoy, and
   Iain S. Duff. All Rights Reserved. Distributed under the GNU LGPL license. 



4. BLAS
   The reference BLAS is a freely-available software package. It is available
   from netlib via anonymous ftp and the World Wide Web. Thus, it can be 
   included in commercial software packages (and has been). We only ask that 
   proper credit be given to the authors. 

   Like all software, it is copyrighted. It is not trademarked, but we do ask
   the following: 

   If you modify the source for these routines we ask that you change the name
   of the routine and comment the changes made to the original. 

   The authors of BLAS  will gladly answer any questions regarding the 
   software. If a modification is done, however, it is the responsibility of
   the person who modified the routine to provide support. 



5. LAPACK
   The complete LAPACK package or individual routines from LAPACK are freely
   available on netlib and can be obtained via the World Wide Web or anonymous
   ftp. 

   The LAPACK homepage can be accessed on the World Wide Web via the URL 
   address:

             http://www.netlib.org/lapack/ 



