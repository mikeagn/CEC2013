Benchmark Functions for CEC'2013 Special Session and Competition on Niching Methods for Multimodal Function Optimization
========================================================================================================================
*Organizers : Xiaodong Li, Andries Engelbrecht, and Michael G. Epitropakis*

*Version    : 1.0.1*

*Developers : Michael G. Epitropakis and Xiaodong Li*

Please refer to:
- http://goanna.cs.rmit.edu.au/~xiaodong/cec13-niching/
- http://goanna.cs.rmit.edu.au/~xiaodong/cec13-niching/competition/


If you are using any material from the current competition please cite the
following work:

X. Li, A. Engelbrecht, and M.G. Epitropakis, ``Benchmark Functions for
CEC'2013 Special Session and Competition on Niching Methods for Multimodal
Function Optimization'', Technical Report, Evolutionary Computation and
Machine Learning Group, RMIT University, Australia, 2013.

--------------------------------------------------------------------------------
- About:
--------------------------------------------------------------------------------

In this folder you can find all necessary source files for the benchmark suite
of the CEC'2013 Special Session and Competition on Niching Methods for
Multimodal Function Optimization. 

The Test suite for the competition is implemented in MATLAB and C++.

--------------------------------------------------------------------------------
- Documentation:
--------------------------------------------------------------------------------

For more information please refer to the Technical Report of the Special
Session/Competition

--------------------------------------------------------------------------------
- Installation instructions:
--------------------------------------------------------------------------------
Please unpack the archive and extract its contents in a folder:

  unzip FILENAME.zip
  cd FILENAME

In the FILENAME folder you will find the following folder structure:

--------------------------------------------------------------------------------
- Directory Structure:
--------------------------------------------------------------------------------

After unpacking the archive file, you should end up with the following
structure:
<pre>

  ./			The MAIN directory, created when unpacking
   |
   +-- matlab		Source code of the benchmark functions in MATLAB
   |   |
   |   +- figs		Figures of the benchmark functions (for validation)
   |   |
   |   +- data		Data files for the benchmark suite 
   |
   +-- c++		Source code of the benchmark functions in C
   |   |
   |   +- plots		Figures of the benchmark function (for validation)
   |   |
   |   +- data		Data files for the benchmark suite 
   |
   +-- java 	Source code of the benchmark functions in JAVA
       |
       +- CEC2013	Project files
</pre>

--------------------------------------------------------------------------------
-- MATLAB: 
--------------------------------------------------------------------------------
The matlab folder contains: 
<pre>
+ niching_func.m
|-- 	The source code of the benchmark functions
|
+ demo_suite.m
|-- 	A demonstration file for using the competition source code 
|
+ count_goptima.m
|-- 	function for counting the number of global optima in a population 
|
+ plots.m
|-- 	MATLAB script for reproducing 1D/2D plots for the benchmark functions
|
+ data/
|-- 	Folder with data files
|
+ figs/
|-- 	Folder with the 1D/2D figures of the benchmark suite
|
+ get_fgoptima.m
|-- 	For each function, get the global's optima fitness value
|
+ get_no_goptima.m
|-- 	For each function, get the amount of known global optima 
|
+ get_rho.m
|-- 	For each function, get the rho value
|
+ get_lb.m
|-- 	For each function, get the lower bounds of the optimization box
|
+ get_ub.m
|-- 	For each function, get the upper bounds of the optimization box
</pre>

--------------------------------------------------------------------------------
-- C++:
--------------------------------------------------------------------------------
The c++ folder contains: 

<pre>
+ main.cpp:
|-- 	Examples on how to use each benchmark functions.
|
+ cec2013.h cec2013.cpp
|-- 	Interface and implementation of the CEC 2013 competition
|
+ cfunction.h cfunction.cpp
|-- 	The implementation of the benchmark function suite.
|
+ plots.cpp
|-- 	Source code to produce data files for plotting the benchmark functions
|
+ plots/
|-- 	It includes plotting data files, and the benchmark functions' figures. 
|
+ plots/plots.m
|-- 	MATLAB script for reproducing 1D/2D plots for the benchmark functions
|
+ Makefile: 
|-- 	A simple Makefile for easy compilation. 
|
+ data/
|-- 	Folder with data files
|
+ rand2.h rand2.c
|-- 	A simple random number generator 
|
+ randfile
|-- 	In Linux just a symbolic link on /dev/urandom 
|
+ sortidx.h
|-- 	Help header file for sorting
</pre>

--------------------------------------------------------------------------------
- Compilation 
--------------------------------------------------------------------------------
The Makefile provides the following three options:
make clean 	: Cleans all unnecessary files.  

make all 	: Compiles the main.cpp file. To see a demonstration of all
		  benchmark functions, execute ./a.out

make plots 	: Compiles the plots.cpp file. To produce the necessary data
		  files execute ./plots. It will produce 12 data files in 
		  the plots/ folder, one for each benchmark function. To 
		  produce the figures just run the plots.m MATLAB script.

--------------------------------------------------------------------------------
-- Java:
--------------------------------------------------------------------------------

The Java folder contains all necessary files to use the competition function
set. Please follow the demonstration file ExampleUsage.java on how to use the
source code. To reproduce the figures from the technical report please execute
the Plots.java file.

--------------------------------------------------------------------------------
-- ACKNOWLEDGEMENT: 
--------------------------------------------------------------------------------

- We really want to thank Dr. Jerry Swan for his fruitful comments and
  suggestions that helped to improve, the design and the implementation of the
  competition in the Java programming language. 
- We also thank Dr. Daniel Molina Cabrera, for finding and correcting a bug in
  the Java source code.

--------------------------------------------------------------------------------
-- LICENSE
--------------------------------------------------------------------------------
This folder should contain a file with the license statement (LICENSE.txt)
--------------------------------------------------------------------------------
Copyright 2013, Michael G. Epitropakis, Xiaodong Li, and Andries Engelbrecht. 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

(This is the so-called Simplified BSD License AKA FreeBSD License)
