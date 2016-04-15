 
ISR Higgs lineshape calculator



Authors:

Stanislaw Jadach

Stanislaw.Jadach@cern.ch

The Henryk Niewodniczański Institute of Nuclear Physics,
Polish Academy of Sciences,
ul. Radzikowskiego 152, 
31-342 Kraków, Poland



Radoslaw A. Kycia

rkycia@pk.edu.pl

The Faculty of Physics, Mathematics and Computer Science
Cracow University of Technology
ul. Warszawska 24,
31-155 Kraków, Poland


Brief description:

Program calculates ISR infuence for the Higgs lineshape and its broadening by beam energy spread.
Description can be found in:

[1] S. Jadach, R.A. Kycia Software for calculations of the Higgs boson lineshape in future lepton colliders; arXiv:1601.02854 [hep-ph]

[2] S. Jadach, R.A. Kycia Lineshape of the Higgs boson in future lepton colliders; Physics Letters B, Volume 755, 10 April 2016, Pages 58–63; doi:10.1016/j.physletb.2016.01.065


Warinigs:

1. The programs were written by means of high standards, however, authors are not responsible for any damages that result from the use of these program. Use them at your own risk.

2. If you have any comments, suggestions, questions please write. We will try to answer them in reasonable time.



Dependence:

- GNU C++ compiler - for compilation
- GNU Make - handling compilation, running visualization. Treat as interpreter.
- CERN ROOT library
- Doxygen  - generation of documentation
- Valgrind - for debugging
- OpenMPI - for MPI programs

Content:
- ./Simple  - contains simplified version of program for educational purposes;
- ./Full    - contains full version of program that creates results published in  arXiv:1509.02406 [hep-ph];
- ./MPI     - contains MPI version of program from ./Full directory;
- ./MPI/sigEPlots - contains program that prepares plots of the cross section dependence on beam energy spread;
- ./MPI/Basic  - contains program that prepare allother plots of ./Full version of the program;


There are following options (type in terminal to see the output):

- make run - compile and run program. 
- make Generate-doc - generates documentation from the code in html and TeX formats using Doxygen.
- make clean - clean the directory from compilation and output files.


License:

This file is part of ISR Higgs lineshape calculator.

    It is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License any later version.

    ISR Higgs lineshape calculator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the ISR Higgs lineshape calculator. If not, see <http://www.gnu.org/licenses/>.


