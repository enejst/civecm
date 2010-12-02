James G. MacKinnon, Alfred Haug, and Leo Michelis, "Numerical
Distribution Functions of Likelihood Ratio Tests for Cointegration",
Journal of Applied Econometrics, Vol. 14, No. 5, 1999, pp. 563-577.

The RePEc link to this paper is wrong. The correct link is:

  http://www3.interscience.wiley.com/journal/66001312/abstract

There are three zip files: mhm-tabs.zip, mhm-progs-unix.zip, and
mhm-progs-win.zip.

The file mhm-tabs.zip contains all the tables needed by the programs that
compute critical values and P values. On a Unix system, these should, if
possible, be put in the directory /usr/local/urcdist. On a DOS, OS/2, or
Windows system, they should be put in the directory \URCDIST. The tables
files were in DOS format (with CR/LF at line ends) when they were zipped. To
unzip them on a Unix system, use "unzip -a mhm-tabs.zip". This will convert
them to Unix format. No conversion should be necessary on a DOS, OS/2, or
Windows system.

The file mhm-progs-unix.zip contains three files of Fortran source code:
lrcdist.f, lrcrouts.f, and johrouts.f. It should be possible to compile them
with any reasonably modern F77 or F90 compiler. The file lrcdist.f contains
a program that is intended to be used interactively. The other two files
contain subroutines that are intended to be called by other programs.

The file mhm-progs-win.zip contains two files, lrcdist.for, which is a
slightly modified version of lrcdist.f, and lrcdist.exe, which is an
executable for Windows that must be run from the command line.
It was compiled using CVF 6.6 on a machine running Windows 2000. The
default location for the .tab files is \urcdist.

Note: The errors that were discovered on 2001-1-08 and 2003-5-05 have
now been corrected in all source files. These errors caused minor
discrepancies in P values and critical values in some cases. The second
error has not been corrected in the Windows executable, because I do not
use Windows and have no access to a Fortran compiler that runs on
Windows.

James G. MacKinnon                       Department of Economics
    phone: 613 533-2293                  Queen's University
      Fax: 613 533-2257                  Kingston, Ontario, Canada
    Email: jgm <AT> qed.econ.queensu.ca  K7L 3N6
Home Page: http://www.econ.queensu.ca/pub/faculty/mackinnon/
