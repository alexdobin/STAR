SGT
===

The files here contain a C++ class for implementing simple Good-Turing
re-estimation, as described by Geoff Sampson in the book Empirical Linguistics
(2001), and on the web at http://www.grsampson.net/RGoodTur.html. The code
here is a C++ adaptation of the published code by Sampson and Gale, with the
bug fix issued in 2000. It is encapsulated as a class to allow it to be
incorporated into other programs. An additional coding change is that the data
can be presented in any order, whereas the original code required the data to
be in ascending order.

Sampson's original code was issued with no restrictions on use. In keeping
with the spirit of this, the code here is issued under an open source licence
which allows essentially unrestricted use.

LICENCE
-------
Copyright (c) David Elworthy 2004.
All rights reserved.

Redistribution and use in source and binary forms for any purpose, with or
without modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions, and the following disclaimer.
 
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the disclaimer that follows 
   these conditions in the documentation and/or other materials 
   provided with the distribution.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact details
---------------
You may contact me at david@friendlymoose.com. I would be happy to hear of any
experiences you have with the code; please feel free to send me updated
versions. The reference site for the code is http://www.friendlymoose.com/.

Files and use
-------------
There are three files:
sgt.h       SGT header file
sgttest.cpp A test and example program

There is no source file, as the SGT class is a template over the observation
type, typically either an int or a double.

Information about using the class is included in the header file. The code has
been tested with g++ version 3.2 on cygwin and Microsoft Visual Studio version
6 on Windows 2000. You can compile and link the test program using g++ using
the command
     g++ -o sgttest sgttest.cpp

For Visual Studio, from the command line, you can compile and link with
     cl -GX sgttest.cpp

Version history
---------------
Initial version released January 2004.
Updated to a better implementation April 2004.
