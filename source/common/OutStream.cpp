// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "common/OutStream.hh"

#include "common/Exceptions.hh"

#include <fstream>
#include <iostream>


OutStream::
OutStream(const std::string& fileName) :
    _isInit(false),
    _fileName(fileName),
    _osptr(&std::cout),
    _ofsptr(new std::ofstream)
{
    if (! _fileName.empty())
    {
        std::ofstream test;
        openFile(_fileName,test);
    }
}



// required for unique_ptr:
OutStream::
~OutStream() {}



void
OutStream::
initStream()
{
    if (! _fileName.empty())
    {
        openFile(_fileName,*_ofsptr);
        _osptr=_ofsptr.get();
    }
    _isInit=true;
}

void
OutStream::
openFile(
    const std::string& filename,
    std::ofstream& ofs)
{
    ofs.open(filename.c_str());
    if (ofs) return;
    std::ostringstream oss;
    oss << "ERROR: Can't open output file: " << filename << "\n";
    BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
}
