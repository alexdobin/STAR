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

#pragma once

#include <iosfwd>
#include <memory>
#include <string>



/// provide an output stream which comes from either a file or a tty default
///
struct OutStream
{
    OutStream(const std::string& fileName);

    ~OutStream();

    std::ostream&
    getStream()
    {
        if (! _isInit) initStream();
        return *_osptr;
    }

private:

    void
    initStream();

    static
    void
    openFile(
        const std::string& filename,
        std::ofstream& ofs);

    bool _isInit;
    std::string _fileName;
    std::ostream* _osptr;
    std::unique_ptr<std::ofstream> _ofsptr;
};
