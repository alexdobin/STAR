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

/**
 ** \brief Encapsulation of the concept of a read pair relative orientation.
 **
 ** Encapsulation of the concept of a read pair relative orientation.
 **
 ** \author Richard Shaw
 **/

#include "common/ReadPairOrient.hh"

#include <iostream>


std::ostream&
operator<<(std::ostream& os, const ReadPairOrient& rpo)
{
    os << PAIR_ORIENT::label(rpo.val());
    return os;
}
