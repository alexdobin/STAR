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

#pragma once

#include "blt_util/blt_types.hh"

#include "blt_util/thirdparty_push.h"

#include "boost/serialization/access.hpp"
#include "boost/serialization/level.hpp"
#include "boost/serialization/nvp.hpp"

#include "blt_util/thirdparty_pop.h"

#include <cassert>
#include <cstring>
#include <iosfwd>
#include <string>


namespace PAIR_ORIENT
{

enum index_t { UNKNOWN, Fm, Fp, Rm, Rp, SIZE };

inline
const char*
label(const index_t i)
{
    switch (i)
    {
    case Fm:
        return "Fm";
    case Fp:
        return "Fp";
    case Rm:
        return "Rm";
    case Rp:
        return "Rp";
    default:
        return "UNKNOWN";
    }
}

inline
index_t
get_index(const pos_t pos1, const bool is_fwd_strand1,
          const pos_t pos2, const bool is_fwd_strand2)
{
    const bool is_read1_left(pos1 < pos2);

    if (is_fwd_strand1 != is_fwd_strand2)
    {
        // special-case very short fragments as innies:
        if (pos1 == pos2) return Rp;

        const bool left_strand(is_read1_left
                               ? is_fwd_strand1
                               : is_fwd_strand2);
        return (left_strand ? Rp : Rm);
    }
    else
    {
        return ((is_read1_left == is_fwd_strand1) ? Fp : Fm);
    }
}

/// inefficient label to id lookup, returns SIZE for unknown string:
inline
index_t
get_index(const char* str)
{
    for (int i(0); i<SIZE; ++i)
    {
        if (0==strcmp(str,label(static_cast<index_t>(i)))) return static_cast<index_t>(i);
    }
    return SIZE;
}
}


/// pair orientation status wrapper:
struct ReadPairOrient
{
    ReadPairOrient()
        : _val(PAIR_ORIENT::UNKNOWN)
    {}

    PAIR_ORIENT::index_t
    val() const
    {
        return _val;
    }

    void
    setVal(const unsigned newVal)
    {
        assert(newVal < PAIR_ORIENT::SIZE);
        _val=static_cast<PAIR_ORIENT::index_t>(newVal);
    }

private:
    PAIR_ORIENT::index_t _val;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned /*version*/)
    {
        std::string strval(PAIR_ORIENT::label(_val));
        ar& boost::serialization::make_nvp("pairOrientation", strval);
        _val = PAIR_ORIENT::get_index(strval.c_str());
    }
};

BOOST_CLASS_IMPLEMENTATION(ReadPairOrient, boost::serialization::object_serializable)


std::ostream&
operator<<(std::ostream& os, const ReadPairOrient& rpo);
