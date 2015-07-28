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

#include "boost/test/unit_test.hpp"

#include "ReadPairOrient.hh"


BOOST_AUTO_TEST_SUITE( test_ReadPairOrient )


BOOST_AUTO_TEST_CASE( test_PairTypes )
{
    // innies:
    {
        const pos_t readAPos(10);
        const bool readAFwd(true);
        const pos_t readBPos(20);
        const bool readBFwd(false);

        // map A->1 B->2
        PAIR_ORIENT::index_t res = PAIR_ORIENT::get_index(readAPos, readAFwd, readBPos, readBFwd);
        BOOST_REQUIRE_EQUAL(res, PAIR_ORIENT::Rp);

        // map A->2 B->1
        PAIR_ORIENT::index_t res2 = PAIR_ORIENT::get_index(readBPos, readBFwd, readAPos, readAFwd);
        BOOST_REQUIRE_EQUAL(res2, PAIR_ORIENT::Rp);
    }

    // outties
    {
        const pos_t readAPos(30);
        const bool readAFwd(true);
        const pos_t readBPos(20);
        const bool readBFwd(false);

        // map A->1 B->2
        PAIR_ORIENT::index_t res = PAIR_ORIENT::get_index(readAPos, readAFwd, readBPos, readBFwd);
        BOOST_REQUIRE_EQUAL(res, PAIR_ORIENT::Rm);

        // map A->2 B->1
        PAIR_ORIENT::index_t res2 = PAIR_ORIENT::get_index(readBPos, readBFwd, readAPos, readAFwd);
        BOOST_REQUIRE_EQUAL(res2, PAIR_ORIENT::Rm);
    }

    // short fragments should resolve to innies:
    {
        const pos_t readAPos(10);
        const bool readAFwd(true);
        const pos_t readBPos(10);
        const bool readBFwd(false);

        // map A->1 B->2
        PAIR_ORIENT::index_t res = PAIR_ORIENT::get_index(readAPos, readAFwd, readBPos, readBFwd);
        BOOST_REQUIRE_EQUAL(res, PAIR_ORIENT::Rp);

        // map A->2 B->1
        PAIR_ORIENT::index_t res2 = PAIR_ORIENT::get_index(readBPos, readBFwd, readAPos, readAFwd);
        BOOST_REQUIRE_EQUAL(res2, PAIR_ORIENT::Rp);
    }

}

BOOST_AUTO_TEST_SUITE_END()

