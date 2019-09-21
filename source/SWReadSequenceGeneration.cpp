//
//  SWReadSequenceGeneration.cpp
//  
//
//  Created by Fahimeh Mirhaj on 6/18/19.
//

#include "SWReadSequenceGeneration.hpp"
#include "SmithWatermanAlignment.h"
#include <cstdlib> //rand()
#include <iostream>
#include <vector>
#include <string>
using namespace std;
vector<uint8> SmithWatermanAlignment::generateReadSeq(vector<uint8> initialSeq, scoreType &trueScore, uint &randomStartIndex, uint &readEndIndex) {
    
    uint errorRate = 10; // percent
    string operations;
    uint error;
    trueScore = 0;
    uint32 length = initialSeq.size();
    vector<uint8> read;
    //randomStartIndex = rand() % (length/3);
    randomStartIndex = 0;
    int randomModificationIndex = 0;
    //readEndIndex = 2 * (length / 3);
    readEndIndex = length;
    for (uint i = randomStartIndex; i < readEndIndex; i++) {
        error = rand() % 100;
        if(error < errorRate) { // do mutation, insertion, deletion randomly
            randomModificationIndex = rand() % 3;
            if(randomModificationIndex == 0) { // mutation
                int randomeGeneratedChar;
                while(randomeGeneratedChar = rand() % 4 == initialSeq[i]) {
                    
                }
                read.push_back(randomeGeneratedChar); // randomly pick 0, 1, 2, 3 for 'A', 'C', 'G', 'T' respectively
                operations += "m";
            }
            else if(randomModificationIndex == 1) { // insertion
                read.push_back(rand() % 4);
                operations += "I";
                i--;
            }
            //randomModificationIndex == 2, deletion, do nothing!
            else {
                operations += "D";
            }
            trueScore--;
        }
        else {
            read.push_back(initialSeq[i]);
            trueScore++;
            operations += "M";
        }
    }
//    cout << "string operations is:" << operations << endl;
    return read;
}
