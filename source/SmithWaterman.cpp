//
//  SmithWaterman.cpp
//  
//
//  Created by Fahimeh Mirhaj on 6/10/19.
//

#include "SmithWaterman.hpp"
#include "GTF.h"
#include <string>
#include <cstdlib>
void GTF::computeScore(vector<uint8> read) {
    scoreType** scoringMatrix = new scoreType*[seqLenType];
    
    for(uint i = 0; i < scoringMatrix[0].size(); i++) {
        scoringMatrix[0][i] = 0;
    }
    for(uint j = 0; j < scoringMatrix.size(); j++) {
        scoringMatrix[j][0] = 0;
    }
    int w = -4;
    for(uint col = 1; col < scoringMatrix[0].size(); col++) {
        for(uint row = 1; row < scoringMatrix.size(); row++) {
            scoringMatrix[row][col] = max(0, max(max(scoring[row-1][col] + w, scoringMatrix[row][col-1] + w), scoringMatrix[row-1][col-1] + similarity(read[row], superTrancript[col])));
        }
    }
    
    // Find max score
    // what if more than 1 max score is present?
    int maxScore = INT_MIN;
    array<uint8, 2> maxScoreIndexes; // required for traceback step
    for(i = 0; i < scoringMatrix.size(); i++) {
        for(j = 0; j < scoringMatrix[i].size(); j++) {
            if(maxScore < scoringMatrix[i][j]) {
                maxScore = scoringMatrix[i][j];
                maxScoreIndexes[0] = i;
                maxScoreIndexes[1] = j;
            }
        }
    }
    printMatrix(scoringMatrix);
    // Traceback and find correspinding sequence
        // check the equality of the current score with score on diagon +/- the match/ mis match score or the gap penalty and trace back to the cell which has the founded score. until you reach a cell with score = 0
}

int similarity(uint8 a, uint8 b) {
    if(a == b) {
        return 5;
    }
    else {
        return -3;
    }
}

vector<uint8> generateReadSeq(vector<uint8> initialSeq, seqLenType length, uint errorRate) {
    int error;
    vector<uint8> read;
    int randomStartIndex = rand() % (1/3 * initialSeq.length());
    int j = 0; // index on read sequence
    int randomModificationIndex = 0;
    for (int i = randomStartIndex; i < (2/3 * length); i++) {
        if(randomModificationIndex == 1) {
            i--;
        }
        error = rand() % 100;
        if(error < errorRate) { // do mutation, insertion, deletion randomly
            randomModificationIndex = rand() % errorRate;
            if(randomModificationIndex == 0) { // mutation
                read[j] = rand() % 4; // randomly pick 0, 1, 2, 3 for 'A', 'C', 'G', 'T' respectively
                j++;
            }
            if(randomModificationIndex == 1) { // insertion
                read[j] = rand() % 4;
                j++;
            }
            // randomModificationIndex == 2, deletion, do nothing!
        }
        else {
            read[j] = initialSeq[i];
            j++;
        }
    }
    return read;
}
void printMatrix(vector<array<scoreType, seqLenType>> scoringMatrix) {
    for(i = 0; i < scoringMatrix.size(); i++) {
        for(j = 0; j < scoringMatrix[i].size(); j++) {
            if(j != scoringMatrix[i].size() - 1) {
                cout << scoringMatrix[i][j] << ", ";
            }
            else {
                cout << scoringMatrix[i][j];
            }
        }
        cout << endl;
    }
}
