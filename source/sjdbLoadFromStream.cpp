#include "sjdbLoadFromStream.h"
void sjdbLoadFromStream(ifstream &sjdbStreamIn, SjdbClass &sjdbLoci) {
    while (sjdbStreamIn.good()) {
        string oneLine,chr1;
        uint u1,u2;
        char str1;
        getline(sjdbStreamIn,oneLine);
        istringstream oneLineStream (oneLine);
        oneLineStream >> chr1 >> u1 >> u2 >> str1;
        if (chr1!="") {
            sjdbLoci.chr.push_back(chr1);
            sjdbLoci.start.push_back(u1);
            sjdbLoci.end.push_back(u2);
            switch (str1) {//convert numbers to symbols
                case '1':
                case '+':
                    str1='+';
                    break;
                case '2':
                case '-':
                    str1='-';
                    break;
                default:
                    str1='.';
            };
            sjdbLoci.str.push_back(str1);
        };
    };
};