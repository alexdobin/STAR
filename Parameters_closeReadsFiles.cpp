#include "Parameters.h"
#include "ErrorWarning.h"
#include <fstream>
#include <sys/stat.h>
void Parameters::closeReadsFiles() {
    for (uint imate=0;imate<readNmates;imate++) {//open readIn files
        if ( inOut->readIn[imate].is_open() ) inOut->readIn[imate].close();
        if (readFilesCommandPID[imate]>0) {
            kill(readFilesCommandPID[imate],SIGKILL);
        };
    };
};