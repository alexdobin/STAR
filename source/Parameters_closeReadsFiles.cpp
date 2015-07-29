#include "Parameters.h"
#include "ErrorWarning.h"
#include <fstream>
#include <sys/stat.h>
#ifdef _WIN32
#include <Windows.h>
#endif

void Parameters::closeReadsFiles() {
    for (uint imate=0;imate<readNmates;imate++) {//open readIn files
        if ( inOut->readIn[imate].is_open() ) inOut->readIn[imate].close();
        if (readFilesCommandPID[imate]>0) {

#ifdef _WIN32
			HANDLE hProcess = OpenProcess(SYNCHRONIZE | PROCESS_TERMINATE, TRUE, readFilesCommandPID[imate]);
			if (hProcess == NULL)
				return;
			BOOL result = TerminateProcess(hProcess, 0);
			CloseHandle(hProcess);
#else
			kill(readFilesCommandPID[imate],SIGKILL);
#endif

        };
    };
};