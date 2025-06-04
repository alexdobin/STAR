#include "Parameters.h"
#include "ErrorWarning.h"
#include <fstream>
#include <sys/stat.h>
#include <sys/wait.h>
void Parameters::closeReadsFiles() {
    for (uint imate=0; imate<readFilesIn.size(); imate++) {//open readIn files
        if ( inOut->readIn[imate].is_open() )
            inOut->readIn[imate].close();
        if (readFilesCommandPID[imate]>0) {
            kill(readFilesCommandPID[imate],SIGTERM);  // Use SIGTERM for graceful shutdown
            int status;
            waitpid(readFilesCommandPID[imate], &status, 0);  // Wait for process to terminate
            
            // Check if process had issues (but don't fail - this is cleanup)
            if (WIFEXITED(status)) {
                int exit_code = WEXITSTATUS(status);
                if (exit_code != 0) {
                    warningMessage("readFilesCommand process for mate " + to_string(imate+1) + 
                                 " exited with non-zero code: " + to_string(exit_code), 
                                 std::cerr, inOut->logMain, *this);
                }
            } else if (WIFSIGNALED(status)) {
                int signal_num = WTERMSIG(status);
                if (signal_num != SIGTERM) {  // Don't warn about our own SIGTERM
                    warningMessage("readFilesCommand process for mate " + to_string(imate+1) + 
                                 " was terminated by signal: " + to_string(signal_num), 
                                 std::cerr, inOut->logMain, *this);
                }
            }
        }
    };
};