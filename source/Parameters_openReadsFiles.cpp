#ifdef _WIN32
#include <Windows.h>
#include <map>
#endif

#include <fstream>
#include <sys/stat.h>

#include "Parameters.h"
#include "ErrorWarning.h"

#ifdef _WIN32

typedef std::map<int ,HANDLE> MapHandles; 

typedef enum 
{
	Read = 0,
	Write = 1
} HandleType;

bool CreateReadWritePipe(MapHandles& handles)
{
	SECURITY_ATTRIBUTES saAttr = {0};

	// Set the bInheritHandle flag so pipe handles are inherited. 
	saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	saAttr.bInheritHandle = TRUE;
	saAttr.lpSecurityDescriptor = NULL;

	HANDLE hChildStd_OUT_Rd = NULL;
	HANDLE hChildStd_OUT_Wr = NULL;

	// Create a pipe for the child process's STDOUT. 
	if (!CreatePipe(&hChildStd_OUT_Rd, &hChildStd_OUT_Wr, &saAttr, 0))
	{
		return false; 
	}

	// Ensure the read handle to the pipe for STDOUT is not inherited.
	if (!SetHandleInformation(hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0))
	{
		CloseHandle(hChildStd_OUT_Rd); 
		CloseHandle(hChildStd_OUT_Wr);
		return false; 

	}
	handles[HandleType::Read] = hChildStd_OUT_Rd;
	handles[HandleType::Write] = hChildStd_OUT_Wr;
	return true; 
}

DWORD CreateChildProcess(HANDLE hPipeWrite, const std::string& filenames, const std::string& command)
{
	std::string cmdLine = "STAR_ReadFile"; 
	cmdLine.append(" "); 
	cmdLine.append(filenames);
	cmdLine.append(" ");
	cmdLine.append(command); 

	char*  szCmdline = new char[cmdLine.size() + 1]; 
	strcpy(szCmdline, cmdLine.c_str()); 

	PROCESS_INFORMATION piProcInfo;
	STARTUPINFO siStartInfo;
	BOOL bSuccess = FALSE;

	// Set up members of the PROCESS_INFORMATION structure. 
	ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));

	// Set up members of the STARTUPINFO structure. 
	// This structure specifies the STDIN and STDOUT handles for redirection.
	ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
	siStartInfo.cb = sizeof(STARTUPINFO);
	siStartInfo.hStdError = hPipeWrite;
	siStartInfo.hStdOutput = hPipeWrite;
	//siStartInfo.hStdInput = g_hChildStd_IN_Rd;
	siStartInfo.dwFlags |= STARTF_USESTDHANDLES;

	// Create the child process. 
	bSuccess = CreateProcess(NULL,
		szCmdline,     // command line 
		NULL,          // process security attributes 
		NULL,          // primary thread security attributes 
		TRUE,          // handles are inherited 
		0,             // creation flags 
		NULL,          // use parent's environment 
		NULL,          // use parent's current directory 
		&siStartInfo,  // STARTUPINFO pointer 
		&piProcInfo);  // receives PROCESS_INFORMATION 

	// If an error occurs, exit the application. 
	if (!bSuccess)
		return -1;
	else
	{
		// Close handles to the child process and its primary thread.
		// Some applications might keep these handles to monitor the status
		// of the child process, for example. 

		CloseHandle(piProcInfo.hProcess);
		CloseHandle(piProcInfo.hThread);
	}
	// This will make sure the read on pipe not block when 
	// writting to pipe finished by child process and it already exited.
	CloseHandle(hPipeWrite); 
	return piProcInfo.dwProcessId; 
}

void Parameters::openReadsFiles() 
{
	string readFilesCommandString("");
	if (readFilesCommand.at(0) == "-") 
	{
		// TODO
		if (readFilesIn.at(0).find(',') < readFilesIn.at(0).size()) 
			readFilesCommandString = "cat";//concatenate multiple files
	}
	else 
	{
		// TODO
		for (uint ii = 0; ii < readFilesCommand.size(); ii++) 
			readFilesCommandString += readFilesCommand.at(ii) + "   ";
	};

	if (readFilesCommandString == "") 
	{
		//read from file
		for (uint ii = 0; ii < readNmates; ii++) 
		{
			//open readIn files
			readFilesCommandPID[ii] = 0;//no command process IDs
			if (inOut->readIn[ii].is_open())
				inOut->readIn[ii].close();
			inOut->readIn[ii].open(readFilesIn.at(ii).c_str()); //try to open the Sequences file right away, exit if failed
			if (inOut->readIn[ii].fail()) 
			{
				ostringstream errOut;
				errOut << "EXITING because of fatal input ERROR: could not open readFilesIn=" << readFilesIn.at(ii) << "\n";
				exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
			};
		};
	}
	else 
	{
		vector<string> readsCommandFileName;
		readFilesNames.resize(readNmates);

		for (uint imate = 0; imate < readNmates; imate++) 
		{
			MapHandles handles; 
			if (!CreateReadWritePipe(handles))
			{
				ostringstream errOut;
				errOut << "EXITING: because of fatal EXECUTION error: Failed CreatePipe\n";
				errOut << errno << ": " << strerror(errno) << "\n";
				errOut << "GetLastError()" << ": " << GetLastError() << "\n";
				exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
				break;
			}

			inOut->logMain << "\n   Input read files for mate " << imate + 1 << ", from input string " << readFilesIn.at(imate) << endl;

			string readFilesInString(readFilesIn.at(imate));
			size_t pos = 0;
			readFilesN = 0;
			do 
			{
				//cycle over multiple files separated by comma
				pos = readFilesInString.find(',');
				string file1 = readFilesInString.substr(0, pos);
				readFilesInString.erase(0, pos + 1);
				readFilesNames.at(imate).push_back(file1);

				++readFilesN;//only increase file count for one mate

			} while (pos != string::npos);

			readFilesCommandPID[imate] = 0;

			ostringstream errOut;
			
			DWORD pid = CreateChildProcess(handles[HandleType::Write], readFilesInString, readFilesCommandString); 
			if (pid > 0)
			{
				if (inOut->readIn[imate].open_pipe_read(handles[HandleType::Read]))
				{
					readFilesCommandPID[imate] = pid;
				}
				else
				{
					errOut << "EXITING: because of fatal EXECUTION error: Failed open_pipe_read\n";
					errOut << errno << ": " << strerror(errno) << "\n";
					exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
					break;
				}

			}
			else
			{
				errOut << "EXITING: because of fatal EXECUTION error: Failed CreateChildProcess\n";
				errOut << errno << ": " << strerror(errno) << "\n";
				errOut << "GetLastError()" << ": " << GetLastError() << "\n";
				exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
				break;
			}
		};
		if (readNmates == 2 && readFilesNames.at(0).size() != readFilesNames.at(1).size()) 
		{
			ostringstream errOut;
			errOut << "EXITING: because of fatal INPUT ERROR: number of input files for mate1: " << readFilesNames.at(0).size() << " is not equal to that for mate2: " << readFilesNames.at(1).size() << "\n";
			errOut << "Make sure that the number of files in --readFilesIn is the same for both mates\n";
			exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
		};

		if (outSAMattrRG.size()>1 && outSAMattrRG.size() != readFilesN) 
		{
			ostringstream errOut;
			errOut << "EXITING: because of fatal INPUT ERROR: number of input read files: " << readFilesN << " does not agree with number of read group RG entries: " << outSAMattrRG.size() << "\n";
			errOut << "Make sure that the number of RG lines in --outSAMattrRGline is equal to either 1, or the number of input read files in --readFilesIn\n";
			exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
		}
		else if (outSAMattrRG.size() == 1) 
		{
			//use the same read group for all files              
			for (uint32 ifile = 1; ifile < readFilesN; ifile++) 
			{
				outSAMattrRG.push_back(outSAMattrRG.at(0));
			};
		};
	};
	readFilesIndex = 0;
};


#else // ~ #ifdef _WIN32

void Parameters::openReadsFiles() {
    string readFilesCommandString("");
    if (readFilesCommand.at(0)=="-") {
        if (readFilesIn.at(0).find(',')<readFilesIn.at(0).size()) readFilesCommandString="cat";//concatenate multiple files
    } else {
        for (uint ii=0; ii<readFilesCommand.size(); ii++) readFilesCommandString+=readFilesCommand.at(ii)+"   ";
    };

    if (readFilesCommandString=="") {//read from file
        for (uint ii=0;ii<readNmates;ii++) {//open readIn files
            readFilesCommandPID[ii]=0;//no command process IDs
            if ( inOut->readIn[ii].is_open() ) inOut->readIn[ii].close();
            inOut->readIn[ii].open(readFilesIn.at(ii).c_str()); //try to open the Sequences file right away, exit if failed
            if (inOut->readIn[ii].fail()) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal input ERROR: could not open readFilesIn=" << readFilesIn.at(ii) <<"\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);                
            };  
        };    
    } else {//create fifo files, execute pre-processing command

         vector<string> readsCommandFileName;

         readFilesNames.resize(readNmates); 

         for (uint imate=0;imate<readNmates;imate++) {//open readIn files
            ostringstream sysCom;
            sysCom << outFileTmp <<"tmp.fifo.read"<<imate+1;
            readFilesInTmp.push_back(sysCom.str());
            remove(readFilesInTmp.at(imate).c_str());   
            mkfifo(readFilesInTmp.at(imate).c_str(), S_IRUSR | S_IWUSR );

            inOut->logMain << "\n   Input read files for mate "<< imate+1 <<", from input string " << readFilesIn.at(imate) <<endl;

            readsCommandFileName.push_back(outFileTmp+"/readsCommand_read" + to_string(imate+1));
			fstream readsCommandFile( readsCommandFileName.at(imate).c_str(), ios::out | std::ios::binary);
            readsCommandFile.close();
            readsCommandFile.open( readsCommandFileName.at(imate).c_str(), ios::in | ios::out);
            //first line in the commands file
            if (sysShell!="-") {//executed via specified shell
                readsCommandFile << "#!" <<sysShell <<"\n";
            };
            readsCommandFile << "exec > \""<<readFilesInTmp.at(imate)<<"\"\n" ; // redirect stdout to temp fifo files
                        
            string readFilesInString(readFilesIn.at(imate));
            size_t pos=0;
            readFilesN=0;
            do {//cycle over multiple files separated by comma
                pos = readFilesInString.find(',');
                string file1 = readFilesInString.substr(0, pos);
                readFilesInString.erase(0, pos + 1);
                readFilesNames.at(imate).push_back(file1);
                
                system(("ls -lL " + file1 + " > "+ outFileTmp+"/readFilesIn.info 2>&1").c_str());
				ifstream readFilesIn_info((outFileTmp+"/readFilesIn.info").c_str(), std::ios::binary);
                inOut->logMain <<readFilesIn_info.rdbuf();

                readsCommandFile << "echo FILE " <<readFilesN << "\n";
                readsCommandFile << readFilesCommandString << "   " <<("\""+file1+"\"") <<"\n";
                ++readFilesN;//only increase file count for one mate

            } while (pos!= string::npos);

            readsCommandFile.flush();
            readsCommandFile.seekg(0,ios::beg);
            inOut->logMain <<"\n   readsCommandsFile:\n"<<readsCommandFile.rdbuf()<<endl;
            readsCommandFile.close();            
            
            chmod(readsCommandFileName.at(imate).c_str(),S_IXUSR | S_IRUSR | S_IWUSR);
            
            readFilesCommandPID[imate]=0;
            
            ostringstream errOut;
			pid_t PID = vfork();
            switch (PID) {
                case -1:
                    errOut << "EXITING: because of fatal EXECUTION error: Failed vforking readFilesCommand\n";
                    errOut << errno << ": " << strerror(errno) << "\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    break;
                
                case 0:
                    //this is the child
                    execlp(readsCommandFileName.at(imate).c_str(), readsCommandFileName.at(imate).c_str(), (char*) NULL); 
                    exit(0);
                    
                default:
                    //this is the father, record PID of the children
					readFilesCommandPID[imate] = PID;
            };
            
//             system((("\""+readsCommandFileName.at(imate)+"\"") + " & ").c_str());
            inOut->readIn[imate].open(readFilesInTmp.at(imate).c_str());                
        };
        if (readNmates==2 && readFilesNames.at(0).size() != readFilesNames.at(1).size()) {
            ostringstream errOut;
            errOut <<"EXITING: because of fatal INPUT ERROR: number of input files for mate1: "<<readFilesNames.at(0).size()  << " is not equal to that for mate2: "<< readFilesNames.at(1).size() <<"\n";
            errOut <<"Make sure that the number of files in --readFilesIn is the same for both mates\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);                
        };

        if (outSAMattrRG.size()>1 && outSAMattrRG.size()!=readFilesN) {
            ostringstream errOut;
            errOut <<"EXITING: because of fatal INPUT ERROR: number of input read files: "<<readFilesN << " does not agree with number of read group RG entries: "<< outSAMattrRG.size() <<"\n";
            errOut <<"Make sure that the number of RG lines in --outSAMattrRGline is equal to either 1, or the number of input read files in --readFilesIn\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);                
        } else if (outSAMattrRG.size()==1) {//use the same read group for all files              
            for (uint32 ifile=1;ifile<readFilesN;ifile++) {
                outSAMattrRG.push_back(outSAMattrRG.at(0));
            };
        };
    };
    readFilesIndex=0;
};

#endif // ~ #ifdef _WIN32
