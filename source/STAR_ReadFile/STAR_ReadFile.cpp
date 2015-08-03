#include <stdio.h>
#include <zlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <Windows.h>

#define BUFSIZE 1024 

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		// no file name, return error.
		return 1;
	}

	// Get current path 
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	std::string modulepath(buffer);
	std::string filedir = modulepath.substr(0, modulepath.find_last_of("\\/"));
	filedir.append("\\");
	
	std::vector<std::string> files;
	std::string filenames = argv[1]; 
	std::string delimiter = ","; 
	
	//get filename(s) from commandline.
	size_t  start = 0, end = 0;
	while (end != std::string::npos)
	{
		end = filenames.find(delimiter, start);

		std::string name = filenames.substr(start, (end == std::string::npos) ? std::string::npos : end - start);
		std::string fullpath = filedir + name; 
		files.push_back(fullpath);

		start = ((end > (std::string::npos - delimiter.size()))
			? std::string::npos : end + delimiter.size());
	}

	int i = 0; 
	for (auto& f : files)
	{
		printf("\nFILE %d", i++); 
		// Read gz file and write to stdout using printf 
		gzFile inFileZ = gzopen(f.c_str(), "rb");
		if (inFileZ == NULL) {
			//printf("Error: Failed to gzopen %s with error %d \n", fileName.c_str(), errno);
			return 1;
		}
		unsigned char unzipBuffer[BUFSIZE];
		unsigned int unzippedBytes = 0;
		while (true)
		{
			unzippedBytes = gzread(inFileZ, unzipBuffer, BUFSIZE);
			if (unzippedBytes > 0)
			{
				for (unsigned int i = 0; i < unzippedBytes; i++)
				{
					// write bytes to stdout
					printf("%s", unzipBuffer);
				}
			}
			else
			{
				break;
			}
		}
		gzclose(inFileZ);
	}
	return 0;
}