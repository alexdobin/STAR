#include <stdio.h>
#include <zlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <Windows.h>
#include <fstream>

#define BUFSIZE 1024


int extractAndOutputFile(const std::string& filepath); 
int outputFile(const std::string& filepath); 

/*	argv[0] - exe name
	argv[1] - filenames separated by ","
	argv[2] - "cat" or "zcat"
			  if "cat", output file to stdout
			  else if "zcat", extract file using zlib and output to stdout.
			  else, return error code 1 
*/
int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		// no file name, return error.
		return 1;
	}

	// Get current path 
	char buffer[MAX_PATH] = {0};
	GetModuleFileName(NULL, buffer, MAX_PATH);
	std::string modulepath(buffer);
	std::string filedir = modulepath.substr(0, modulepath.find_last_of("\\/"));
	filedir.append("\\");
	

	std::vector<std::string> files;
	
	// argv[1] is file name(s)
	std::string filenames = argv[1]; 
	std::string delimiter = ","; 
	
	//get filename(s) separated by ","
	size_t  start = 0, end = 0;
	while (end != std::string::npos)
	{
		end = filenames.find(delimiter, start);

		std::string name = filenames.substr(start, (end == std::string::npos) ? std::string::npos : end - start);
		std::string fullpath = filedir + name; 
		files.push_back(fullpath);

		start = ((end > (std::string::npos - delimiter.size()))	? std::string::npos : end + delimiter.size());
	}

	int i = 0; 
	std::string command = argv[2]; 
	for (auto& f : files)
	{
		std::cout<< "FILE " << i++ << std::endl;

		if (command == "cat")
		{
			outputFile(f); 
		}
		else if (command == "zcat")
		{
			extractAndOutputFile(f); 
		}
	}
	return 1;
}

int extractAndOutputFile(const std::string& filepath)
{
	// Read gz file and write to stdout using printf 
	gzFile inFileZ = gzopen(filepath.c_str(), "rb");
	if (inFileZ == NULL) {
		// failed to open gz file, return failed code 1
		return 1;
	}
	
	unsigned int unzippedBytes = 0;
	int i = 0; 
	while (true)
	{
		unsigned char unzipBuffer[BUFSIZE + 1] = { 0 };
		unzippedBytes = gzread(inFileZ, unzipBuffer, BUFSIZE);
		if (unzippedBytes > 0)
		{
			// write bytes to stdout
			unzipBuffer[BUFSIZE] = '\0'; 
			std::cout<< unzipBuffer;
		}
		else
		{
			// Write End Marker
			std::cout << "--END--" << std::endl;
			break;
		}
	}
	gzclose(inFileZ);
	return 0; 
}
int outputFile(const std::string& filepath)
{
	std::ifstream file(filepath);
	if (!file.good())
	{
		// failed to open file, retrurn failed code 1
		return 1;
	}
	char buffer[BUFSIZE] = {0};
	while (!file.eof())
	{
		file.read(buffer, BUFSIZE);
		// write to stdout.
		printf(buffer);
	}
	// Write End Marker
	printf("%s","--END--");
	file.close();
	// return success code 0
	return 0;
}