#include <string>
#include <boost/filesystem/operations.hpp>
#include <Windows.h>
#include "DirFunctions.h"

// Get current path 
void GetCurrentPath(std::string& path)
{
	char buffer[MAX_PATH] = { 0 };
	GetModuleFileName(NULL, buffer, MAX_PATH);
	std::string modulepath(buffer);
	path = modulepath.substr(0, modulepath.find_last_of("\\/"));
	path.append("\\");
}

// Remove directory recursively.
void removeDir(const std::string& dirName)
{
	boost::filesystem::remove_all(dirName);
}

// Create Directory
bool createDir(const std::string& dirName)
{
	return boost::filesystem::create_directory(dirName);
}

int mkdir(const std::string& dir, mode_t permission)
{
	// permission not used in Windows.
	return createDir(dir.c_str()) ? 0 : -1;
}
