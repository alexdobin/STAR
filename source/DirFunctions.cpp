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
	// Get current path and create full path
	std::string fullPath;
	GetCurrentPath(fullPath);
	fullPath.append(dirName);
	boost::filesystem::path path(fullPath); 
	boost::filesystem::remove_all(path);
}

// Create Directory
bool createDir(const std::string& dirName)
{
	// Get current path and create full path
	std::string fullPath;
	GetCurrentPath(fullPath);
	fullPath.append(dirName);

	return boost::filesystem::create_directory(fullPath);
}
