#ifndef DIRFUNCTIONS_H
#define DIRFUNCTIONS_H

#include "CrossPlatform.h"

void GetCurrentPath(std::string& path);
void removeDir(const std::string& dirName);
bool createDir(const std::string& dirName);
int mkdir(const std::string& dir, mode_t mode); 

#endif // DIRFUNCTIONS_H