#ifndef IFSTREAMREADIN_H

#include <fstream>
#include <string>
#include <Windows.h>

class IfstreamReadIn
{
private:
	// private fields
	std::ifstream* pStream; 
	HANDLE hPipe; 

public:
	// Constructors 
	IfstreamReadIn(); // default one
	IfstreamReadIn(std::string file); 
	IfstreamReadIn(HANDLE hPipe); 

	// operators 
	std::ifstream& operator >> (const std::string& st); 

	// public methods.
	void open(const char* filename); 
	bool is_open() const; 
	bool fail() const; 
	void close(); 
	int peek(); 
	bool good() const;
	std::streamsize gcount() const;
	std::istream& getline(char* s, std::streamsize n);
	std::istream& ignore(std::streamsize n = 1, int delim = EOF);
			
};


#define IFSTREAMREADIN_H
#endif // IFSTREAMREADIN_H