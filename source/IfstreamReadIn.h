#ifndef IFSTREAMREADIN_H
#define IFSTREAMREADIN_H

#include <io.h>
#include <fcntl.h>
#include <fstream>
#include <string>
#include <Windows.h>

class IfstreamReadIn
{
private:
	// private fields
	std::ifstream* _pStream; 

public:
	// Constructors 
	IfstreamReadIn(); // default one
	IfstreamReadIn(const std::string& file);

	//destructor
	~IfstreamReadIn();

	// operators 
	std::ifstream& operator >> (char* s); 
	std::ifstream& operator >> (std::string s);
	std::ifstream& operator >> (int n);

	// public methods.
	void open(const char* filename);
	bool open_pipe_read(HANDLE hPipeRead); 
	bool is_open() const; 
	bool fail() const; 
	void close(); 
	int peek(); 
	bool good() const;
	std::streamsize gcount() const;
	std::istream& getline(char* s, std::streamsize n);
	std::istream& ignore(std::streamsize n = 1, int delim = EOF);
			
};


#endif // IFSTREAMREADIN_H