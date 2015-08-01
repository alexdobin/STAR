#include "IfstreamReadIn.h"

IfstreamReadIn::IfstreamReadIn()
{

}
IfstreamReadIn::IfstreamReadIn(std::string file)
{

}
IfstreamReadIn::IfstreamReadIn(HANDLE hPipe)
{

}

std::ifstream& IfstreamReadIn::operator >> (const std::string& st)
{ 
	// TODO : validation.
	return (*pStream);
}

void IfstreamReadIn::open(const char* filename)
{

}

bool IfstreamReadIn::is_open() const
{
	//TODO valiation.
	return pStream->is_open(); 
}

bool IfstreamReadIn::fail() const
{
	// TODO validation.
	return pStream->fail(); 
}

void IfstreamReadIn::close()
{
	// TODO validation.
	pStream->close();
}

int IfstreamReadIn::peek()
{
	// TODO : validation.
	return pStream->peek(); 
}

bool IfstreamReadIn::good() const
{
	// TODO : validation.
	return pStream->good();
}

std::streamsize IfstreamReadIn::gcount() const
{
	// TODO : validation.
	return pStream->gcount();
}

std::istream& IfstreamReadIn::getline(char* s, std::streamsize n)
{
	// TODO : validation.
	return pStream->getline(s,n);
}

std::istream& IfstreamReadIn::ignore(std::streamsize n, int delim)
{
	// TODO : validation.
	return pStream->ignore(n,delim);
}