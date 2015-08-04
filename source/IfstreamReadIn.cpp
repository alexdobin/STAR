
#include "IfstreamReadIn.h"

IfstreamReadIn::IfstreamReadIn() : _pStream(nullptr)
{

}
IfstreamReadIn::IfstreamReadIn(const std::string& file)
{

}

IfstreamReadIn::~IfstreamReadIn() 
{
	if (_pStream )
	{
		if (_pStream->is_open())
			_pStream->close(); 
		delete _pStream; 
	}
 
}

std::ifstream& IfstreamReadIn::operator >> (char* s)
{ 
	if (_pStream && _pStream->is_open())
	{
		(*_pStream) >> s;
		return *_pStream;
	}
	
}

std::ifstream& IfstreamReadIn::operator >> (std::string s)
{
	if (_pStream && _pStream->is_open())
	{
		(*_pStream) >> s;
		return *_pStream;
	}
}

std::ifstream& IfstreamReadIn::operator >> (int n)
{
	if (_pStream && _pStream->is_open())
	{
		(*_pStream) >> n;
		return *_pStream;
	}
}

void IfstreamReadIn::open(const char* filename)
{
	_pStream = new std::ifstream(filename); 
}

bool  IfstreamReadIn::open_pipe_read(HANDLE hPipeRead)
{
	if (hPipeRead != INVALID_HANDLE_VALUE) 
	{
		int fd = _open_osfhandle((intptr_t)hPipeRead, _O_RDONLY);
		if (fd != -1)
		{
			FILE* f = _fdopen(fd, "r");
			if (f != 0)
			{
				_pStream = new std::ifstream(f);
			}
		}
		else
		{
			_close(fd); // Also calls CloseHandle()
			return false; 
		}
	}
	else
	{
		CloseHandle(hPipeRead);
		return false; 
	}
	return true; 
}

bool IfstreamReadIn::is_open() const
{
	if (_pStream)
		return _pStream->is_open();
	else
		return false; 
}

bool IfstreamReadIn::fail() const
{
	if (_pStream && _pStream->is_open())
		return _pStream->fail();
	else
		return false; 
}

void IfstreamReadIn::close()
{
	if (_pStream && _pStream->is_open())
		_pStream->close();
}

int IfstreamReadIn::peek()
{
	if (_pStream && _pStream->is_open())
		return _pStream->peek(); 
}

bool IfstreamReadIn::good() const
{
	if (_pStream && _pStream->is_open())
		return _pStream->good();
	else
		return false; 
}

std::streamsize IfstreamReadIn::gcount() const
{
	if (_pStream && _pStream->is_open())
		return _pStream->gcount();
}

std::istream& IfstreamReadIn::getline(char* s, std::streamsize n)
{
	if (_pStream && _pStream->is_open())
		return _pStream->getline(s,n);
}

std::istream& IfstreamReadIn::ignore(std::streamsize n, int delim)
{
	if (_pStream && _pStream->is_open())
		return _pStream->ignore(n,delim);
}