#include "Transcript.h"

#include <cstring>
#include <iostream>
#include <new>

int main()
{
    alignas(Transcript) unsigned char storage[sizeof(Transcript)];
    memset(storage, 0xbe, sizeof(storage));

    Transcript *transcript=new (storage) Transcript;
    bool initialized=!transcript->sjYes && !transcript->primaryFlag;
    transcript->~Transcript();

    if (!initialized)
    {
        std::cerr << "Transcript boolean state was not initialized\n";
        return 1;
    };
    return 0;
}
