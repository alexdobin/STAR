#include <string>
#include <cstring>
//#define _XOPEN_SOURCE 500
#include <ftw.h>
#include <unistd.h>

int removeFileOrDir(const char *fpath,const struct stat *sb, int typeflag, struct FTW *ftwbuf) {
    
    {//to avoid unused variable warning
        (void) sb;
        (void) ftwbuf;
    };
    
    if (typeflag==FTW_F) {//file
        remove(fpath);
    } else if (typeflag==FTW_DP) {//dir
        rmdir(fpath);
    } else {//something went wrong, stop the removal
        return -1;
    };
    return 0;
};


void sysRemoveDir(std::string dirName) {//remove directory and all its contents
    int nftwFlag=FTW_DEPTH;
    nftw(dirName.c_str(), removeFileOrDir, 100, nftwFlag);
};
