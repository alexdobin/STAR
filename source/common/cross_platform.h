#ifndef HTSLIB_CROSS_PLATFORM_H
#define HTSLIB_CROSS_PLATFORM_H

#ifdef _WIN32

#ifndef __cplusplus
#define inline __inline
#endif

#define __func__ __FUNCTION__
#define strcasecmp _stricmp 
#define strncasecmp _strnicmp 
#define alloca _alloca
#define PATH_MAX _MAX_PATH

#endif //_WIN32

#endif //HTSLIB_CROSS_PLATFORM_H