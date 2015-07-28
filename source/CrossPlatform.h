#ifdef _WIN32

#define sleep(x) Sleep(x)
typedef __int64 int64;
typedef int key_t;

#endif