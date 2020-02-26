// Utility for maintaining separate copies of an object, each of which is used in a different thread in OpenMP

#ifndef THREAD_H
#define THREAD_H

#include <omp.h>
#include <vector>

template <typename T>
class OMPWrapper
{
public:
    static const int MaxThreads = 16;
    OMPWrapper() : t(MaxThreads) { }
    OMPWrapper(const T& value) : t(MaxThreads, value) { }
    T& operator  * () { return t[omp_get_thread_num()]; }
    T* operator -> () { return &t[omp_get_thread_num()]; }
    T& operator [] (int i) { return t[i]; }
    auto begin() { return t.begin(); }
    auto end()   { return t.end(); }
    auto size()  { return t.size(); }

private:
    std::vector<T> t;
};

/*

could add bounds checking above:
    if (omp_get_thread_num() >= MaxThreads)
        throw logic_error("Thread number exceeds MaxThreads.");


*/
#endif
