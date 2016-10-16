#ifndef write_log_h
#define write_log_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>

using namespace std;

//method declarations
template <class T>
void write_log(T log_input);

#include "write_log.hpp"


#endif // write_log_h