/*
 *  stl.h
 *  seqtools
 *
 *  Created by John Mu on 3/4/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 * 
 *  john.mu@ieee.org
 *
 */


/* The MIT License

   Copyright (c) 2012 John C. Mu.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#ifndef STL_H

#define STL_H

// #define DEBUG 1
// #define DEBUG_SSE 1
// #define DEBUG_KMER 1
// #define DEBUG_PAIR 1
// #define NO_OUTPUT 1
// #define TIMER 1

#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include <map>
#include <sstream>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <assert.h>
#include <ctime>
#include <cctype>
#include <climits>
#include <stdint.h>
#include <bitset>
#include <queue>
#include <deque>
#include <list>
#include <queue>
#include <zlib.h>



// linux specific stuff
#include <sys/time.h>
#include <pthread.h>

// oh well
using namespace std;

namespace global {
    const float mm_error_prob = 0.02;
    const int MAGIC_INT = 314;
    const int HASH_VERSION = 6;
}

namespace qual {
    const int SANGER = 0;
    const int ILL_13 = 1;
    const int ILL_15 = 2;
    const int SOLEXA = 3;

    // define the lookup table here
    // -5 ... 62 (thanks to MATLAB!)
    const char solexa2phred32_table[68] ={1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7,
        8, 9, 10, 10, 11, 12, 13, 14, 15, 16, 17, 18,
        19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
        43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
        55, 56, 57, 58, 59, 60, 61, 62};

}


namespace index_type {
    const int ALL_KMER = 0;
    const int SAMPLE   = 1;
    const int SPACED   = 2;
    const int SAMPLE_AND_SPACED = 3;
    const int ALL_KMER_SPACED = 4;
}


namespace aln_type{
    const char MM_ALIGN = 'M';
    const char NW_ALIGN = 'W';
    const char SW_ALIGN = 'S';
    const char UNIQUE   = 'U';
    const char REPEAT   = 'R';
}



#endif
