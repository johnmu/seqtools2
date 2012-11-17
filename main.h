/*
 *  main.h
 *  seqtools
 *
 *  Created by John Mu on 3/3/11.
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

#ifndef MAIN_H

#define MAIN_H



#include "stl.h"
#include "general_utils.h"
#include "fast.h"


using namespace std;

const string PROG_NAME = "seqtools2";
const string BUILD = "0.1-r17";


class twoint {
public:
    uint64_t a;
    int64_t b;

    twoint(const uint64_t &a, const int64_t &b) {
        this->a = a;
        this->b = b;
    }

    twoint() {
        a = 0;
        b = 0;
    }

    twoint(const twoint &x) {
        this->a = x.a;
        this->b = x.b;
    }

    bool operator<(const twoint& y) const {
        return a < y.a;
    }

    void set_val(const uint64_t &a, const int64_t &b) {
        this->a = a;
        this->b = b;
    }
};

struct b_store {
    string code;
    string name;
    int count;
    ofstream* outfile;


    bool operator<(const b_store & a) const {
        return code < a.code;
    }

};


struct b_store_pair {
    string code;
    string name;
    int count;
    ofstream* outfile1;
    ofstream* outfile2;

    bool operator<(const b_store_pair & a) const {
        return code < a.code;
    }

};


struct chr_idx_t {
    int count;
    string name;

    bool operator<(const chr_idx_t & a) const {
        return count > a.count;
    }

};


void print_usage_and_exit();

int merge_sam(vector<string> params);
int sample_reads(vector<string> params);
int fastq_filter(vector<string> params);
int fq_len_filter(vector<string> params);
int qual_convert(vector<string> params);
int fastq_stats(vector<string> params);
int fastq_quals(vector<string> params);
int sam_quals(vector<string> params);
int replace_sam_quals(vector<string> params);
int quant_sam_quals(vector<string> params);
int depth_stats(vector<string> params);
int sam_stats(vector<string> params);
int fastq_trim(vector<string> params);
int barcode_split(vector<string> params);
int pair_split(vector<string> params);
int raw2fasta(vector<string> params);

int time_test(vector<string> params);


int subseq(vector<string> params);

#endif



