/*
 *  fast.h
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


#ifndef FAST_H

#define FAST_H

#include "stl.h"
#include "general_utils.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


struct fast_t{
    string name;
    string seq;
    string qual;
    bool qual_exist;
    //int worst_qual;
};

class fast{
private:
    string filename;
    gzFile fp;
    FILE* f;
	kseq_t *seq;
    int lines_read;
    int l;
    bool file_open;
    bool trim_slash;


public:
    fast(){
        file_open = false;
        trim_slash = true;
        l = 0;
        filename = "";
        lines_read = 0;
    }

    fast(string filename,bool trim_slash){
        this->trim_slash = trim_slash;
        open(filename);

    }

    ~fast(){ // desctuctor ... let's not memory leak
        if(file_open && is_open()){
            file_open = false;
            kseq_destroy(seq);
            gzclose(fp);
        }
    }

    fast_t get_next_record(){
        fast_t output = {"","","",100};

        if (!is_open()) {
            output.qual = "File not open!";
            return output;
        }

        bool line_gotten = false;

        while(!line_gotten){

            l = kseq_read(seq);

            if (l < 0) {
                return output;
            }

            if (l == 0) {
                continue;
            }

            if (seq->qual.l == 0) {  // fasta file
                output.name = seq->name.s;
                output.seq = seq->seq.s;
                output.qual = "*";
                output.qual_exist = false;
                //output.worst_qual = 62;
                line_gotten = true;

                if(trim_slash){
                        size_t slash_pos = output.name.find_last_of("/ \t");
                        output.name = output.name.substr(0,slash_pos);
                }


            }else { // fastq file... we don't read BAM files yet
                output.name = seq->name.s;
                output.seq = seq->seq.s;
                output.qual = seq->qual.s;
                output.qual_exist = true;
                //output.worst_qual = 62;
                line_gotten = true;

                if(trim_slash){
                        size_t slash_pos = output.name.find_last_of("/ \t");
                        output.name = output.name.substr(0,slash_pos);
                }


            }

        }

        return output;
    }

    bool open(string filename){
        file_open = true;
        l = 0;
        lines_read = 0;
        this->filename = filename;
        f = fopen(filename.c_str(), "r");
        fp = gzdopen(fileno(f), "r");

        seq = kseq_init(fp);



        return is_open();
    }

    z_off_t get_gzoffset(){
        //return gzoffset(seq->f->f);
        return 0;
    }

    bool eof(){
        return l<0;
    }

    bool is_open(){
        return seq->f != 0;
    }


    void close(){
        if(file_open && is_open()){
            file_open = false;
            kseq_destroy(seq);
            gzclose(fp);
        }
    }


    string get_name(){
        return filename;
    }

    int num_read(){
        return lines_read;
    }
};




#endif



