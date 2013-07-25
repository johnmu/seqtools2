/*
 *  main.cpp
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

#include "main.h"

//const double LOG_2 = log(2);

int main(int argc, char * const argv[]) {

#if defined DEBUG || defined DEBUG_PAIR
    cerr << "DEBUG MODE!... :(" << '\n';
#endif


    cerr << "**" + PROG_NAME << '\n';
    cerr << "Version: " << BUILD << "\n";
    cerr << "John C. Mu" << '\n';
    cerr << "http://www.stanford.edu/group/wonglab" << "\n";
    cerr << "\n";


    srand((unsigned int) time(NULL));

    string mode = "";
    string temp = "";
    int error_num = 0;
    vector<string> params;

    if (argc < 2) {
        print_usage_and_exit();
    }

    mode = argv[1];

    for (int i = 2; i < argc; i++) {
        temp = argv[i];
        params.push_back(temp);
    }

    if (mode == "raw2fasta") {

        error_num = raw2fasta(params);
    } else if (mode == "subseq") {

        error_num = subseq(params);
    } else if (mode == "qual_convert") {

        error_num = qual_convert(params);
    } else if (mode == "fastq_trim") {

        error_num = fastq_trim(params);
    } else if (mode == "fastq_pico_trim") {

        error_num = fastq_pico_trim(params);
    } else if (mode == "merge_sam") {

        error_num = merge_sam(params);
    } else if (mode == "sample_reads") {

        error_num = sample_reads(params);
    } else if (mode == "fastq_filter") {

        error_num = fastq_filter(params);
    } else if (mode == "fq_len_filter") {

        error_num = fq_len_filter(params);
    } else if (mode == "fastq_stats") {

        error_num = fastq_stats(params);
    } else if (mode == "fastq_quals") {

        error_num = fastq_quals(params);
    } else if (mode == "sam_quals") {

        error_num = sam_quals(params);
    } else if (mode == "replace_sam_quals") {

        error_num = replace_sam_quals(params);
    } else if (mode == "replace_sam_quals_diff") {

        error_num = replace_sam_quals_diff(params);
    } else if (mode == "quant_sam_quals") {

        error_num = quant_sam_quals(params);
    } else if (mode == "depth_stats") {

        error_num = depth_stats(params);
    } else if (mode == "sam_stats") {

        error_num = sam_stats(params);
    } else if (mode == "barcode_split") {

        error_num = barcode_split(params);
    } else if (mode == "pair_split") {

        error_num = pair_split(params);
    } else if (mode == "select_haps") {

        error_num = select_haps(params);
    } else if (mode == "vcf_to_hap") {

        error_num = vcf_to_hap(params);
    } else {
        print_usage_and_exit();

    }

    return error_num;
}

inline vector<b_store_pair>::iterator search_barcode(vector<b_store_pair> &vec,string barcode, int num_mm){
    vector<b_store_pair>::iterator out = vec.end();

    for (vector<b_store_pair>::iterator i = vec.begin();i!= vec.end();i++){

        int mm = 0;

        if(barcode.length() != i->code.length()){
            cerr << "Warning: Barcode length wrong, " << barcode << '\n';
        }

        for(int j = 0;j<(int)barcode.length();j++){
            if(barcode[j] != i->code[j]){
                mm++;
                if(mm>num_mm) break;
            }
        }

        if(mm < num_mm){
            if(out == vec.end()){
                out = i;
            }else{
                return vec.end();
            }
        }

    }


    return out;
}


int pair_split(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " pair_split <barcode_file> <file1> <file2>\n"
            + "     barcode_file - tab-delimited text file (barcode, filename_prefix)\n"
            + "     file1   - The first pair to split\n"
            + "     file2   - The second pair to split\n"
            + "Split a file into multiple files by bar code, FASTQ can be gziped";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    string barcode_filename = params[0];
    string pair1 = params[1];
    string pair2 = params[2];


    ifstream infile1;
    ifstream infile2;


    infile1.open(pair1.c_str(), ios::in);

    if (!infile1.is_open()) {
        cerr << "Error: cannot open input file1: " << pair1 << '\n';
        return 3;
    }

    infile1.close();

    infile2.open(pair2.c_str(), ios::in);

    if (!infile2.is_open()) {
        cerr << "Error: cannot open input file1: " << pair2 << '\n';
        return 3;
    }

    infile2.close();


    ifstream barcodefile;

    barcodefile.open(barcode_filename.c_str(), ios::in);

    if (!barcodefile.is_open()) {
        cerr << "Error: cannot open barcode file: " << barcode_filename << '\n';
        return 3;
    }


    // read in the barcode list


    int barcode_len = 0;
    vector<b_store_pair> barcode_vec;

    int num_unmatched = 0;
    ofstream unmatched_out1;
    string temp = "unmatched_"+pair1;
    unmatched_out1.open(temp.c_str(),ios::out);
    if(!unmatched_out1.is_open()){
        cerr << "ERROR: could not write, " << temp << '\n';
        return 2;
    }
    ofstream unmatched_out2;
    temp = "unmatched_"+pair2;
    unmatched_out2.open(temp.c_str(),ios::out);
    if(!unmatched_out2.is_open()){
        cerr << "ERROR: could not write, " << temp << '\n';
        return 2;
    }


    string line;

    while(!barcodefile.eof()){
        getline(barcodefile,line);
        trim2(line);

        if(line.length() == 0){
            continue;
        }

        vector<string> line_list = split(line);

        if(line_list.size() != 2){
            cerr << "ERROR: Bad barcode file format: " << line << '\n';
            return 1;
        }

        if(barcode_len == 0){
            barcode_len = line_list[0].length();
        }else if(barcode_len != (int)line_list[0].length()){
            cerr << "ERROR: Barcode not all same length\n";
            return 2;

        }

        b_store_pair contents;

        contents.code  = line_list[0];
        contents.name  = line_list[1];
        contents.count = 0;
        string temp_name = line_list[1] + "_" + line_list[0] + "_" + pair1;
        contents.outfile1 = new ofstream(temp_name.c_str(),ios::out);

        if(!contents.outfile1->is_open()){
            cerr << "ERROR: cannot create file: " << temp_name << '\n';
            return 2;
        }

        temp_name = line_list[1] + "_" + line_list[0] + "_" + pair2;
        contents.outfile2 = new ofstream(temp_name.c_str(),ios::out);

        if(!contents.outfile2->is_open()){
            cerr << "ERROR: cannot create file: " << temp_name << '\n';
            return 2;
        }

        barcode_vec.push_back(contents);



    }

    barcodefile.close();

    if(barcode_vec.size() == 0){
        cerr << "ERROR: empty barcode file\n";
        return 1;
    }

    sort(barcode_vec.begin(),barcode_vec.end());


    fast in1(pair1, false);
    fast in2(pair2, false);

    while (!in1.eof()) {
        fast_t line1 = in1.get_next_record();
        fast_t line2 = in2.get_next_record();

        if (line1.seq.length() == 0) {
            continue;
        }
        if (line2.seq.length() == 0) {
            continue;
        }

        //int len1 = (int) line1.seq.length();
        //int len2 = (int) line2.seq.length();


        string barcode1 = "";
        string barcode2 = "";

        size_t loc = line1.name.rfind(':');
        if (loc != string::npos) {
            barcode1 = line1.name.substr(loc + 1, barcode_len);
        } else {
            cerr << "ERROR: no : found, " << line1.name << '\n';
            return 2;
        }

        loc = line2.name.rfind(':');
        if (loc != string::npos) {
            barcode2 = line2.name.substr(loc + 1, barcode_len);
        } else {
            cerr << "ERROR: no : found, " << line2.name << '\n';
            return 2;
        }

        vector<b_store_pair>::iterator code_loc = barcode_vec.end();
        if(barcode1 == barcode2){
            code_loc = search_barcode(barcode_vec,barcode1, 1);
        }

        if(!(code_loc == barcode_vec.end())){

            if (line1.qual_exist) {
                *(code_loc->outfile1) << "@" << line1.name << '\n';
                *(code_loc->outfile1) << line1.seq << '\n';
                *(code_loc->outfile1) << "+" << '\n';
                *(code_loc->outfile1) << line1.qual << '\n';
            }else{
                *(code_loc->outfile1) << ">" << line1.name << '\n';
                *(code_loc->outfile1) << line1.seq << '\n';
            }

            if (line2.qual_exist) {
                *(code_loc->outfile2) << "@" << line2.name << '\n';
                *(code_loc->outfile2) << line2.seq << '\n';
                *(code_loc->outfile2) << "+" << '\n';
                *(code_loc->outfile2) << line2.qual << '\n';
            }else{
                *(code_loc->outfile2) << ">" << line2.name << '\n';
                *(code_loc->outfile2) << line2.seq << '\n';
            }

            code_loc->count++;
        }else{
            if (line1.qual_exist) {
                unmatched_out1 << "@" << line1.name << '\n';
                unmatched_out1 << line1.seq << '\n';
                unmatched_out1 << "+" << '\n';
                unmatched_out1 << line1.qual << '\n';


            } else {
                unmatched_out1 << ">" << line1.name << '\n';
                unmatched_out1 << line1.seq << '\n';

            }

            if (line2.qual_exist) {
                unmatched_out2 << "@" << line2.name << '\n';
                unmatched_out2 << line2.seq << '\n';
                unmatched_out2 << "+" << '\n';
                unmatched_out2 << line2.qual << '\n';
            }else{
                unmatched_out2 << ">" << line2.name << '\n';
                unmatched_out2 << line2.seq << '\n';
            }
            num_unmatched++;
        }



    }

    if(!in2.eof()){
        cerr << "Warning... files different lengths :S \n";
    }

    unmatched_out1.close();
    unmatched_out2.close();

    // output stats
    for(vector<b_store_pair>::iterator it = barcode_vec.begin();it != barcode_vec.end();it++){
        cout << it->name << "," << it->code << "," << it->count << '\n';
        it->outfile1->close();
        it->outfile2->close();
    }
    cout << "Unmatched,Unmatched," << num_unmatched << '\n';




    return 0;
}

int barcode_split(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " barcode_split <barcode_file> <file> <mode>\n"
            + "     barcode_file - tab-delimited text file (barcode, filename_prefix)\n"
            + "     file   - The file to split\n"
            + "     mode   - 2 = barcode in name[SAM], 1=barcode in name[FASTQ], 0=unstripped[FASTQ]\n"
            + "Split a file into multiple files by bar code, FASTQ can be gziped";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    string barcode_filename = params[0];
    string filename = params[1];
    int mode = strTo<int>(params[2]);

    if(mode < 0 || mode > 2){
        cerr << "ERROR: mode wrong \n";
        return 3;
    }

    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    ifstream barcodefile;

    barcodefile.open(barcode_filename.c_str(), ios::in);

    if (!barcodefile.is_open()) {
        cerr << "Error: cannot open barcode file: " << barcode_filename << '\n';
        return 3;
    }


    // read in the barcode list


    int barcode_len = 0;
    vector<b_store> barcode_vec;

    int num_unmatched = 0;
    ofstream unmatched_out;
    string temp = "unmatched_"+filename;
    unmatched_out.open(temp.c_str(),ios::out);
    if(!unmatched_out.is_open()){
        cerr << "ERROR: could not write, " << temp << '\n';
        return 2;
    }


    string line;

    while(!barcodefile.eof()){
        getline(barcodefile,line);
        trim2(line);

        if(line.length() == 0){
            continue;
        }

        vector<string> line_list = split(line);

        if(line_list.size() != 2){
            cerr << "ERROR: Bad barcode file format: " << line << '\n';
            return 1;
        }

        if(barcode_len == 0){
            barcode_len = line_list[0].length();
        }else if(barcode_len != (int)line_list[0].length()){
            cerr << "ERROR: Barcode not all same length\n";
            return 2;

        }

        b_store contents;

        contents.code  = line_list[0];
        contents.name  = line_list[1];
        contents.count = 0;
        string temp_name = line_list[1] + "_" + line_list[0] + "_" + filename;
        contents.outfile = new ofstream(temp_name.c_str(),ios::out);

        if(!contents.outfile->is_open()){
            cerr << "ERROR: cannot create file: " << temp_name << '\n';
            return 2;
        }

        barcode_vec.push_back(contents);



    }

    barcodefile.close();

    if(barcode_vec.size() == 0){
        cerr << "ERROR: empty barcode file\n";
        return 1;
    }

    sort(barcode_vec.begin(),barcode_vec.end());



    if (mode == 0 || mode == 1) {
        // FASTQ

        fast in(filename, false);

        while (!in.eof()) {
            fast_t line = in.get_next_record();

            if (line.seq.length() == 0) {
                continue;
            }

            int len = (int) line.seq.length();

            if (mode == 0 && len <= barcode_len) {
                cerr << "Warning: len <= barcode_len, " << line.name << '\n';
            } else {
                b_store search;

                if (mode == 0) {
                    search.code = line.seq.substr(0, barcode_len);
                } else if (mode == 1) {
                    size_t loc = line.name.find('#');
                    if (loc != string::npos) {
                        search.code = line.name.substr(loc + 1, barcode_len);
                    } else {
                        cerr << "ERROR: no # found, " << line.name << '\n';
                        return 2;
                    }

                }

                string seq;
                string qual;

                if (mode == 0){
                    seq = line.seq.substr(barcode_len);
                    if (line.qual_exist) qual = line.qual.substr(barcode_len);
                }else if(mode == 1){
                    seq = line.seq;
                    if (line.qual_exist) qual = line.qual;
                }

                vector<b_store>::iterator it = lower_bound(barcode_vec.begin(), barcode_vec.end(), search);

                if (it != barcode_vec.end() && it->code == search.code) {
                    if (line.qual_exist) {
                        *(it->outfile) << "@" << line.name << '\n';
                        *(it->outfile) << seq << '\n';
                        *(it->outfile) << "+" << '\n';
                        *(it->outfile) << qual << '\n';
                    } else {
                        *(it->outfile) << ">" << line.name << '\n';
                        *(it->outfile) << seq << '\n';
                    }

                    it->count++;

                } else {
                    if (line.qual_exist) {
                        unmatched_out << "@" << line.name << '\n';
                        unmatched_out << seq << '\n';
                        unmatched_out << "+" << '\n';
                        unmatched_out << qual << '\n';
                    } else {
                        unmatched_out << ">" << line.name << '\n';
                        unmatched_out << seq << '\n';
                    }

                    num_unmatched++;
                }


            }

        }
    }else if(mode == 2){
        // SAM

        ifstream in(filename.c_str(), ios::in);
        string line;

        while (!in.eof()) {
            getline(in,line);

            trim2(line);

            if (line.length() == 0 || line[0] == '@') {
                continue;
            }

            b_store search;

            size_t loc = line.find('#');
            if (loc != string::npos) {
                search.code = line.substr(loc + 1, barcode_len);
            } else {
                cerr << "ERROR: no # found, " << line << '\n';
                return 2;
            }


            vector<b_store>::iterator it = lower_bound(barcode_vec.begin(), barcode_vec.end(), search);

            if (it != barcode_vec.end() && it->code == search.code) {
                *(it->outfile) <<  line << '\n';
                it->count++;

            } else {
                unmatched_out << line << '\n';

                num_unmatched++;
            }




        }

    }

    unmatched_out.close();

    // output stats
    for(vector<b_store>::iterator it = barcode_vec.begin();it != barcode_vec.end();it++){
        cout << it->name << "," << it->code << "," << it->count << '\n';
        it->outfile->close();
    }
    cout << "Unmatched,Unmatched," << num_unmatched << '\n';




    return 0;
}



int fastq_pico_trim(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " fastq_pico_trim <FASTA/Q_file>\n"
            + "    FASTA/Q_file     -- \n"
            + "Trim a FASTA/Q file to remove the picoplex adaptor";

    if (params.size() != 1) {
        cerr << usage_text << endl;
        return 3;
    }

    string filename = params[0];

    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    fast in(filename, false);

    while (!in.eof()) {
        fast_t line = in.get_next_record();

        if ((int)line.seq.length() == 0) {
            continue;
        }

        

        int start_idx = trim_phase_seq_primer(line.seq);

        
        if(line.qual_exist) {

            if (start_idx >= 0) {
                if (start_idx < (int) line.seq.length() - 1) {
                    line.seq = line.seq.substr(start_idx);
                    line.qual = line.qual.substr(start_idx);
                } else {
                    line.seq = "N";
                    line.qual = "#";
                }
            }

            cout << '@' << line.name << '\n';
            cout << line.seq << '\n';
            cout << "+\n";
            cout << line.qual << '\n';
        } else {
            
            if(start_idx >= 0) {
                if (start_idx < (int) line.seq.length() - 1) {
                    line.seq = line.seq.substr(start_idx);
                } else {
                    line.seq = "N";
                }
            }

            cout << '>' << line.name << '\n';
            cout << line.seq << '\n';;
        }

    }


    return 0;
}


int fastq_trim(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " fastq_trim <min_len> <Phred_score> <FASTQ_file>\n"
            + "    min_len        -- Minimum length of read after trimming\n"
            + "    Phred_score    -- Phred score threshold (reads must be Sanger quals)\n"
            + "    FASTQ_file     -- Pred score (reads must be Sanger quals)\n"
            + "Trim a FASTQ file BWA style, reads trimmed to zero are replaced by a single N";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    int min_len = strTo<int>(params[0]);
    int qual    = strTo<int>(params[1]);
    string filename = params[2];

    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    fast in(filename, false);

    while (!in.eof()) {
        fast_t line = in.get_next_record();

        if (line.seq.length() == 0) {
            continue;
        }

        if (!line.qual_exist){
            cerr << "Error: not FASTQ file: " << line.name << '\n';
            return 2;
        }

        int new_len = bwa_trim_len(line.qual, qual, min_len);

        if(new_len > 0){

            line.seq = line.seq.substr(0,new_len);
            line.qual = line.qual.substr(0,new_len);
        }else{
            line.seq = "N";
            line.qual = "#";
        }

        cout << '@' << line.name << '\n';
        cout << line.seq << '\n';
        cout << "+\n";
        cout << line.qual << '\n';


    }


    return 0;
}



int depth_stats(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " depth_stats <genome_file> <depth_file>\n"
            + "Get stats for a depth file";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

    string genome_fai_file = params[0] + ".fai";
    string depth_file  = params[1];

    ifstream infile;
    int64_t chr_len = 0;
    string chr_name = "";

    infile.open(genome_fai_file.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open genome fai file: " << genome_fai_file << '\n';
        return 3;
    }else{
        string temp = "";

        getline(infile,temp);

        trim2(temp);

        if(temp.length() != 0){
            vector<string> ll = split(temp);

            if(ll.size() == 5){
                chr_name = ll[0];
                chr_len = strTo<int64_t>(ll[1]);
            }else{
                cerr << "Error: Bad fai file: " << temp << '\n';
                return 3;
            }

        }


    }



    infile.close();

    infile.open(depth_file.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open depth file: " << depth_file << '\n';
        return 3;
    }

    // stats to record
    int64_t num_pos = 0; // number of positions covered
    int64_t num_pos_10 = 0; // number of positions covered with depth greater or equal to 10

    vector<int64_t> contig_store;
    vector<int64_t> gap_store;

    vector<int> depth_100;

    // temp variables
    string line;
    int avg_depth_100 = 0;
    int64_t pos_counter = 0;
    double depth_sum = 0.0;

    int64_t prev_pos_lim = 0;
    int64_t prev_pos = 0;

    while(!infile.eof()){
        getline(infile,line);

        trim2(line);

        if(line.length() == 0) continue;

        vector<string> ll = split(line);

        if(ll.size() != 3) continue;

        if(ll[0] != chr_name) continue;

        int64_t pos   = strTo<int64_t>(ll[1]);
        int     depth = strTo<int>(ll[2]);

        while(pos_counter < pos){
            pos_counter++;

            if(pos_counter % 100 == 0){
                depth_100.push_back(avg_depth_100);
                avg_depth_100 = 0;
            }
        }

        avg_depth_100 += depth;
        depth_sum += depth;

        num_pos++;
        if(depth >= 10)num_pos_10++;

        if(prev_pos==0){
            prev_pos_lim = pos;
            prev_pos = pos;
        }

        if(pos == prev_pos + 1){
            prev_pos = pos;
        }else{
            // we have hit a gap!!
            int64_t contig_size = prev_pos - prev_pos_lim + 1;
            int64_t gap_size = pos - prev_pos - 1;

            contig_store.push_back(contig_size);
            gap_store.push_back(gap_size);


            prev_pos = pos;
            prev_pos_lim = pos;
        }

    }

    while(pos_counter < chr_len) {
        pos_counter++;

        if (pos_counter % 100 == 0) {
            depth_100.push_back(avg_depth_100);
            avg_depth_100 = 0;
        }
    }


    // add last one
    contig_store.push_back(prev_pos - prev_pos_lim + 1);


    infile.close();

    sort(contig_store.begin(),contig_store.end());
    sort(gap_store.begin(),gap_store.end());

    // output all the good stuff
    cout << "Total % of " << chr_name << " covered: ";
    output_percentage(cout, chr_len, num_pos);
    cout << '\n';

    cout << "Total % of " << chr_name << " covered (depth >= 10): ";
    output_percentage(cout, chr_len, num_pos_10);
    cout << '\n';

    cout << "Average depth: " << (depth_sum/chr_len) << '\n';

    cout << "Largest contig: " << contig_store.back() << '\n';
    
    if(gap_store.size() > 0) {
        cout << "Largest gap:    " << gap_store.back() << '\n';


        for (int i = 1; i < min(gap_store.size() - 1, 10); i++) {
            cout << toStr<int>(i+1) << "th Largest gap:    " << *(gap_store.end() - i - 1) << '\n';
        }
        // compute average gap length
        int64_t avg_len = 0;
        for (vector<int64_t>::iterator it = gap_store.begin(); it != gap_store.end(); it++) {
            avg_len += *it;
        }

        cout << "Average gap:    " << (avg_len - gap_store.back()) / ((double) gap_store.size() - 1.0) << '\n';

        cout << "Median gap:    " << *(gap_store.end()-(gap_store.size() / 2)) << '\n';

        cout << "--- Contig lengths -- \n";

        int64_t length_count[9] = {1ll, 10ll, 100ll, 1000ll, 10000ll, 100000ll, 1000000ll, 10000000ll, 100000000ll};
        int64_t mult = 10ll;
        for (int i = 0; i < 9; i++) {
            int64_t lower = length_count[i];
            int64_t upper = ((mult * length_count[i]) - 1);
            cout << lower << " - " << upper << ":";

            vector<int64_t>::iterator start = lower_bound(contig_store.begin(), contig_store.end(), lower);
            vector<int64_t>::iterator end = lower_bound(contig_store.begin(), contig_store.end(), upper);

            int counter = 0;
            for (vector<int64_t>::iterator it = start; it != end; it++) {
                if (*it <= upper && *it >= lower) {
                    counter++;
                }
            }

            cout << counter << '\n';

        }

        cout << "--- Gap lengths -- \n";


        for (int i = 0; i < 9; i++) {
            int64_t lower = length_count[i];
            int64_t upper = ((mult * length_count[i]) - 1);
            cout << lower << " - " << upper << ":";

            vector<int64_t>::iterator start = lower_bound(gap_store.begin(), gap_store.end(), lower);
            vector<int64_t>::iterator end = lower_bound(gap_store.begin(), gap_store.end(), upper);

            int counter = 0;
            for (vector<int64_t>::iterator it = start; it != end; it++) {
                if (*it <= upper && *it >= lower) {
                    counter++;
                }
            }

            cout << counter << '\n';

        }

        for (vector<int>::iterator it = depth_100.begin(); it != depth_100.end(); it++) {
            cerr << (*it) / 100.0 << '\n';
        }

    }
    
    return 0;
}

struct base_count_t{
        int64_t A;
        int64_t T;
        int64_t C;
        int64_t G;
        int64_t N;

        base_count_t(){
            A = 0;
            C = 0;
            T = 0;
            G = 0;
            N = 0;
        }
    };


int fastq_stats(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " fastq_stats <FASTQ_file>\n"
            + "Get stats for a FASTQ file";




    if (params.size() != 1) {
        cerr << usage_text << endl;
        return 3;
    }

    string filename = params[0];

    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    fast in(filename, false);

    map<int, int> read_len;
    map<int, int>::iterator read_len_it;
    int total_read_count = 0;

    vector<int64_t> qual_sum; // sum of the base qualities
    vector<base_count_t> base_sum; // sum of the base qualities

    while (!in.eof()) {
        fast_t line = in.get_next_record();

        if (line.seq.length() == 0) {
            continue;
        }

        int len = (int) line.seq.length();

        read_len_it = read_len.find(len);

        if (read_len_it == read_len.end()) {
            read_len.insert(pair<int, int>(len, 1));
        } else {
            read_len_it->second++;
        }


        if (line.qual_exist) {
            // count the quality scores
            int qual_len = (int) line.qual.size();

            while ((int) qual_sum.size() < qual_len) {
                qual_sum.push_back(0);
            }

            for (int i = 0; i < qual_len; i++) {
                qual_sum[i] += ((int64_t) line.qual[i] - 33);
            }
        }



        while ((int) base_sum.size() < len) {
            base_sum.push_back(base_count_t());
        }

        for (int i = 0; i < len; i++) {
            switch(line.seq[i]) {
                case 'A':
                    base_sum[i].A++;
                    break;
                case 'T':
                    base_sum[i].T++;
                    break;
                case 'C':
                    base_sum[i].C++;
                    break;
                case 'G':
                    base_sum[i].G++;
                    break;
                default:
                    base_sum[i].N++;

            }
        }


        total_read_count++;
    }

    int comu_read_count = 0;

    for (read_len_it = read_len.begin(); read_len_it != read_len.end(); read_len_it++) {
        comu_read_count += read_len_it->second;
        cout << read_len_it->first << "," << read_len_it->second << ",";
        output_percentage(cout, total_read_count, comu_read_count);
        cout << '\n';
    }

    cout << "Avg base quality (base_idx,avg_qual):" << '\n';

    for(int i = 0;i<(int)qual_sum.size();i++){
        cout << (i+1) << ',' << ((double)qual_sum[i])/total_read_count << '\n';
    }


    cout << "Base ratio (base_idx,A,T,C,G,N):" << '\n';

    for(int i = 0;i<(int)base_sum.size();i++){

        cout << (i+1) << ',';
        output_percentage(cout, total_read_count, base_sum[i].A);
        cout << ',';
        output_percentage(cout, total_read_count, base_sum[i].T);
        cout << ',';
        output_percentage(cout, total_read_count, base_sum[i].C);
        cout << ',';
        output_percentage(cout, total_read_count, base_sum[i].G);
        cout << ',';
        output_percentage(cout, total_read_count, base_sum[i].N);
        cout << '\n';
    }

    return 0;
}


int fastq_quals(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " fastq_quals <FASTQ_file> <offset>\n"
            + "Offset typically 33 (illumina 1.8+,sanger) or 64 (illumina 1.3+)\n"
            + "Get quals from a FASTQ file normalized";


    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

    string filename = params[0];
    int offset = strTo<int>(params[1]);

    
    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    fast in(filename, false);


    while (!in.eof()) {
        fast_t line = in.get_next_record();

        if (line.seq.length() == 0) {
            continue;
        }

        if(!line.qual_exist){
            cerr << "Error: no qualities..." << '\n';
            return 1;
        }
        
        int len = (int) line.seq.length();
        
        if(len != (int)line.qual.length() || (int)line.qual.length() == 0){
            cerr << "Error: mismatched quality len: " << line.name << '\n';
            continue;
        }
        
        cout << (int)(line.qual[0] - offset);
        
        for (int i = 1;i < (int)line.qual.length(); i++){
            cout << '\t' << (int)(line.qual[i] - offset) ;
        }
        
        cout << "\n";
                

    }


    return 0;
}


int sam_quals(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " sam_quals <SAM_file> <offset>\n"
            + "Offset typically 33 (illumina 1.8+,sanger) or 64 (illumina 1.3+)\n"
            + "Get quals from a SAM file normalized and orientation corrected";


    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

    string filename = params[0];
    int offset = strTo<int>(params[1]);

    
    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    infile.open(filename.c_str(), ios::in);
    string line = "";

    while (!infile.eof()) {
        getline(infile,line);
        
        trim2(line);
        
        if(line.length() < 1){
            continue;            
        }
        
        if(line[0] == '@'){
            continue;
        }

        vector<string> line_list = split(line);
        
        if(line_list.size() < 10){
            cerr << "Bad line: " << line << '\n';
            continue;
        }
        
        // get flag
        int flag = strTo<int>(line_list[1]);
        
        
        // get quality value
        string quals = line_list[10];
        
        // invert if necessary
        bool invert = false;
        
        if(!((flag & 4) && (flag & 8))){
            // at least on read is mapped
            if(flag & 64){
                // first in pair
                invert = flag & 16;
            }else if(flag & 128){
                // second in pair
                invert = (!(flag & 16));
            }
            
        }
        
        if(invert){
            reverse(quals.begin(),quals.end());
        }
        
        // print it out
        
        cout << (int)(quals[0] - offset);
        
        for (int i = 1;i < (int)quals.length(); i++){
            cout << '\t' << (int)(quals[i] - offset) ;
        }
        
        cout << "\n";
                

    }


    return 0;
}


int replace_sam_quals(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " replace_sam_quals [-e] <SAM_file> <offset> <quals_file>\n"
            + "    -e  --  Take log\n"
            + "Offset typically 33 (illumina 1.8+,sanger) or 64 (illumina 1.3+)\n"
            + "replace quals in SAM file with new ones!";


    if (params.size() < 3 || params.size() > 4) {
        cerr << usage_text << endl;
        return 3;
    }
    
    bool take_log = false;
    int param_offset = 0;
    
    if(params.size() == 4){
        if(params[0] == "-e"){
            take_log = true;
            param_offset = 1;
        }else{
            cerr << usage_text << endl;
            return 3;
        }
    }

    string filename = params[0 + param_offset];
    int offset = strTo<int>(params[1 + param_offset]);
    string quals_filename = params[2 + param_offset];
    
    ifstream infile;

    infile.open(quals_filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open quals file: " << quals_filename << '\n';
        return 3;
    }

    infile.close();
    

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    
    ifstream infile2;
    infile2.open(quals_filename.c_str(), ios::in);
    string line = "";
    string line2 = "";
    
    infile.open(filename.c_str(), ios::in);
    

    while (!infile.eof() && !infile2.eof()) {
        getline(infile,line);
        
        trim2(line);
        
        if(line.length() < 1){
            continue;            
        }
        
        if(line[0] == '@'){
            cout << line << '\n';
            continue;
        }

        vector<string> line_list = split(line);
        
        if(line_list.size() < 10){
            cerr << "Bad line: " << line << '\n';
            cout << line << '\n';
            continue;
        }
        
        // read the quals
        getline(infile2,line2);
        trim2(line2);
        
        if(line2.length() < 1){
            cerr << "Bad quals file..." << '\n';
            cerr << line << '\n';
            cerr << line2 << '\n';
            return 1;
        }
        
        vector<string> line_list2 = split(line2);
        
        string new_quals = "";
        for(vector<string>::iterator ll = line_list2.begin();
                ll != line_list2.end();ll++){
            
            int val = 0;
            if(take_log){
                val = (int)round(-10 * log(strTo<double>(*ll)));
            }else{
                val = (int)round(strTo<double>(*ll));
            }
            
            if(val <= 0 || val >= 255){
                cerr << "ERROR: bad quals("<< val <<"|" << *ll <<"): " << line2 << '\n';
            }
            
            if(val < 2){
                val = 2;
            }else if(val > 40){
                val = 40;
            }
            
            val += offset;
            
            
            
            
            
            new_quals.push_back((char)val);
            
        }
        
        // get flag
        int flag = strTo<int>(line_list[1]);
        
       
        
        // invert if necessary
        bool invert = false;
        
        if(!((flag & 4) && (flag & 8))){
            // at least on read is mapped
            if(flag & 64){
                // first in pair
                invert = flag & 16;
            }else if(flag & 128){
                // second in pair
                invert = (!(flag & 16));
            }
            
        }
        
        if(invert){
            reverse(new_quals.begin(),new_quals.end());
        }
        
        // print it out
        
        line_list[10] = new_quals;
        
        cout << line_list[0];
        
        for (int i = 1;i < (int)line_list.size(); i++){
            cout << '\t' << line_list[i] ;
        }
        
        cout << "\n";
                

    }
    
    infile.close();
    infile2.close();


    return 0;
}

// temp thing for annoying thing
int replace_sam_quals_diff(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " replace_sam_quals_diff [-e] <SAM_file> <window_size> <offset> <quals_file>\n"
            + "    -e    -- de-Exponentiate first\n"
            + "replace quals in SAM file with new ones, but diff...!";


    if (params.size() < 4 || params.size() > 5) {
        cerr << usage_text << endl;
        return 3;
    }
    
    bool take_log = false;
    int offset = 0;
    
    if(params.size() == 5){
        if(params[0] == "-e"){
            take_log = true;
            offset = 1;
        }else{
            cerr << usage_text << endl;
            return 3;
        }
    }

    string filename = params[0 + offset];
    int wsize  = strTo<int>(params[1 + offset]);
    int qual_offset  = strTo<int>(params[2 + offset]);
    string quals_filename = params[3 + offset];
    
    ifstream infile;

    infile.open(quals_filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open quals file: " << quals_filename << '\n';
        return 3;
    }

    infile.close();
    

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    
    ifstream infile2;
    infile2.open(quals_filename.c_str(), ios::in);
    string line = "";
    string line2 = "";
    
    infile.open(filename.c_str(), ios::in);
    

    while (!infile.eof() && !infile2.eof()) {
        getline(infile,line);
        
        trim2(line);
        
        if(line.length() < 1){
            continue;            
        }
        
        if(line[0] == '@'){
            cout << line << '\n';
            continue;
        }

        vector<string> line_list = split(line);
        
        if(line_list.size() < 10){
            cerr << "Bad line: " << line << '\n';
            cout << line << '\n';
            continue;
        }
        
        // read the quals
        getline(infile2,line2);
        trim2(line2);
        
        if(line2.length() < 1){
            cerr << "Bad quals file..." << '\n';
            cerr << line << '\n';
            cerr << line2 << '\n';
            return 1;
        }
        
        vector<string> line_list2 = split(line2);
        
        string new_quals = "";
        for(vector<string>::iterator ll = line_list2.begin();
                ll != line_list2.end();ll++){

            int val = 0;
            if(take_log){
                val = (int)round(-10.0*log(strTo<double>(*ll)));
            }else{
                val = (int)round(strTo<double>(*ll));
            }
            new_quals.push_back((char)val);
        }
        
        // get flag
        int flag = strTo<int>(line_list[1]);
        
       
        
        // invert if necessary
        bool invert = false;
        
        if(!((flag & 4) && (flag & 8))){
            // at least on read is mapped
            if(flag & 64){
                // first in pair
                invert = flag & 16;
            }else if(flag & 128){
                // second in pair
                invert = (!(flag & 16));
            }
            
        }
        
        if(invert){
            reverse(new_quals.begin(),new_quals.end());
        }
        
        // do the diff
        int llen = (int)line_list[10].size();
        int num_wind = llen/wsize;
        
        bool last_wind = false;
        if(num_wind * wsize != llen){
            last_wind = true;
        }
        
        for(int i = 0;i<num_wind;i++){
            
            int offset = i*wsize;
            int noffset = i*(wsize-1);
            for(int j = 1;j<wsize;j++){
                int new_val = line_list[10][offset+j-1]+new_quals[noffset+j-1];
                
                if(new_val-qual_offset <0){
                    new_val = qual_offset;
                }
                
                if(new_val-qual_offset >40){
                    new_val = qual_offset+40;
                }

                
                line_list[10][offset+j] = new_val;
            }
        }
        
        if(last_wind){
            int old_wsize = wsize;
            wsize = llen - num_wind * old_wsize;
            
            int offset = num_wind*old_wsize;
            int noffset = num_wind*(old_wsize-1);
            for(int j = 1;j<wsize;j++){
                
                int new_val = line_list[10][offset+j-1]+new_quals[noffset+j-1];
                
                if(new_val-qual_offset <0){
                    new_val = qual_offset;
                }
                
                if(new_val-qual_offset >40){
                    new_val = qual_offset+40;
                }

                line_list[10][offset+j] = new_val;
            }
        }
        
        cout << line_list[0];
        
        for (int i = 1;i < (int)line_list.size(); i++){
            cout << '\t' << line_list[i] ;
        }
        
        cout << "\n";
                

    }
    
    infile.close();
    infile2.close();


    return 0;
}



int quant_sam_quals(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " quant_sam_quals <SAM_file> <offset> <quant_levels>\n"
            + "Offset typically 33 (illumina 1.8+,sanger) or 64 (illumina 1.3+)\n"
            + "replace quals in SAM file with new quantized ones!";


    const int LOW = 4;
    const int HIGH = 40;
    
    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    string filename = params[0];
    int offset = strTo<int>(params[1]);
    int quant_levels = strTo<int>(params[2]);
    
    if(quant_levels < 2 || quant_levels >= 15){
        cerr << "Error: too few/many quant levels\n";
        return 2;
    }
    
    vector<int> levels(quant_levels,0);
    
    levels[0] = LOW;
    levels[levels.size()-1] = HIGH;

    if (quant_levels >= 3) {
        
        int gap = (int)(((double)HIGH-LOW)/((double)quant_levels-1));
        
        int len = (int) levels.size() - 1;
        for (int i = 1; i < len; i++) {
            levels[i] = levels[i-1]+gap;
        }
    }
    
    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();


    string line = "";
    
    infile.open(filename.c_str(), ios::in);
    

    while (!infile.eof()) {
        getline(infile,line);
        
        trim2(line);
        
        if((int)line.length() < 1){
            continue;            
        }
        
        if(line[0] == '@'){
            cout << line << '\n';
            continue;
        }

        vector<string> line_list = split(line);
        
        if(line_list.size() < 10){
            cerr << "Bad line: " << line << '\n';
            cout << line << '\n';
            continue;
        }
        

        // quantize
        string new_quals = line_list[10];
       
        int len = (int)new_quals.length();
        for (int i = 0;i<len;i++){
            int qual = new_quals[i]-offset;
            
            
            int diff = abs(qual-levels[0]);
            int val = 0;
            for(int j = 1;j<quant_levels;j++){
                int new_diff = abs(qual-levels[j]);
                
                if(new_diff >= diff){
                    val = j-1;
                    break;
                }else{
                    val = j;
                    diff = new_diff;
                }
            }
            
            new_quals[i] = levels[val]+offset;
            
        }
       
        
        // print it out
        
        line_list[10] = new_quals;
        
        cout << line_list[0];
        
        for (int i = 1;i < (int)line_list.size(); i++){
            cout << '\t' << line_list[i] ;
        }
        
        cout << "\n";
                

    }
    
    infile.close();

    return 0;
}


int sam_stats(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " sam_stats <SAM_file>\n"
            + "Get stats for a SAM file";

    if (params.size() != 1) {
        cerr << usage_text << endl;
        return 3;
    }

    string filename = params[0];

    ifstream infile;

    infile.open(filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << filename << '\n';
        return 3;
    }

    infile.close();

    cout << "\n\n============\n";
    cout << "Filename:" << filename << '\n';

    ifstream in(filename.c_str(), ios::in);
    string line;

    // some stats to remember
    int num_reads = 0;
    int num_reads_aligned = 0;
    map<int,int> mapq_count; // count for each MAPQ
    map<int,int> edit_count; // count for each edit distance
    map<int,int> gap_ext_count; // count for each gap_extend
    map<int,int> mm_count; // count for each mismatch


    map<string,int> chr_count;    // count aligned to each chromosome
    map<string,int> chr_count_20; // count aligned to each chromosome (MAPQ>20)


    while (!in.eof()) {
        getline(in,line);

        trim2(line);

        if (line.length() == 0) {
            continue;
        }

        //cerr << line[0] << "\n";

        if(line[0] == '@') continue;

        vector<string> line_list = split(line);
        map<int,int>::iterator it;




        if(line_list.size() < 5){
            continue;
        }

        //cerr << line_list[5] << ':';
        //cerr << line_list[5][0] << ':';

        if(line_list[5][0] != '*'){
            num_reads_aligned++;

            //cerr << num_reads_aligned << '\n';


            int mapq = strTo<int>(line_list[4]);

            it = mapq_count.find(mapq);

            if(it == mapq_count.end()){
                mapq_count.insert(pair<int,int>(mapq,1));
            }else{
                it->second++;
            }

            for(int i = 11;i<(int)line_list.size();i++){
                if(line_list[i].size() >= 6){
                    if(line_list[i][0] == 'N' && line_list[i][1] == 'M'){
                        int edit = strTo<int>(line_list[i].substr(5));

                        it = edit_count.find(edit);

                        if (it == edit_count.end()) {
                            edit_count.insert(pair<int, int>(edit, 1));
                        } else {
                            it->second++;
                        }
                    }else if(line_list[i][0] == 'X' && line_list[i][1] == 'G'){
                        int ext = strTo<int>(line_list[i].substr(5));

                        it = gap_ext_count.find(ext);

                        if (it == gap_ext_count.end()) {
                            gap_ext_count.insert(pair<int, int>(ext, 1));
                        } else {
                            it->second++;
                        }
                    }else if(line_list[i][0] == 'X' && line_list[i][1] == 'M'){
                        int mm = strTo<int>(line_list[i].substr(5));

                        it = mm_count.find(mm);

                        if (it == mm_count.end()) {
                            mm_count.insert(pair<int, int>(mm, 1));
                        } else {
                            it->second++;
                        }
                    }
                }
            }

            map<string,int>::iterator it2;


            it2 = chr_count.find(line_list[2]);

            if(it2 == chr_count.end()){
                chr_count.insert(pair<string,int>(line_list[2],1));
            }else{
                it2->second++;
            }


            if (mapq >= 20) {
                it2 = chr_count_20.find(line_list[2]);

                if (it2 == chr_count_20.end()) {
                    chr_count_20.insert(pair<string, int>(line_list[2], 1));
                } else {
                    it2->second++;
                }
            }

        }





        num_reads++;
    }


    // output everything

    cout << "Total Reads:" << num_reads << '\n';
    cout << "Total Reads Aligned:" << num_reads_aligned << '\n';
    output_percentage(cout, num_reads, num_reads_aligned);
    cout << '\n';

    vector<chr_idx_t> temp_vec;

    for(map<string,int>::iterator it = chr_count.begin();it != chr_count.end();it++){

        chr_idx_t contents;

        contents.name = it->first;
        contents.count = it->second;

        temp_vec.push_back(contents);
    }

    sort(temp_vec.begin(),temp_vec.end());

    cout << "Chromosome count (chr_name,count): " << '\n';
    cout << temp_vec.size() << '\n';
    for(vector<chr_idx_t>::iterator it = temp_vec.begin();it != temp_vec.end();it++){
        cout << it->name << ',' << it->count << '\n';
    }


    temp_vec.clear();
    for(map<string,int>::iterator it = chr_count_20.begin();it != chr_count_20.end();it++){

        chr_idx_t contents;

        contents.name = it->first;
        contents.count = it->second;

        temp_vec.push_back(contents);
    }

    sort(temp_vec.begin(),temp_vec.end());

    cout << "MAPQ 20 Chromosome count[MAPQ>=20] (chr_name,count): " << '\n';
    cout << temp_vec.size() << '\n';
    for(vector<chr_idx_t>::iterator it = temp_vec.begin();it != temp_vec.end();it++){
        cout << it->name << ',' << it->count << '\n';
    }


    cout << "MAPQ (MAPQ,count):" << '\n';
    cout << mapq_count.size() << '\n';
    for(map<int,int>::iterator it = mapq_count.begin();it != mapq_count.end();it++){
        cout << it->first << ',' << it->second << '\n';
    }

    cout << "Edit distance (edits,count):" << '\n';
    cout << edit_count.size() << '\n';
    for(map<int,int>::iterator it = edit_count.begin();it != edit_count.end();it++){
        cout << it->first << ',' << it->second << '\n';
    }

    cout << "Gap extension length (ext_len,count):" << '\n';
    cout <<  gap_ext_count.size() << '\n';
    for(map<int,int>::iterator it = gap_ext_count.begin();it != gap_ext_count.end();it++){
        cout << it->first << ',' << it->second << '\n';
    }

    cout << "Mismatch number (mismatch,count):" << '\n';
    cout << mm_count.size() << '\n';
    for(map<int,int>::iterator it = mm_count.begin();it != mm_count.end();it++){
        cout << it->first << ',' << it->second << '\n';
    }

    return 0;
}

int fastq_filter(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " fastq_filter <phred_score> <percent> <output_prefix> <first_pair> [<second_pair>]\n"
            + "    phred_score    -- Phred score (reads must be Sanger quals)\n"
            + "    percent        -- Percent of bases with the specified phred score\n"
            + "    output_prefix  -- Output will be output_prefix.fq or output_prefix_1.fq,output_prefix_2.fq\n"
            + "    first_pair     -- First pair or single-end FASTQ file\n"
            + "    second_pair    -- Second pair of FASTQ file [optional]\n"
            + "For paired-end reads only one end needs to satisfy this criteria";

    if (params.size() < 4 || params.size() > 5) {
        cerr << usage_text << endl;
        return 3;
    }

    int num_pair = 1;
    if (params.size() == 5) {
        num_pair = 2;
    }


    string temp;

    temp = params[0];

    int phred_score = strTo<int>(temp);

    if (phred_score < 0) {
        cerr << "Error: Invalid Phred score: " << phred_score << '\n';
        return 3;
    }


    temp = params[1];

    double percent = strTo<double>(temp);

    if (percent < 0 || percent > 100) {
        cerr << "Error: Invalid percent: " << percent << '\n';
        return 3;
    }

    string output_prefix = params[2];

    string filename[2];

    for (int i = 0; i < num_pair; i++) {
        filename[i] = params[i + 3];

        ifstream infile;

        infile.open(filename[i].c_str(), ios::in);

        if (!infile.is_open()) {
            cerr << "Error: cannot open: " << filename[i] << '\n';
            return 1;
        }

        infile.close();
    }

    ofstream out[2];
    ofstream out_reject[2];

    if (num_pair == 1) {
        temp = output_prefix + ".fq";

        out[0].open(temp.c_str(), ios::out);

        temp = output_prefix + "_reject.fq";

        out_reject[0].open(temp.c_str(), ios::out);

    } else {
        temp = output_prefix + "_1.fq";

        out[0].open(temp.c_str(), ios::out);

        temp = output_prefix + "_2.fq";

        out[1].open(temp.c_str(), ios::out);


        temp = output_prefix + "_reject_1.fq";

        out_reject[0].open(temp.c_str(), ios::out);

        temp = output_prefix + "_reject_2.fq";

        out_reject[1].open(temp.c_str(), ios::out);
    }


    // count the number of reads in the file

    fast in[2];

    for (int i = 0; i < num_pair; i++) {
        in[i].open(filename[i]);
    }

    fast_t line[2];

    while (!in[0].eof()) {
        for (int i = 0; i < num_pair; i++) {
            line[i] = in[i].get_next_record();

            if (!line[i].qual_exist) {
                cerr << "Error: Input must be a FASTQ file" << '\n';

                return 2;
            }
        }

        if (line[0].seq.length() == 0) {

            if (num_pair == 2) {
                if (line[1].seq.length() != 0) {
                    cerr << "Error: FASTQ files are different length" << '\n';

                    return 1;
                }
            }

            continue;
        }

        if (num_pair == 2) {
            if (line[1].seq.length() == 0) {
                cerr << "Error: FASTQ files are different length" << '\n';

                return 1;
            }
        }


        bool passed = false; // if read passes criteria then we output

        for (int i = 0; i < num_pair; i++) {
            string::iterator sit = line[i].qual.begin();

            int count_total = line[i].qual.length();
            int count_passed = 0;

            while (sit != line[i].qual.end()) {
                int qual = (int) *sit - 33;

                if (qual < 0) {
                    cerr << "Error invalid quality score: " << line[i].name << '\n';

                    return 1;
                }

                if (qual >= phred_score) {
                    count_passed++;
                }

                sit++;

            }

            if ((100.0 * (double) count_passed / (double) count_total) >= percent) {
                passed = true;
            }

        }


        if (passed) {
            for (int i = 0; i < num_pair; i++) {
                // output the reads

                out[i] << '@' << line[i].name;

                if (num_pair == 2) {
                    out[i] << '/' << (i + 1);
                }

                out[i] << '\n';

                out[i] << line[i].seq << '\n';

                out[i] << '+' << '\n';

                out[i] << line[i].qual << '\n';


            }
        } else {
            for (int i = 0; i < num_pair; i++) {
                // output the reads

                out_reject[i] << '@' << line[i].name;

                if (num_pair == 2) {
                    out_reject[i] << '/' << (i + 1);
                }

                out_reject[i] << '\n';

                out_reject[i] << line[i].seq << '\n';

                out_reject[i] << '+' << '\n';

                out_reject[i] << line[i].qual << '\n';


            }
        }

    }

    for (int i = 0; i < num_pair; i++) {
        in[i].close();

        out[i].close();
        out_reject[i].close();
    }



    return 0;
}


int fq_len_filter(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " fq_len_filter <min_len> <output_prefix> <first_pair> [<second_pair>]\n"
            + "    min_len        -- Pairs with one end shorter than this will be removed\n"
            + "    output_prefix  -- Output will be output_prefix.fq or output_prefix_1.fq,output_prefix_2.fq\n"
            + "    first_pair     -- First pair or single-end FASTQ file\n"
            + "    second_pair    -- Second pair of FASTQ file [optional]\n"
            + "For paired-end reads only one end needs to satisfy this criteria";

    if (params.size() < 3 || params.size() > 4) {
        cerr << usage_text << endl;
        return 3;
    }

    int num_pair = 1;
    if (params.size() == 4) {
        num_pair = 2;
    }


    string temp;

    temp = params[0];

    int min_len = strTo<int>(temp);

    if (min_len < 0) {
        cerr << "Error: Invalid min_len: " << min_len << '\n';
        return 3;
    }


    string output_prefix = params[1];

    string filename[2];

    for (int i = 0; i < num_pair; i++) {
        filename[i] = params[i + 2];

        ifstream infile;

        infile.open(filename[i].c_str(), ios::in);

        if (!infile.is_open()) {
            cerr << "Error: cannot open: " << filename[i] << '\n';
            return 1;
        }

        infile.close();
    }

    ofstream out[2];
    ofstream out_reject[2];

    if (num_pair == 1) {
        temp = output_prefix + ".fq";

        out[0].open(temp.c_str(), ios::out);

        temp = output_prefix + "_reject.fq";

        out_reject[0].open(temp.c_str(), ios::out);

    } else {
        temp = output_prefix + "_1.fq";

        out[0].open(temp.c_str(), ios::out);

        temp = output_prefix + "_2.fq";

        out[1].open(temp.c_str(), ios::out);


        temp = output_prefix + "_reject_1.fq";

        out_reject[0].open(temp.c_str(), ios::out);

        temp = output_prefix + "_reject_2.fq";

        out_reject[1].open(temp.c_str(), ios::out);
    }


    // count the number of reads in the file

    fast in[2];

    for (int i = 0; i < num_pair; i++) {
        in[i].open(filename[i]);
    }

    fast_t line[2];

    while (!in[0].eof()) {
        for (int i = 0; i < num_pair; i++) {
            line[i] = in[i].get_next_record();

            if (!line[i].qual_exist) {
                cerr << "Error: Input must be a FASTQ file" << '\n';

                return 2;
            }
        }

        if (line[0].seq.length() == 0) {

            if (num_pair == 2) {
                if (line[1].seq.length() != 0) {
                    cerr << "Error: FASTQ files are different length" << '\n';

                    return 1;
                }
            }

            continue;
        }

        if (num_pair == 2) {
            if (line[1].seq.length() == 0) {
                cerr << "Error: FASTQ files are different length" << '\n';

                return 1;
            }
        }


        bool passed = true; // if read passes criteria then we output

        for (int i = 0; i < num_pair; i++) {
            if((int)line[i].seq.length() < min_len){
                passed = false;
            }

        }


        if (passed) {
            for (int i = 0; i < num_pair; i++) {
                // output the reads

                out[i] << '@' << line[i].name;

                if (num_pair == 2) {
                    out[i] << '/' << (i + 1);
                }

                out[i] << '\n';

                out[i] << line[i].seq << '\n';

                out[i] << '+' << '\n';

                out[i] << line[i].qual << '\n';


            }
        } else {
            for (int i = 0; i < num_pair; i++) {
                // output the reads

                out_reject[i] << '@' << line[i].name;

                if (num_pair == 2) {
                    out_reject[i] << '/' << (i + 1);
                }

                out_reject[i] << '\n';

                out_reject[i] << line[i].seq << '\n';

                out_reject[i] << '+' << '\n';

                out_reject[i] << line[i].qual << '\n';


            }
        }

    }

    for (int i = 0; i < num_pair; i++) {
        in[i].close();

        out[i].close();
        out_reject[i].close();
    }



    return 0;
}


int sample_reads(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " sample_reads <num_reads> <output_prefix> <first_pair> [<second_pair>]\n"
            + "    num_reads      -- Number of reads to extract\n"
            + "    output_prefix  -- Output will be output_prefix.fq or output_prefix_1.fq,output_prefix_2.fq\n"
            + "    first_pair     -- First pair or single-end FASTQ/A file\n"
            + "    second_pair    -- Second pair of FASTQ/A file [optional]\n"
            + "Uniformly randomly extract reads from a FASTQ file";

    if (params.size() < 3 || params.size() > 4) {
        cerr << usage_text << endl;
        return 3;
    }

    int num_pair = 1;
    if (params.size() == 4) {
        num_pair = 2;
    }


    string temp;

    temp = params[0];

    int num_reads = strTo<int>(temp);

    if (num_reads < 0) {
        cerr << "Error: Invalid number of reads: " << num_reads << '\n';
        return 3;
    }

    string output_prefix = params[1];

    string filename[2];

    for (int i = 0; i < num_pair; i++) {
        filename[i] = params[i + 2];

        ifstream infile;

        infile.open(filename[i].c_str(), ios::in);

        if (!infile.is_open()) {
            cerr << "Error: cannot open: " << filename[i] << '\n';
            return 1;
        }

        infile.close();
    }

    ofstream out[2];

    if (num_pair == 1) {
        temp = output_prefix + ".fq";

        out[0].open(temp.c_str(), ios::out);
    } else {
        temp = output_prefix + "_1.fq";

        out[0].open(temp.c_str(), ios::out);

        temp = output_prefix + "_2.fq";

        out[1].open(temp.c_str(), ios::out);
    }


    // count the number of reads in the file

    fast in[2];


    for (int i = 0; i < num_pair; i++) {
        in[i].open(filename[i]);
    }


    fast_t line;


    vector<int64_t> indexes;
    MT_random randgen;

    int64_t k = 0;

    while (!in[0].eof()) {

        line = in[0].get_next_record();

        if (line.seq.length() == 0) {
            continue;
        }

        k = k + 1;

        if (k <= num_reads) {
            indexes.push_back(k);
        } else {
            if (randgen.genrand_bern(num_reads / (double) k)) {
                // we replace a random current sample
                int idx = randgen.genrand_int_range(0, num_reads - 1);

                indexes[idx] = k;
            }
        }


    }

    sort(indexes.begin(), indexes.end());

    in[0].close();

    in[0].open(filename[0]);

    // pick out the indexes of all the num_reads
    if (num_reads >= k) {
        cerr << "Error: the number of reads in the file is not greater than num_reads" << '\n';
        return 3;
    }

    // now iterate through both files and output the reads

    for (int i = 0; i < num_pair; i++) {

        int64_t curr_idx = 0;

        int64_t k = 0;

        while (!in[i].eof()) {

            line = in[i].get_next_record();

            if (line.seq.length() == 0) {
                continue;
            }

            k = k + 1;

            if (k == indexes[curr_idx]) {

                if (curr_idx < num_reads) {
                    curr_idx++;
                }

                // do the outputting

                if (line.qual_exist) {
                    if (num_pair == 2) {
                        out[i] << "@" << line.name << "/" << (i + 1) << '\n';
                    } else {
                        out[i] << "@" << line.name << '\n';
                    }
                    out[i] << line.seq << '\n';
                    out[i] << "+\n";
                    out[i] << line.qual << '\n';
                } else {
                    if (num_pair == 2) {
                        out[i] << ">" << line.name << "/" << (i + 1) << '\n';
                    } else {
                        out[i] << ">" << line.name << '\n';
                    }
                    out[i] << line.seq << '\n';
                }
            }


        }


        in[i].close();
        out[i].close();
    }

    return 0;

}

int merge_sam(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " merge_sam <sam_file> <other_sam_file>\n"
            + "    sam_file       -- First SAM file, the headers will be preserved\n"
            + "    other_sam_file -- Second SAM file\n"
            + "Merge two SAM files, probably generated from hybrid mode alignment";

    if (params.size() != 2) {
        cerr << usage_text << endl;
        return 3;
    }

    string sam1_name = params[0];
    string sam2_name = params[1];


    ifstream sam1;

    sam1.open(sam1_name.c_str(), ios::in);

    if (!sam1.is_open()) {
        cerr << "Error cannot open: " << sam1_name << '\n';
        return 1;
    }

    ifstream sam2;

    sam2.open(sam2_name.c_str(), ios::in);

    if (!sam2.is_open()) {
        cerr << "Error cannot open: " << sam2_name << '\n';
        return 1;
    }


    // write out first SAM file

    string line;

    while (!sam1.eof()) {
        getline(sam1, line);

        trim2(line);

        if (line.length() == 0) {
            continue;
        }

        cout << line << '\n';
    }


    while (!sam2.eof()) {
        getline(sam2, line);

        trim2(line);

        if (line.length() == 0) {
            continue;
        }

        if (line[0] == '@') {
            continue;
        }

        cout << line << '\n';
    }

    return 0;
}

int qual_convert(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " qual_convert <original_format> <new_format> <input_file>\n"
            + "-- Format --\n"
            + "    S - Sanger or Illumina 1.8+ format [Phred + 33]\n"
            + "    X - Solexa format [Solexa + 64]\n"
            + "    I - Illumina 1.3+ or 1.5+ format [Phred + 64]\n"
            + "-- Params --\n"
            + "    seq       -- Read sequence\n"
            + "    kmer_size -- k-mer size to analyze\n"
            + "Covert the quality format of FASTQ files";

    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    string orig = params[0];
    string newf = params[1];
    string in_filename = params[2];

    cerr << "Original format: ";
    if (orig == "S") {
        cerr << "Sanger\n";
    } else if (orig == "X") {
        cerr << "Solexa\n";
    } else if (orig == "I") {
        cerr << "Illumina 1.3+\n";
    } else {
        cerr << "Error: Wrong original format!" << '\n';
        return 3;
    }


    cerr << "New format: ";
    if (newf == "S") {
        cerr << "Sanger\n";
    } else if (newf == "X") {
        cerr << "Error: Cannot convert to Solexa format\n";
        return 3;
    } else if (newf == "I") {
        cerr << "Illumina 1.3+\n";
    } else {
        cerr << "Error: Wrong new format!" << '\n';
        return 3;
    }

    if (orig == newf) {
        cerr << "Error: Original and New format are the same" << '\n';
        return 3;
    }

    ifstream infile;

    infile.open(in_filename.c_str(), ios::in);

    if (!infile.is_open()) {
        cerr << "Error: cannot open input file: " << in_filename << '\n';
        return 3;
    }

    infile.close();


    fast in(in_filename, false);

    while (!in.eof()) {
        fast_t line = in.get_next_record();

        if (line.seq.length() == 0) {
            continue;
        }

        if (!line.qual_exist) {
            cerr << "Error: File is not a FASTQ file" << '\n';
            return 3;
        }

        // Convert to Sanger format
        if (orig[0] == 'X') {
            if (solexatosanger(line.qual) != 0) {
                cerr << "Error: line " << line.name << "," << line.qual << '\n';
                return 3;
            }

        } else if (orig[0] == 'I') {
            if (ill_13tosanger(line.qual) != 0) {
                cerr << "Error: line " << line.name << "," << line.qual << '\n';
                return 3;
            }
        }// Convert from Sanger to new format
        else if (newf[0] == 'I') {
            if (sangertoill_13(line.qual) != 0) {
                cerr << "Error: line " << line.name << "," << line.qual << '\n';
                return 3;
            }
        }

        // write out the line
        cout << "@" << line.name << '\n';
        cout << line.seq << '\n';
        cout << "+\n";
        cout << line.qual << '\n';

    }

    return 0;
}

int raw2fasta(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " raw2fasta <raw_file>"
            + "    raw_file  -- Input text file, one line per read"
            + "Converts raw reads into FASTA file, Outputs to stdout";


    int read_index = 1;
    string line;



    if (params.size() != 1) {
        cerr << usage_text << endl;
        return 3;
    }

    string raw_filename = params[0];

    ifstream raw_file;

    raw_file.open(raw_filename.c_str());
    if (!raw_file.is_open()) {
        cerr << "ERROR: Cannot open RAW file, " << raw_filename << endl;
        return 1;
    }

    while (!raw_file.eof()) {
        getline(raw_file, line);

        trim2(line);

        if (line.length() == 0) {
            continue;
        }

        cout << '>' << read_index << '\n';
        cout << line << '\n';

        read_index++;

    }

    return 0;
}

int subseq(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " subseq <chr_file> <start> <end>\n"
            + "    chr_file -- FASTA file containing only one sequence\n"
            + "    start    -- one based start position (inclusive)\n"
            + "    end      -- one based end position (inclusive)\n"
            + "Extract a sub-sequence from a FASTA file, must be only 1 sequence in the FASTA file";

    string genome_filename;
    ifstream genome_file;

    string chr_line = "";
    string chr_name;

    string line;
    string temp;



    if (params.size() != 3) {
        cerr << usage_text << endl;
        return 3;
    }

    genome_filename = params[0];

    int64_t a = strTo<int64_t > (params[1]);
    int64_t b = strTo<int64_t > (params[2]);

    if (b < a) {
        cerr << "ERROR: End before Start!";
        return 2;
    }

    // load the genome into memory

    genome_file.open(genome_filename.c_str(), ios::in);

    if (!genome_file.is_open()) {
        cout << "ERROR: I'm sorry, the genome file cannot be opened. " << endl;
        exit(1);
    }

    getline(genome_file, chr_name);
    chr_name.erase(0, 1);
    cout << "Chromosone Name: " << chr_name << endl;

    while (!genome_file.eof()) {
        getline(genome_file, line);

        trim2(line);

        if (line.length() == 0)
            continue;

        chr_line.append(line);
    }

    cout << "Chromosone Length: " << chr_line.length() << endl;

    cout << "Parameters: [" << a << "," << b << "]" << endl;

    genome_file.close();


    // Display the results
    // NOTE: Need to finish this

    temp = chr_line.substr(a - 1, b - a + 1);
    transform(temp.begin(), temp.end(), temp.begin(), (int(*)(int))toupper); // make uppercase

    cout << temp << endl;
    compleseq(temp);
    cout << temp << endl;


    return 0;

}


int select_haps(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " select_haps <hap_file> <sample_file> <legend_file> <test_size> [seed]\n"
            + "    hap_file           -- File with the haplotypes\n"
            + "    sample_file        -- The sample file associated with the hap file\n"
            + "    legend_file        -- The legend file associated with the hap file\n"
            + "    test_size          -- The size of test set (each test sample needs 2 individuals)\n"
            + "    seed               -- Optional seed for random number generator\n"
            + "Randomly select some haplotypes and generate test and train set";



    if (params.size() < 4 || params.size() > 5) {
        cerr << usage_text << endl;
        return 3;
    }

    string hap_filename = params[0];
    string sample_filename = params[1];
    string legend_filename = params[2];
    int number_test = strTo<int>(params[3]);
    
    if(number_test <= 0){
        cerr << "Must have at least one test\n";
        return 2;
    }
    
    uint64_t seed = 0;
    if(params.size() == 5){
        seed = strTo<int>(params[3]);
    }
    
    MT_random randgen(seed);

    { // Check each of the files can be opened
        ifstream hap_file(hap_filename.c_str(), ios::in);
        if (!hap_file.is_open()) {
            cerr << "Cannot open: " << hap_filename << '\n';
            return 1;
        }
        hap_file.close();


        ifstream sample_file(sample_filename.c_str(), ios::in);
        if (!sample_file.is_open()) {
            cerr << "Cannot open: " << sample_filename << '\n';
            return 1;
        }
        sample_file.close();
    }


    // read the haplotypes into a matrix
    // also count the number of individuals
    vector < vector<char> > all_haps;
    {
        ifstream hap_file(hap_filename.c_str(), ios::in);

        string temp;
        
        while(getline(hap_file,temp)){
            trim2(temp);
            
            if(temp.length() == 0){
                continue;
            }
            all_haps.push_back(vector<char>());
            vector < vector<char> >::iterator it = all_haps.end();
            it--;
            for (int i = 0;i<(int)temp.length();i=i+4){
                switch(temp[i]){
                    case '0':
                        it->push_back(temp[i]);
                        break;
                    case '1':
                        it->push_back(temp[i]);
                        break;
                    case '-':
                        it->push_back(temp[i]);
                        break;
                    default:
                        cerr << "Error bad hap: " << temp << '\n';
                }
            }
        }
        hap_file.close();
    }
    
    // check if we have enough individuals
    int num_ind = (int)all_haps[0].size();
    int num_snp = (int)all_haps.size();
    
    if(num_ind < number_test*2){
        cerr << "Not enough individuals for testing: " << num_ind << '\n';
        return 2;
    }

    // read in the sample file
    
    vector<string> all_samples;
    {
        ifstream sample_file(sample_filename.c_str(), ios::in);
        string temp;
        while(getline(sample_file,temp)){
            trim2(temp);
            
            if(temp.length() == 0){
                continue;
            }
            all_samples.push_back(temp);
        }
        
        sample_file.close();
    }
    
    // read in the legend file
    vector<vector<string> > all_legend;
    all_legend.reserve(num_snp);
    {
        ifstream legend_file(legend_filename.c_str(), ios::in);
        string temp;
        while(getline(legend_file,temp)){
            trim2(temp);
            
            if(temp.length() == 0){
                continue;
            }
            
            vector<string> ll = split(temp);
            
            all_legend.push_back(ll);
        }
        
        legend_file.close();
    }
    
    
    
    // sample without replacement the test set indexes
    
    int num_test_ind = number_test*2;
    vector<int> indexes;

    {
        for (int k = 0;k<num_ind;k++) {

            if (k <= num_test_ind) {
                indexes.push_back(k);
            } else {
                if (randgen.genrand_bern(num_test_ind / (double) k)) {
                    // we replace a random current sample
                    int idx = randgen.genrand_int_range(0, num_test_ind - 1);
                    indexes[idx] = k;
                }
            }
        }
    }
    sort(indexes.begin(),indexes.end());
    
    
    // output the training set files
    {
        string temp = "train_"+hap_filename;
        
        ofstream hap_file(temp.c_str(), ios::out);
        
        for(int i = 0;i<num_snp;i++) {

            bool first_one = true;
            for (int k = 0; k < num_ind; k++) {

                if (binary_search(indexes.begin(), indexes.end(), k)) {
                    continue;
                }

                char h = all_haps[i][k];
                
                if(first_one) {
                    hap_file << h << ' ' << h;
                    first_one = false;
                } else {
                    hap_file << ' ' << h << ' ' << h;
                }
                
            }
            
            hap_file << '\n';
        }
        
        hap_file.close();
        
        temp = "train_"+sample_filename;
        
        ofstream sample_file(temp.c_str(), ios::out);

        for(int k = 0; k < num_ind; k++){
            
            if(binary_search(indexes.begin(),indexes.end(),k)){
                continue;
            }
            
            sample_file << all_samples[k] << '\n';
        }
        sample_file.close();
    }
    // output the true test samples
    vector<int> perm_indexes = indexes;
    random_shuffle(perm_indexes.begin(),perm_indexes.end());
    
    {
        string temp = "true_"+hap_filename;
        
        ofstream hap_file(temp.c_str(), ios::out);
        

        for (int i = 0; i < num_snp; i++) {
            bool first_one = true;
            for (int k = 0; k < num_ind; k++) {

                if (!binary_search(indexes.begin(), indexes.end(), k)) {
                    continue;
                }

                char h = all_haps[i][k];

                if (first_one) {
                    hap_file << h << ' ' << h;
                    first_one = false;
                } else {
                    hap_file << ' ' << h << ' ' << h;
                }

            }
            hap_file << '\n';
        }
        
        hap_file.close();
        
        temp = "true_"+sample_filename;
        
        ofstream sample_file(temp.c_str(), ios::out);
        for(int k = 0; k < num_ind; k++){
            
            if(!binary_search(indexes.begin(),indexes.end(),k)){
                continue;
            }
            
            sample_file << all_samples[k] << '\n';
        }
        sample_file.close();
    }
    
    // output the test samples for input into the program
    // need to output in gen/sample format
    {
        int test_idx = 0;
        int test_pair[2] = {0,0};
        
        vector<string> ll = split(hap_filename,'.');
        string prefix_name = hap_filename;
        if(ll.size() == 2){
            prefix_name = ll[0];
        }
        
        string temp = "test_"+prefix_name+".gen";

        ofstream gen_file(temp.c_str(), ios::out);

        for (int i = 0; i < num_snp; i++) {
            bool first_one = true;
            
            // output SNP details
            // always output chr24...
            gen_file << "24 ";
            
            gen_file << all_legend[i][0] << ' ' << all_legend[i][1] << ' ';
            gen_file << all_legend[i][2] << ' ' << all_legend[i][3] << ' ';
            

            for (int k = 0; k < num_ind; k++) {

                if (!binary_search(indexes.begin(), indexes.end(), k)) {
                    continue;
                }

                int p = test_idx % 2;

                test_pair[p] = k;

                if (p == 1) {
                    char h1 = all_haps[i][test_pair[0]];
                    char h2 = all_haps[i][test_pair[1]];
                    
                    if(first_one){
                        first_one = false;
                    }else{
                        gen_file << ' ';
                    }
                    if(h1 != h2){
                        gen_file << "0 1 0";
                    }else if(h1 == '0'){
                        gen_file << "1 0 0";
                    }else if(h1 == '1'){
                        gen_file << "0 0 1";
                    }else{
                        cerr << "strange hap: " << h1 << '\n';
                    }
                    
                }

                test_idx++;
            }
            
            gen_file << '\n';
        }
        gen_file.close();
        
        test_idx = 0;
        
        temp = "test_"+sample_filename; // this sample file is fake
        
        ofstream sample_file(temp.c_str(), ios::out);
        for(int k = 0; k < num_ind; k++){
            
            if(!binary_search(indexes.begin(),indexes.end(),k)){
                continue;
            }
            
            int p = test_idx % 2;
            
            if (p == 1) {
                sample_file << all_samples[k] << '\n';
            }
            test_idx++;
        }
        sample_file.close();
    }
    
    return 0;

}



int vcf_to_hap(vector<string> params) {
    string usage_text = "Usage: " + PROG_NAME + " vcf_to_hap <vcf_file>\n"
            + "    vcf_file -- Input VCF file\n"
            + "Convert a VCF file to haps, sample, legend";

    string vcf_filename;


    if (params.size() != 1) {
        cerr << usage_text << endl;
        return 3;
    }

    vcf_filename = params[0];
    
    {
        ifstream vcf_file(vcf_filename.c_str(),ios::in);
        if(!vcf_file.is_open()){
            cerr << "Could not open VCF file\n";
            return 1;
        }
        vcf_file.close();
    }
    
    vector<string> ind_list;
    int num_snp = 0;
    
    {
        string temp = "";
        ifstream vcf_file(vcf_filename.c_str(),ios::in);
        
        while(getline(vcf_file,temp)){
            trim2(temp);
            if(temp.length()<2){
                continue;
            }
            if(temp[0]=='#' && temp[1] == '#'){
                continue;
            }
            if(temp[0] == '#'){
                // we read in the individuals now
                vector<string> ll = split(temp);
                for(int i=9;i< (int)ll.size();i++){
                    ind_list.push_back(ll[i]);
                }
                
                continue;
            }
            
            num_snp++;
            
            if(num_snp % 10000 == 0){
                cerr << "Reading SNP: " << num_snp << '\n';
            }
        }
        
        vcf_file.close();
    }
    
    int num_ind = (int)ind_list.size();
    
    cerr << "Number of SNPs: " << num_snp << '\n';
    cerr << "Number of individuals: " << num_ind << '\n';
    cerr << "-------------\n";
    
    vector<vector<vector<int> > > snp_list;
    vector<vector<string> > lib_list;


    {  // initialise snp_list
        
        
        snp_list.reserve(num_snp);
        
    }
    
    { // read in the SNPs
        vector<int> empty;
        empty.push_back(-1);
        empty.push_back(-1);
        
        string temp = "";
        ifstream vcf_file(vcf_filename.c_str(),ios::in);
        
        int snp_idx = 0;
        
        int count = 0;
        
        while(getline(vcf_file,temp)){
            count++;
            if (count % 10000 == 0){
                cerr << "Reading line: " << count << '\n';
            }

            trim2(temp);
            if(temp.length()<2){
                continue;
            }
            if(temp[0]=='#'){
                continue;
            }
            
            vector<string> ll = split(temp);
            
            if(ll[4].length()>1){
                continue; // multi-allellic loci
            }
            
            vector<string> contents;
            contents.push_back(ll[0]); // chr_name
            contents.push_back(ll[1]); // chr_loc
            
            if(ll[2][0]=='.'){
                ll[2] = "SNP_"+toStr<int>(snp_idx);
            }
            contents.push_back(ll[2]); // snp_id
            contents.push_back(ll[3]); // allele1
            contents.push_back(ll[4]); // allele2
            
            cerr << "---" << ll[1] << "\n";
            
            snp_list.push_back(vector<vector<int> >());
            vector<vector<vector<int> > >::iterator it = snp_list.end();
            it--;
            it->reserve(num_ind);
            for (int j = 0; j < num_ind; j++) {
                it->push_back(empty);
            }
            
            for(int i = 9;i<(int)ll.size();i++){
                string snpl = ll[i];
                
                if(!(snpl[0] == '.')){
                    snp_list[snp_idx][i-9][0] = snpl[0]-48;
                    snp_list[snp_idx][i-9][1] = snpl[2]-48;

                    if (i == 9){
                        cerr << "SNP string: " << snpl << '\n';
                        cerr << (int) snp_list[snp_idx][i - 9][0] << ","
                            << snp_list[snp_idx][i - 9][1] << '\n';
                        if(snpl[0] != snpl[2]){
                            cerr << "HETT\n";
                        }
                    }
                    
                }
                

            }
            
            lib_list.push_back(contents);
            
            snp_idx++;
            
            if(num_snp % 10000 == 0){
                cerr << "Reading SNP: " << num_snp << " : " << ll[2] << '\n';
            }
        }
        
        vcf_file.close();
        
        num_snp = snp_idx;
    }
    
    cerr << "Total SNPs read: " << num_snp << '\n';
    
    vector<double> het_ratio(num_ind,0.0);
    vector<double> mis_ratio(num_ind,0.0);
    
    {// determine males and females
        for(int i = 0;i < num_snp; i++) {
            for (int j = 0; j < num_ind; j++) {
                if (snp_list[i][j][0] >= 0) {
                    if (snp_list[i][j][0] != snp_list[i][j][1]) {
                        het_ratio[j]++;
                    }
                } else {
                    mis_ratio[j]++;
                }
            }
        }
        
        for(int j = 0;j<num_ind;j++){
            if((num_snp - mis_ratio[j]) == 0){
                het_ratio[j] = 1;
            }else{
                het_ratio[j] = het_ratio[j]/(num_snp-mis_ratio[j]);
            }
            mis_ratio[j] = mis_ratio[j]/num_snp;
        }
        
        // output the ratio
        string temp = vcf_filename + "_ind_stats.txt";
        ofstream stats_file(temp.c_str(),ios::out);
        for(int j = 0;j<num_ind;j++){
            stats_file << het_ratio[j] << '\t' << mis_ratio[j] << '\n';
        }
        stats_file.close();
                
    }
    
    vector<char> bad_ind(num_ind,0);

    {
        for(int j = 0;j<num_ind;j++){
            if(het_ratio[j]>0.01){
                bad_ind[j] = 1;
                cerr << "Female Individual: " << j << '\n';
            }
            if(mis_ratio[j] > 0.05){
                bad_ind[j] = 1;
                cerr << "Bad Missing Individual: " << j << "\n";
            }
        }    
    }
    
    vector<double> snp_het_ratio(num_snp,0.0);
    vector<double> snp_mis_ratio(num_snp,0.0);
    int num_good_ind = 0;
    
    {
        for(int i = 0;i < num_snp; i++) {
            for (int j = 0; j < num_ind; j++) {
                
                if(bad_ind[j] > 0){
                    continue;
                }
                
                if (snp_list[i][j][0] >= 0) {
                    if (snp_list[i][j][0] != snp_list[i][j][1]) {
                        snp_het_ratio[i]++;
                    }
                } else {
                    snp_mis_ratio[i]++;
                }
            }
        }


        for (int j = 0; j < num_ind; j++) {
            if (bad_ind[j] == 0) {
                num_good_ind++;
            }
        }

        cerr << "Number of good individuals: " << num_good_ind << "/" << num_ind << '\n';


        for (int i = 0; i < num_snp; i++) {
            if (num_good_ind - snp_mis_ratio[i] > 0) {
                snp_het_ratio[i] = snp_het_ratio[i]/(num_good_ind - snp_mis_ratio[i]);
            } else {
                snp_het_ratio[i] = 1;
            }
            
            snp_mis_ratio[i] = snp_mis_ratio[i] / num_good_ind;
        }
        
        // output the ratio
        string temp = vcf_filename + "_snp_stats.txt";
        ofstream stats_file(temp.c_str(),ios::out);
        for(int i = 0; i < num_snp; i++){
            stats_file << snp_het_ratio[i] << '\t' << snp_mis_ratio[i] << '\n';
        }
        stats_file.close();
        
    }
    
    vector<char> bad_snp(num_snp,0);
    int num_bad_snp = 0;
    {
        for(int i = 0;i<num_snp;i++){
            if(snp_het_ratio[i] > 0.05){
                bad_snp[i] = 1;
            }
            if(snp_mis_ratio[i] > 0.05){
                bad_snp[i] = 1;
            }
            
            if(bad_snp[i] > 0){
                num_bad_snp++;
            }
        }    
        
        cerr << "Number of good SNPs: " << (num_snp-num_bad_snp) << "/" << num_snp << '\n';
    }
    
    { // write hap file
        string temp = vcf_filename + ".hap";
        ofstream out_file(temp.c_str(),ios::out);
        for(int i = 0;i<num_snp;i++){
            if(bad_snp[i] > 0){
                continue;
            }
            
            for(int j = 0; j < num_ind; j++){
                if(bad_ind[j] > 0){
                    continue;
                }
                
                if(j!=0){
                    out_file << ' ';
                }
                if (snp_list[i][j][0] >= 0) {
                    if (not snp_list[i][j][0] == snp_list[i][j][1]) {
                        out_file << "- -";
                    } else {
                        out_file << snp_list[i][j][0] << ' ' << snp_list[i][j][0];
                    }
                } else {
                    out_file << "- -";
                }
            }
            out_file << '\n';
        }
        out_file.close();
    }
    
    { // write SAMPLE file
        string temp = vcf_filename + ".sample";
        ofstream out_file(temp.c_str(),ios::out);
        for(int j = 0; j < num_ind; j++){
            if(bad_ind[j] > 0){
                continue;
            }
            out_file << ind_list[j] << " CEU EUR 1\n";
        }
        out_file.close();
    }

    { // write LEGNED file
        string temp = vcf_filename + ".legend";
        ofstream out_file(temp.c_str(), ios::out);
        for(int i = 0;i<num_snp;i++){
            if(bad_snp[i] > 0){
                continue;
            }
            out_file << lib_list[i][2] << ' ' << lib_list[i][1] << ' ' 
                    << lib_list[i][3] << ' ' << lib_list[i][4] << '\n';
        }
        out_file.close();
    }
    
    
    return 0;

}



void print_usage_and_exit() {
    cerr << "Usage: " + PROG_NAME + " <option>" << "\n";
    cerr << "Options:" << "\n";
    cerr << "-== Useful Tools ==-" << '\n';
    cerr << "  raw2fasta      -- Convert RAW reads to FASTA format" << "\n";
    cerr << "  subseq         -- Extract a region from a single chromosome FASTA file" << "\n";
    cerr << "  qual_convert      -- Convert the quality of FASTQ files" << "\n";
    cerr << "  fastq_trim        -- Trim BWA style [not implemented yet]" << "\n";
    cerr << "  fastq_pico_trim   -- Trim the picoplex adaptor" << "\n";
    cerr << "  fastq_filter      -- Only output reads with a percentage of bases above a certain quality" << "\n";
    cerr << "  fq_len_filter     -- Reject reads less than some length" << "\n";
    cerr << "  fastq_stats       -- Get some statistics from a FASTQ file" << "\n";
    cerr << "  fastq_quals       -- Get some qualities from a FASTQ file" << "\n";
    cerr << "  sam_quals         -- Get some qualities from a SAM file" << "\n";
    cerr << "  replace_sam_quals -- Replace qualities in a SAM file" << "\n";
    cerr << "  quant_sam_quals   -- Quantize some qualities from a SAM file" << "\n";
    cerr << "  depth_stats       -- Get some statistics from a samtools depth file" << "\n";
    cerr << "  sam_stats         -- Get some statistics from a SAM file" << "\n";
    cerr << "  barcode_split     -- Do barcode splitting on FASTA/Q files" << "\n";
    cerr << "  pair_split        -- Do barcode splitting on paired FASTA/Q files" << "\n";
    cerr << "  merge_sam         -- Concatenate two SAM files, not BAM" << "\n";
    cerr << "  sample_reads      -- Randomly select some reads from a FASTQ file" << "\n";
    exit(2);
}

