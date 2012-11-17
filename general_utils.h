/*
 *  general_utils.h
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

#ifndef GENERAL_UTILS_H

#define GENERAL_UTILS_H

#include "stl.h"
#include "fast.h"



// SAM flags
namespace sam {

    const uint16_t MULTI_FRAG = 0x1; // i.e. paired end in sequencing
    const uint16_t PROPER_ALN = 0x2;
    const uint16_t UNMAPPED = 0x4;
    const uint16_t NEXT_UMAPPED = 0x8; // same as pair unmapped

    const uint16_t REVERSE_COMP = 0x10;
    const uint16_t NEXT_REVERSE_COMP = 0x20; // same as pair reversed
    const uint16_t FIRST_FRAG = 0x40; // first pair
    const uint16_t LAST_FRAG = 0x80; // second  pair

    const uint16_t SECONDARY_ALN = 0x100; // unused most of the time
    const uint16_t FAIL_QUAL = 0x200;
    const uint16_t DUPLICATE = 0x400;

}

void trim2(string& str);
void rtrim(string& str);
void ltrim(string& str);
vector<string> split(const string &s, char delim);
vector<string> split(const string &s);

template <typename T> string toStr(T tmp) {
    ostringstream out;
    out << tmp;
    return out.str();
}

template <typename T> T strTo(string tmp) {
    T output;
    istringstream in(tmp);
    in >> output;
    return output;
}

struct upper {

    int operator()(int c) {
        return std::toupper((unsigned char) c);
    }

};


// http://www.johndcook.com/cpp_phi.html

inline double phi(double x) {
    // constants
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x) / sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5 * (1.0 + sign * y);
}


// from knuth

inline uint64_t combination(uint32_t n, uint32_t k) {
    if (k > n) {
        return 0;
    }

    uint64_t r = 1;

    for (uint32_t d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }

    return r;
}

inline uint64_t power(unsigned int n, unsigned int m) {
    uint64_t output = 1;

    if (n == 0 || n == 1) {
        return 1u;
    }

    for (; m > 0; m--) {
        output *= n;
    }

    return output;
}



// http://www.johndcook.com/blog/2008/04/24/how-to-calculate-binomial-probabilities/

inline double binom_p(double p, double q, int m, int n) {
    double temp = lgamma(m + n + 1.0);
    temp -= lgamma(n + 1.0) + lgamma(m + 1.0);
    temp += m * log(p) + n * log(q);
    return exp(temp);
}

inline double binom_prob(uint32_t x, uint32_t n, double p) {
    return binom_p(p, 1 - p, x, n - x);
}


template <typename T> bool get_bit(T a, int i) {
    return a & (((T) 1) << i);
}

inline void set_bit(int &a, int i) {
    a |= 1u << i;
}

inline void print_bits(ostream &o, const vector<int> &a) {

    for (vector<int>::const_reverse_iterator rit = a.rbegin(); rit != a.rend(); rit++) {

        for (int i = ((int) (sizeof (int) *8) - 1); i >= 0; i--) {
            o << get_bit<int>(*rit, i);
        }
    }
}

template <typename T> void print_bits(ostream &o, T a) {

    for (int i = ((int) (sizeof (T)*8) - 1); i >= 0; i--) {
        o << get_bit<T > (a, i);
    }

}




// MT code taken from
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
// Ported to C++ by John Mu

class MT_random {
private:
    static const uint64_t NN = 312;
    static const uint64_t MM = 156;
    static const uint64_t MATRIX_A = 0xB5026F5AA96619E9ULL;
    static const uint64_t UM = 0xFFFFFFFF80000000ULL; /* Most significant 33 bits */
    static const uint64_t LM = 0x7FFFFFFFULL; /* Least significant 31 bits */




    /* The array for the state vector */
    uint64_t mt[NN];
    /* mti==NN+1 means mt[NN] is not initialized */
    uint64_t mti;

    void init_all() {
        mti = NN + 1;
    }

public:

    MT_random() {

        init_all();
        //init_genrand64(time(NULL));
        init_genrand64(0);
    }

    MT_random(unsigned long long seed) {

        init_all();
        init_genrand64(seed);
    }

    /* initializes mt[NN] with a seed */
    void init_genrand64(unsigned long long seed) {
        mt[0] = seed;
        for (mti = 1; mti < NN; mti++)
            mt[mti] = (6364136223846793005ULL * (mt[mti - 1] ^ (mt[mti - 1] >> 62)) + mti);
    }

    /* initialize by an array with array-length */
    /* init_key is the array for initializing keys */

    /* key_length is its length */
    void init_by_array64(unsigned long long init_key[],
            unsigned long long key_length) {
        unsigned long long i, j, k;
        init_genrand64(19650218ULL);
        i = 1;
        j = 0;
        k = (NN > key_length ? NN : key_length);
        for (; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 3935559000370003845ULL))
                    + init_key[j] + j; /* non linear */
            i++;
            j++;
            if (i >= NN) {
                mt[0] = mt[NN - 1];
                i = 1;
            }
            if (j >= key_length) j = 0;
        }
        for (k = NN - 1; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 2862933555777941757ULL))
                    - i; /* non linear */
            i++;
            if (i >= NN) {
                mt[0] = mt[NN - 1];
                i = 1;
            }
        }

        mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
    }

    /* generates a random number on [0, 2^64-1]-interval */
    unsigned long long genrand64_int64(void) {
        uint64_t i;
        uint64_t x;
        static uint64_t mag01[2] = {0ULL, MATRIX_A};

        if (mti >= NN) { /* generate NN words at one time */

            /* if init_genrand64() has not been called, */
            /* a default initial seed is used     */
            if (mti == NN + 1)
                init_genrand64(5489ULL);

            for (i = 0; i < NN - MM; i++) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] = mt[i + MM] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            for (; i < NN - 1; i++) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] = mt[i + (MM - NN)] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            x = (mt[NN - 1] & UM) | (mt[0] & LM);
            mt[NN - 1] = mt[MM - 1] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];

            mti = 0;
        }

        x = mt[mti++];

        x ^= (x >> 29) & 0x5555555555555555ULL;
        x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
        x ^= (x << 37) & 0xFFF7EEE000000000ULL;
        x ^= (x >> 43);

        return x;
    }

    /* generates a random number on [0, 2^63-1]-interval */
    long long genrand64_int63(void) {
        return (long long) (genrand64_int64() >> 1);
    }

    /* generates a random number on [0,1]-real-interval */
    double genrand64_real1(void) {
        return (genrand64_int64() >> 11) * (1.0 / 9007199254740991.0);
    }

    /* generates a random number on [0,1)-real-interval */
    double genrand64_real2(void) {
        return (genrand64_int64() >> 11) * (1.0 / 9007199254740992.0);
    }

    /* generates a random number on (0,1)-real-interval */
    double genrand64_real3(void) {
        return ((genrand64_int64() >> 12) + 0.5) * (1.0 / 4503599627370496.0);
    }

    // this is not the "correct" way... but good enough

    double genrand_norm(double mean, double sd) {
        return (genrand_norm() * sd) +mean;
    }

    double genrand_norm() {
        double sum = 0;

        for (int i = 0; i < 48; i++) {
            sum = sum + (double) (rand() + 1.0) / (RAND_MAX + 1.0);
        }

        // only for 48
        return (sum - 24) / 2;
    }

    double genrand_exp(double mean) {
        return -mean * log(1 - genrand64_real2());
    }

    int64_t genrand_binom(uint64_t n, double p) {

        if (p < 0 || p > 1) {
            return -1;
        }

        int64_t total = 0;
        for (uint64_t i = 0; i < n; i++) {
            total += (genrand64_real1() < p ? 1 : 0);
        }

        return total;
    }

    // generate random int inclusive of range [start,end]

    uint64_t genrand_int_range(uint64_t start, uint64_t end) {

        if(start == end){
            return start;
        }

        return (uint64_t) ((genrand64_real2()*(end - start + 1)) + start);
    }

    // bernouli p is true

    bool genrand_bern(double p) {
        return genrand64_real2() <= p;
    }



};



// This class is specific to linux :(

class mu_timer {
private:
    struct timeval start_tv;
    bool enabled;
public:

    mu_timer() {
        gettimeofday(&start_tv, NULL);
        enabled = true;
    }

    void set_enabled(bool e){
        enabled = e;
    }

    void reset() {
        if (enabled) {
            gettimeofday(&start_tv, NULL);
        }
    }

    double elapsed_time() {
        if (enabled) {
            struct timeval tv;
            gettimeofday(&tv, NULL);

            return (double) (tv.tv_sec - start_tv.tv_sec) + (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
        } else {
            return 0;
        }
    }

    void print_elapsed_time(ostream &o, string s) {
        if (enabled) {
            o << s << " time: " << elapsed_time() << "s." << '\n';
        }
    }
};

class param_reader {
private:
    map<string, string> params_map;
    string single_char; // string of params which don't have an argument
    string single_char_exist; // string of all params
    set<string> multi_set; // mutichar ones with out argument
    set<string> multi_set_exist; // mutichar ones

public:

    static const char BLANK = '\n';

    param_reader(string single_char, vector<string> multi_char, string single_char_exist, vector<string> multi_char_exist) {
        this->single_char = single_char;
        this->single_char_exist = single_char_exist;

        for (int i = 0; i < (int) multi_char.size(); i++) {
            multi_set.insert(multi_char[i]);
        }
        for (int i = 0; i < (int) multi_char_exist.size(); i++) {
            multi_set_exist.insert(multi_char_exist[i]);
        }
    }


    // adds BLANK for blank parametrs

    int add_params(const vector<string> &params) {
        vector<string>::const_iterator it = params.begin();
        int count = 1;

        while (it != params.end()) {
            if (it->length() > 1 && (*it)[0] == '-') {


                if ((*it)[1] != '-') { // single char mode

                    // check for single char args
                    for (int i = 1; i < (int) it->length() - 1; i++) {
                        if (single_char_exist.find((*it)[i]) == string::npos) {
                            cerr << "Option doesn't exist: " << (*it)[i] << '\n';
                            return count;
                        }

                        if (single_char.find((*it)[i]) != string::npos) {
                            pair<string, string> contents;
                            contents.first = (*it)[i];
                            contents.second = BLANK;

                            params_map.insert(contents);

                        } else {
                            cerr << "Option needs Argument: " << (*it)[i] << '\n';
                            return count;
                        }
                    }

                    pair<string, string> contents;
                    contents.first = (*it)[it->length() - 1];

                    // the last one can be an argumented parameter
                    if (single_char_exist.find((*it)[it->length() - 1]) == string::npos) {
                        cerr << "Option doesn't exist: " << contents.first << '\n';
                        return count;
                    }


                    if (single_char.find((*it)[it->length() - 1]) != string::npos) {

                        contents.second = BLANK;

                        params_map.insert(contents);

                    } else {
                        it++;

                        if (it != params.end()) {
                            contents.second = *it;

                            params_map.insert(contents);
                        } else {
                            cerr << "Option needs Argument: " << contents.first << '\n';
                            return count;
                        }

                        count++;
                    }


                } else { // multi char mode
                    pair<string, string> contents;
                    contents.first = it->substr(2);

                    if (multi_set_exist.find(contents.first) == multi_set_exist.end()) {
                        cerr << "Option doesn't exist: " << contents.first << '\n';
                        return count;
                    }

                    if (multi_set.find(contents.first) != multi_set.end()) {
                        contents.second = BLANK;

                        params_map.insert(contents);

                    } else {

                        it++;

                        if (it != params.end()) {
                            contents.second = *it;

                            params_map.insert(contents);
                        } else {
                            cerr << "Option needs Argument: " << contents.first << '\n';
                            return count;
                        }

                        count++;

                    }

                }


            } else {
                return count;
            }

            it++;
            count++;
        }



        return 0;
    }

    string get_val(string key) {
        string output = "";

        map<string, string>::iterator it = params_map.find(key);

        if (it != params_map.end()) {
            output = it->second;
        }

        return output;
    }

    void clear() {
        params_map.clear();
    }

};




///////////
// Inlines

inline char dna_upper(char c) {
    return (char) c & (0xDFu);
}

inline int min(int a, int b) {
    return ((a) < (b) ? (a) : (b));
}

inline int max(int a, int b) {
    return ((a) > (b) ? (a) : (b));
}

inline int max3(int a, int b, int c) {
    return max((a), max(b, c));
}


inline string reverse_string(const string &s) {
    string output;
    string::const_reverse_iterator rit = s.rbegin();
    while (rit != s.rend()) {
        output += *rit;
        rit++;
    }

    return output;
}


inline uint64_t factorial(unsigned int n) {
    uint64_t output = 1;
    if (n == 0 || n == 1) {
        // do nothing
    } else {
        for (; n > 1; n--) {
            output = output * n;
        }
    }
    return output;
}


// k!*stirling2(n,k)

inline int balls(unsigned int n, unsigned int k) {
    int output = 0;

    for (unsigned int j = 0; j <= k; j++) {
        if (j % 2 == 0) {
            output += (int) combination(k, j)*(int) power((k - j), n);
        } else {
            output -= (int) combination(k, j)*(int) power((k - j), n);
        }
    }


    return output;
}



inline void compleseq(string &seq) {
    string temp = seq;
    string::iterator it = seq.begin();
    string::reverse_iterator rit = temp.rbegin();
    while (it != seq.end()) {
        switch (*rit) {
            case 'A':
                *it = 'T';
                break;
            case 'C':
                *it = 'G';
                break;
            case 'G':
                *it = 'C';
                break;
            case 'T':
                *it = 'A';
                break;
            case 'a':
                *it = 't';
                break;
            case 'c':
                *it = 'g';
                break;
            case 'g':
                *it = 'c';
                break;
            case 't':
                *it = 'a';
                break;
            default:
                *it = *rit;
                break;
        }
        it++;
        rit++;
    }
}




vector<string> list_dir(string path);
bool mkdir(string dir_name);
bool program_exists(string name);



class seen_t {
private:
    uint32_t* locs; // hash table
    vector<uint32_t> locs_added;

public:
    static const uint32_t hash_size = 1048576;
    static const uint32_t shift_bits = 32 - 20;
    static const uint32_t magic = 2654435761u;

    //static const uint32_t mask = 0xFFFFF;

    seen_t() {
        locs = new uint32_t[hash_size];
        memset(locs, hash_size, sizeof (uint32_t));
        locs_added.reserve(5000);
    }

    ~seen_t() {
        clear();
        delete [] locs;
    }

    // true if inserted

    bool insert(uint32_t a) {
        uint32_t hash_val = (a * magic) >> shift_bits;


        if (locs[hash_val] == 0) {
            locs[hash_val] = a;
            locs_added.push_back(hash_val);
            return true;
        } else if (locs[hash_val] == a) {
            return false;
        } else {
            hash_val++;

            uint32_t num_seen = 1;
            while (num_seen < hash_size) {
                uint32_t hash_val_temp = hash_val;
                if (hash_val >= hash_size) {
                    hash_val_temp = hash_val_temp - hash_size;
                }


                if (locs[hash_val_temp] == 0) {
                    locs[hash_val_temp] = a;
                    locs_added.push_back(hash_val_temp);
                    return true;
                } else if (locs[hash_val_temp] == a) {
                    return false;
                }

                hash_val++;
                num_seen++;

            }

            return true; // hash is full :( ... should not happen
        }
    }

    // true if found

    bool find(uint32_t a) {
        uint32_t hash_val = (a * magic) >> shift_bits;


        if (locs[hash_val] == 0) {
            return false;
        } else if (locs[hash_val] == a) {
            return true;
        } else {
            hash_val++;

            uint32_t num_seen = 1;
            while (num_seen < hash_size) {
                uint32_t hash_val_temp = hash_val;
                if (hash_val >= hash_size) {
                    hash_val_temp = hash_val_temp - hash_size;
                }


                if (locs[hash_val_temp] == 0) {
                    return false;
                } else if (locs[hash_val_temp] == a) {
                    return true;
                }

                hash_val++;
                num_seen++;

            }

            return true; // hash is full :( ... should not happen
        }
    }

    void clear() {
        vector<uint32_t>::iterator added_it = locs_added.begin();
        while (added_it != locs_added.end()) {

            locs[*added_it] = 0;

            added_it++;
        }

        locs_added.clear();
    }

};



////////

// convert illumina 1.3+ format to Sanger format

inline int ill_13tosanger(string &q) {
    string::iterator sit = q.begin();

    while (sit != q.end()) {
        int new_val = (int) (*sit) - 31;
        if (new_val > 126) {
            return 1; // error
        }

        *sit = (char) new_val;

        sit++;
    }

    return 0;
}

// convert Sanger format to illumina 1.3+ format

inline int sangertoill_13(string &q) {
    string::iterator sit = q.begin();

    while (sit != q.end()) {
        int new_val = (int) (*sit) + 31;
        if (new_val < 33) {
            return 1; // error
        }

        *sit = (char) new_val;

        sit++;
    }

    return 0;
}


// convert Solexa format to Sanger format

inline int solexatosanger(string &q) {
    string::iterator sit = q.begin();

    while (sit != q.end()) {
        int new_val = (int) (*sit);
        if (new_val < 0 || new_val > 67) {
            return 1; // error
        }

        *sit = (char) qual::solexa2phred32_table[new_val];

        sit++;
    }

    return 0;
}

// qual is phred
// code from BWA
inline int bwa_trim_len(string &qual_str, int qual, int min_len) {
    int s = 0;
    int l = 0;
    int max = 0;
    int line_len = (int)qual_str.length();
    int max_l = line_len - 1;

    if (qual < 1 || line_len == 0) {
        return 0;
    }

    for (l = line_len - 1; l >= min_len - 1; --l) {
        s += qual - (qual_str[l] - 33);
        if (s < 0) break;
        if (s > max) {
            max = s;
            max_l = l;
        }
    }
    return max_l + 1;
}


/////




void output_percentage(ostream &out, uint64_t total, uint64_t num);


#endif




