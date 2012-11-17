/*
 *  general_utils.cpp
 *  seqtools
 *
 *  Created by John Mu on 3/5/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 */

#include "general_utils.h"

bool program_exists(string name)
{
	int out = 1;

	string cmd = "which " + name + " > /dev/null 2> /dev/null";
	out = system(cmd.c_str());


	return (out == 0);
}

bool mkdir(string dir_name)
{
	int out = 1;

	string cmd = "mkdir " + dir_name + " 2> /dev/null";
	out = system(cmd.c_str());


	return (out == 0);
}

vector<string> list_dir(string path)
{
	vector<string> result;
	string cmd;
	string temp;

	//int out;
	FILE *fp;
	char line[1024];

	cmd = "ls " + path;

	fp = popen(cmd.c_str(), "r");

	if (fp == NULL) {
		cout << "ERROR: failed to open pipe" << endl;

		exit(2);
	}else {
		while ( fgets( line, 1024, fp))
		{
			temp = line;
			trim2(temp);
			if (temp.length() > 0) {
				result.push_back(temp);
			}
		}
	}
	pclose(fp);



	return result;
}



// Source: http://mlawire.blogspot.com/2009/07/c-whitespace-trimming-functions.html
void ltrim(string& str)
{
	string::size_type pos = 0;
	while (pos < str.size() && (isspace(str[pos]))) pos++;
	str.erase(0, pos);
}
void rtrim(string& str)
{
	string::size_type pos = str.size();
	while (pos > 0 && (isspace(str[pos - 1]))) pos--;
	str.erase(pos);
}
void trim2(string& str)
{
	ltrim(str);
	rtrim(str);
}

vector<string> split(const string &s, char delim)
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// Split using all white space as delims
vector<string> split(const string &s)
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(ss >> item) {
        elems.push_back(item);
    }
    return elems;
}





string str_toupper(const string &s){
    string output = "";
    for (uint64_t i = 0; i<s.length(); i++) {
        output += toupper(s[i]);
    }

    return output;
}



void output_percentage(ostream &out, uint64_t total, uint64_t num)
{

    if(total == num){
        out << "(100.0%)";
    }else if(num == 0){
        out << "(0.0%)";
    }else{
        out << " (" << setprecision(6) << (num*100.0/(double)total) << resetiosflags (ios_base::fixed) << "%)";
    }
}






