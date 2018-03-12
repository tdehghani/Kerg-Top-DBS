#pragma once
#include <vector>
#include <iostream>
using namespace std;
class strand
{
private:
	int strand_id;
	int strand_begin;
	int strand_end;
	vector<char> strand_seq;
public:
	strand();
	void set_strand(int, int, int, vector <char>);
	~strand();
	/*int get_strand_id();
	int get_strand_begin();
	int get_strand_end();
	vector<char> & get_strand_seq();
	void strand::set_strand_id(int id);
	void strand::set_strand_begin(int begin);
	void strand::set_strand_end(int end);
	void strand::set_strand_seq(vector<char> seq);*/
};

