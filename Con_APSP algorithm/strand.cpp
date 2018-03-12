/*
*   Top-DBS algorithm for Beta-sheet Topology Prediction Problem
*   Copyright (C) 2018 ،Toktam Dehghani
* 	Knowledge Engineering Research Group (KERG), Department of  Computer Engineering,
*   Faculty of Engineering, Ferdowsi University of Mashhad, Iran. https://kerg.um.ac.ir/
*   Toktam Dehghani, email: dehghani.toktam@mail.um.ac.ir
*
* main.c - This file is part of Top-DBS
*
* Top-DBS is a free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundation; either version 3 of the License,
* or (at your option) any later version.
*
* Top-DBS is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
*/

#include "strand.h"


strand::strand()
{}
strand::~strand()
{
	strand_seq.clear();

}
void strand::set_strand(int id, int begin, int end, vector<char> seq)
{
	strand_id = id;
	strand_begin = begin;
	strand_end = end;
	strand_seq = seq;
}
//int strand::get_strand_id()
//{
//	return strand_id;
//}
//int strand::get_strand_begin()
//{
//	return strand_begin;
//}
//int strand::get_strand_end()
//{
//	return strand_end;
//}
//vector<char> & strand::get_strand_seq()
//{
//	return strand_seq;
//}
//
//void strand::set_strand_id(int id)
//{
//	strand_id = id;
//}
//void strand::set_strand_begin(int begin)
//{
//	strand_begin = begin;
//}
//void strand::set_strand_end(int end)
//{
//	strand_end = end;
//}
//void strand::set_strand_seq(vector<char> seq)
//{
//	strand_seq = seq;
//}