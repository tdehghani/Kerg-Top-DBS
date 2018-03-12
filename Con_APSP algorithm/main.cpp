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
#define _CRT_SECURE_NO_DEPRECATE
#pragma comment(linker, "/STACK:8000000")
#pragma comment(linker, "/HEAP:8000000")
#include<stdio.h>
#include<time.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include<conio.h>
#include <math.h>
#include <Windows.h>
#include <direct.h>
#include <sstream>
#include "strand.h"
# include "find_max_beta_topology.h"

 
#define max_num_strand = 30;// maximume number of strand in each protein

int read_int_from_file(char * str, int n, FILE * ifp)
{
	if (fgets(str, n, ifp) == NULL) {
		str[0] = 'e';
	}
	int length = (int)strlen(str);
	if ((str[0] != 'e') && (str[length - 1] == '\n')) {
		str[length - 1] = '\0';
	}
	int temp1 = 0;
	temp1 = atoi(str);
	return temp1;
}
double string_to_double(const std::string& s)
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}
double read_score_from_file(char * str, int n, FILE * ifp)
{
	if (fgets(str, n, ifp) == NULL) {
		str[0] = 'e';
	}
	int length = (int)strlen(str);
	if ((str[0] != 'e') && (str[length - 1] == '\n')) {
		str[length - 1] = '\0';
	}
	double out = string_to_double(str);
	return out;
}
string read_string_from_file(char * str, int n, FILE * ifp)
{
	if (fgets(str, n, ifp) == NULL) {
		str[0] = 'e';
	}
	int length = (int)strlen(str);
	if ((str[0] != 'e') && (str[length - 1] == '\n')) {
		str[length - 1] = '\0';
	}
	return(str);
}
void discard_filename_extension(char * filename_input) {
	int filename_length;
	int i;
	filename_length = (int)strlen(filename_input);
	for (i = 0; i < filename_length; i++) {
		if (filename_input[i] == '.')
			break;
	}
	filename_input[i] = '\0';
}
void intial_clique(clique_info &avg_clique)
{
	avg_clique.run_time = 0;
};

void main()
{
	int max_length = 3000;//for length of amino asid sequence	
	/********************************************************************************************************************************/
	/********************************************************************************************************************************/
	printf("********************************************************************************************************\n");
	printf("********************************************************************************************************\n");
	printf("********************************************************************************************************\n");
	printf(" 	Welcome to Top-DBS method: Protein Beta-Sheet Topology Prediction              \n");
	printf(" 	Knowledge Engineering Research Group(KERG), Department of  Computer Engineering \n");
	printf(" 	Faculty of Engineering, Ferdowsi University of Mashhad, Iran.https://kerg.um.ac.ir/\n");
	printf(" 	Toktam Dehghani, email : dehghani.toktam@mail.um.ac.ir\n");
	printf("********************************************************************************************************\n");
	printf("	Please input protein PDB ID and Chain (as example: 1nega): ");
	char * pdb_id = new char[200];
	scanf("%5c",pdb_id);
	pdb_id[(int)strlen (pdb_id) - 1] = '\0';
	//strcpy_s(pdb_id, max_length, "1goua");
	/********************************************************************************************************************************/
	/********************************************************************************************************************************/
	int number_of_beta_strands;//number_of_strands
	FILE * input_file;// input file= strands pairwise alignments scores
	//FILE* ifp_list;
	//FILE* ofp_list;
	//FILE* run_time;
	char * filename_input = new char[200]; // input filename address
	char* fold_number_str = new char[3];
	char * fold_path_str = new char[1000];
	char * result_fold_path_str = new char[1000];// address of out put file
	char * result_file_path_str = new char[1000];// address of out put file
	char * result_file2_path_str = new char[1000];// address of list of out put file
	//char * path_floyed_matrix = new char[1000];
	//char * fold_list_str = new char[1000];
	char * clique_fold_path_str = new char[1000];
	char * protein_information_path_str = new char[1000];
	//char * run_time_file_path= new char[1000];
	
	char * s = new char[max_length];
	vector<strand > strand_list;// each strand info
	vector<vector<double> > max_alignment_scores;
	vector<vector<double> > max_alignment_prob;
	float runtime_floyed = 0, num_running_folyed = 0;
	double  num_running_clique = 0;
	clique_info avg_clique;
	intial_clique(avg_clique);	
	
	// read input data :strand pairwise alignmnets
	strcpy_s(fold_path_str, max_length, "example\\");
	strcat_s(fold_path_str, max_length, pdb_id);
	strcat_s(fold_path_str, max_length, "\\");
	strcat_s(fold_path_str, max_length, pdb_id);
	strcat_s(fold_path_str, max_length, ".txt");
	input_file = fopen(fold_path_str, "r");
	if (input_file == NULL)
		{
		printf("******************************************************  End  ******************************************\n");

		printf("	Error in the openning of input files \n");
		printf("	For more information about file formats, please, refer to README.txt \n");
		}
	else
	{
		printf("	Steps:\n");
		printf("	1.	Read beta-strand pairwise alignment matrix. \n");
		max_alignment_scores.clear();
		max_alignment_prob.clear();
		read_string_from_file(s, max_length, input_file);
		number_of_beta_strands = 0;
		for (int i = 0; i < (int)strlen(s); i++)
		{
			if (s[i] == ' ')
				number_of_beta_strands++;
			if (s[i] == '\n')
				break;
		}
		number_of_beta_strands = number_of_beta_strands / (2);

		//read input file
		max_alignment_scores.resize((number_of_beta_strands), vector<double>(number_of_beta_strands * 2));
		double max = 0;
		for (int i = 0; i < number_of_beta_strands; i++)
		{
			int j = 0;
			double x = 0;
			string str;
			for (int z = 0; z < (int)strlen(s); z++)
			{
				if (s[z] == ' ')
				{
					x = string_to_double(str);
					max_alignment_scores[i][j] = x;
					if (x > max)
						max = x;
					str.clear();
					j++;
				}
				else
					str.push_back(s[z]);
			}
			read_string_from_file(s, max_length, input_file);
		}
		max_alignment_prob.resize((number_of_beta_strands), vector<double>(number_of_beta_strands * 2));
		for (int i = 0; i < (number_of_beta_strands); i++)
		{
			for (int j = 0; j < (2 * number_of_beta_strands); j++)
			{
				if ((i < number_of_beta_strands) && (j < number_of_beta_strands))
					max_alignment_prob[i][j] = max_alignment_scores[i][j];
				if ((i < number_of_beta_strands) && (j >= number_of_beta_strands))
					max_alignment_prob[i][j] = max_alignment_scores[i][j];
				if (i == j || (i == (j - number_of_beta_strands)) || (j == (i - number_of_beta_strands)))
					max_alignment_prob[i][j] = 0;
			}
		}
		//out put files
		strcpy_s(clique_fold_path_str, max_length, "example\\");
		strcat_s(clique_fold_path_str, max_length, pdb_id);
		strcat_s(clique_fold_path_str, max_length, "\\");
		strcat_s(clique_fold_path_str, max_length, pdb_id);
		strcat_s(clique_fold_path_str, max_length, ".out");
		clique_info my_clique;
		//Input matrices
		vector<vector<vector<double>>> SCORE_matrix;//all-pirs shortest path score matrix output of Con-APSP algorithm
		vector<vector<vector<vector<int>>>> PATH_matrix;//all-pirs shortest paths matrix output of Con-APSP algorithm
		//Con-APSP algorithm
		printf("	2.	Generate a set of candidate beta-sheets (Con-APSP). \n");
		Con_APSP(SCORE_matrix, PATH_matrix, max_alignment_prob, 2 * number_of_beta_strands, strand_list, result_file_path_str);//find shortest paths in the beta strand graph
		//Top-DBS algorithm
		printf("	3.	Extract a subset of maximum weight disjoint beta-sheets (Top-DBS). \n");
	    my_clique = clique_algorithm(SCORE_matrix, PATH_matrix, 2 * number_of_beta_strands, number_of_beta_strands, clique_fold_path_str);
		//Run time
		avg_clique.run_time = my_clique.run_time;
        fclose(input_file);
	
	char *path = NULL;
	path = _getcwd(NULL, 0);
	printf("******************************************************  End  ******************************************\n");

	printf("	Top-DBS algorithms predicted successfully the beta-sheet topology\n");
	printf(" 	of protein: %s",pdb_id);	printf(" with %d beta-strands\n", number_of_beta_strands);
	printf(" 	in %f seconds.\n", avg_clique.run_time);
	printf(" 	The results can be found in the current directory:\n");
	printf(" 	%\\example\\%s\\%s.out\n", pdb_id,pdb_id);
	printf("********************************************************************************************************\n");
	printf("********************************************************************************************************\n");
}
	_getch();
}//main