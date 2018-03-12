#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include<conio.h>
#include <math.h>
#include <sstream>
#include "strand.h"
#include <time.h>
#include <algorithm>    
/* functions */
using namespace std;
struct paths
{
	vector<int> path;
	double score;
};
	
struct clique_info
{
	/*int number_edge;
	int number_node;
	double density;
	int tree_width;*/
	float run_time;
	/*int prune_number_edge;
	int prune_number_node;	
	int prune_tree_width;*/
};

void normlized();
int no_similar_strand(int i, int k, int w, int j, int step1, int step2);
//void print_sheet(int i, int j, int step);
vector<int> merge_paths(int i, int j, int k, int step1, int step2);
int no_Complement_strands(int i, int k, int w, int j, int step1, int step2);
void func_write_to_file_all_floyedwarshall_sheets(char * result_fold_path_str, int step, vector<strand > strand_list);
double Con_APSP(vector<vector<vector<double>>> &SCORE_matrix, vector<vector<vector<vector<int>>>> &PATH_matrix, vector<vector<double> > weights_matrix, int num, vector<strand > strand_list, char *result_fold_path_str);// , char *path_floyed_matrix);
	void func_P_or_AP(int x, int y, int num_row);
void func_print_sheet(int max_row, int max_col, int step);
void func_write_to_file(FILE* ofp_conformation, vector<strand > strand_list);
int lookfor(int index1, int index2, vector<int> ch, int k);
//clique
//void func_print_sheet_clique(vector<int> PATH, int max_row, int max_col, int step, int num_row, char *clique_fold_path_str);
void func_P_or_AP_clique(int x, int y, int num_row, vector<vector<int>> &pair_strands, vector<char> &interaction_type);
clique_info clique_algorithm(vector<vector<vector<double>>> SCORE_matrix, vector<vector<vector<vector<int>>>> PATH_matrix, int num_row, int num_of_sheet, char *clique_fold_path_str);