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

#include "find_max_beta_topology.h"
using namespace std;
int num_row;
vector<vector<vector<double>>> A_matrix;//  APSP Scores
vector<vector<vector<int>>> P_matrix;//  APSP Paths_length
vector<vector<vector<vector<int>>>> Path_matrix;//  APSP Paths
vector<vector<int>> pair_strands;
vector<char> interaction_type;
double Infinite = 99999;
vector<vector<vector<double>>> A2_matrix;
bool sortByScore(paths &lhs, paths &rhs) { return lhs.score < rhs.score; }
int lookfor(int index1, int index2, vector<int> ch, int k)
{
	int stop = -1;
	for (int i = k; i >= 0 && (stop == -1); i--)
	{
		int key = 0;
		for (int j = 0; j < P_matrix[index1][index2][i]; j++)
		{
			int index = Path_matrix[index1][index2][i][j];
			if (index<num_row / 2)
				if (ch[index]>0)
					key = 1;
			else
				if (ch[index - (num_row / 2)]>0)
					key = 1;
		}
		if (key == 0)
			stop = i;
	}
	return(stop);
}
void quickSort(int left, int right) {
	int i = left, j = right;
	vector<double> tmp;
	double pivot = A2_matrix[((left + right) / 2) / (num_row)][((left + right) / 2) % (num_row)][2];
	while (i <= j) {
		while (A2_matrix[i / (num_row)][i % (num_row)][2] < pivot)
			i++;
		while (A2_matrix[j / (num_row)][j % (num_row)][2] > pivot)
			j--;
		if (i <= j) {
			tmp = A2_matrix[i / (num_row)][i % (num_row)];
			A2_matrix[i / (num_row)][i % (num_row)] = A2_matrix[j / (num_row)][j % (num_row)];
			A2_matrix[j / (num_row)][j % (num_row)] = tmp;
			i++;
			j--;
		}
	};
	if (left < j)
		quickSort(left, j);
	if (i < right)
		quickSort(i, right);
}
vector<int> merge_paths(int i, int j, int k, int step1, int step2)
{
	vector<int> out;
	if ((P_matrix[i][k][step1] == 0) && (P_matrix[k][j][step2] == 0))
			out.push_back(k);
	else
	{
		if (P_matrix[i][k][step1] == 0)
		{
			out.push_back(k);
			for (int z = 0; z < P_matrix[k][j][step2]; z++)
				out.push_back(Path_matrix[k][j][step2][z]);
		}
		else
		{
			if (P_matrix[k][j][step2] == 0)
			{
				out = Path_matrix[i][k][step1];
				out.push_back(k);
			}
			else
			{
				out = Path_matrix[i][k][step1];
				out.push_back(k);
				for (int z = 0; z < P_matrix[k][j][step2]; z++)
					out.push_back(Path_matrix[k][j][step2][z]);
			}
		}
	}
	return(out);
}
int no_similar_strand(int i, int k, int w, int j, int step1, int step2)// find simillar strands in two paths strands i to k compare to strands w to j
{
	int out = 0;
	vector<int> check;
	check.resize(num_row);
	if ((i == w) || (i == k) || (w == j) || (j == k))
		out = 1;
	for (int p = 0; p < num_row; p++)
		check[p] = 0;
	if (check[i] == 1)
		out = 1;
	else
		check[i] = 1;
	if (check[k] == 1)
		out = 1;
	else
		check[k] = 1;
	if (check[j] == 1)
		out = 1;
	else
		check[j] = 1;
	if (w != k)
	{
		if (check[w] == 1)
			out = 1;
		else
			check[w] = 1;
	}
	if (out == 0)
	{
		for (int z = 0; z < P_matrix[i][k][step1]; z++)
		{
			int qq = Path_matrix[i][k][step1][z];
			if (check[qq] == 1)
				out = 1;
			else
				check[qq] = 1;
		}
	}
	if (out == 0)
	{
		for (int y = 0; y < P_matrix[w][j][step2]; y++)
		{
			int qq = Path_matrix[w][j][step2][y];
			if (check[qq] == 1)
				out = 1;
			else
				check[qq] = 1;
		}
	}
	return(out);
}
int no_Complement_strands(int i, int k, int w, int j, int step1, int step2)// find complement strands in two paths
{
	int out = 0;
	vector<int> check;
	check.resize(num_row / 2);
	if ((i == w - (num_row / 2)) || (w == i - (num_row / 2)))
		out = 1;
	if ((i == k - (num_row / 2)) || (k == i - (num_row / 2)))
		out = 1;
	if ((j == w - (num_row / 2)) || (j == i - (num_row / 2)))
		out = 1;
	if ((j == k - (num_row / 2)) || (k == j - (num_row / 2)))
		out = 1;
	for (int p = 0; p < num_row / 2; p++)
		check[p] = 0;
	int qq;
	qq = i;
	if (qq >= (num_row / 2))
		qq = qq - (num_row / 2);
	if (check[qq] == 1)
		out = 1;
	else
		check[qq] = 1;
	qq = k;
	if (qq >= (num_row / 2))
		qq = qq - (num_row / 2);
	if (check[qq] == 1)
		out = 1;
	else
		check[qq] = 1;
	qq = j;
	if (qq >= (num_row / 2))
		qq = qq - (num_row / 2);
	if (check[qq] == 1)
		out = 1;
	else
		check[qq] = 1;
	if (w != k)
	{
		qq = w;
		if (qq >= (num_row / 2))
			qq = qq - (num_row / 2);
		if (check[qq] == 1)
			out = 1;
		else
			check[qq] = 1;
	}
	if (out == 0)
	{
		for (int z = 0; z < P_matrix[i][k][step1]; z++)
		{
			qq = Path_matrix[i][k][step1][z];
			if (qq >= (num_row / 2))
				qq = qq - (num_row / 2);
			if (check[qq] == 1)
				out = 1;
			else
				check[qq] = 1;
		}
	}
	if (out == 0)
	{
		for (int y = 0; y < P_matrix[w][j][step2]; y++)
		{
			qq = Path_matrix[w][j][step2][y];
			if (qq >= (num_row / 2))
				qq = qq - (num_row / 2);
			if (check[qq] == 1)
				out = 1;
			else
				check[qq] = 1;
		}
	}
	return(out);
}
void func_write_to_file_all_floyedwarshall_sheets(char * result_fold_path_str, int step, vector<strand > strand_list)// prints all paths
{
	FILE* ofp_conformation;// = fopen(result_fold_path_str, "w");	
	ofp_conformation = fopen(result_fold_path_str, "w");
	A2_matrix.clear();
	A2_matrix.resize((num_row), vector<vector<double>>((num_row), vector<double>(3)));
	int num_of_sheets = 0;
	for (int i = 0; i < (num_row); i++)
	{
		for (int j = 0; j < (num_row); j++)
		{
			A2_matrix[i][j][0] = i;
			A2_matrix[i][j][1] = j;
			A2_matrix[i][j][2] = A_matrix[i][j][step];
			if (A2_matrix[i][j][2] != (-1 * Infinite))
				num_of_sheets++;
		}
	}
	quickSort(0, (((num_row)*(num_row)) - 1));
	int count = num_of_sheets;
	double percent = 1;
	fprintf(ofp_conformation, "Num of strands: ");
	fprintf(ofp_conformation, "%d\n", (num_row / 2));
	fprintf(ofp_conformation, "Num of sheets in the pool: ");
	fprintf(ofp_conformation, "%d\n", num_of_sheets);
	for (int i = (num_row)-1; i >= 0; i--)
	{
		for (int j = (num_row)-1; j >= 0; j--)
		{
			if ((A2_matrix[i][j][2] != (-1 * Infinite)))
			{

				fprintf(ofp_conformation, "id: %d\n", count);
				//fprintf(ofp_conformation, "score: \n");
				//fprintf(ofp_conformation, "%f\n", A2_matrix[i][j][2]);
				/*fprintf(ofp_conformation, "ranked: \n");
				percent = ((double)percent / (double)num_of_sheets) * 100;
				fprintf(ofp_conformation, "%f\n", percent);*/
				interaction_type.clear();
				pair_strands.clear();
				func_print_sheet(A2_matrix[i][j][0], A2_matrix[i][j][1], step);
				func_write_to_file(ofp_conformation, strand_list);
				fprintf(ofp_conformation, "*\n");
				count -= 1;
			}
		}
	}
	fclose(ofp_conformation);
}
void normlized()
{
	for (int k = 0; k < num_row + 1; k++)
	{
		for (int i = 0; i < num_row; i++)
		{
			for (int j = 0; j < num_row; j++)
			{
				if (A_matrix[i][j][k] == Infinite)
						A_matrix[i][j][k] = A_matrix[i][j][k] * -1;
				else
					if (A_matrix[i][j][k] != 0)
						A_matrix[i][j][k] = (-1)*A_matrix[i][j][k]; // return to positive values
			}
		}
	}
}
void func_P_or_AP(int x, int y, int num_row)
{
	vector<int> new_pair;
	char new_interaction;
	if (x < (num_row / 2) && y < (num_row / 2))
	{
		new_pair.push_back(x);
		new_pair.push_back(y);
		new_interaction = 'P';
	}
	if (x >= (num_row / 2) && y >= (num_row / 2))
	{
		new_pair.push_back(x - (num_row / 2));
		new_pair.push_back(y - (num_row / 2));
		new_interaction = 'P';
	}
	if (x < (num_row / 2) && y >= (num_row / 2))
	{
		new_pair.push_back(x);
		new_pair.push_back(y - (num_row / 2));
		new_interaction = 'A';
	}
	if (x >= (num_row / 2) && y< (num_row / 2))
	{
		new_pair.push_back(x - (num_row / 2));
		new_pair.push_back(y);
		new_interaction = 'A';
	}
	interaction_type.push_back(new_interaction);
	pair_strands.push_back(new_pair);
}
void func_print_sheet(int max_row, int max_col, int step)// , vector<vector<double>> A_matrix, vector<vector<int>> P_matrix, vector<vector<vector<int>>> Path_matrix, int num_row)// clique algorithm find maximum sheets )
{
	if (P_matrix[max_row][max_col][step] != 0)
	{
		for (int k = 1; k < P_matrix[max_row][max_col][step]; k++)
			func_P_or_AP(Path_matrix[max_row][max_col][step][k - 1], Path_matrix[max_row][max_col][step][k], num_row);
	}
	else
		if (max_row != max_col)
			func_P_or_AP(max_row, max_col, num_row);
}
void func_write_to_file(FILE* ofp_conformation, vector<strand > strand_list)
{
	fprintf(ofp_conformation, "Num of pairs: %d\n", interaction_type.size());
	for (int i = 0; i < interaction_type.size(); i++)
	{
		fprintf(ofp_conformation, "%d-", pair_strands[i][0] + 1);
		fprintf(ofp_conformation, "%d-", pair_strands[i][1] + 1);
		fprintf(ofp_conformation, "%c ", interaction_type[i]);
	}
}
void initialization(vector<vector<double> > weights_matrix, int num, vector<strand > strand_list)
{
	num_row = num;
	A_matrix.resize(num_row, vector<vector<double>>(num_row, vector<double>(num_row + 1)));
	P_matrix.resize(num_row, vector<vector<int>>(num_row, vector<int>(num_row + 1)));// Matrix for storing length of the path
	Path_matrix.resize(num_row, vector<vector<vector<int>>>(num_row, vector<vector<int>>(num_row + 1, vector<int>(num_row)))); // Matrix for storing best path(include all nodesin the path) between pairs in each step

	for (int i = 0; i < num_row / 2; i++)
	{
		for (int j = 0; j < num_row; j++)
		{
			if ((i != j) && (i + (num_row / 2) != j) && (j + (num_row / 2) != i))
			{
				A_matrix[i][j][0] = (-1) * weights_matrix[i][j];// change from finding longest path problem to shortest path
				if (weights_matrix[i][j] == 0)
					A_matrix[i][j][0] = Infinite;
			}
			else
				A_matrix[i][j][0] = Infinite;
			P_matrix[i][j][0] = 0;
			for (int k = 0; k < num_row; k++)
				Path_matrix[i][j][0][k] = -1;
		}
	}
	for (int i = num_row / 2; i < num_row; i++)
	{
		for (int j = 0; j < num_row / 2; j++)
		{
			if ((i != j) && (i - (num_row / 2) != j) && (j + (num_row / 2) != i))
			{
				A_matrix[i][j][0] = (-1) * weights_matrix[j][i];

				if (weights_matrix[j][i] == 0)
					A_matrix[i][j][0] = Infinite;
			}
			else
				A_matrix[i][j][0] = Infinite;
			P_matrix[i][j][0] = 0;
			for (int k = 0; k < num_row; k++)
				Path_matrix[i][j][0][k] = -1;
		}
		for (int j = num_row / 2; j < num_row; j++)
		{
			if ((i != j) && (i - (num_row / 2) != j) && (j - (num_row / 2) != i))
			{
				A_matrix[i][j][0] = (-1) * weights_matrix[i - (num_row / 2)][j - (num_row / 2)];
				if (weights_matrix[i - (num_row / 2)][j - (num_row / 2)] == 0)
					A_matrix[i][j][0] = Infinite;
			}
			else
				A_matrix[i][j][0] = Infinite;
			P_matrix[i][j][0] = 0;
			for (int k = 0; k < num_row; k++)
				Path_matrix[i][j][0][k] = -1;
		}
	}
}
vector<int> which_to_be_omitted(int i, int j, int k)
{
	vector<int> omit; vector<int> check;
	omit.resize(num_row / 2);
	check.resize(num_row / 2);
	vector<vector<int>> position;
	position.resize((num_row / 2), vector <int>(2));
	for (int q = 0; q < num_row / 2; q++)
	{
		check[q] = 0;
		position[q][0] = -1;
		position[q][1] = -1;
	}
	if (i < (num_row / 2))
	{
		check[i] = 1;
		position[i][0] = -1;
	}
	else
	{
		check[i - (num_row / 2)] = 1;
		position[i - (num_row / 2)][0] = -1;
	}
	if (j < (num_row / 2))
	{
		check[j] = 1;
		position[j][1] = P_matrix[k][j][k];
	}
	else
	{
		check[j - (num_row / 2)] = 1;
		position[j - (num_row / 2)][1] = P_matrix[k][j][k];
	}
	if (k < (num_row / 2))
	{
		check[k] = 1;
		position[k][1] = -1;
		position[k][0] = P_matrix[i][k][k];
	}
	else
	{
		check[k - (num_row / 2)] = 1;
		position[k - (num_row / 2)][1] = -1;
		position[k - (num_row / 2)][0] = P_matrix[i][k][k];
	}
	for (int q = 0; q < P_matrix[i][k][k]; q++)
	{
		if (Path_matrix[i][k][k][q] < (num_row / 2))
		{
			check[Path_matrix[i][k][k][q]] += 1;
			position[Path_matrix[i][k][k][q]][0] = q;
		}
		else
		{
			check[Path_matrix[i][k][k][q] - (num_row / 2)] += 1;
			position[Path_matrix[i][k][k][q] - (num_row / 2)][0] = q;

		}
	}
	for (int q = 0; q < P_matrix[k][j][k]; q++)
	{
		if (Path_matrix[k][j][k][q] < (num_row / 2))
		{
			check[Path_matrix[k][j][k][q]] += 1;
			position[Path_matrix[k][j][k][q]][1] = q;
		}
		else
		{
			check[Path_matrix[k][j][k][q] - (num_row / 2)] += 1;
			position[Path_matrix[k][j][k][q] - (num_row / 2)][1] = q;

		}
	}
	//select which commen elements is omitted
	double cs1 = 0, cs2 = 0;
	int index1 = -1, index2 = -1, index3 = -1;
	vector<int>  path1, path2;
	for (int q = 0; q < num_row / 2; q++)
	{
		omit[q] = -1;
	}
	for (int q = 0; q < num_row / 2; q++)
	{
		path1.clear();
		path2.clear();
		if (check[q] > 1)
		{
			if ((q == i) || (q == i - (num_row / 2)))
				omit[q] = 1;
			else
				if (q == j || (q == j - (num_row / 2)))
					omit[q] = 0;
				else
					if (q == k || (q == k - (num_row / 2)))
					{
						//cout << "error";
					}
					else
					{
						if (position[q][0] == 0)
						{

							if (P_matrix[i][k][k] == 1)
							{
								cs1 = A_matrix[i][k][0];
								//path1.push_back(-1);
							}
							else
							{
								index2 = Path_matrix[i][k][k][1];
								cs1 = A_matrix[i][index2][0];// +A_matrix[index2][k][k];//A_matrix[index1][index2][0] + A_matrix[index2][index3][0];
								path1.push_back(index2);
								for (int gg = 1; gg < Path_matrix[i][k][k].size() - 1; gg++)
								{
									cs1 += A_matrix[Path_matrix[i][k][k][gg]][Path_matrix[i][k][k][gg + 1]][0];
									path1.push_back(Path_matrix[i][k][k][gg + 1]);
								}
								cs1 += A_matrix[Path_matrix[i][k][k][Path_matrix[i][k][k].size() - 1]][k][0];
							}
						}
						else
						{
							if (position[q][0] == P_matrix[i][k][k] - 1)
							{
								index1 = Path_matrix[i][k][k][P_matrix[i][k][k] - 1];
								index3 = k;
								cs1 = A_matrix[i][Path_matrix[i][k][k][0]][0];
								path1.push_back(Path_matrix[i][k][k][0]);
								for (int gg = 0; gg < Path_matrix[i][k][k].size() - 2; gg++)
								{
									cs1 += A_matrix[Path_matrix[i][k][k][gg]][Path_matrix[i][k][k][gg + 1]][0];
									path1.push_back(Path_matrix[i][k][k][gg + 1]);
								}
								cs1 += A_matrix[Path_matrix[i][k][k][P_matrix[i][k][k] - 2]][k][0];
							}
							else
							{
								index1 = Path_matrix[i][k][k][position[q][0] - 1];
								index3 = Path_matrix[i][k][k][position[q][0] + 1];
								cs1 = A_matrix[index1][index3][0] + A_matrix[i][Path_matrix[i][k][k][0]][0];
								path1.push_back(Path_matrix[i][k][k][0]);
								for (int gg = 0; gg < position[q][0] - 2; gg++)
								{
									cs1 += A_matrix[Path_matrix[i][k][k][gg]][Path_matrix[i][k][k][gg + 1]][0];
									path1.push_back(Path_matrix[i][k][k][gg + 1]);
								}
								path1.push_back(index3);
								for (int gg = position[q][0] + 1; gg < Path_matrix[i][k][k].size() - 1; gg++)
								{
									cs1 += A_matrix[Path_matrix[i][k][k][gg]][Path_matrix[i][k][k][gg + 1]][0];
									path1.push_back(Path_matrix[i][k][k][gg + 1]);
								}
								cs1 += A_matrix[Path_matrix[i][k][k][Path_matrix[i][k][k].size() - 1]][k][0];
							}
						}
						if (position[q][1] == 0)
						{
							if (P_matrix[k][j][k] == 1)
							{
								cs2 = A_matrix[k][j][0];//A_matrix[index1][index2][0] + A_matrix[index2][index3][0];
								//path2.push_back(-1);
							}
							else
							{
								index3 = Path_matrix[k][j][k][1];
								cs2 = A_matrix[k][index3][0];
								path2.push_back(index3);
								for (int gg = 1; gg < Path_matrix[k][j][k].size() - 1; gg++)
								{
									cs2 += A_matrix[Path_matrix[k][j][k][gg]][Path_matrix[k][j][k][gg + 1]][0];
									path2.push_back(Path_matrix[k][j][k][gg + 1]);
								}
								cs2 += A_matrix[Path_matrix[k][j][k][Path_matrix[k][j][k].size() - 1]][j][0];
							}
						}
						else
						{
							if (position[q][1] == P_matrix[k][j][k] - 1)
							{
								index1 = Path_matrix[k][j][k][P_matrix[k][j][k] - 1];
								cs2 = A_matrix[k][Path_matrix[k][j][k][0]][0];
								path2.push_back(Path_matrix[k][j][k][0]);
								for (int gg = 0; gg < Path_matrix[k][j][k].size() - 2; gg++)
								{
									cs2 += A_matrix[Path_matrix[k][j][k][gg]][Path_matrix[k][j][k][gg + 1]][0];
									path2.push_back(Path_matrix[k][j][k][gg + 1]);
								}
								cs2 += A_matrix[Path_matrix[k][j][k][P_matrix[k][j][k] - 2]][j][0];
							}
							else
							{
								index1 = Path_matrix[k][j][k][position[q][1] - 1];
								index2 = Path_matrix[k][j][k][position[q][1]];
								index3 = Path_matrix[k][j][k][position[q][1] + 1];
								cs2 = A_matrix[index1][index3][0] + A_matrix[k][Path_matrix[k][j][k][0]][0];
								path2.push_back(Path_matrix[k][j][k][0]);
								for (int gg = 0; gg < position[q][1] - 2; gg++)
								{
									cs2 += A_matrix[Path_matrix[k][j][k][gg]][Path_matrix[k][j][k][gg + 1]][0];
									path2.push_back(Path_matrix[k][j][k][gg + 1]);
								}
								path2.push_back(index2);
								for (int gg = position[q][1] + 1; gg < Path_matrix[k][j][k].size() - 1; gg++)
								{
									cs2 += A_matrix[Path_matrix[k][j][k][gg]][Path_matrix[k][j][k][gg + 1]][0];
									path2.push_back(Path_matrix[k][j][k][gg + 1]);
								}
								cs2 += A_matrix[Path_matrix[k][j][k][Path_matrix[k][j][k].size() - 1]][j][0];
							}
						}
						int length_path1 = path1.size(), length_path2 = path2.size();
						for (int qq = 0; qq < path1.size(); qq++)
						{
							if (path1[qq] == j || (path1[qq] > (num_row / 2) && (path1[qq] == j - (num_row / 2))))
							{
								length_path1--;
								if (path1.size() == 1)
								{
									cs1 = A_matrix[i][k][0];
								}
								else
								{
									if (qq == 0)
									{
										cs1 = cs1 - (A_matrix[i][path1[0]][0] + A_matrix[path1[0]][path1[1]][0]) + A_matrix[i][path1[1]][0];
									}
									else
									{
										if (path1.size() - 1 == qq)
											cs1 = cs1 - (A_matrix[path1[qq]][k][0] + A_matrix[path1[qq - 1]][path1[qq]][0]) + A_matrix[path1[qq - 1]][k][0];
										else
											cs1 = cs1 - (A_matrix[path1[qq - 1]][path1[qq]][0] + A_matrix[path1[qq]][path1[qq + 1]][0]) + A_matrix[path1[qq - 1]][path1[qq + 1]][0];
									}
								}
							}
						}

						for (int qq = 0; qq < path2.size(); qq++)
						{
							if (path2[qq] == i || (path2[qq] > (num_row / 2) && path2[qq] == i - (num_row / 2)))
							{
								length_path2--;
								if (path2.size() == 1)
								{
									cs2 = A_matrix[k][j][0];
								}
								else
								{
									if (qq == 0)
									{
										cs2 = cs2 - (A_matrix[k][path2[0]][0] + A_matrix[path2[0]][path2[1]][0]) + A_matrix[k][path2[1]][0];
									}
									else
									{
										if (path2.size() - 1 == qq)
											cs2 = cs2 - (A_matrix[path2[qq]][j][0] + A_matrix[path2[qq - 1]][path2[qq]][0]) + A_matrix[path2[qq - 1]][j][0];
										else
											cs2 = cs2 - (A_matrix[path2[qq - 1]][path2[qq]][0] + A_matrix[path2[qq]][path2[qq + 1]][0]) + A_matrix[path2[qq - 1]][path2[qq + 1]][0];
									}
								}
							}
						}
						cs1 += A_matrix[k][j][k] ;
						cs2 += A_matrix[i][k][k];
						//avg: cs1 = cs1 / ((double)(length_path1+1));
						//avg: cs2 = cs2 / ((double)(length_path2+1));
						if (cs1 < cs2)
							omit[q] = 0;
						else
							omit[q] = 1;
					}
		}
	}
	return(omit);
}
void delete_common_nodes(vector<int> omit,int i, int j, int k)
{
	vector<int> check;
	double cs1 = 0, cs2 = 0;
	int index1 = -1, index2 = -1, index3 = -1;
	int acceptable = 0;
	check.resize(num_row / 2);
	if (i < (num_row / 2))
		check[i] = 1;
	else
		check[i - (num_row / 2)] = 1;
	if (j < (num_row / 2))
		check[j] = 1;
	else
		check[j - (num_row / 2)] = 1;
	if (k < (num_row / 2))
		check[k] = 1;
	else
		check[k - (num_row / 2)] = 1;
	index1 = i; index2 = -1; cs1 = 0;
	vector<int> new_path;
	new_path.clear();
	for (int q = 0; q < P_matrix[i][k][k]; q++)
	{
		index3 = Path_matrix[i][k][k][q];
		if (index3 >= (num_row / 2))
			index3 -= (num_row / 2);
		if (omit[index3] != 0)
		{
			if (check[index3] == 0)
			{
				index2 = Path_matrix[i][k][k][q];
				if (A_matrix[index1][index2][0] == Infinite)
					acceptable = -1;
				cs1 += A_matrix[index1][index2][0];
				if (index2 < (num_row / 2))
				{
					if (check[index2] > 0)
					{//cout << "error";
					}

					else
						check[index2] += 1;
				}
				else
				{
					if (check[index2 - (num_row / 2)] > 0)
					{//cout << "error";
					}
					else
						check[index2 - (num_row / 2)] += 1;
				}
				new_path.push_back(index2);
				index1 = index2;
			}
		}
		else
		{
			if (q < P_matrix[i][k][k] - 1 && (P_matrix[i][k][k] > 1))
			{
				index3 = Path_matrix[i][k][k][q + 1];
				if (index3 >= (num_row / 2))
					index3 -= (num_row / 2);
				if (omit[index3] != 0)
				{
					index2 = Path_matrix[i][k][k][q + 1];
					if (A_matrix[index1][index2][0] == Infinite)
						acceptable = -1;
					cs1 += A_matrix[index1][index2][0];
					if (index2 < (num_row / 2))
					{
						if (check[index2] > 0)
						{//cout << "error";
						}
						else
							check[index2] += 1;
					}
					else
					{
						if (check[index2 - (num_row / 2)] > 0)
						{//cout << "error";
						}
						else
							check[index2 - (num_row / 2)] += 1;
					}
					new_path.push_back(index2);
					index1 = index2;
				}
			}
		}
	}
	index2 = k;
	if (A_matrix[index1][index2][0] == Infinite)
		acceptable = -1;
	cs1 += A_matrix[index1][index2][0];
	index1 = index2;
	new_path.push_back(k);
	for (int q = 0; q < P_matrix[k][j][k]; q++)
	{
		index3 = Path_matrix[k][j][k][q];
		if (index3 >= (num_row / 2))
			index3 -= (num_row / 2);
		if (omit[index3] != 1)
		{
			if (check[index3] == 0)
			{
				index2 = Path_matrix[k][j][k][q];
				if (A_matrix[index1][index2][0] == Infinite)
					acceptable = -1;
				cs1 += A_matrix[index1][index2][0];
				if (index2 < (num_row / 2))
				{
					if (check[index2] > 0)
					{//cout << "error";
					}
					else
						check[index2] += 1;
				}
				else
				{
					if (check[index2 - (num_row / 2)] > 0)
					{//cout << "error";
					}
					else
						check[index2 - (num_row / 2)] += 1;
				}
				new_path.push_back(index2);
				index1 = index2;
			}
		}
		else
		{
			if (q < P_matrix[k][j][k] - 1 && (P_matrix[k][j][k] > 1))
			{
				index3 = Path_matrix[k][j][k][q + 1];
				if (index3 >= (num_row / 2))
					index3 -= (num_row / 2);
				if (omit[index3] != 1)
				{
					index2 = Path_matrix[k][j][k][q + 1];
					if (A_matrix[index1][index2][0] == Infinite)
						acceptable = -1;
					cs1 += A_matrix[index1][index2][0];
					if (index2 < (num_row / 2))
					{
						if (check[index2] > 0)
						{//cout << "error";
						}
						else
							check[index2] += 1;
					}
					else
					{
						if (check[index2 - (num_row / 2)] > 0)
						{//cout << "error";
						}
						else
							check[index2 - (num_row / 2)] += 1;
					}
					new_path.push_back(index2);
					index1 = index2;
				}
			}
		}
	}
	index2 = j;
	if (A_matrix[index1][index2][0] == Infinite)
		acceptable = -1;
	cs1 += A_matrix[index1][index2][0];
	for (int gg = 0; gg < num_row / 2; gg++)
		if (check[gg] > 1)
		{//cout << "error";
		}
	/*if (acceptable == -1)
		cout << "error";*/
	if (cs1 <= A_matrix[i][j][k] && (acceptable == 0))
	{
		//avg: A_matrix[i][j][k + 1] = cs1/((double)(new_path.size()+1));
		A_matrix[i][j][k + 1] = cs1;
		P_matrix[i][j][k + 1] = new_path.size();
		Path_matrix[i][j][k + 1] = new_path;
	}
}
double Con_APSP(vector<vector<vector<double>>> &SCORE_matrix, vector<vector<vector<vector<int>>>> &PATH_matrix, vector<vector<double> > weights_matrix, int num, vector<strand > strand_list, char *result_fold_path_str)// , char *path_floyed_matrix)
{
	double score = 0;
	char * floyed_result_fold_path_str = new char[1000];
	strcpy(floyed_result_fold_path_str, result_fold_path_str);
	FILE* ofp;
	clock_t begin, end;
	double time_spent = 0;
	initialization(weights_matrix, num, strand_list);
	begin = clock();
	for (int k = 0; k < num_row; k++)
	{
		for (int i = 0; i < num_row; i++)
		{
			for (int j = 0; j < num_row; j++)
			{
				A_matrix[i][j][k + 1] = A_matrix[i][j][k];
				P_matrix[i][j][k + 1] = P_matrix[i][j][k];
				Path_matrix[i][j][k + 1] = Path_matrix[i][j][k];
			}
		}
		for (int i = 0; i < (num_row); i++)
		{
			for (int j = 0; j < num_row; j++)
			{
				if ((i != j) && (j != k) && (k != i) && (A_matrix[k][j][k] != Infinite) && (A_matrix[i][k][k] != Infinite) && (P_matrix[i][j][k] < (num_row / 2)) && (i != j + (num_row / 2)) && (j != k + (num_row / 2)) && (k != i + (num_row / 2)) && (j != i + (num_row / 2)) && (k != j + (num_row / 2)) && (i != k + (num_row / 2)))
				{
					if ((A_matrix[i][k][k]  + A_matrix[k][j][k] ) <= A_matrix[i][j][k])
				  //avg: if (((((A_matrix[i][k][k]*(P_matrix[i][k][k]+1)) + (A_matrix[k][j][k] * (P_matrix[k][j][k]+1)))/(double)(P_matrix[i][k][k]+ P_matrix[k][j][k]+2))) <= A_matrix[i][j][k] )
					{
						if ((no_similar_strand(i, k, k, j, k, k) == 0) && (no_Complement_strands(i, k, k, j, k, k) == 0))// check all except equality of k,k
						{
							//avg: A_matrix[i][j][k + 1] = ((((A_matrix[i][k][k] * (P_matrix[i][k][k]+1)) + (A_matrix[k][j][k] * (P_matrix[k][j][k]+1))) / (double)(P_matrix[i][k][k] + P_matrix[k][j][k]+2)));
								A_matrix[i][j][k + 1] =  (A_matrix[i][k][k] + A_matrix[k][j][k]);// / (A_matrix[i][k] * A_matrix[k][j]);
							P_matrix[i][j][k + 1] = P_matrix[i][k][k] + P_matrix[k][j][k] + 1;
							vector<int> temp = merge_paths(i, j, k, k, k);
							for (int w = 0; w < P_matrix[i][j][k]; w++)
								Path_matrix[i][j][k + 1][w] = -1;
							Path_matrix[i][j][k + 1] = temp;
						}
						else
						{
							vector<int> omit;
							omit= which_to_be_omitted(i, j, k);	//select which commen elements is omitted						
							delete_common_nodes(omit, i, j, k);	//omit commen elements

						}//else similar 
					}//	if ((i != j) && (j != k) && (k != i) && (A_matrix[k][j][k] != Infinite)&& (A_matrix[i][j][k] != Infinite) && (A_matrix[i][k][k] != Infinite) && (P_matrix[i][j][k] < (num_row / 2)))
				}//for (int j = 0; j < num_row; j++)
			}
		}//for (int i = 0; i < num_row; i++)		
	}//for (int k = 0; k < num_row; k++)
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	normlized();
	
	double max_score = 0;
	
	for (int i = 0; i < num_row; i++)
	{
		for (int j = 0; j < num_row; j++)
		{
			int num_non_double = 0;
			double dd = -99999;
			//insert start and end of paths in to matrix path
			for (int k = 0; k < num_row + 1; k++)//copy entire path
			{
				//	if (((k == num_row)&&((i < (num_row / 2))) )||(k==0))//omit 	Symmetry
				if (((k <= num_row))) //|| (k == 0))

				{
					if ((dd < A_matrix[i][j][k]) || (num_row == k))
					{
						num_non_double++;
						dd = A_matrix[i][j][k];
						vector<int> p;
						p.push_back(i);
						for (int w = 0; w < P_matrix[i][j][k]; w++)
						{
							p.push_back(Path_matrix[i][j][k][w]);
						}
						p.push_back(j);
						Path_matrix[i][j][k] = p;
						P_matrix[i][j][k] = p.size();
						if (max_score < A_matrix[i][j][k])
							max_score = A_matrix[i][j][k];
					}
					else//delete repitations
					{
						Path_matrix[i][j][k].clear();
						P_matrix[i][j][k] = 0;
						A_matrix[i][j][k] = (-1 * Infinite);

					}

				}
				else //delete intermidate paths
				{
					Path_matrix[i][j][k].clear();
					P_matrix[i][j][k] = 0;
					A_matrix[i][j][k] = (-1 * Infinite);
				}
			}

			//find weak interactions in the paths at last stepymmetry
			
			vector<paths> new_paths_array;
			if (P_matrix[i][j][num_row] > 2 && A_matrix[i][j][num_row] != (-1 * Infinite))
			{
				//less than avg at each sheet
				vector<bool> week_interactions;
				week_interactions.clear();
				week_interactions.resize(P_matrix[i][j][num_row] - 1);
				for (int k = 0; k < P_matrix[i][j][num_row] - 1; k++)
					if (A_matrix[Path_matrix[i][j][num_row][k]][Path_matrix[i][j][num_row][k + 1]][0] < (double)A_matrix[i][j][num_row] / (double)(P_matrix[i][j][num_row]))
						week_interactions[k] = true;
				paths new_paths;
				for (int k = 0; k < P_matrix[i][j][num_row] - 1; k++)
				{
					if (week_interactions[k] == true)
					{
						if ((0 < k) && (k < P_matrix[i][j][num_row] - 2))
						{
							vector<int> p;
							if (k > 1)
							{
								for (int w = 0; w < k + 1; w++)
									p.push_back(Path_matrix[i][j][num_row][w]);
								new_paths.path=p;
								new_paths.score = 0;
								new_paths_array.push_back(new_paths);
							}
							p.clear();
							if (k < P_matrix[i][j][num_row] - 3)
							{
								for (int w = k + 1; w < P_matrix[i][j][num_row]; w++)
									p.push_back(Path_matrix[i][j][num_row][w]);
								new_paths.path=p; 
								new_paths.score = 0;
								new_paths_array.push_back(new_paths);
							}
						}
						else
						{
							vector<int> p;
							if (k == 0)
							{
								for (int w = 1; w < P_matrix[i][j][num_row]; w++)
									p.push_back(Path_matrix[i][j][num_row][w]);
							}
							else
							{
								for (int w = 0; w < P_matrix[i][j][num_row] - 1; w++)
									p.push_back(Path_matrix[i][j][num_row][w]);

							}
							new_paths.path = p;
							new_paths.score = 0;
							new_paths_array.push_back(new_paths);
						}
					}
				}
				for (int k = 0;  k < new_paths_array.size(); k++)
				{
					new_paths_array[k].score = 0;
					for(int q=1;q<new_paths_array[k].path.size()-1;q++)
						new_paths_array[k].score+= A_matrix[new_paths_array[k].path[q-1]][new_paths_array[k].path[q]][0];
				}
				sort(new_paths_array.begin(), new_paths_array.end(), sortByScore);
				reverse(new_paths_array.begin(), new_paths_array.end());		

				//replace new_paths with intermidiate paths
				if (new_paths_array.size() + num_non_double > num_row)
				{
				//	cout << "error in new paths";
					new_paths_array._Pop_back_n( (new_paths_array.size() + num_non_double)-num_row-1);
				}
				int qq = 0;
				for (int k = 0; k < num_row && qq < (new_paths_array.size()); k++)
				{
					if (A_matrix[i][j][k] == (-1 * Infinite))
					{
						
						Path_matrix[i][j][k] = new_paths_array[qq].path;
						P_matrix[i][j][k] = new_paths_array[qq].path.size();
						A_matrix[i][j][k] = new_paths_array[qq].score;
						qq++;
					}
				}
			}
		}
	}
	SCORE_matrix = A_matrix;
	PATH_matrix = Path_matrix;
	//fclose(ofp_floyed_matrix);
	//strcat(floyed_result_fold_path_str, ".txt");
	//func_write_to_file_all_floyedwarshall_sheets(floyed_result_fold_path_str, num_row, strand_list);
	
	// free memory
	for (int i = 0; i < A_matrix.size(); i++)
		A_matrix[i].clear();
	A_matrix.clear();
	for (int i = 0; i < P_matrix.size(); i++)
		P_matrix[i].clear();
	P_matrix.clear();
	for (int i = 0; i < Path_matrix.size(); i++)
		Path_matrix[i].clear();
	Path_matrix.clear();
	interaction_type.clear();
	for (int i = 0; i < pair_strands.size(); i++)
		pair_strands[i].clear();
	pair_strands.clear();
	return(time_spent);
	
}
