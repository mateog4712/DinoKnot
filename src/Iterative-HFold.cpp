// Iterative HFold files
#include "iterative-hfold.hh"
#include "hotspot.hh"
#include "Result.hh"
#include "cmdline.hh"
#include "W_final.hh"
#include "h_externs.hh"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <getopt.h>


std::string remove_structure_intersection(std::string restricted, std::string structure){
	cand_pos_t length = structure.length();
	for(cand_pos_t i=0; i< length; ++i){
		if(restricted[i] == '(' || restricted[i] == ')') structure[i] = '.';
		
		if (structure[i] == '[') structure[i] = '(';
		
		if (structure[i] == ']') structure[i] = ')';

		if(restricted[i] == 'x') structure[i] = 'x';
	}
	return structure;
}
// /**
//  * @brief returns a vector of pairs which represent the start and end indices for each disjoint substructure in the structure
//  * 
//  * @param CL_ Candidate list
//  * @return total number of candidates
//  */
void find_disjoint_substructure(std::string structure, std::vector< std::pair<int,int> > &pair_vector){
	cand_pos_t n = structure.length();
	cand_pos_t count = 0;
	cand_pos_t i = 0;
	for(cand_pos_t k=0; k<n;++k){
		if(structure[k] == '('){
			if(count == 0) i = k;
			count++;

		}else if(structure[k] == ')'){
			count--;
			if(count == 0){
				std::pair <int,int> ij_pair (i,k);
				pair_vector.push_back(ij_pair);
			}
		}
	}
}
/**
 * @brief Fills the pair array
 * p_table will contain the index of each base pair
 * X or x tells the program the base cannot pair and . sets it as unpaired but can pair
 * @param structure Input structure
 * @param p_table Restricted array
 */
void detect_pairs(const std::string &structure, std::vector<cand_pos_t> &p_table){
	cand_pos_t i, j, length = structure.length();
	std::vector<cand_pos_t>  pairs;
	pairs.push_back(length);

	for (i=length-1; i >=0; --i){
		if ((structure[i] == 'x') || (structure[i] == 'X'))
			p_table[i] = -1;
		else if (structure[i] == '.')
			p_table[i] = -2;
		if (structure[i] == ')'){
			pairs.push_back(i);
		}
		if (structure[i] == '('){
			j = pairs[pairs.size()-1];
			pairs.erase(pairs.end()-1);
			p_table[i] = j;
			p_table[j] = i;
		}
	}
	pairs.pop_back();
	if (pairs.size() != 0) std::cout << pairs[0] << std::endl << structure << std::endl;
	if (pairs.size() != 0)
	{
		fprintf (stderr, "The given structure is not valid: more left parentheses than right parentheses: \n");
		exit (1);
	}
}

cand_pos_t paired_structure(cand_pos_t i, cand_pos_t j, std::vector<cand_pos_t> &pair_index){
	cand_pos_t n = pair_index.size();
	return (i >= 0 && j < n && (pair_index[i] == j));
}
/**
 * @brief Takes the input constraint structure as the base output structure. Finds the stacking bases, bulges of size 1, internal loops of size 1x1, 2x1, and 1x2 on the pk_free structure
 * which have these substructures forming on the base pairs and adds them to the output structure.
 * 
 * @param restricted input constraint structure used as the base output structure
 * @param pkfree_structure structure post pseudoknot-free base pair filling.
*/
std::string obtainRelaxedStems(std::string restricted, std::string pkfree_structure){
	cand_pos_t n = restricted.length();

	//Gresult <- G1
	std::string relaxed = restricted;

	cand_pos_t i = 0;
	cand_pos_t j = 0;
	
	std::vector<cand_pos_t> G1_pair;
	std::vector<cand_pos_t> G2_pair;
	G1_pair.resize(n,-2);
	G2_pair.resize(n,-2);
	detect_pairs(restricted,G1_pair);
	detect_pairs(pkfree_structure,G2_pair);

	
	for(int k=0;k<n;++k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (restricted[i] != pkfree_structure[i])){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i-1,j+1,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include bulges of size 1
					}else if(paired_structure(i-2,j+1,G1_pair) || paired_structure(i-1,j+2,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i-2,j+2,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if( paired_structure(i-3,j+2,G1_pair) || paired_structure(i-2,j+3,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}
	for(int k=n-1;k>=0;--k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (restricted[i] != pkfree_structure[i])){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i+1,j-1,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
						
					//include bulges of size 1
					}else if(paired_structure(i+1,j-2,G1_pair) || paired_structure(i+2,j-1,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i+2,j-2,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if(paired_structure(i+2,j-3,G1_pair) || paired_structure(i+3,j-2,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}
	return relaxed;
}

std::string remove_x(std::string structure){
	for(char &c: structure) {if(c == 'x') c = '.';}
	return structure; 
}

std::string hfold(std::string seq,std::string res, energy_t &energy, bool pk_free, bool pk_only, vrna_param_s *params1, vrna_param_s *params2){
	sparse_tree tree(res,res.length());
	W_final min_fold(seq,res, pk_free, pk_only,params1,params2);
	energy = min_fold.hfold_interacting(tree);
    std::string structure = min_fold.structure;
    return structure;
}

std::string method2(std::string &seq, std::string &restricted, energy_t &method2_energy, vrna_param_s *params1, vrna_param_s *params2){

	std::string pk_only_output = hfold(seq,restricted,method2_energy,false,true,params1,params2);
	std::string pk_free_removed = remove_structure_intersection(restricted,pk_only_output);
	std::string no_x_restricted = remove_x(restricted);

	if(pk_only_output != no_x_restricted) return hfold(seq,pk_free_removed,method2_energy,false,false,params1,params2);
	else return pk_only_output;
}

std::string Iterative_HFold_interacting (std::string seq,std::string res, double &en, vrna_param_s *params1, vrna_param_s *params2, int &method_chosen, bool hard){

	// Iterate through all hotspots or the single given input structure
	energy_t method1_energy = INF;
	energy_t method2_energy = INF;
	energy_t method3_energy = INF;
	energy_t method4_energy = INF;
	energy_t final_en = INF;
	std::string final_structure;

	//Method1
	std::string method1_structure = hfold(seq,res,method1_energy,false,false,params1,params2);
	// std::string method1_structure = hfold(seq,res,method1_energy,false,false,2);
	if(method1_energy < final_en){
	final_en = method1_energy;
	final_structure=method1_structure;
	method_chosen = 1;
	}

	if(!hard){
		//Method2
		std::string method2_structure = method2(seq,res,method2_energy,params1,params2);
		if(method2_energy < final_en){
		final_en = method2_energy;
		final_structure=method2_structure;
		method_chosen = 2;
		}

		//Method3
		std::string pk_free = hfold(seq,res,method3_energy,true,false,params1,params2);
		std::string relaxed = obtainRelaxedStems(res,pk_free);
		cand_pos_t n = res.length();
		for(cand_pos_t i =0; i< n;++i) if(res[i] == 'x') relaxed[i] = 'x';
		std::string method3_structure = method2(seq,relaxed,method3_energy,params1,params2);
		if(method3_energy < final_en){
			final_en = method3_energy;
			final_structure=method3_structure;
		method_chosen = 3;
		}

		//Method4
		std::vector< std::pair<int,int> > disjoint_substructure_index;
		find_disjoint_substructure(res,disjoint_substructure_index);
		std::string disjoint_structure = res;
		for(auto current_substructure_index : disjoint_substructure_index){
			cand_pos_t i = current_substructure_index.first;
			cand_pos_t j = current_substructure_index.second;
			energy_t energy = INF;

			std::string subsequence = seq.substr(i,j-i+1);
			std::string substructure = res.substr(i,j-i+1);

			std::string pk_free = hfold(subsequence,substructure,energy,true,false,params1,params2);
			std::string relaxed = obtainRelaxedStems(substructure,pk_free);
			cand_pos_t sub_n = substructure.length();
			for(int i =0; i< sub_n;++i) if(substructure[i] == 'x') relaxed[i] = 'x';
			disjoint_structure.replace(i,j-i+1,relaxed);
		}
		std::string method4_structure = method2(seq,disjoint_structure,method4_energy,params1,params2);
		if(method4_energy < final_en){
			final_en = method4_energy;
			final_structure=method4_structure;
		method_chosen = 4;
		}
	}

	double final_energy = final_en/100.0;
	en = final_energy;

    return final_structure;
}
