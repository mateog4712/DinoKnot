#include "W_final.hh"
#include "h_struct.hh"
#include "h_externs.hh"
#include "common.hh"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#define debug 0

// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(std::string seq,std::string res,bool pk_free, bool pk_only, vrna_param_s *params1, vrna_param_s *params2)
{
	seq_ = seq;
	this->res = res;
	this->n = seq.length();
	params_ = params1;
	params2_ = params2;
	make_pair_matrix();
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);
	this->pk_free = pk_free;
	this->pk_only = pk_only;
	cand_pos_t total_length = ((n+1) *(n+2))/2;
	W.resize(total_length,0);
	space_allocation();
}


W_final::~W_final()
{
	delete WMB;
	delete V;
	delete [] f;
	// free(params_);
	free(S_);
	free(S1_);
}

// Hosna June 20th, 2007
// allocates space for WMB object and V_final
void W_final::space_allocation(){

	// From simfold
	f = new minimum_fold [n+1];

    V = new s_energy_matrix (seq_, n,S_,S1_,params_,params2_);
	structure = std::string (n+1,'.');

	// Hosna: June 20th 2007
    WMB = new pseudo_loop (seq_,res,V,S_,S1_,params_,params2_);

}
// m1 = V(i,j), m2= WMB(i,j), m3 = W(i,k-1) + V(k,j), m4 = W(i,k-1) + WMB(k,j), m5 = W(i,j-1)
void W_final::compute_W(cand_pos_t i, cand_pos_t j, sparse_tree tree){
	energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5=INF;
	cand_pos_t ij = V->index[i] + j - i;
	const bool cross_model = is_cross_model(i,j);
	if(cross_model){
		energy_t en1 = E_ext_Stem(V->get_energy(i,j),V->get_energy(i+1,j),V->get_energy(i,j-1),V->get_energy(i+1,j-1),S_,params_,i,j,n,tree.tree);
		energy_t en2 = E_ext_Stem(V->get_energy(i,j),V->get_energy(i+1,j),V->get_energy(i,j-1),V->get_energy(i+1,j-1),S_,params2_,i,j,n,tree.tree);
		m1 = emodel_energy_function(i,j,en1,en2);
		for(cand_pos_t k = i+1;k<j-TURN;++k){
			en1 = E_ext_Stem(V->get_energy(k,j),V->get_energy(k+1,j),V->get_energy(k,j-1),V->get_energy(k+1,j-1),S_,params_,k,j,n,tree.tree);
			en2 = E_ext_Stem(V->get_energy(k,j),V->get_energy(k+1,j),V->get_energy(k,j-1),V->get_energy(k+1,j-1),S_,params_,k,j,n,tree.tree);
			m3 = std::min(m3,get_energy(i,k-1) + emodel_energy_function(i,j,en1,en2));
			m4 = std::min(m4,get_energy(i,k-1) + WMB->get_WMB(k,j) + PS_penalty);
		}
	}
	else{
		if(j<linker_pos){
			m1 = E_ext_Stem(V->get_energy(i,j),V->get_energy(i+1,j),V->get_energy(i,j-1),V->get_energy(i+1,j-1),S_,params_,i,j,n,tree.tree);
			for(cand_pos_t k = i+1;k<j-TURN;++k){
				m3 = std::min(m3,get_energy(i,k-1) + E_ext_Stem(V->get_energy(k,j),V->get_energy(k+1,j),V->get_energy(k,j-1),V->get_energy(k+1,j-1),S_,params_,k,j,n,tree.tree));
				m4 = std::min(m4,get_energy(i,k-1) + WMB->get_WMB(k,j) + PS_penalty);
			}
		}
		else if(i>linker_pos_right){
			m1 = E_ext_Stem(V->get_energy(i,j),V->get_energy(i+1,j),V->get_energy(i,j-1),V->get_energy(i+1,j-1),S_,params2_,i,j,n,tree.tree);
			for(cand_pos_t k = i+1;k<j-TURN;++k){
				m3 = std::min(m3,get_energy(i,k-1) + E_ext_Stem(V->get_energy(k,j),V->get_energy(k+1,j),V->get_energy(k,j-1),V->get_energy(k+1,j-1),S_,params2_,k,j,n,tree.tree));
				m4 = std::min(m4,get_energy(i,k-1) + WMB->get_WMB(k,j) + PS_penalty);
			}
		}
	}
	m2 = WMB->get_WMB(i,j) + PS_penalty;
	if(tree.tree[j].pair<0) m5 = get_energy(i,j-1);
	W[ij] = std::min({m1,m2,m3,m4,m5});
}

energy_t W_final::hfold_interacting(sparse_tree &tree){

	for (cand_pos_t i = n; i >=1; --i)
	{	
		for (cand_pos_t j =i; j<=n; ++j)
		{
			const bool evaluate = tree.weakly_closed(i,j);
			const pair_type ptype_closing = pair[S_[i]][S_[j]];
			const bool restricted = tree.tree[i].pair == -1 || tree.tree[j].pair == -1;
			const bool paired = (tree.tree[i].pair == j && tree.tree[j].pair == i);

			const bool pkonly = (!pk_only || paired);

			if(ptype_closing> 0 && evaluate && !restricted && pkonly) V->compute_energy_restricted_emodel (i,j,tree);

			if(!pk_free) WMB->compute_energies_emodel(i,j,tree);

			V->compute_WMv_WMp_emodel(i,j,WMB->get_WMB(i,j),tree.tree);

			V->compute_energy_WM_restricted_emodel(i,j,tree,WMB->WMB); // if i>linkerpos || j < linker_pos_right?
			compute_W(i,j,tree);
			if(i<linker_pos && j>linker_pos_right) V->compute_VMprime(i,j,tree,WMB->WMB);
		}

	}
    energy_t energy = get_energy(1,n);
    // backtrack
    // first add (1,n) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 1;
    stack_interval->j = n;
    stack_interval->energy = get_energy(1,n);
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack_restricted_emodel (cur_interval,tree);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }
	this->structure = structure.substr(1,n);
    return energy;

}

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * Until the changes to fres, I am adding +1 to the ptype closing and Si and Sj's to make them match - Mateo 2024
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t W_final::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n, std::vector<Node> &tree){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[i]][S[j]];
	
    if ((tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j && tree[j].pair == i)) {
				en = vij; // i j
				if (en != INF) {
					if (params->model_details.dangles == 2){
						base_type si1 = (i>1 && S[i-1]!=0) ? S[i-1] : -1;
                		base_type sj1 = (j<n && S[j+1]!=0) ? S[j+1] : -1;
                        en += vrna_E_ext_stem(tt, si1, sj1, params);
					}
                    else{
                        en += vrna_E_ext_stem(tt, -1, -1, params);
					}

                    e = MIN2(e, en);
					
				}

	}
	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
		
        if (((tree[i+1].pair <-1 && tree[j].pair <-1) || (tree[i+1].pair == j)) && tree[i].pair<0) {
            en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

            if (en != INF) {

                base_type si1 = S[i];
                en += vrna_E_ext_stem(tt, si1, -1, params);
            }

            e = MIN2(e,en);
        }
        tt  = pair[S[i]][S[j-1]];
        if (((tree[i].pair <-1 && tree[j-1].pair <-1) || (tree[i].pair == j-1)) && tree[j].pair<0) {
            en = (j-1-i>TURN) ? vij1 : INF; // i j-1
            if (en != INF) {

                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, -1, sj1, params);
            }
            e = MIN2(e,en);

        }
        tt  = pair[S[i+1]][S[j-1]];
        if (((tree[i+1].pair <-1 && tree[j-1].pair <-1) || (tree[i+1].pair == j-1)) && tree[i].pair < 0 && tree[j].pair<0) {
            en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1

            if (en != INF) {

                base_type si1 = S[i];
                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, si1, sj1, params);
            }
            e = MIN2(e,en);
        }
	}
	return e;
}

void W_final::backtrack_restricted_emodel(seq_interval *cur_interval, sparse_tree &tree){
    char type;
	switch (cur_interval->type){
		case LOOP:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			if (i >= j) return;
			f[i].pair = j;
			f[j].pair = i;

			structure[i] = '(';
			structure[j] = ')';		

			type = V->get_type (i,j);
			if(debug) printf("At %d at %d and %c\n",cur_interval->i,cur_interval->j,type);
			switch (type){
				case HAIRP:
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
				{
					f[i].type = INTER;
					f[j].type = INTER;

					int skip = 0;
					if(is_cross_model(i,j)){
						skip = linker_length;
					}
					// detect the other closing pair
					cand_pos_t best_ip=j, best_jp=i;
					energy_t min = INF;
					cand_pos_t max_ip = std::min(j-TURN-2,i+MAXLOOP+1+skip);
					for (cand_pos_t k = i+1; k <= max_ip; ++k)
					{
						if (tree.up[k-1]>=(k-i-1)){
							cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i-skip) - MAXLOOP-2;
							for (cand_pos_t l = j-1; l >= min_l; --l)
							{
								
								if(tree.up[j-1]>=(j-l-1)){
							
									energy_t tmp = emodel_energy_function(i,j,V->compute_int_emodel(i,j,k,l,params_),V->compute_int_emodel(i,j,k,l,params2_));
									if( is_cross_model(i,j) && !(is_cross_model(k,l))) tmp += start_hybrid_penalty;
									if (tmp < min)
									{
										min = tmp;
										best_ip = k;
										best_jp = l;
									}
								}
							}
						}
					}

					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						fprintf(stderr,"NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					bool cross_model = is_cross_model(i,j);
					int best_k = -1, best_row = -1;
					energy_t tmp= INF, min = INF;
					pair_type tt  = pair[S_[j]][S_[i]];
					if(cross_model){
						for (cand_pos_t k = i+2; k <= j-3; k++){
							tmp = V->get_energy_VMprime(i+1,k-1) + V->get_energy_WM(k,j-1) + emodel_energy_function(i,j,E_MLstem(tt, -1, -1, params_),E_MLstem(tt, -1, -1, params2_)) + params_->MLclosing;
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 9;
							}

							tmp = V->get_energy_VMprime(i+2,k-1) + V->get_energy_WM(k,j-1) + emodel_energy_function(i,j,E_MLstem(tt,-1,S_[i+1],params_),E_MLstem(tt,-1,S_[i+1],params2_)) + params_->MLclosing;
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 10;
							}

							tmp = V->get_energy_VMprime(i+1,k-1) + V->get_energy_WM(k,j-2) + emodel_energy_function(i,j,E_MLstem(tt,S_[j-1],-1,params_),E_MLstem(tt,S_[j-1],-1,params2_)) + params_->MLclosing;
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 11;
							}

							tmp = V->get_energy_VMprime(i+2,k-1) + V->get_energy_WM(k,j-2) + emodel_energy_function(i,j,E_MLstem(tt,S_[j-1],S_[i+1],params_),E_MLstem(tt,S_[j-1],S_[i+1],params2_)) + params_->MLclosing;
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 12;
							}
						}
					} else {

					
						for (cand_pos_t k = i+1; k <= j-1; k++){

							if(seq_[k] == 'X') continue;
							

							tmp = V->get_energy_WM (i+1,k-1) + std::min(V->get_energy_WMv(k, j-1),V->get_energy_WMp(k, j-1));
							tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_) + params_->MLclosing,E_MLstem(pair[S_[j]][S_[i]],-1,-1,params2_) + params2_->MLclosing);
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 1;
							}
							// TODO:
							// Hosna, May 1st, 2012
							// do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
							if (tree.tree[i+1].pair <= -1){
								tmp = V->get_energy_WM (i+2,k-1) + std::min(V->get_energy_WMv(k, j-1),V->get_energy_WMp(k, j-1));
								tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase,E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params2_) + params2_->MLclosing + params2_->MLbase);
								if (tmp < min){
									min = tmp;
									best_k = k;
									best_row = 2;
								}
							}
							if (tree.tree[j-1].pair <= -1){
								tmp = V->get_energy_WM (i+1,k-1) + std::min(V->get_energy_WMv(k, j-2),V->get_energy_WMp(k, j-2));
								tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params2_) + params2_->MLclosing + params2_->MLbase);
								if (tmp < min){
									min = tmp;
									best_k = k;
									best_row = 3;
								}
							}
							if (tree.tree[i+1].pair <= -1 && tree.tree[j-1].pair <= -1){
								tmp = V->get_energy_WM (i+2,k-1) + std::min(V->get_energy_WMv(k, j-2),V->get_energy_WMp(k, j-2));
								tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params2_) + params2_->MLclosing + 2*params2_->MLbase);
								
								if (tmp < min){
									min = tmp;
									best_k = k;
									best_row = 4;
								}
							}

							tmp = static_cast<energy_t>((k-i-1)*params_->MLbase + V->get_energy_WMp(k,j-1));
							tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_) + params_->MLclosing,E_MLstem(pair[S_[j]][S_[i]],-1,-1,params2_) + params2_->MLclosing);
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 5;
							}
							// TODO:
							// Hosna, May 1st, 2012
							// do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
							if (tree.tree[i+1].pair <= -1){
								if((k-(i+1)-1) >=0){ 
									tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-1);
									tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase,E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params2_) + params2_->MLclosing + params2_->MLbase);
								}
								if (tmp < min){
									min = tmp;
									best_k = k;
									best_row = 6;
								}
							}
							if (tree.tree[j-1].pair <= -1){
								tmp = static_cast<energy_t>((k-i-1)*params_->MLbase) + V->get_energy_WMp(k,j-2);
								tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params2_) + params2_->MLclosing + params2_->MLbase);
								if (tmp < min){
									min = tmp;
									best_k = k;
									best_row = 7;
								}
							}
							if (tree.tree[i+1].pair <= -1 && tree.tree[j-1].pair <= -1){
								if((k-(i+1)-1) >=0){
									tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-2);
									tmp += emodel_energy_function(i,j,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase,E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params2_) + params2_->MLclosing + 2*params2_->MLbase);
								}
								if (tmp < min){
									min = tmp;
									best_k = k;
									best_row = 8;
								}
							}						
						}
					}
					switch (best_row)
					  {
					  case 1:
						insert_node (i+1, best_k-1, M_WM);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 2:
						insert_node (i+2, best_k-1, M_WM);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 3:
		             	// printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k-1,best_k,j-2);
						insert_node (i+1, best_k-1, M_WM);
						insert_node (best_k, j-2, M_WM);
						break;
					  case 4:
		             	// printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
						insert_node (i+2, best_k-1, M_WM);
						insert_node (best_k, j-2, M_WM);
						break;
					  case 5:
						insert_node (best_k, j-1, M_WM);
						break;
					  case 6:
						insert_node (best_k, j-1, M_WM);
						break;
					  case 7:
						insert_node (best_k, j-2, M_WM);
						break;
					  case 8:
						insert_node (best_k, j-2, M_WM);
						break;
					  case 9:
						insert_node (i+1, best_k-1, M_VMp);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 10:
						insert_node (i+2, best_k-1, M_VMp);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 11:
						insert_node (i+1, best_k-1, M_VMp);
						insert_node (best_k, j-2, M_WM);
						break;
					  case 12:
						insert_node (i+2, best_k-1, M_VMp);
						insert_node (best_k, j-2, M_WM);
						break;
					  }
				}
					break;
			}
		}
			break;
		case FREE:
		{	
			if(debug) printf("At %d and %d and %c\n",cur_interval->i,cur_interval->j,FREE);
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;

			if (j==1) return;

			energy_t min = INF, tmp = INF, acc = INF, energy_kj = INF;
			cand_pos_t best_row = -1, best_k = -1;

			// this case is for j unpaired, so I have to check that.
			if (tree.tree[j].pair <= -1){
				tmp = get_energy(i,j-1);
				if (tmp < min){
					min = tmp;
					best_row = 0;
				}
			}
			for (cand_pos_t k=i; k<j-3; ++k)    // no TURN
			{
				if (seq_[k-1] == 'X' || seq_[j-1] == 'X') continue; // should this be k+1 and j+1?
				// Don't need to make sure k and j don't have to pair with something else
				//  it's INF, done in fold_sequence_restricted
				acc = get_energy(i,k-1);
				energy_kj = V->get_energy(k,j);

				if (energy_kj < INF){	
					if(params_->model_details.dangles == 2){
						base_type sk1 = i>1 ? S_[k-1] : -1;
						base_type sj1 = j<n ? S_[j+1] : -1;
						tmp = energy_kj + emodel_energy_function(k,j,E_ExtLoop(pair[S_[k]][S_[j]],sk1,sj1,params_),E_ExtLoop(pair[S_[k]][S_[j]],sk1,sj1,params2_)) + acc;
					} else 
						tmp = energy_kj + emodel_energy_function(k,j,E_ExtLoop(pair[S_[k]][S_[j]],-1,-1,params_),E_ExtLoop(pair[S_[k]][S_[j]],-1,-1,params2_)) + acc; 
					if (tmp < min){
					min = tmp;
					best_k = k;
					best_row = 1;
					}
					
				}
				if(params_->model_details.dangles ==1){
					if (tree.tree[k].pair <= -1){
						energy_kj = V->get_energy(k+1,j);
						if (energy_kj < INF){
							tmp = energy_kj + emodel_energy_function(k+1,j,E_ExtLoop(pair[S_[k+1]][S_[j]],S_[k],-1,params_),E_ExtLoop(pair[S_[k+1]][S_[j]],S_[k],-1,params2_)) + acc;
							
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 2;
							}
							
						}
					}
					if (tree.tree[j].pair <= -1){
						energy_kj = V->get_energy(k,j-1);
						if (energy_kj < INF){
							tmp = energy_kj + emodel_energy_function(k,j-1,E_ExtLoop(pair[S_[k]][S_[j-1]],-1,S_[j],params_),E_ExtLoop(pair[S_[k]][S_[j-1]],-1,S_[j],params2_)) + acc;
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 3;
							}
						}
					}
					if (tree.tree[k].pair <= -1 && tree.tree[j].pair <= -1){
						energy_kj = V->get_energy(k+1,j-1);
						if (energy_kj < INF){
							tmp = energy_kj + emodel_energy_function(k+1,j-1,E_ExtLoop(pair[S_[k+1]][S_[j-1]],S_[k],S_[j],params_),E_ExtLoop(pair[S_[k+1]][S_[j-1]],S_[k],S_[j],params2_)) + acc;
							if (tmp < min){
								min = tmp;
								best_k = k;
								best_row = 4;
							}
						}
					}
				}
			}
			// Hosna June 30, 2007
			// The following would not take care of when
			// we have some unpaired bases before the start of the WMB
			for (cand_pos_t k=i; k<j-TURN; ++k)
			{
				// Hosna: July 9, 2007
				// We only chop W to W + WMB when the bases before WMB are free
				if (k == 1 || (tree.weakly_closed(i,k-1) && tree.weakly_closed(k,j))){

					acc = get_energy(i,k-1);

					energy_kj = WMB->get_WMB(k,j);

					if (energy_kj < INF)
					{
						tmp = energy_kj + PS_penalty + acc;

						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
							best_k = k;
						}
					}
				}
			}
			switch (best_row)
			{
				case 0:
					insert_node (i, j-1, FREE); break;
				case 1:
					insert_node (best_k, j, LOOP);
					if (best_k-1 > 1)     
						insert_node (i, best_k-1, FREE);
					break;
				case 2:
					insert_node (best_k+1, j, LOOP);
					if (best_k-1 >= 1)
						insert_node (i, best_k-1, FREE);
					break;
				case 3:
					insert_node (best_k, j-1, LOOP);
					if (best_k-1 > 1)
						insert_node (i, best_k-1, FREE);
					break;
				case 4:
					insert_node (best_k+1, j-1, LOOP);
					if (best_k-1 >= 1)
						insert_node (i, best_k-1, FREE);
					break;
				case 5:
					insert_node (best_k, j, P_WMB);
					if (best_k-1 > 1)
						insert_node (i, best_k-1, FREE);
					break;
			}
		}
			break;
		case M_WM:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			if(debug) printf("M_WM at %d and %d\n",i,j);
			energy_t min = INF;
			cand_pos_t best_k = j, best_row = -1;
			  
			for (cand_pos_t k=i; k <= j-TURN-1; k++){	
				energy_t m1 = INF,m2 = INF;
				bool can_pair = tree.up[k-1] >= (k-(i));
				if(can_pair) m1 = emodel_energy_function(i,j,static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMv (k, j),static_cast<energy_t>((k-i)*params2_->MLbase) + V->get_energy_WMv (k, j));
				if (m1 < min){
					min = m1;
					best_k = k;
					best_row = 1;
				}
				if(can_pair) m2 = emodel_energy_function(i,j,static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMp(k, j),static_cast<energy_t>((k-i)*params2_->MLbase) + V->get_energy_WMp(k, j));
				if (m2 < min){
					min = m2;
					best_k = k;
					best_row = 2;
				}
				energy_t m3 = V->get_energy_WM (i, k-1) + V->get_energy_WMv (k, j);
				if (m3 < min){
					min = m3;
					best_k = k;
					best_row = 3;
				}
				energy_t m4 = V->get_energy_WM (i, k-1) + V->get_energy_WMp(k, j);
				if (m4 < min){
					min = m4;
					best_k = k;
					best_row = 4;
				}
			}
			if(tree.tree[j].pair<0){
				if(V->get_energy_WM(i,j-1)< min){
					min = V->get_energy_WM(i,j-1);
					best_row = 5;
				}
			}
			if(seq_[j] == 'X'){
				best_row = 5;
			}
			  
			switch (best_row){
				case 1: insert_node (best_k, j, M_WMv); break;
				case 2: insert_node (best_k, j, M_WMp); break;
				case 3:
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k, j, M_WMv);
				break;
				case 4:
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k, j, M_WMp);
				break;
				case 5:
					insert_node (i,j-1,M_WM); break;
			}
		}
			break;
		case M_WMv:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			if(debug) printf("M_WMv at %d and %d\n",i,j);
			energy_t min = INF;
			cand_pos_t best_row;
			cand_pos_t si = S_[i];
			cand_pos_t sj = S_[j];
			cand_pos_t si1 = (i>1) ? S_[i-1] : -1;
			cand_pos_t sj1 = (j<n) ? S_[j+1] : -1;
			pair_type tt = pair[S_[i]][S_[j]];

			min = V->get_energy(i,j);
			min += ((params_->model_details.dangles == 2) ? emodel_energy_function(i,j,E_MLstem(tt,si1,sj1,params_),E_MLstem(tt,si1,sj1,params2_)) : emodel_energy_function(i,j,E_MLstem(tt,-1,-1,params_),E_MLstem(tt,-1,-1,params2_)));

			best_row = 1;
			if(params_->model_details.dangles == 1){
				if(tree.tree[i].pair<0){
					tt = pair[S_[i+1]][S_[j]];
					energy_t tmp = V->get_energy(i+1,j) + emodel_energy_function(i,j,E_MLstem(tt,si,-1,params_),E_MLstem(tt,si,-1,params2_));
					if(tmp<min){
						min = tmp;
						best_row = 2;
					}
				}
				if(tree.tree[j].pair<0){
					tt = pair[S_[i]][S_[j-1]];
					energy_t tmp = V->get_energy(i,j-1) + emodel_energy_function(i,j,E_MLstem(tt,-1,sj,params_),E_MLstem(tt,-1,sj,params2_));
					if(tmp<min){
						min = tmp;
						best_row = 3;
					}
				}
				if(tree.tree[i].pair<0 && tree.tree[j].pair<0){
					tt = pair[S_[i+1]][S_[j-1]];
					energy_t tmp = V->get_energy(i+1,j-1) + emodel_energy_function(i,j,E_MLstem(tt,si,sj,params_),E_MLstem(tt,si,sj,params2_));
					if(tmp<min){
						min = tmp;
						best_row = 4;
					}
				}
			}
			if(tree.tree[j].pair<0){
				if(V->get_energy_WMv(i,j-1)< min){
					min = V->get_energy_WMv(i,j-1);
					best_row = 5;
				}
			}
			if(seq_[j] == 'X'){
				best_row = 5;
			}
			switch (best_row){
				case 1: insert_node (i, j, LOOP); break;
				case 2: insert_node (i+1, j, LOOP); break;
				case 3: insert_node (i, j-1, LOOP); break;
				case 4: insert_node (i+1, j-1, LOOP); break;
				case 5: insert_node (i, j-1, M_WMv); break;
			}
		}
		break;
		case M_WMp:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			if(debug) printf("M_WMp at %d and %d\n",i,j);
			energy_t min = INF;
			int best_row = -1;

			min = WMB->get_WMB(i,j) + PSM_penalty + b_penalty;
			best_row = 1;
			if(tree.tree[j].pair<0){
				if(V->get_energy_WMp(i,j-1)< min){
					min = V->get_energy_WMp(i,j-1);
					best_row = 2;
				}
			}
			if(seq_[j] == 'X'){
				best_row = 2;
			}
			switch (best_row){
				case 1: insert_node (i, j, P_WMB); break;
				case 2: insert_node (i, j-1, M_WMp); break;
			}
		}
		break;
		case M_VMp:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			if(debug) printf("M_VMp at %d and %d\n",i,j);
			energy_t min = INF;
			cand_pos_t best_k = j, best_row = -1;
			for (cand_pos_t k = i+4; k <= j-3; ++k){
				cand_pos_t sk = S_[k];
				cand_pos_t sj = S_[j];
				cand_pos_t sk1 = (i>1) ? S_[k-1] : -1;
				cand_pos_t sj1 = (j<n) ? S_[j+1] : -1;
				pair_type tt = pair[S_[k]][S_[j]];

				energy_t m1 = V->get_energy_WM(i,k-1) + V->get_energy(k,j);
				m1 += (params_->model_details.dangles == 2) ? emodel_energy_function(k,j,E_MLstem(tt,sk1,sj1,params_),E_MLstem(tt,sk1,sj1,params2_)) : emodel_energy_function(k,j,E_MLstem(tt,-1,-1,params_),E_MLstem(tt,-1,-1,params2_));
				if(m1 < min){
					min = m1;
					best_k = k;
					best_row = 1;
				}

				if(params_->model_details.dangles == 1){
					if(tree.tree[k].pair<0){
						tt = pair[S_[k+1]][S_[j]];
						energy_t tmp = V->get_energy_WM(i,k-1) + V->get_energy(k+1,j) + emodel_energy_function(k+1,j,E_MLstem(tt,sk,-1,params_),E_MLstem(tt,sk,-1,params2_));
						if(tmp<min){
							min = tmp;
							best_row = 2;
						}
					}
					if(tree.tree[j].pair<0){
						tt = pair[S_[k]][S_[j-1]];
						energy_t tmp = V->get_energy_WM(i,k-1) + V->get_energy(k,j-1) + emodel_energy_function(k,j-1,E_MLstem(tt,-1,sj,params_),E_MLstem(tt,-1,sj,params2_));
						if(tmp<min){
							min = tmp;
							best_row = 3;
						}
					}
					if(tree.tree[k].pair<0 && tree.tree[j].pair<0){
						tt = pair[S_[k+1]][S_[j-1]];
						energy_t tmp = V->get_energy_WM(i,k-1) + V->get_energy(k+1,j-1) + emodel_energy_function(k+1,j-1,E_MLstem(tt,sk,sj,params_),E_MLstem(tt,sk,sj,params2_));
						if(tmp<min){
							min = tmp;
							best_row = 4;
						}
					}
				}
				energy_t m2 = std::min(m2,V->get_energy_WM(i,k-1) + WMB->get_WMB(k,j) + PSM_penalty + b_penalty);
				if(m2 < min){
					min = m2;
					best_k = k;
					best_row = 5;
				}
			}
			energy_t m3 = V->get_energy_VMprime(i,j-1) + params2_->MLbase;
			if (m3 < min){
				min = m3;
				best_row = 6;
			}


			switch (best_row){
				case 1: 
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k, j, LOOP);
					break;
				case 2: 
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k+1, j, LOOP);
					break;
				case 3: 
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k, j-1, LOOP);
					break;
				case 4: 
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k+1, j-1, LOOP);
					break;
				case 5: 
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k, j, P_WMB);
					break;
				case 6: insert_node (i, j-1, M_VMp); break;
			}
		}
		break;
		case P_WMB:
		case P_WMBP:
		case P_WMBW:
		case P_VP:
		case P_VPR:
		case P_VPL:
		case P_WI:
		case P_BE:
		case P_WIP:
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

void W_final::insert_node (int i, int j, char type)
  // insert at the beginning
{
    seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}


//Mateo 13 Sept 2023
//return number of bases in between the two inclusive index
int distance(int left, int right){
    return (right-left-1);
}

//Mateo 13 Sept 2023
//given a initial hotspot which is a hairpin loop, keep trying to add a arc to form a larger stack
void expand_hotspot(s_energy_matrix *V, Hotspot &hotspot, int n){
    //calculation for the hairpin that is already in there
    V->compute_hotspot_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index(),0);


    //try to expand by adding a arc right beside the current out most arc
    while(hotspot.get_left_outer_index()-1 >= 1 && hotspot.get_right_outer_index()+1 <= n){
		base_type sim1 = V->S_[hotspot.get_left_outer_index()-1];
		base_type sjp1 = V->S_[hotspot.get_right_outer_index()+1];
		pair_type ptype_closing = pair[sim1][sjp1];
        if(ptype_closing>0){
            hotspot.move_left_outer_index();
            hotspot.move_right_outer_index();
            hotspot.increment_size();
            V->compute_hotspot_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index(),1);
        }else{
            break;
        }
    }
	base_type i = hotspot.get_left_outer_index();
	base_type j = hotspot.get_right_outer_index();
	pair_type tt = pair[V->S_[i]][V->S_[j]];
	base_type si1 = i>1 ? V->S_[i-1] : -1;
	base_type sj1 = j<=n ? V->S_[j+1] : -1;
	energy_t dangle_penalty = vrna_E_ext_stem(tt, si1, sj1, V->params_);


    double energy = V->get_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index());


    energy = (energy + dangle_penalty) / 100;

    hotspot.set_energy(energy);
    return;
}

//Mateo 13 Sept 2023
//look for every possible hairpin loop, and try to add a arc to form a larger stack with at least min_stack_size bases
void get_hotspots(std::string seq,std::vector<Hotspot> &hotspot_list,int max_hotspot, vrna_param_s *params){
    
	int n = seq.length();
	s_energy_matrix *V;
	make_pair_matrix();
	short *S_ = encode_sequence(seq.c_str(),0);
	short *S1_ = encode_sequence(seq.c_str(),1);
	V = new s_energy_matrix (seq,n,S_,S1_,params);
    int min_bp_distance = 3;
    int min_stack_size = 3; //the hotspot must be a stack of size >= 3
    // Hotspot current_hotspot;
    //start at min_stack_size-1 and go outward to try to add more arcs to form bigger stack because we cannot expand more than min_stack_size from there anyway
    for(int i = min_stack_size; i <= n; i++){
        for(int j = i; j <= n; j++){
			int ptype_closing = pair[V->S_[i]][V->S_[j]];
            if(ptype_closing>0 && distance(i,j) >= min_bp_distance){
                // current_hotspot = new Hotspot(i,j,nb_nucleotides);
				
                Hotspot current_hotspot(i,j,n);

                expand_hotspot(V,current_hotspot,n);


                if(current_hotspot.get_size() < min_stack_size || current_hotspot.is_invalid_energy()){

                }else{
                    
                    current_hotspot.set_structure();
                    hotspot_list.push_back(current_hotspot);

                }
            }
        }
    }
    //make sure we only keep top 20 hotspot with lowest energy
    std::sort(hotspot_list.begin(), hotspot_list.end(),compare_hotspot_ptr);
    while((int) hotspot_list.size() > max_hotspot){
        hotspot_list.pop_back();
    }

    //if no hotspot found, add all _ as restricted
    if((int) hotspot_list.size() == 0){
        Hotspot hotspot(1,n,n+1);
        hotspot.set_default_structure();
        hotspot_list.push_back(hotspot);
    }
	delete V;
	free(S_);
	free(S1_);

    return;
}

bool compare_hotspot_ptr(Hotspot &a, Hotspot &b) { 
    return (a.get_energy() < b.get_energy()); 
}