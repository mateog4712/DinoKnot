#ifndef ITERATIVE_HFOLD_H_
#define ITERATIVE_HFOLD_H_

#include "W_final.hh"
#include <string>

int is_invalid_restriction(char* restricted_structure, char* current_structure);

std::string Iterative_HFold(std::string seq,std::string res, double &en);

std::string Iterative_HFold_interacting (std::string seq,std::string res, double &en, vrna_param_s *params1, vrna_param_s *params2);



#endif