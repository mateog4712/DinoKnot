#ifndef ITERATIVE_HFOLD_H_
#define ITERATIVE_HFOLD_H_

#include "W_final.hh"
#include <string>

std::string Iterative_HFold_interacting (std::string seq,std::string res, double &en, vrna_param_s *params1, vrna_param_s *params2, int &method_chosen, bool hard);



#endif