#include "common.hh"
#include "W_final.hh"
#include "h_externs.hh"




//kevin 10 Aug 2017
//return 1 if i and j goes across linker
int is_cross_model(cand_pos_t i, cand_pos_t j){
	if((i < linker_pos) && (j > linker_pos+linker_length-1) ){
		return 1;
	}
	return 0;
}

//kevin 31 july Used to calculate the average energy of two models
int emodel_energy_function (cand_pos_t i, cand_pos_t j, energy_t e1, energy_t e2){

    energy_t energy = INF;
    int index_of_last_linker_position = linker_pos+linker_length-1;

    if(j <= linker_pos){ //left of linker
        energy = e1;
    }else if((i >= index_of_last_linker_position)){ //right of linker
        energy = e2;
    }else if((i < linker_pos) && (j > index_of_last_linker_position)){ //cross linker
        energy = (e1 + e2)/2;
    }else { //when i and j both X
        energy = INF;

    }

    return  energy;
}