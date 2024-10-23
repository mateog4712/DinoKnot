// Dinoknot files
#include "dinoknot.hh"
#include "iterative-hfold.hh"
#include "hotspot.hh"
#include "Result.hh"
#include "cmdline.hh"
#include "W_final.hh"
#include "h_globals.hh"
// a simple driver for HFold
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <getopt.h>

#define RNA 0
#define DNA 1
#define PMO 2 

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence1 or sequence2 is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
		std::cout << "Sequence contains character " << c << " that is not G,C,A,U, or T." << std::endl;
		exit(EXIT_FAILURE);
    }
  }
}

void validateStructure(std::string sequence, std::string structure){
	if(structure.length() != sequence.length()){
		std::cout << " The length of the sequence and corresponding structure must have the same length" << std::endl;
		exit(EXIT_FAILURE);
	}

	//check if any characters are not ._x()
	for(char c : structure) {
		if (!(c == '.' || c == '_' || c == 'x' || c == '(' || c == ')')){
			std::cout << "Structure must only contain ._(): " << c << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

void get_input(std::string file, std::string &sequence1, std::string &sequence2, std::string &structure1, std::string &structure2 ){
	if(!exists(file)){
		std::cout << "Input file does not exist" << std::endl;
		exit(EXIT_FAILURE);
	}
	// Partial IUPAC notation
	std::string bases = "ACGTUWSMKRY";
	std::string s = "._()";
	std::ifstream in(file.c_str());
	bool secondseq = false;
	bool secondstruct = false;
	std::string str;
	while(getline(in,str)){
		if(bases.find(str[0])+1){
			if(secondseq) sequence2 = str;
			else {
				sequence1 = str; 
				secondseq = true;
			}
		}
		if(s.find(str[0])+1){
			
			if(secondstruct) structure2 = str;
			else {
				structure1 = str; 
				secondstruct = true;
			}
		}
	}

	in.close();
}

double get_START_HYBRID_PENALTY(int type1, int type2){
	if(type1 == type2){ //if both model are the same
		if(type1 == RNA){ //if both are RNA
			return 22.96551130344778;
		}else if(type1 == DNA){ //if both are DNA
			return 34.14979695798525;
		}else if(type1 == PMO){
			fprintf(stderr, "ERROR: model cannot be both PMO\n");
			exit(1);
		}
	}
	return 166.0; //when 2 different model old: 58.4511432
}

int main (int argc, char *argv[]) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	int model_1_Type = type_1 <= 2 ? type_1 : 0;
	int model_2_Type = type_2 <= 2 ? type_2 : 0;
	std::string inputSequence1 = "";
	std::string inputSequence2 = "";
	std::string inputStructure1 = "";
	std::string inputStructure2 = "";

	std::string inputFile = args_info.input_given ? input_file : "";
	if(args_info.input_given) get_input(inputFile,inputSequence1,inputSequence2,inputStructure1,inputStructure2);

	if(args_info.sequence1_given) inputSequence1 =  sequence_1;
	if(args_info.sequence2_given) inputSequence2 =  sequence_2;
	if(args_info.structure1_given) inputStructure1 = structure_1;
	if(args_info.structure2_given) inputStructure2 = structure_2;
	validateSequence(inputSequence1);
	validateSequence(inputSequence2);
	if(inputStructure1 != "") validateStructure(inputSequence1,inputStructure1);
	if(inputStructure2 != "") validateStructure(inputSequence2,inputStructure2);
	
	std::string seq = inputSequence1 + "XXXXX" + inputSequence2;

				
	std::string outputDir = args_info.dir_given ? output_dir : "";
	std::string outputFile = args_info.output_given ? output_file : "";
	std::string hotspotDir = args_info.h_only_given ? hotspot_dir : "";

	int max_hotspot = args_info.h_num_given ? hotspot_num : 20;
	int number_of_suboptimal_structure = args_info.subopt_given ? subopt : 400;

	bool hotspot_only = args_info.h_only_given;

	int dangle = args_info.dangles_given ? dangle_model : 2;

	start_hybrid_penalty = args_info.pen_given ? hybrid_pen : get_START_HYBRID_PENALTY(model_1_Type,model_2_Type);

	linker_pos = inputSequence1.length();

	//                                                   Energy Model Portion 
//-----------------------------------------------------------------------------------------------------------
	vrna_param_s *params1;
	vrna_param_s *params2;
	if(args_info.parameter1_given){
		std::string file = parameter1;
		if(file!=""){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
		}
		params1 = scale_parameters();

	}
	else{
		if(type_1 == 0){
			std::string file = "src/params/rna_turner2004.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params1 = scale_parameters();
		}
		else{
			std::string file = "src/params/dna_matthews2004.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params1 = scale_parameters();
		}
	}
	if(args_info.parameter2_given){
		std::string file = parameter2;
		if(file!=""){
			vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
		}
		params2 = scale_parameters();
	}
	else{
		if(type_2 == 0){
			std::string file = "src/params/rna_turner2004.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params2 = scale_parameters();
		}
		else{
			std::string file = "src/params/dna_matthews2004.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params2 = scale_parameters();
		}
	}
	params1->model_details.dangles = dangle;
	params2->model_details.dangles = dangle;
	cmdline_parser_free(&args_info);
//--------------------------------------------------------------------------------------------------------------------------

	std::vector<Hotspot> hotspot_list1;
	std::vector<Hotspot> hotspot_list2;
	// Hotspot* hotspot;
	if(inputStructure1 != ""){
		Hotspot hotspot(1,inputStructure1.length(),inputStructure1.length()+1);
		hotspot.set_structure(inputStructure1);
		hotspot_list1.push_back(hotspot);
	}else {
		get_hotspots(inputSequence1, hotspot_list1,max_hotspot,params1);
	}

	if(inputStructure2 != ""){
		Hotspot hotspot(1,inputStructure2.length(),inputStructure1.length()+1);
		hotspot.set_structure(inputStructure2);
		hotspot_list2.push_back(hotspot);
		
	} else {
		get_hotspots(inputSequence2, hotspot_list2,max_hotspot,params2);
	}

	if(hotspot_only){
		std::ofstream out(hotspotDir.c_str());
		if(!exists(hotspot_dir)){
			std::cout << "Input File does not exist!" << std::endl;
			exit (EXIT_FAILURE);
    	}
		for(int i =0; i < hotspot_list1.size(); i++){
			out << "Seq1_hotspot_" << i << ": " << hotspot_list1[i].get_structure() << "(" << hotspot_list1[i].get_energy() << ")" << std::endl;
		}
		out << "---------------" << std::endl;
		for(int j = 0; j < hotspot_list2.size(); j++){
			out << "Seq2_hotspot_" << j << ": " << hotspot_list2[j].get_structure() << "(" << hotspot_list2[j].get_energy() << ")" << std::endl;
		}
		out.close();

		free(params1);
		free(params2);
	}
	else {

		std::vector<Result> result_list;
		
		cand_pos_t n = seq.length();
		for(int i =0; i < hotspot_list1.size(); i++){
			for(int j = 0; j < hotspot_list2.size(); j++){
				
				double final_energy = 0;
				int method_chosen = 1;
				std::string restricted = hotspot_list1[i].get_structure() + "xxxxx" + hotspot_list2[j].get_structure();
				
				std::string structure = Iterative_HFold_interacting(seq,restricted,final_energy,params1,params2,method_chosen);
				
				Result result(seq,restricted,hotspot_list1[i].get_energy()+hotspot_list2[i].get_energy(),structure,final_energy,method_chosen);
				result_list.push_back(result);
				std::cout << restricted << std::endl;
				std::cout << structure << " (" << final_energy << ")" << std::endl;
			}
		}
	
		Result::Result_comp result_comp;
		std::sort(result_list.begin(), result_list.end(),result_comp );
		free(params1);
		free(params2);


	// 	//kevin 5 oct 2017
		int number_of_output = 1;
	// 	// //printf("number_of_suboptimal_structure: %d\n",number_of_suboptimal_structure);
		if(number_of_suboptimal_structure != 1){
			number_of_output = std::min( (int) result_list.size(),number_of_suboptimal_structure);
		}

	// 	//Mateo 7/19/2023
	// 	//output to file
	// 	if(outputFile != ""){
	// 		std::ofstream out(output_file);
	// 		if(!exists(output_file)){
	// 			std::cout << "file is not valid" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		out << "Seq:          " << seq << std::endl;
	// 		for (int i=0; i < number_of_output; i++) {
	// 			out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
	// 			out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
			
	// 		}
	// 		out.close();
	// 	}
	// 	else if(outputDir != ""){
	// 		// Mateo 2023
	// 		if(exists(outputDir)){
	// 			if(outputDir[outputDir.length()] != '/') outputDir += '/';
	// 			for (int i=0; i < number_of_output; ++i) {
	// 				std::string path_to_file = outputDir + "output_" + std::to_string(i) + ".txt";
	// 				std::ofstream out(path_to_file);
	// 				out << "Seq:          " << seq << std::endl;
	// 				out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
	// 				out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;  
	// 				out.close();
	// 			}
	// 		}
	// 		else{
	// 			std::cout << "Not a valid output directory" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}
	// 	} else{
	// 		// Mateo 2023
	// 		std::cout << "Seq:          " << seq << std::endl;
	// 		for (int i=0; i < number_of_output; i++) {
	// 			std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
	// 			std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
	// 		}
	// 	}
	}
	return 0;
}