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
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <getopt.h>
#include <tuple>

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

void validateStructure(std::string &seq, std::string &structure){
	int n = structure.length();
	std::vector<int> pairs;
	for(int j = 0; j<n;++j){
		if(structure[j] == '(') pairs.push_back(j);
		if(structure[j] == ')'){
			if(pairs.empty()){
				std::cout << "Incorrect input: More left parentheses than right" << std::endl;
				exit(0);
			}
			else {
				int i = pairs.back();
				pairs.pop_back();
				if(seq[i] == 'A' && seq[j] == 'U'){}
				else if (seq[i] == 'C' && seq[j] == 'G'){}
				else if ((seq[i] == 'G' && seq[j] == 'C') || (seq[i] == 'G' && seq[j] == 'U')){}
				else if ((seq[i] == 'U' && seq[j] == 'G') || (seq[i] == 'U' && seq[j] == 'A')){}
				else{
					std::cout << "Incorrect input: " << seq[i] << " does not pair with " << seq[j] << std::endl;
					exit(0);
				}
			}
		}
	}
	if(!pairs.empty()){
		std::cout << "Incorrect input: More left parentheses than right" << std::endl;
		exit(0);
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

void seqtoRNA(std::string &sequence){
    for (char &c : sequence) {
      	if (c == 'T') c = 'U';
    }
}

void can_pair(char a, char b){
	if(!((a == 'A' && b == 'U') || (a == 'C' && b == 'G') || (a == 'G' && b == 'C') || (a == 'G' && b == 'U') || (a == 'U' && b == 'A') || (a == 'U' && b== 'G'))){
		printf("Error, not a valid pair: %c and %c\n",a,b);
		exit(0);
	}
}

void load_base_pairs(std::string file, std::vector< std::tuple<cand_pos_t,cand_pos_t> > &pairs){
	if(!exists(file)) return;
	std::ifstream in (file);
	std::string str;
	while(getline(in,str)){
		std::istringstream iss(str);
		cand_pos_t index1; 
		iss >> index1;
		cand_pos_t index2;
		iss >> index2;
		pairs.push_back(std::make_tuple(index1,index2));
	}
}

int main (int argc, char *argv[]) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	int model_1_Type = args_info.t1_arg;
	int model_2_Type = args_info.t2_arg;

	std::string inputSequence1;
	std::string inputSequence2;
	std::string inputStructure1;
	std::string inputStructure2;
	if(args_info.input_file_given) get_input(args_info.input_file_arg,inputSequence1,inputSequence2,inputStructure1,inputStructure2);

	inputSequence1 = (args_info.s1_given) ? args_info.s1_arg : "";
	inputSequence2 = (args_info.s2_given) ? args_info.s2_arg : "";
	if(model_1_Type == 0) seqtoRNA(inputSequence1);
	if(model_2_Type == 0) seqtoRNA(inputSequence2);

	validateSequence(inputSequence1);
	validateSequence(inputSequence2);
	cand_pos_t n1 = inputSequence1.length();
	cand_pos_t n2 = inputSequence2.length();

	inputStructure1 = (args_info.r1_given) ? args_info.r1_arg : std::string(n1,'.');
	inputStructure2 = (args_info.r2_given) ? args_info.r2_arg : std::string(n2,'.');
	if(inputStructure1 != "") validateStructure(inputSequence1,inputStructure1);
	if(inputStructure2 != "") validateStructure(inputSequence2,inputStructure2);

	int max_hotspot = args_info.hotspot_num_given ? args_info.hotspot_num_arg : 20;
	int number_of_suboptimal_structure = args_info.opt_given ? args_info.opt_arg : std::pow(max_hotspot,2);


	bool micro = args_info.micro_given;

	bool hard = args_info.hard_given;

	start_hybrid_penalty = args_info.pen_given ? args_info.pen_arg : lrint(get_START_HYBRID_PENALTY(model_1_Type,model_2_Type));

	linker_pos = inputSequence1.length()+1;
	linker_pos_right = inputSequence1.length()+5;
	//                                                   Energy Model Portion 
//-----------------------------------------------------------------------------------------------------------
	vrna_param_s *params1;
	vrna_param_s *params2;
	if(args_info.paramFile1_given){
		std::string file = args_info.paramFile1_arg;
		if(file!=""){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
		}
		params1 = scale_parameters();

	}
	else{
		if(model_1_Type == 0){
			std::string file = "params/rna_DirksPierce09.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params1 = scale_parameters();
		}
		else{
			std::string file = "params/dna_Matthews04.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params1 = scale_parameters();
		}
	}
	if(args_info.paramFile2_given){
		std::string file = args_info.paramFile2_arg;
		if(file!=""){
			vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
		}
		params2 = scale_parameters();
	}
	else{
		if(model_2_Type == 0){
			std::string file = "params/rna_DirksPierce09.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params2 = scale_parameters();
		}
		else{
			std::string file = "params/dna_Matthews04.par";
			if(file!=""){
				vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
			}
			params2 = scale_parameters();
		}
	}
	params1->model_details.dangles = args_info.dangles_arg;
	params2->model_details.dangles = args_info.dangles_arg;
//--------------------------------------------------------------------------------------------------------------------------
	if(micro) args_info.r1_given = true;

	std::vector<std::tuple<cand_pos_t,cand_pos_t> > pairs;
	if(args_info.basePairFile_given) load_base_pairs(args_info.basePairFile_arg,pairs);
	if(!pairs.empty()){
		args_info.r1_given = true; 
		args_info.r2_given = true;
		inputStructure1 = std::string(n1,'.');
		inputStructure2 = std::string(n2,'.');
		int npairs = pairs.size();
		for(cand_pos_t i = 0; i<npairs; ++i){
			cand_pos_t k = std::get<0>(pairs[i]);
			cand_pos_t l = std::get<1>(pairs[i]);
			inputStructure1[(k-1)] = '(';
			inputStructure2[(l-1)] = ')';
			can_pair(inputSequence1[(k-1)], inputSequence2[(l-1)]);
		}
				
	}
	std::vector<Hotspot> hotspot_list1;
	std::vector<Hotspot> hotspot_list2;
	
	if(args_info.r1_given){
		Hotspot hotspot(1,n1,n1+1);
		hotspot.set_structure(inputStructure1);
		hotspot_list1.push_back(hotspot);
	}else {
		get_hotspots(inputSequence1, hotspot_list1,max_hotspot,params1);
	}

	if(args_info.r2_given){
		Hotspot hotspot(1,n2,n2+1);
		hotspot.set_structure(inputStructure2);
		hotspot_list2.push_back(hotspot);
		
	} else {
		get_hotspots(inputSequence2, hotspot_list2,max_hotspot,params2);
	}

	// Generate full sequence after reversal of sequence 1 has occurred
	std::string seq = inputSequence1 + "XXXXX" + inputSequence2;

	if(args_info.hotspot_only_given){
		if(!exists(args_info.hotspot_only_arg)){
			std::cout << "Input File does not exist!" << std::endl;
			exit (EXIT_FAILURE);
    	}
		std::ofstream out(args_info.hotspot_only_arg);
		cand_pos_t size1 = hotspot_list1.size();
		cand_pos_t size2 = hotspot_list2.size();
		for(cand_pos_t i =0; i < size1; i++){
			out << "Seq1_hotspot_" << i << ": " << hotspot_list1[i].get_structure() << "(" << hotspot_list1[i].get_energy() << ")" << std::endl;
		}
		out << "---------------" << std::endl;
		for(cand_pos_t j = 0; j < size2; j++){
			out << "Seq2_hotspot_" << j << ": " << hotspot_list2[j].get_structure() << "(" << hotspot_list2[j].get_energy() << ")" << std::endl;
		}
		out.close();

		free(params1);
		free(params2);
	}
	else {

		std::vector<Result> result_list;
		
		cand_pos_t n = seq.length();
		cand_pos_t size1 = hotspot_list1.size();
		cand_pos_t size2 = hotspot_list2.size();

		for(int i =0; i < size1; i++){
			for(int j = 0; j < size2; j++){
				
				double final_energy = 0;
				int method_chosen = 1;
				std::string restricted = hotspot_list1[i].get_structure() + "xxxxx" + hotspot_list2[j].get_structure();
				
				std::string structure = Iterative_HFold_interacting(seq,restricted,final_energy,params1,params2,method_chosen,hard);

				Result result(seq,restricted,hotspot_list1[i].get_energy()+hotspot_list2[i].get_energy(),structure,final_energy,method_chosen);
				result_list.push_back(result);
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

		if(args_info.varna_given && exists(args_info.varna_arg)){
			std::string varna = args_info.varna_arg;
			for(cand_pos_t i = 0; i < number_of_output; ++i){
				std::string command = "java -cp " +  varna +  " fr.orsay.lri.varna.applications.VARNAcmd -algorithm line -resolution 15.0 -basesStyle1 \"fill=##0000FF\" -basesStyle2 \"fill=#0000FF\" -basesStyle3 \"fill=#FFFF00\" -applyBasesStyle1on \"1-" + std::to_string(linker_pos-1) + "\" -applyBasesStyle2on \"" +  std::to_string(linker_pos) + "-" +  std::to_string(linker_pos+linker_length) + "\" -applyBasesStyle3on \"" +  std::to_string(linker_pos+linker_length+1) + "-" +  std::to_string(n) + "\" -sequenceDBN \"" + seq + "\" -structureDBN \"" + result_list[i].get_final_structure() + "\"" + " -o \"varna/file" + std::to_string(i) + ".png\"";
				system(command.c_str());
			}
		}

		//Mateo 7/19/2023
		//output to file
		if(args_info.output_file_given){
			if(!exists(args_info.output_file_arg)){
				std::cout << "file is not valid" << std::endl;
				exit(EXIT_FAILURE);
			}
			std::ofstream out(args_info.output_file_arg);

			out << "Seq:          " << seq << std::endl;
			out << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
			out << "Result_" << 0 << ":     " << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
			for (int i=1; i < number_of_output; i++) {
				if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
				out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
				out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
			}
			out.close();
		}
		else if(args_info.dir_given){
			// Mateo 2023
			if(exists(args_info.dir_arg)){
				std::string dir = args_info.dir_arg;
				if(dir[dir.length()] != '/') dir += '/';
				std::string path_to_file = dir + "output_" + std::to_string(0) + ".txt";
				std::ofstream out(path_to_file);
				out << "Seq:          " << seq << std::endl;
				out << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
				out << "Result_" << 0 << ":     " << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;  
				out.close();
				for (int i=1; i < number_of_output; ++i) {
					if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
					std::string path_to_file = dir + "output_" + std::to_string(i) + ".txt";
					std::ofstream out(path_to_file);
					out << "Seq:          " << seq << std::endl;
					out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
					out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;  
					out.close();
				}
			}
			else{
				std::cout << "Not a valid output directory" << std::endl;
				exit(EXIT_FAILURE);
			}
		} else{
			// Mateo 2023
			std::cout << "Seq:          " << seq << std::endl;
			std::cout << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
			std::cout << "Result_" << 0 << ":     " << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
			for (int i=1; i < number_of_output; i++) {
				if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
				std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
				std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
			}
		}
	}
	cmdline_parser_free(&args_info);
	return 0;
}