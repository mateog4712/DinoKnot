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
	if(!exists(base_pair_file)) return;
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

	int model_1_Type = type_1;
	int model_2_Type = type_2;

	std::string inputFile = args_info.input_given ? input_file : "";

	std::string inputSequence1;
	std::string inputSequence2;
	std::string inputStructure1;
	std::string inputStructure2;
	if(args_info.input_given) get_input(inputFile,inputSequence1,inputSequence2,inputStructure1,inputStructure2);

	inputSequence1 = (args_info.sequence1_given) ? sequence_1 : "";
	inputSequence2 = (args_info.sequence2_given) ? sequence_2 : "";
	if(model_1_Type == 0) seqtoRNA(inputSequence1);
	if(model_2_Type == 0) seqtoRNA(inputSequence2);

	validateSequence(inputSequence1);
	validateSequence(inputSequence2);
	cand_pos_t n1 = inputSequence1.length();
	cand_pos_t n2 = inputSequence2.length();

	inputStructure1 = (args_info.structure1_given) ? structure_1 : std::string(n1,'.');
	inputStructure2 = (args_info.structure2_given) ? structure_2 : std::string(n2,'.');
	if(inputStructure1 != "") validateStructure(inputSequence1,inputStructure1);
	if(inputStructure2 != "") validateStructure(inputSequence2,inputStructure2);
	
	

				
	std::string outputDir = args_info.dir_given ? output_dir : "";
	std::string outputFile = args_info.output_given ? output_file : "";
	std::string hotspotDir = args_info.h_only_given ? hotspot_dir : "";
	std::string varnaFile = args_info.varna_given ? varna : "";
	std::string basepairFile = args_info.basePairFile_given ? base_pair_file : "";
	
	std::vector<std::tuple<cand_pos_t,cand_pos_t> > pairs;
	load_base_pairs(base_pair_file,pairs);

	int max_hotspot = args_info.h_num_given ? hotspot_num : 20;
	int number_of_suboptimal_structure = args_info.subopt_given ? subopt : 100;

	bool hotspot_only = args_info.h_only_given;

	bool micro = args_info.micro_given;

	bool hard = args_info.hard_given;

	int dangle = args_info.dangles_given ? dangle_model : 1;

	start_hybrid_penalty = args_info.pen_given ? hybrid_pen : get_START_HYBRID_PENALTY(model_1_Type,model_2_Type);

	linker_pos = inputSequence1.length()+1;

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
	if(args_info.parameter2_given){
		std::string file = parameter2;
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
	params1->model_details.dangles = dangle;
	params2->model_details.dangles = dangle;
//--------------------------------------------------------------------------------------------------------------------------
	if(micro) args_info.structure1_given = true;
	if(!pairs.empty()){
		args_info.structure1_given = true; 
		args_info.structure2_given = true;
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
	
	if(args_info.structure1_given){
		Hotspot hotspot(1,n1,n1+1);
		hotspot.set_structure(inputStructure1);
		hotspot_list1.push_back(hotspot);
	}else {
		get_hotspots(inputSequence1, hotspot_list1,max_hotspot,params1);
	}

	if(args_info.structure2_given){
		Hotspot hotspot(1,n2,n2+1);
		hotspot.set_structure(inputStructure2);
		hotspot_list2.push_back(hotspot);
		
	} else {
		get_hotspots(inputSequence2, hotspot_list2,max_hotspot,params2);
	}

	cmdline_parser_free(&args_info);

	// Generate full sequence after reversal of sequence 1 has occurred
	std::string seq = inputSequence1 + "XXXXX" + inputSequence2;

	if(hotspot_only){
		std::ofstream out(hotspotDir.c_str());
		if(!exists(hotspot_dir)){
			std::cout << "Input File does not exist!" << std::endl;
			exit (EXIT_FAILURE);
    	}
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

		if(varnaFile != "" && exists(varnaFile)){
			for(cand_pos_t i = 0; i < number_of_output; ++i){
				std::string command = "java -cp " +  varnaFile +  " fr.orsay.lri.varna.applications.VARNAcmd -algorithm line -resolution 15.0 -basesStyle1 \"fill=##0000FF\" -basesStyle2 \"fill=#0000FF\" -basesStyle3 \"fill=#FFFF00\" -applyBasesStyle1on \"1-" + std::to_string(linker_pos-1) + "\" -applyBasesStyle2on \"" +  std::to_string(linker_pos) + "-" +  std::to_string(linker_pos+linker_length) + "\" -applyBasesStyle3on \"" +  std::to_string(linker_pos+linker_length+1) + "-" +  std::to_string(n) + "\" -sequenceDBN \"" + seq + "\" -structureDBN \"" + result_list[i].get_final_structure() + "\"" + " -o \"varna/file" + std::to_string(i) + ".png\"";
				system(command.c_str());
			}
		}

		//Mateo 7/19/2023
		//output to file
		if(outputFile != ""){
			std::ofstream out(output_file);
			if(!exists(output_file)){
				std::cout << "file is not valid" << std::endl;
				exit(EXIT_FAILURE);
			}

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
		else if(outputDir != ""){
			// Mateo 2023
			if(exists(outputDir)){
				if(outputDir[outputDir.length()] != '/') outputDir += '/';
				std::string path_to_file = outputDir + "output_" + std::to_string(0) + ".txt";
				std::ofstream out(path_to_file);
				out << "Seq:          " << seq << std::endl;
				out << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
				out << "Result_" << 0 << ":     " << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;  
				out.close();
				for (int i=1; i < number_of_output; ++i) {
					if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
					std::string path_to_file = outputDir + "output_" + std::to_string(i) + ".txt";
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
	return 0;
}