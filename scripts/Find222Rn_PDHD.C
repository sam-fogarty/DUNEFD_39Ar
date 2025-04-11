
#include "tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>


// get list of filepaths from a text file
std::vector<std::string> readFilePaths(const std::string& fileName) {
    std::ifstream file(fileName);
    std::vector<std::string> filePaths;
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << fileName << std::endl;
        return filePaths;
    }

    while (std::getline(file, line)) {
        filePaths.push_back(line);
    }

    file.close();
    return filePaths;
}

void Find222Rn_PDHD(){
    // Try to find 222Rn events from low energy charge clusters
 
    // get list of ART ROOT files to analyze
    string filelist = "28850_filepaths.txt";
    std::vector<std::string> filenames = readFilePaths(filelist);

    // get list of ROOT ntuple files with TPC clusters
    string run = "028850";
    string singlehit_runlist="/pnfs/dune/scratch/users/lavaut/03684/1/run28850.list";
    bool verbose = true;
    // constants
    const double dt    = 16e-3; // DTS sampling, us
    const double dtcrp = 0.512; // CRO sampling, us
    const double t_window = 2000.0;
    const double bubble_radius = 1.0;
    size_t start;
    size_t end; 

    std::string output_filename = "Rn222_PDHD_test_"+run+".root";
    
    //double r_charge_z_1 = 0.0;
    //double r_charge_z_2 = 0.0;
    //double r_charge_y_1 = 0.0;
    //double r_charge_y_2 = 0.0;
    //double r_charge_t_1 = 0.0;
    //double r_charge_t_2 = 0.0;
    //int r_charge_apa = 0;
    //int r_charge_nplanes_1 = 0;
    //int r_charge_nplanes_2 = 0;
    //double r_charge_energy_1 = 0.0;
    //double r_charge_energy_2 = 0.0;

    std::string identifier;
    int nfiles = 680;//filenames.size()
    for (int file_index = 0; file_index < nfiles; file_index++) {
        cout << "File " << file_index << "/" << nfiles << " processed" << endl;
        // this block finds the corresponding TPC clusters file
        std::ifstream singlehit_files(singlehit_runlist);
        start = filenames[file_index].find("raw_");
        end = filenames[file_index].find("_datawriter");
        if (start != std::string::npos && end != std::string::npos && start < end) {
            identifier = filenames[file_index].substr(start + 4, end - (start + 4));
            } else {
                std::cout << "Error getting charge cluster file!" << std::endl;
                continue;
                   }
        std::string singlehit_file = "";
        std::string file_line;
        while (std::getline(singlehit_files, file_line)) { 
            if (file_line.find(identifier) != std::string::npos) {
                singlehit_file = file_line;
                break;
                }
            }
        if (singlehit_file == ""){
            cout << "Cannot locate singlehit clusters file." << endl;
            continue;
        } else {
	    if (verbose){
            cout << "Using TPC clusters file: " << singlehit_file << endl;
            }
	}
     double r_charge_z_1 = 0.0;
             double r_charge_z_2 = 0.0;
             double r_charge_y_1 = 0.0;
             double r_charge_y_2 = 0.0;
             double r_charge_t_1 = 0.0;
             double r_charge_t_2 = 0.0;
             int r_charge_apa = 0;
             int r_charge_nplanes_1 = 0;
             int r_charge_nplanes_2 = 0;
             double r_charge_energy_1 = 0.0;
             double r_charge_energy_2 = 0.0;


        //TFile *outputFile = new TFile(output_filename.c_str(), "UPDATE");
        //TTree *results_tree = (TTree*) outputFile->Get("results");
        //if (!results_tree){
        //results_tree = new TTree("results", "Low energy charge-light matching results");
        //results_tree->Branch("z_1", &r_charge_z_1, "z_1/D");
        //results_tree->Branch("z_2", &r_charge_z_2, "z_2/D");
        //results_tree->Branch("y_1", &r_charge_y_1, "y_1/D");
        //results_tree->Branch("y_2", &r_charge_y_2, "y_2/D");
        //results_tree->Branch("t_1", &r_charge_t_1, "t_1/D");
        //results_tree->Branch("t_2", &r_charge_t_2, "t_2/D");
        //results_tree->Branch("apa", &r_charge_apa, "apa/I");
        //results_tree->Branch("nplanes_1", &r_charge_nplanes_1, "nplanes_1/I");
        //results_tree->Branch("nplanes_2", &r_charge_nplanes_2, "nplanes_2/I");
        //results_tree->Branch("E_1", &r_charge_energy_1, "E_1/D");
        //results_tree->Branch("E_2", &r_charge_energy_2, "E_2/D"); 
        //} else {
        //results_tree->SetBranchAddress("z_1", &r_charge_z_1);
        //results_tree->SetBranchAddress("z_2", &r_charge_z_2);
        //results_tree->SetBranchAddress("y_1", &r_charge_y_1);
        //results_tree->SetBranchAddress("y_2", &r_charge_y_2);
        //results_tree->SetBranchAddress("t_1", &r_charge_t_1);
        //results_tree->SetBranchAddress("t_2", &r_charge_t_2);
        //results_tree->SetBranchAddress("apa", &r_charge_apa);
        //results_tree->SetBranchAddress("nplanes_1", &r_charge_nplanes_1);
        //results_tree->SetBranchAddress("nplanes_2", &r_charge_nplanes_2);
        //results_tree->SetBranchAddress("E_1", &r_charge_energy_1);
        //results_tree->SetBranchAddress("E_2", &r_charge_energy_2);        
        //}
        // collect TPC clusters data from file
        TFile *clusters_file = TFile::Open(singlehit_file.c_str());
        if (!clusters_file || clusters_file->IsZombie()) {
            std::cerr << "Error: Cannot open ROOT file." << std::endl;
        }

        TDirectory *dir = (TDirectory*)clusters_file->Get("ana");
        TTree *clusters_tree = nullptr;
        dir->GetObject("ClusterTree", clusters_tree);
         
        std::vector<float>* clusters_Z = nullptr;
        clusters_tree->SetBranchAddress("Z", &clusters_Z);
        std::vector<float>* clusters_Y = nullptr; 
        clusters_tree->SetBranchAddress("Y", &clusters_Y);
        std::vector<double>* clusters_PeakTime = nullptr;
        clusters_tree->SetBranchAddress("PeakTime", &clusters_PeakTime);
        std::vector<int>* clusters_NOFTB_pre = nullptr;
        std::vector<int> clusters_NOFTB;
        clusters_tree->SetBranchAddress("NearOrFarToTheBeam", &clusters_NOFTB_pre);
        std::vector<float>* clusters_E = nullptr;
        clusters_tree->SetBranchAddress("EnergyCollection", &clusters_E);
        Double_t ChargeT0;
        clusters_tree->SetBranchAddress("CRP_T0", &ChargeT0);
        UInt_t eventID;
        clusters_tree->SetBranchAddress("eventID", &eventID);
        std::vector<int>* NumberOfPlane0 = nullptr;
        clusters_tree->SetBranchAddress("NumberOfPlane0", &NumberOfPlane0);
        std::vector<int>* NumberOfPlane1 = nullptr;
        clusters_tree->SetBranchAddress("NumberOfPlane1", &NumberOfPlane1);
        std::vector<int>* NumberOfCollection = nullptr;
        clusters_tree->SetBranchAddress("NumberOfCollection", &NumberOfCollection);
        Long64_t nEntries = clusters_tree->GetEntries();
        if (verbose){cout << "Total events in file = " << nEntries << endl;}

        // setup vectors to store TPC clusters data
        std::vector<double> charge_t;
        std::vector<int> charge_nplanes;
        std::vector<int> charge_apa;
        std::vector<float> charge_E;
        std::vector<float> charge_Z;
        std::vector<float> charge_Y;
        std::vector<int> charge_PeakTime;
        int lastLength = 0;
        double lastT0 = 0;
        std::vector<double> ChargeT0s;
        int vector_index = 0;
        for (Long64_t i = 0; i < nEntries; ++i) {
            clusters_tree->GetEntry(i);
            // combine results from each event
            vector_index = 0;
            for (int j = lastLength; j < clusters_NOFTB_pre->size(); j++){
	        //charge_event_indices.push_back(charge_event_index);
                clusters_NOFTB.push_back((*clusters_NOFTB_pre)[j]);
		//all_charge_t.push_back(ChargeT0*dt + (*clusters_PeakTime)[vector_index]*dtcrp);
                charge_t.push_back(ChargeT0*dt + (*clusters_PeakTime)[vector_index]*dtcrp);
                charge_E.push_back((*clusters_E)[vector_index]);
		charge_PeakTime.push_back((*clusters_PeakTime)[vector_index]);
                charge_Z.push_back((*clusters_Z)[vector_index]);
                charge_Y.push_back((*clusters_Y)[vector_index]);
                if (((*NumberOfPlane0)[vector_index]) && ((*NumberOfPlane1)[vector_index])) {
			charge_nplanes.push_back(3);
			}
		else {
			charge_nplanes.push_back(2);
			}
                // assign APA number to clusters, 
                // account for incorrectly formatted NearOrFarToTheBeam branch
		if (((*clusters_Z)[vector_index] < 230) && (clusters_NOFTB[vector_index] == 1)){
			charge_apa.push_back(1);
			}
		else if (((*clusters_Z)[vector_index] < 230) && (clusters_NOFTB[vector_index] == -1)){
			charge_apa.push_back(3);
			}
		else if (((*clusters_Z)[vector_index] > 230) && (clusters_NOFTB[vector_index] == 1)){
			charge_apa.push_back(2);
			}
                else if (((*clusters_Z)[vector_index] > 230) && (clusters_NOFTB[vector_index] == -1)){
                        charge_apa.push_back(4);
                        }
		else {
			charge_apa.push_back(-1);
			}
                vector_index++;
                }
        lastLength = clusters_NOFTB_pre->size();
        }
        clusters_file->Close(); 
        int total_matches = 0; 
        for (int i = 0; i < charge_t.size(); i++){
	     int nmatches = 0;
	     std::vector<int> matched_indices;
             for (int j = i+1; j < charge_t.size(); j++){
                 if ((charge_t[j] < charge_t[i] + t_window) &&
		     (charge_t[j] > charge_t[i]) &&
		     (sqrt(pow(charge_Z[i]-charge_Z[j], 2)+pow(charge_Y[i]-charge_Y[j], 2)) < bubble_radius) &&
 		     (charge_apa[i] == charge_apa[j])){
			nmatches++;
			matched_indices.push_back(j);
			}  
      	     } // end j loop
             if (nmatches != 1){continue;}
      
             double r_charge_z_1 = 0.0;
             double r_charge_z_2 = 0.0;
             double r_charge_y_1 = 0.0;
             double r_charge_y_2 = 0.0;
             double r_charge_DelT = 0.0;
             int r_charge_apa = 0;
             int r_charge_nplanes_1 = 0;
    	     int r_charge_nplanes_2 = 0;
             double r_charge_energy_1 = 0.0;
             double r_charge_energy_2 = 0.0;

             TFile *outputFile = new TFile(output_filename.c_str(), "UPDATE");
             TTree *results_tree = (TTree*) outputFile->Get("results");
             if (!results_tree){
              results_tree = new TTree("results", "Low energy charge-light matching results");
              results_tree->Branch("z_1", &r_charge_z_1, "z_1/D");
              results_tree->Branch("z_2", &r_charge_z_2, "z_2/D");
              results_tree->Branch("y_1", &r_charge_y_1, "y_1/D");
              results_tree->Branch("y_2", &r_charge_y_2, "y_2/D");
              results_tree->Branch("DelT", &r_charge_DelT, "DelT/D");
              results_tree->Branch("apa", &r_charge_apa, "apa/I");
        results_tree->Branch("nplanes_1", &r_charge_nplanes_1, "nplanes_1/I");
        results_tree->Branch("nplanes_2", &r_charge_nplanes_2, "nplanes_2/I");
        results_tree->Branch("E_1", &r_charge_energy_1, "E_1/D");
        results_tree->Branch("E_2", &r_charge_energy_2, "E_2/D");
        } else {
        results_tree->SetBranchAddress("z_1", &r_charge_z_1);
        results_tree->SetBranchAddress("z_2", &r_charge_z_2);
        results_tree->SetBranchAddress("y_1", &r_charge_y_1);
        results_tree->SetBranchAddress("y_2", &r_charge_y_2);
        results_tree->SetBranchAddress("DelT", &r_charge_DelT);
        results_tree->SetBranchAddress("apa", &r_charge_apa);
        results_tree->SetBranchAddress("nplanes_1", &r_charge_nplanes_1);
        results_tree->SetBranchAddress("nplanes_2", &r_charge_nplanes_2);
        results_tree->SetBranchAddress("E_1", &r_charge_energy_1);
        results_tree->SetBranchAddress("E_2", &r_charge_energy_2);
        }

             total_matches++;
	     int k = matched_indices[0];
             r_charge_z_1 = charge_Z[i];
             r_charge_z_2 = charge_Z[k];
             r_charge_y_1 = charge_Y[i];
	     r_charge_y_2 = charge_Y[k];
             r_charge_DelT = charge_t[i] - charge_t[k];
             r_charge_apa = charge_apa[i];
             r_charge_nplanes_1 = charge_nplanes[i];
             r_charge_nplanes_2 = charge_nplanes[k];
             r_charge_energy_1 = charge_E[i];
             r_charge_energy_2 = charge_E[k];
             results_tree->Fill();
             if (total_matches){results_tree->Write("", TObject::kOverwrite);}
             outputFile->Close();
             delete outputFile;
             //delete results_tree;
             //if (nmatches > 0){cout << "Number of clusters matched to cluster = " << nmatches << endl;}
        } // end i loop
        //results_tree->Write("", TObject::kOverwrite);
        //if (total_matches){results_tree->Write("", TObject::kOverwrite);}
        //outputFile->Close();
        //delete outputFile;
    }
}
