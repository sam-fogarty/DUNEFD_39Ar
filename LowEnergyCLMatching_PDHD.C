R__ADD_INCLUDE_PATH("gallery/Event.h")

#include "tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>

using namespace art;

struct pdschannelmap {
    int APA;
    int TPC;
    double X;
    double Y;
    double Z;
};

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

void leclm(){
std::unordered_map<int, pdschannelmap> PDSChannelMap;
std::ifstream channelmap_file("PDHD_PDS_ChannelMap.csv");
if (!channelmap_file.is_open()) {
        std::cerr << "Error opening channel map file!" << std::endl;
    }
std::string line;
bool firstLine = true;
while (std::getline(channelmap_file, line)) {
        if (firstLine) {firstLine = false; continue;}
        std::stringstream ss(line);
        std::string token;
        int channelNumber;
        int APA;
        int TPC;
        double X, Y, Z, Z_det;

        std::getline(ss, token, ',');
        channelNumber = std::stoi(token);
	std::getline(ss, token, ',');
        APA = std::stoi(token);
	std::getline(ss, token, ',');
        TPC = std::stoi(token);
        std::getline(ss, token, ',');
        X = std::stod(token);
        std::getline(ss, token, ',');
        Y = std::stod(token);
        std::getline(ss, token, ',');
        Z_det = std::stod(token);
        std::getline(ss, token, ',');
        Z = std::stod(token);
        PDSChannelMap[channelNumber] = { APA, TPC, X, Y, Z };
}
channelmap_file.close();

string filelist = "28850_filepaths_test.txt";
std::vector<std::string> filenames = readFilePaths(filelist);
//std::string line;

string filelist_root = "28850_filepaths_root.txt";
std::vector<std::string> filenames_root = readFilePaths(filelist_root);

string run = "028850";
string singlehit_runlist="/pnfs/dune/scratch/users/lavaut/03684/1/run28850.list";
//std::ifstream singlehit_files(singlehit_runlist);

string optag = "pdhddaphne:daq";
InputTag rawop_tag(optag);

double dt    = 16e-3; // us
double dtcrp = 0.512;
int selftrig_length = 1024;
int evt = 0;
int nhit_limit=1;

// clustering parameters
const double T_window = 0.1; // usec
const double outer_T_window = 15*dt;
const double space_window = 600; // cm
const int min_samples_T = 1;
const int min_samples_X = 1;

// matching paramters
const double lower_t_window = 5000;
const double upper_t_window = 5000;
const double z_window = 60;
const double y_window = 15;

int baseline_ticks = 50; // how many ticks at start of wvfm to avg for baseline

int chFound = 0;
int chNotFound = 0;
size_t start;
size_t end;
//std::vector<std::vector<double>> pdsTimestamp;
std::vector<double> pdsTimestamp;
std::vector<int> pdsAPA;
std::vector<int> pdsTPC;
std::vector<double> pdsX;
std::vector<double> pdsY;
std::vector<double> pdsZ;
std::vector<int> pdsChannel;
std::vector<double> pdsAmplitude;

std::vector<double> all_pds_t;
std::vector<double> all_charge_t;
std::vector<int> pds_event_indices;
std::vector<int> charge_event_indices;
int pds_event_index = 0;
int charge_event_index = 0;
std::string identifier;
int nfiles = 10;//filenames.size()
for (int file_index = 0; file_index < nfiles; file_index++) {
   
   cout << " " << endl;
   cout << "Processing file " << file_index+1 << "/" << filenames.size() << endl;
   std::ifstream singlehit_files(singlehit_runlist);
   start = filenames[file_index].find("raw_");
   end = filenames[file_index].find("_datawriter");
   //cout << filenames[file_index] << endl;
   if (start != std::string::npos && end != std::string::npos && start < end) {
         // extract the substring between "raw_" and "_datawriter"
         identifier = filenames[file_index].substr(start + 4, end - (start + 4));
         //cout << identifier << endl;
         } else {
             std::cout << "Error getting charge cluster file.!" << std::endl;
             continue;
         }
   std::string singlehit_file = "";
   std::string file_line;
   while (std::getline(singlehit_files, file_line)) { 
         //if (line.find(identifier) != std::string::npos){cout<<"yes"<< endl;} else{cout<<"no"<<endl;}
         if (file_line.find(identifier) != std::string::npos) {
             singlehit_file = file_line;
             break;
             }
         //else{cout << "Cannot find identifier in filenames" << endl;}
      }
   //cout << singlehit_file << endl;
   if (singlehit_file == ""){
              cout << "Cannot locate singlehit clusters file." << endl;
              continue;
   } else {
              cout << "Using singlehit clusters file: " << singlehit_file << endl;
   }
   
   std::string filename_root = filenames_root[file_index];
   cout << filename_root << endl;
   TFile *pds_root_file = TFile::Open(filename_root.c_str());
   if (!pds_root_file){cout << "ROOT file not found" << endl;}
   TTree *pds_hits_tree = nullptr;
   pds_root_file->GetObject("hits", pds_hits_tree);
   if (!pds_hits_tree){cout << "No hits dataset found" << endl;}
   std::vector<float> pds_z_root_vec;
   std::vector<float> pds_y_root_vec;
   std::vector<Long64_t> pds_t_root_vec;
   std::vector<Short_t> pds_apa_root_vec;
   float pds_z_root;
   float pds_y_root;
   Long64_t pds_t_root;
   Short_t pds_apa_root;

   pds_hits_tree->SetBranchAddress("Z", &pds_z_root);
   pds_hits_tree->SetBranchAddress("Y", &pds_y_root);
   pds_hits_tree->SetBranchAddress("T", &pds_t_root);
   pds_hits_tree->SetBranchAddress("APA", &pds_apa_root);
   Long64_t nEntries_root = pds_hits_tree->GetEntries();
   for (Long64_t i = 0; i < nEntries_root; ++i) {
        pds_hits_tree->GetEntry(i);
        pds_z_root_vec.push_back(pds_z_root); 
        pds_y_root_vec.push_back(pds_y_root);
	pds_t_root_vec.push_back(pds_t_root);
	pds_apa_root_vec.push_back(pds_apa_root);  
   }
   pds_root_file->Close();
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
   //clusters_tree->Print();
   std::vector<double>* clusters_PeakTime = nullptr;
   clusters_tree->SetBranchAddress("PeakTime", &clusters_PeakTime);
   std::vector<int>* clusters_NOFTB_pre = nullptr;
   std::vector<int> clusters_NOFTB;
   clusters_tree->SetBranchAddress("NearOrFarToTheBeam", &clusters_NOFTB_pre);
   std::vector<float>* clusters_E = nullptr;
   clusters_tree->SetBranchAddress("EnergyCollection", &clusters_E);
   double ChargeT0;
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
   cout << "Total events in file = " << nEntries << endl;

   std::vector<double> charge_t;
   std::vector<int> charge_nplanes;
   std::vector<int> charge_apa;
   std::vector<float> charge_E;
   std::vector<float> charge_Z;
   std::vector<float> charge_Y;
   int lastLength = 0;
   double lastT0 = 0;
   std::vector<double> ChargeT0s;
   for (Long64_t i = 0; i < nEntries; ++i) {
        clusters_tree->GetEntry(i);
        //cout << ChargeT0*dt << endl;
        //cout << lastT0*dt << endl;
        //cout << (ChargeT0 - lastT0)*dt << endl;
	lastT0 = ChargeT0;
        ChargeT0s.push_back(ChargeT0*dt);
        //cout << " " << endl;
        // combine results from each event
        int vector_index = 0;
        for (int j = lastLength; j < clusters_NOFTB_pre->size(); j++){
		charge_event_indices.push_back(charge_event_index);
		clusters_NOFTB.push_back((*clusters_NOFTB_pre)[j]);
		all_charge_t.push_back(ChargeT0*dt + (*clusters_PeakTime)[vector_index]*dtcrp);
                charge_t.push_back(ChargeT0*dt + (*clusters_PeakTime)[vector_index]*dtcrp);
                charge_E.push_back((*clusters_E)[vector_index]);
                charge_Z.push_back((*clusters_Z)[vector_index]);
                charge_Y.push_back((*clusters_Y)[vector_index]);
                if (((*NumberOfPlane0)[vector_index]) && ((*NumberOfPlane1)[vector_index])) {
			charge_nplanes.push_back(3);
			}
		else {
			charge_nplanes.push_back(2);
			}
                // account for incorrect formatted NearOrFarToTheBeam branch
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
        charge_event_index++;
        lastLength = clusters_NOFTB_pre->size();
     }
   clusters_file->Close();    
   // get PDS data
   std::vector<int> ch_sofar;
   std::vector<int> ch_sofar2;
   std::vector<std::string> filename = {filenames[file_index]};
   int evt_index = 0;
   for (gallery::Event ev(filename); !ev.atEnd(); ev.next()) {
       pds_event_index++;
       auto const& rawops = *ev.getValidHandle<vector<raw::OpDetWaveform>>(rawop_tag);
       //cout << rawops[0].TimeStamp()*dt - ChargeT0s[evt_index] << endl;
       evt_index++;
       for (int i = 0; i < rawops.size(); i++){
            auto const& ADC_Count = rawops[i].Waveform();
            double baseline = 0;
            for (int j = 0; j < baseline_ticks; j++) {
		baseline = baseline + ADC_Count[j];	
		}
            baseline = baseline / baseline_ticks;
            double amplitude = baseline - std::abs(*std::min_element(ADC_Count.begin(), ADC_Count.end()));
	    if (ADC_Count.size() > selftrig_length) {continue;} // ignore streaming wvfms
            double ts = rawops[i].TimeStamp() * dt;
            all_pds_t.push_back(ts);
	    pds_event_indices.push_back(pds_event_index);
            int ch = rawops[i].ChannelNumber(); 
            //cout << rawops[i].TimeStamp()*dt - rawops[i+1].TimeStamp()*dt << endl;
            // look up position of waveform's channel, and APA and TPC for channel
            if (PDSChannelMap.find(ch) != PDSChannelMap.end()) {
                chFound++;
                pdschannelmap values = PDSChannelMap[ch];
                pdsAPA.push_back(values.APA);
		pdsTPC.push_back(values.TPC);
		pdsX.push_back(values.X);
		pdsY.push_back(values.Y);
		pdsZ.push_back(values.Z);
                pdsChannel.push_back(ch);
		pdsTimestamp.push_back(ts);
                pdsAmplitude.push_back(amplitude);
            } else {
                chNotFound++;
                continue;
            }

       }
       evt++;    
   }
   std::set<int> unique_set(pdsChannel.begin(), pdsChannel.end());
   std::vector<int> unique_vec(unique_set.begin(), unique_set.end());
   //for (int u = 0; u < unique_vec.size(); u++){cout<<unique_vec[u] << endl;}
   std::map<int, int> count_map;

    // Count occurrences of each element
         for (int val : ch_sofar) {
                 count_map[val]++;
                     }
    
    //                     // Print each unique value and its count
                             //for (const auto& pair : count_map) {
                             //        std::cout << "Value: " << pair.first << ", Count: " << pair.second << "\n";
                             //            }

   // do clustering of PDS hits
   // 1D
   double last_t = pdsTimestamp[0];
   std::vector<std::vector<int>> all_cluster_indices;
   std::vector<int> cluster_indices;
   cluster_indices.push_back(0);
   for (int i = 1; i < pdsTimestamp.size(); i++){
	double current_t = pdsTimestamp[i];
        if (current_t <= last_t + T_window) {cluster_indices.push_back(i);}
	else {
		all_cluster_indices.push_back(cluster_indices);
		cluster_indices.clear();
		cluster_indices.push_back(i);
	     }
        last_t = current_t;
	}
   cout << "dt clusters found = " << all_cluster_indices.size() << endl;
   
   std::vector<bool> test_this_index(pdsX.size(), true);
   std::vector<std::vector<int>> hit_keep_indices;
   std::vector<int> hit_cluster_indices;
   int cluster_index = 0;
   for (int i = 0; i < all_cluster_indices.size(); i++){
        std::vector<int> cluster_hit_indices;
        for (int j = 0; j < all_cluster_indices[i].size(); j++){
          cluster_hit_indices.push_back(all_cluster_indices[i][j]);
          int nhit = 1;
	  for (int k = j+1; k < all_cluster_indices[i].size(); k++)
            {
             if ((pdsZ[all_cluster_indices[i][k]] < pdsZ[all_cluster_indices[i][j]] + space_window/2)
		&& (pdsZ[all_cluster_indices[i][k]] > pdsZ[all_cluster_indices[i][j]] - space_window/2) 
		&& (pdsY[all_cluster_indices[i][k]] < pdsY[all_cluster_indices[i][j]] + space_window/2) 
		&& (pdsY[all_cluster_indices[i][k]] > pdsY[all_cluster_indices[i][j]] - space_window/2)
                && (pdsX[all_cluster_indices[i][k]] < pdsX[all_cluster_indices[i][j]] + space_window/2)
                && (pdsX[all_cluster_indices[i][k]] > pdsX[all_cluster_indices[i][j]] - space_window/2))
	            {
		
		     nhit++;
                     cluster_hit_indices.push_back(all_cluster_indices[i][k]);
                     if (nhit > nhit_limit) {break;}
                
		    }
            }
            if (nhit <= nhit_limit){
		hit_keep_indices.push_back(cluster_hit_indices);
		if (!(std::find(ch_sofar2.begin(), ch_sofar2.end(), pdsChannel[cluster_hit_indices[0]]) != ch_sofar2.end()))
		{
			ch_sofar2.push_back(pdsChannel[cluster_hit_indices[0]]);
			//cout << pdsChannel[cluster_hit_indices[0]] << endl;	
		}
		}
            cluster_hit_indices.clear();	
	    nhit = 0;	 
	}
   }
   pdsTimestamp.clear();
   pdsZ.clear();
   pdsY.clear();
   pdsAPA.clear();
   hit_keep_indices.clear();
   for (int j = 0; j < pds_z_root_vec.size(); j++){
	pdsTimestamp.push_back(static_cast<double>(pds_t_root_vec[j]));
        pdsZ.push_back(pds_z_root_vec[j]);
        pdsY.push_back(pds_y_root_vec[j]);
        pdsAPA.push_back(pds_apa_root_vec[j]);
        std::vector<int> hit_keep_i;
        hit_keep_i.push_back(j);
        hit_keep_indices.push_back(hit_keep_i);
	}
   cout << "Total number of pds triggers = " << pdsTimestamp.size() << endl; 
   cout << "# of isolated pds hits found = " << hit_keep_indices.size() << endl;
   int totalMatches = 0; 
   int totalTwoViewMatches = 0;
   int totalThreeViewMatches = 0;

   std::vector<int> charge_indices;
   std::vector<int> light_indices;
   int cluster_match_limit = 1;
   for (int j = 0; j < hit_keep_indices.size(); j++){
        int hit_index = hit_keep_indices[j][0]; // only single hit clusters for now 
        int total_cluster_matches = 0;
        for (int k = 0; k < charge_t.size(); k++){
		//cout << pdsTimestamp[hit_index] - charge_t[k] << endl;
		//cout << "charge t = " << charge_t[k] << endl;
                //cout << charge_t[k] - pdsTimestamp[hit_index] << endl;
                if ((charge_t[k] < pdsTimestamp[hit_index] + upper_t_window) 
			&& (charge_t[k] > pdsTimestamp[hit_index]-lower_t_window)
			&& (charge_Z[k] < pdsZ[hit_index] + z_window/2)
			&& (charge_Z[k] > pdsZ[hit_index] - z_window/2)
			&& (charge_Y[k] < pdsY[hit_index] + y_window/2)
			&& (charge_Y[k] > pdsY[hit_index] - y_window/2)) {
				//cout << pdsTimestamp[hit_index] - charge_t[k] << endl;
				total_cluster_matches++;
				if (total_cluster_matches > cluster_match_limit){break;}
                                totalMatches++;	
     				charge_indices.push_back(k);
				light_indices.push_back(hit_index);
 				if (charge_nplanes[k] == 2){totalTwoViewMatches++;}
				else {totalThreeViewMatches++;}
			}
		}
	}
   cout << "Total number of charge clusters = " << charge_t.size() << endl;
   cout << "# of charge clusters matched to pds hits = " << totalMatches << endl;
   //cout << "Fraction of charge clusters with a match = " << totalMatches/charge_t.size() << endl;
   cout << "# of two-view charge clusters matched to pds hits = " << totalTwoViewMatches << endl;
   cout << "# of three-view charge clusters matched to pds hits = " << totalThreeViewMatches << endl;
  
   std::string output_filename = "pdhd_"+run+"_LowECLMatching_test"+".root"; 
   TFile *outputFile = new TFile(output_filename.c_str(), "UPDATE");
   TTree *results_tree = (TTree*) outputFile->Get("results");
   
   double r_pds_z = 0.0;
   double r_pds_y = 0.0;
   double r_pds_t = 0.0;
   double r_pds_amplitude = 0.0;
   float r_charge_z = 0.0;
   float r_charge_y = 0.0;
   double r_charge_t = 0.0;
   double r_charge_energy = 0.0;
   int r_charge_nplanes = 0;
   double r_charge_x = 0.0;
   double r_dt = 0.0;
   int r_apa = 0;
   if (results_tree) {
         results_tree->SetBranchAddress("pds_z", &r_pds_z);
         results_tree->SetBranchAddress("pds_y", &r_pds_y);
         results_tree->SetBranchAddress("pds_t", &r_pds_t);
         results_tree->SetBranchAddress("pds_amplitude", &r_pds_amplitude);
         results_tree->SetBranchAddress("charge_z", &r_charge_z);
         results_tree->SetBranchAddress("charge_y", &r_charge_y);
         results_tree->SetBranchAddress("charge_t", &r_charge_t);
         results_tree->SetBranchAddress("charge_energy", &r_charge_energy);
         results_tree->SetBranchAddress("charge_nplanes", &r_charge_nplanes);
         results_tree->SetBranchAddress("charge_x", &r_charge_x);
         results_tree->SetBranchAddress("dt", &r_dt);
         results_tree->SetBranchAddress("APA", &r_apa);
       }
   else {
         results_tree = new TTree("results", "Low energy charge-light matching results");      
   	 results_tree->Branch("pds_z", &r_pds_z, "pds_z/D");
         results_tree->Branch("pds_y", &r_pds_y, "pds_y/D");
         results_tree->Branch("pds_t", &r_pds_t, "pds_t/D");
         //results_tree->Branch("pds_amplitude", &r_pds_amplitude, "pds_amplitude/D");
         results_tree->Branch("charge_z", &r_charge_z, "charge_z/F");
         results_tree->Branch("charge_y", &r_charge_y, "charge_y/F");
         results_tree->Branch("charge_t", &r_charge_t, "charge_t/D");
         results_tree->Branch("charge_energy", &r_charge_energy, "charge_energy/D");
         results_tree->Branch("charge_nplanes", &r_charge_nplanes, "charge_nplanes/I");
         results_tree->Branch("charge_x", &r_charge_x, "charge_x/D");
         results_tree->Branch("dt", &r_dt, "dt/D");
         results_tree->Branch("APA", &r_apa, "APA/I");
      }
   for (int j = 0; j < light_indices.size(); j++){
	r_pds_z = pdsZ[light_indices[j]];
	r_pds_y = pdsY[light_indices[j]];
        r_pds_t = pdsTimestamp[light_indices[j]];
	r_pds_amplitude = pdsAmplitude[light_indices[j]];
	r_charge_z = charge_Z[charge_indices[j]];
	r_charge_y = charge_Y[charge_indices[j]];
	r_charge_energy = charge_E[charge_indices[j]];
	r_charge_t = charge_t[charge_indices[j]];
	r_charge_nplanes = charge_nplanes[charge_indices[j]];
        r_apa = pdsAPA[light_indices[j]];
        r_dt = charge_t[charge_indices[j]] - pdsTimestamp[light_indices[j]];
        if (r_apa >= 3){
		r_charge_x = 353.0 - r_dt*0.16;
             }
        else {
		r_charge_x = -353.2 + r_dt*0.16;
	}
        results_tree->Fill();
	}
   results_tree->Write("", TObject::kOverwrite);
   outputFile->Close();
   delete outputFile;
   
   pdsAPA.clear();
   pdsTPC.clear();
   pdsAmplitude.clear();
   pdsX.clear();
   pdsY.clear();
   pdsZ.clear();
   pdsChannel.clear();
   pdsTimestamp.clear();
   
}
double r_all_pds_t = 0.0;
double r_all_charge_t = 0.0;
int r_pds_event_index = 0;
int r_charge_event_index = 0;

std::string output_filename2 = "all_timestamps.root";
//TTree *light = (TTree*) outputFile->Get("results");
TFile *outputFile2 = new TFile(output_filename2.c_str(), "UPDATE");
TTree *light_tree = (TTree*) outputFile2->Get("light_tree");
TTree *charge_tree = (TTree*) outputFile2->Get("charge_tree");
if (!light_tree){
	light_tree = new TTree("light", "Low energy charge-light matching results");
}
if (!charge_tree){
        charge_tree = new TTree("charge", "Low energy charge-light matching results");
}
light_tree->Branch("t", &r_all_pds_t, "t/D");
light_tree->Branch("eventID", &r_pds_event_index, "eventID/I");

charge_tree->Branch("t", &r_all_charge_t, "t/D");
charge_tree->Branch("eventID", &r_charge_event_index, "eventID/I");

for (int i = 0; i < all_pds_t.size(); i++){
	r_all_pds_t = all_pds_t[i];
	r_pds_event_index = pds_event_indices[i];	
	light_tree->Fill();
	}
for (int i = 0; i < all_charge_t.size(); i++){
        r_all_charge_t = all_charge_t[i];
        r_charge_event_index = charge_event_indices[i];
        charge_tree->Fill();
        }
light_tree->Write("", TObject::kOverwrite);
charge_tree->Write("", TObject::kOverwrite);
outputFile2->Close();

cout << "Channels found = " << chFound << endl;
cout << "Channels not found = " << chNotFound << endl;
}
