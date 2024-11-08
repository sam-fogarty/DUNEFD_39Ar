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

void LowEnergyCLMatching_PDHD(){
    //// Script for matching charge and light for low energy TPC clusters
    //// It requires 1) a runlist file pointing to the ART ROOT files containing PDS data
    ////             2) a runlist file pointing to ROOT ntuple files with TPC clusters
    //// Also the ART ROOT filepaths ideally need to be xrootd compatible, and user would
    ////    setup xrootd, kx509, and Gallery.
    //// Run `kx509;voms-proxy-init -noregen -rfc -voms dune:/dune/Role=Analysis`, and then
    //// precede the root command with `LD_PRELOAD=$XROOTD_LIB/libXrdPosixPreload.so`.

    // setup PDS channel map to lookup channel positions
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
        int channelNumber, APA, TPC;
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

    // get list of ART ROOT files to analyze
    string filelist = "28850_filepaths.txt";
    std::vector<std::string> filenames = readFilePaths(filelist);

    // get list of ROOT ntuple files with TPC clusters
    string run = "028850";
    string singlehit_runlist="/pnfs/dune/scratch/users/lavaut/03684/1/run28850.list";
    //string singlehit_runlist="run28850.list";
    string optag = "pdhddaphne:daq";
    InputTag rawop_tag(optag);
    bool verbose = true;
    // constants
    const double dt    = 16e-3; // DTS sampling, us
    const double dtcrp = 0.512; // CRO sampling, us
    const int selftrig_length = 1024;
    const int nhit_limit=2;

    // clustering parameters
    const double T_window = 0.1; // window to group pds hits in time, us
    const double z_space_window = 600; // window in Z and Y to check for isolated hits, cm
    const double y_space_window = 600;
    // matching paramters
    const double lower_t_window = 0;
    const double upper_t_window = 600;
    const double z_window = 55.0; // acceptance window in Z around PDS hit for TPC clusters
    const double y_window = 15.0; // acceptance window in Y around PDS hit for TPC clusters
    const double z_veto = 10.0;
    const double y_veto = 10.0;
    // misc
    const int baseline_ticks = 50; // how many ticks at start of wvfm to avg for baseline
    int chFound = 0;
    int chNotFound = 0;
    size_t start;
    size_t end;
   
    //std::vector<double> all_pds_t;
    //std::vector<double> all_charge_t;
    //std::vector<int> pds_event_indices;
    //std::vector<int> charge_event_indices;

    int pds_event_index = 0;
    int charge_event_index = 0;
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
	std::vector<std::vector<int>> charge_event_indices;
        for (Long64_t i = 0; i < nEntries; ++i) {
            clusters_tree->GetEntry(i);
            //cout << "Charge t0 = " <<  std::fixed << std::setprecision(18) << ChargeT0*dtcrp + upper_t_window << endl;
	    lastT0 = ChargeT0;
            ChargeT0s.push_back(ChargeT0*dt);
            
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
		std::vector<int> charge_event_start_stop = {lastLength, static_cast<int>(charge_Z.size())-1};
		charge_event_indices.push_back(charge_event_start_stop);
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
            std::vector<double> pdsTimestamp;
	    std::vector<int> pdsAPA;
    	    std::vector<int> pdsTPC;
   	    std::vector<double> pdsX;
            std::vector<double> pdsY;
            std::vector<double> pdsZ;
    	    std::vector<int> pdsChannel;
            std::vector<double> pdsAmplitude;

            pds_event_index++;
            auto const& rawops = *ev.getValidHandle<vector<raw::OpDetWaveform>>(rawop_tag);
            for (int i = 0; i < rawops.size(); i++){
                auto const& ADC_Count = rawops[i].Waveform();
                double baseline = 0;
                for (int j = 0; j < baseline_ticks; j++) {
		    baseline = baseline + ADC_Count[j];	
		}
                baseline = baseline / baseline_ticks;
                double amplitude = baseline - std::abs(*std::min_element(ADC_Count.begin(), ADC_Count.end()));
	        if (ADC_Count.size() > selftrig_length) {continue;} // ignore streaming wvfms
                //cout << std::fixed << std::setprecision(18) << rawops[i].TimeStamp() << endl;
                double ts = rawops[i].TimeStamp() * dt;
                int ch = rawops[i].ChannelNumber(); 
    
                // look up position of waveform's channel, and APA and TPC for channel 
                if (PDSChannelMap.find(ch) != PDSChannelMap.end()) {
                    chFound++;
                    pdschannelmap values = PDSChannelMap[ch];
                    //cout << "ch = " << ch << " , Y = " << values.Y << " , Z = " << values.Z << " , APA = " << values.APA << endl;
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

           } // end rawops loop
        //} // end gallery loop
     
        // do clustering of PDS hits
        // 1D
        double last_t = pdsTimestamp[0];
        std::vector<std::vector<int>> all_cluster_indices;
        std::vector<int> cluster_indices;
	
	int nhit_limit_dt = 5;
        for (int i = 0; i < pdsTimestamp.size(); i++){
	    int nhit = 1;
	    cluster_indices.push_back(i);
	    for (int j = 0; j < pdsTimestamp.size(); j++){
 	 	if (j == i){continue;}
	        if (abs(pdsTimestamp[j] - pdsTimestamp[i]) < T_window){
			nhit++;
			cluster_indices.push_back(j);
			}
		if (nhit > nhit_limit_dt){break;}
		}
            if (nhit <= nhit_limit_dt){all_cluster_indices.push_back(cluster_indices);}
            cluster_indices.clear();
	}
	
        if (verbose){cout << "dt clusters found = " << all_cluster_indices.size() << endl;}
        // group clusters in 2D, look for isolated single hits
        std::vector<bool> test_this_index(pdsX.size(), true);
        std::vector<std::vector<int>> hit_keep_indices;
        std::vector<int> hit_cluster_indices;
        int cluster_index = 0;
        int total_nhit1 = 0;
        int total_nhit2 = 0;
        for (int i = 0; i < all_cluster_indices.size(); i++){
             std::vector<int> cluster_hit_indices;
             for (int j = 0; j < all_cluster_indices[i].size(); j++){
                 cluster_hit_indices.push_back(all_cluster_indices[i][j]);
                 int nhit = 1;
	         for (int k = j+1; k < all_cluster_indices[i].size(); k++){
                     if ((abs(pdsZ[all_cluster_indices[i][k]] - pdsZ[all_cluster_indices[i][j]]) < z_space_window/2)
		     && (abs(pdsY[all_cluster_indices[i][k]] - pdsY[all_cluster_indices[i][j]]) < y_space_window/2) 
                     && (abs(pdsX[all_cluster_indices[i][k]] - pdsX[all_cluster_indices[i][j]]) < z_space_window/2))
	            {
		        nhit++;
                        cluster_hit_indices.push_back(all_cluster_indices[i][k]);
			if (nhit > nhit_limit) {break;}
		    }
		} // end k loop
                if (nhit == 1){
		    hit_keep_indices.push_back(cluster_hit_indices);
		    total_nhit1++;
		}
                cluster_hit_indices.clear();	
	        nhit = 0;	 
	     } // end j loop
        } // end i loop
	if (verbose){
        cout << "Total number of pds triggers = " << pdsTimestamp.size() << endl; 
        cout << "# of isolated pds hits found = " << total_nhit1 << endl;
        cout << "# of double isolated pds hits found = " << total_nhit2 << endl;
        }
	int totalMatches = 0; 
        int totalTwoViewMatches = 0;
        int totalThreeViewMatches = 0;
        
        double upper_z_window = z_window/2;
        double lower_z_window = z_window/2;
        double upper_y_window = y_window/2;
        double lower_y_window = y_window/2;
 
        // match charge and light
        std::vector<int> charge_indices;
        std::vector<int> light_indices;
        int cluster_match_limit = 1;
        for (int j = 0; j < hit_keep_indices.size(); j++){
            int hit_index = hit_keep_indices[j][0]; // only single hit clusters for now 
            int total_cluster_matches = 0;
            for (int k = charge_event_indices[evt_index][0]; k < charge_event_indices[evt_index][1]; k++){
		 if ((pdsChannel[hit_index] == 89) || (pdsChannel[hit_index] == 9)) {
			upper_z_window = 40;
			lower_z_window = z_window/2;
			upper_y_window = y_window/2;
			lower_y_window = 32;
		 } else {
			upper_z_window = z_window/2;
                        lower_z_window = z_window/2;
                        upper_y_window = y_window/2;
                        lower_y_window = y_window/2;
		 }
		 bool is_matched = (charge_t[k] < pdsTimestamp[hit_index] + upper_t_window) 
                     && (charge_t[k] > pdsTimestamp[hit_index]-lower_t_window)
                     && (charge_Z[k] < pdsZ[hit_index] + upper_z_window)
                     && (charge_Z[k] > pdsZ[hit_index] - lower_z_window)
                     && (charge_Y[k] < pdsY[hit_index] + upper_y_window)
                     && (charge_Y[k] > pdsY[hit_index] - lower_y_window)
                     && (charge_apa[k] == pdsAPA[hit_index]);
                 if (is_matched) {
			for (int l = charge_event_indices[evt_index][0]; l < charge_event_indices[evt_index][1]; l++){
			   bool rightSide = ((charge_Z[l] < pdsZ[hit_index]+upper_z_window+z_veto) &&
			      (charge_Z[l] > pdsZ[hit_index]+upper_z_window) && 
			      (charge_Y[l] < pdsY[hit_index]+upper_y_window+y_veto) &&
			      (charge_Y[l] > pdsY[hit_index]-lower_y_window-y_veto));
                           bool leftSide = ((charge_Z[l] > pdsZ[hit_index]-lower_z_window-z_veto) &&
                              (charge_Z[l] < pdsZ[hit_index]-lower_z_window) &&
                              (charge_Y[l] < pdsY[hit_index]+upper_y_window+y_veto) &&
                              (charge_Y[l] > pdsY[hit_index]-lower_y_window-y_veto));
	                   bool topSide = ((charge_Z[l] > pdsZ[hit_index]-lower_z_window-z_veto) &&
                              (charge_Z[l] < pdsZ[hit_index]+upper_z_window+z_veto) &&
                              (charge_Y[l] < pdsY[hit_index]+upper_y_window+y_veto) &&
                              (charge_Y[l] > pdsY[hit_index]+upper_y_window));
                           bool bottomSide = ((charge_Z[l] > pdsZ[hit_index]-lower_z_window-z_veto) &&
                              (charge_Z[l] < pdsZ[hit_index]+upper_z_window+z_veto) &&
                              (charge_Y[l] > pdsY[hit_index]-lower_y_window-y_veto) &&
                              (charge_Y[l] < pdsY[hit_index]-lower_y_window));
		           if (rightSide || leftSide || topSide || bottomSide){is_matched=false;}
			}
			if (!is_matched){break;}
 			total_cluster_matches++;
		        if (total_cluster_matches > cluster_match_limit){break;}
                            totalMatches++;
			    if (hit_keep_indices[j].size() == 1){
     			        charge_indices.push_back(k);
			        light_indices.push_back(hit_index);
			    }
 			    if (charge_nplanes[k] == 2){totalTwoViewMatches++;}
				else {totalThreeViewMatches++;}
			    }
		    } // end k loop  
	} // end j loop
        if (verbose){
	cout << "Total number of charge clusters = " << charge_t.size() << endl;
        cout << "# of charge clusters matched to pds hits = " << totalMatches << endl;
        cout << "# of two-view charge clusters matched to pds hits = " << totalTwoViewMatches << endl;
        cout << "# of three-view charge clusters matched to pds hits = " << totalThreeViewMatches << endl;
        cout << "# of single PDS hit coincidences = " << charge_indices.size() << endl;
	cout << "# of double PDS hit coincidences = " << charge_indices.size() << endl;
        }
	// save results to ROOT file
        std::string output_filename = "pdhd_"+run+"_LowECLMatching_test2"+".root"; 
        TFile *outputFile = new TFile(output_filename.c_str(), "UPDATE");
        TTree *results_tree = (TTree*) outputFile->Get("results");
        TTree *pds_tree = (TTree*) outputFile->Get("isolated_pds"); 
        TTree *charge_tree = (TTree*) outputFile->Get("isolated_charge");
        double r_all_pds_z = 0.0;
	double r_all_pds_y = 0.0;
        double r_all_pds_x = 0.0;
	double r_all_pds_t = 0.0;
        int r_all_pds_ch = 0;
        double r_all_pds_amplitude = 0.0;
        int r_all_pds_apa = 0;

	double r_all_charge_z = 0.0;
	double r_all_charge_y = 0.0;
	double r_all_charge_t = 0.0;
        double r_all_charge_energy = 0.0;	
        int r_all_charge_apa = 0;
	int r_all_charge_nplanes = 0;
        
        if (charge_tree){
           charge_tree->SetBranchAddress("z", &r_all_charge_z);
           charge_tree->SetBranchAddress("y", &r_all_charge_y);
           charge_tree->SetBranchAddress("t", &r_all_charge_t);
           charge_tree->SetBranchAddress("apa", &r_all_charge_apa);
	   charge_tree->SetBranchAddress("nplanes", &r_all_charge_nplanes);
           charge_tree->SetBranchAddress("E", &r_all_charge_energy);
        } else {
           charge_tree = new TTree("isolated_charge", "Low energy charge-light matching results");
           charge_tree->Branch("z", &r_all_charge_z, "z/D");
           charge_tree->Branch("y", &r_all_charge_y, "y/D");
           charge_tree->Branch("t", &r_all_charge_t, "t/D");
           charge_tree->Branch("apa", &r_all_charge_apa, "apa/I");
	   charge_tree->Branch("nplanes", &r_all_charge_nplanes, "nplanes/I");
           charge_tree->Branch("E", &r_all_charge_energy, "E/D");
        }
	if (pds_tree){
	   pds_tree->SetBranchAddress("z", &r_all_pds_z);
           pds_tree->SetBranchAddress("y", &r_all_pds_y);
           pds_tree->SetBranchAddress("x", &r_all_pds_x);
           pds_tree->SetBranchAddress("t", &r_all_pds_t);
	   pds_tree->SetBranchAddress("ch", &r_all_pds_ch);
	   pds_tree->SetBranchAddress("apa", &r_all_pds_apa);
           pds_tree->SetBranchAddress("amplitude", &r_all_pds_amplitude);
	} else {
	   pds_tree = new TTree("isolated_pds", "Low energy charge-light matching results");
           pds_tree->Branch("z", &r_all_pds_z, "z/D");
           pds_tree->Branch("y", &r_all_pds_y, "y/D");
           pds_tree->Branch("x", &r_all_pds_x, "x/D");
           pds_tree->Branch("t", &r_all_pds_t, "t/D");
           pds_tree->Branch("ch", &r_all_pds_ch, "ch/I");
	   pds_tree->Branch("apa", &r_all_pds_apa, "apa/I");
           pds_tree->Branch("amplitude", &r_all_pds_amplitude, "amplitude/D");
	}
        for (int j = 0; j < charge_Z.size(); j++){
	   r_all_charge_z = charge_Z[j];
	   r_all_charge_y = charge_Y[j];
	   r_all_charge_t = charge_t[j];
	   r_all_charge_apa = charge_apa[j];
	   r_all_charge_nplanes = charge_nplanes[j];	
	   r_all_charge_energy = charge_E[j];
           charge_tree->Fill();
	}
	charge_tree->Write("", TObject::kOverwrite);
        
        for (int j = 0; j < hit_keep_indices.size(); j++){
	   r_all_pds_z = pdsZ[hit_keep_indices[j][0]];
	   r_all_pds_y = pdsY[hit_keep_indices[j][0]];
	   r_all_pds_x = pdsX[hit_keep_indices[j][0]];
           r_all_pds_ch = pdsChannel[hit_keep_indices[j][0]];
	   r_all_pds_t = pdsTimestamp[hit_keep_indices[j][0]];
	   r_all_pds_apa = pdsAPA[hit_keep_indices[j][0]];
	   r_all_pds_amplitude = pdsAmplitude[hit_keep_indices[j][0]];
	   pds_tree->Fill();
	}
        pds_tree->Write("", TObject::kOverwrite);
        int r_pds_ch = 0;
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
	int r_charge_peaktime = 0;
        if (results_tree) {
	    results_tree->SetBranchAddress("charge_peaktime", &r_charge_peaktime);
            results_tree->SetBranchAddress("pds_ch", &r_pds_ch);
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
       } else {
           results_tree = new TTree("results", "Low energy charge-light matching results");      
   	   results_tree->Branch("pds_ch", &r_pds_ch, "pds_ch/I");
	   results_tree->Branch("charge_peaktime", &r_charge_peaktime, "charge_peaktime/I");
	   results_tree->Branch("pds_z", &r_pds_z, "pds_z/D");
           results_tree->Branch("pds_y", &r_pds_y, "pds_y/D");
           results_tree->Branch("pds_t", &r_pds_t, "pds_t/D");
           results_tree->Branch("pds_amplitude", &r_pds_amplitude, "pds_amplitude/D");
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
          r_pds_ch = pdsChannel[light_indices[j]];
          r_pds_z = pdsZ[light_indices[j]];
          r_pds_y = pdsY[light_indices[j]];
          r_charge_peaktime = charge_PeakTime[charge_indices[j]];
          r_pds_t = static_cast<double>(pdsTimestamp[light_indices[j]]);
	  r_pds_amplitude = pdsAmplitude[light_indices[j]];
	  r_charge_z = charge_Z[charge_indices[j]];
	  r_charge_y = charge_Y[charge_indices[j]];
	  r_charge_energy = charge_E[charge_indices[j]];
	  r_charge_t = static_cast<double>(charge_t[charge_indices[j]]);
	  r_charge_nplanes = charge_nplanes[charge_indices[j]];
          r_apa = pdsAPA[light_indices[j]];
          r_dt = static_cast<double>(charge_t[charge_indices[j]] - pdsTimestamp[light_indices[j]]);
          //cout << "dt = " << r_dt << endl;
          if (r_apa >= 3){
	      r_charge_x = 353.0 - r_dt*0.16;
          } else {
	      r_charge_x = -353.2 + r_dt*0.16;
	  }
          results_tree->Fill();
	  }
       results_tree->Write("", TObject::kOverwrite);
       
       outputFile->Close();
       delete outputFile;
   
       // clear vectors for next file set
       pdsAPA.clear();
       pdsTPC.clear();
       pdsAmplitude.clear();
       pdsX.clear();
       pdsY.clear();
       pdsZ.clear();
       pdsChannel.clear();
       pdsTimestamp.clear();
      evt_index++;
      } // end gallery loop 
    }
    //cout << "Channels found = " << chFound << endl;
    //cout << "Channels not found = " << chNotFound << endl;
}
