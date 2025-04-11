R__ADD_INCLUDE_PATH("gallery/Event.h")

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <stdexcept>
#include <regex>
#include <filesystem>

using namespace art;

void handleError(const std::string& message) {
    std::cerr << message << std::endl;
    throw std::runtime_error(message);
}

struct pdschannelmap {
    int APA;
    int TPC;
    double X;
    double Y;
    double Z;
};

struct chargeData {
    std::vector<double> charge_Z;
    std::vector<double> charge_Y;
    std::vector<double> charge_T;
    std::vector<int> charge_apa;
    std::vector<int> charge_nplanes;
    std::vector<double> charge_E;
    std::vector<int> charge_PeakTime;
    std::vector<double> ChargeT0s;
    std::vector<std::vector<int>> charge_event_indices;
};

struct pdsData {
    std::vector<double> pdsZ;
    std::vector<double> pdsY;
    std::vector<double> pdsX;
    std::vector<int> pdsChannel;
    std::vector<double> pdsTimestamp;
    std::vector<int> pdsAPA;
    std::vector<int> pdsTPC;
    std::vector<double> pdsAmplitude;
    std::vector<double> pdsBaseline;
    std::vector<std::vector<double>> pdsWaveform;
};

std::string getIthLine(const std::string& filename, int i) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return "";
    }

    std::string line;
    for (int currentLine = 0; currentLine <= i; ++currentLine) {
        if (!std::getline(file, line)) {
            std::cerr << "Line " << i << " not found in file." << std::endl;
            return "";
        }
    }

    return line;
}

// get list of filepaths from a text file
std::vector<std::string> readFilePaths(const std::string& fileName) {
    std::ifstream file(fileName);
    std::vector<std::string> filePaths;
    std::string line;

    if (!file.is_open()) {
        handleError("Unable to open filelist file: " + fileName);
    }

    while (std::getline(file, line)) {
        filePaths.push_back(line);
    }

    file.close();
    return filePaths;
}

std::string findRecoFileIdentifier(const std::string& singlehit_file) {
    //// Find identifier for SingleHit file that can be used to find corresponding ART ROOT file
    // Input: Path to txt file with paths to SingleHit ROOT ntuple files
    size_t startPos = singlehit_file.find("singleHit_output_");
    if (startPos == std::string::npos) {
        handleError("Invalid root file name format. Expected to find 'singleHit_output_' in the filename.");
    }
    namespace fs = std::filesystem;
    fs::path rootPath(singlehit_file);
    std::string logFilename = "log_singleHit" + singlehit_file.substr(startPos + 16);
    fs::path logPath = rootPath.parent_path() / logFilename;
    logFilename = logPath.string();
    size_t pos = logFilename.rfind(".root");
    if (pos != std::string::npos) {
        logFilename.replace(pos, 5, ".log");
    } else {
        handleError("Invalid root file name format. Expected to find '.root' at the end of the filename.");
    } 
    std::cout << "log filename: " << logFilename << std::endl;
    std::ifstream logFile(logFilename);
    if (!logFile.is_open()) { 
        handleError("Unable to open log file: " + logFilename);
    }
        
    std::string line;
    std::regex pattern(R"(Opened input file \"(.*?)\")"); // Regex to match the required line
    std::smatch match;
    std::string orig_filename = "";
    while (std::getline(logFile, line)){
        if (std::regex_search(line, match, pattern)){
            orig_filename = match[1].str();
        }
    }
    if (orig_filename == ""){
        handleError("Could not locate filename in log file, exiting");
    }
    auto start = orig_filename.find("raw_");
    auto end = orig_filename.find("_datawriter");
    if (start != std::string::npos && end != std::string::npos && start < end) {
        return orig_filename.substr(start + 4, end - (start + 4));
    } else {
        handleError("Error getting file!");
    }
    return "";
}

std::string findRecoFilename(const std::string& inputFilename, const std::string& textFilePath) {
    // find ART ROOT filename from list of ART ROOT files that corresponds to SingleHit results ntuple
    size_t lastSlash_1 = inputFilename.find_last_of('/');
    std::string singlehit_filename = (lastSlash_1 != std::string::npos) ? inputFilename.substr(lastSlash_1 + 1) : inputFilename;
    std::string searchPattern = singlehit_filename.substr(17);
    std::ifstream inputFile(textFilePath);
    if (!inputFile.is_open()) {
        handleError("Could not find text file with reco filepaths!");
        return "";
    }
    std::string line;
    while (std::getline(inputFile, line)) {
        size_t lastSlash_2 = line.find_last_of("/\\");
        std::string filenameInPath = (lastSlash_2 == std::string::npos) ? line : line.substr(lastSlash_2 + 1);
        
        if (filenameInPath == searchPattern) {
            inputFile.close();
            std::cout << "Using ART ROOT file: " << line << std::endl;
            return line;
        }
    }
    inputFile.close();
    handleError("Could not locate reco file!");
    return ""; 
}

std::unordered_map<int, pdschannelmap> loadPDSChannelMap(const std::string& channelmap_path) {
    std::unordered_map<int, pdschannelmap> PDSChannelMap;
    std::ifstream channelmap_file(channelmap_path);

    if (!channelmap_file.is_open()) {
        handleError("Error opening channel map file!");
    }

    std::string line;
    bool firstLine = true;

    while (std::getline(channelmap_file, line)) {
        if (firstLine) {
            firstLine = false;
            continue;  // Skip the header line
        }

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
    return PDSChannelMap;
}

chargeData LoadClustersData(const std::string& singlehit_file){
    constexpr double DTS_sampling           = 16e-3; // DTS sampling, us
    constexpr double CRO_sampling           = 0.512; // CRO sampling, us

    // collect TPC clusters data from file
    constexpr int CRO_time_offset = -250; // to fix the larsoft bug
    constexpr const char* ClusterTreeName = "ClusterTree";
    constexpr const char* ClusterBranchName = "r2";
    TFile *clusters_file = TFile::Open(singlehit_file.c_str());
    if (!clusters_file || clusters_file->IsZombie()) {
        handleError("Error: Cannot open ROOT file with clusters data.");
    }
    std::string cluster_tree_name = std::string(ClusterTreeName);
    std::string cluster_branch_name = std::string(ClusterBranchName);
    //TDirectory *dir = (TDirectory*)clusters_file->Get(cluster_branch_name);
    TDirectory *dir = clusters_file->GetDirectory(cluster_branch_name.c_str());
    TTree *clusters_tree = nullptr;
    dir->GetObject(cluster_tree_name.c_str(), clusters_tree);
        
    std::vector<double>* clusters_Z = nullptr;
    clusters_tree->SetBranchAddress("Z", &clusters_Z);
    std::vector<double>* clusters_Y = nullptr; 
    clusters_tree->SetBranchAddress("Y", &clusters_Y);
    std::vector<double>* clusters_PeakTime = nullptr;
    clusters_tree->SetBranchAddress("PeakTime", &clusters_PeakTime);
    std::vector<int>* clusters_NOFTB_pre = nullptr;
    std::vector<int> clusters_NOFTB;
    clusters_tree->SetBranchAddress("NearOrFarToTheBeam", &clusters_NOFTB_pre);
    std::vector<double>* clusters_E = nullptr;
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
    std::cout << "Total events in clusters file = " << nEntries << std::endl;

    // setup vectors to store TPC clusters data
    chargeData charge_data;
    int lastLength = 0;
    double lastT0 = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        clusters_tree->GetEntry(i);
        lastT0 = ChargeT0;
        charge_data.ChargeT0s.push_back(ChargeT0 * DTS_sampling);
        
        // combine results from each event
        int vector_index = 0;
        for (int j = 0; j < clusters_NOFTB_pre->size(); j++){
            if (!((*NumberOfPlane0)[vector_index] && (*NumberOfPlane1)[vector_index])){
                charge_data.charge_nplanes.push_back(2);
            } else {
                charge_data.charge_nplanes.push_back(3);
            }
            clusters_NOFTB.push_back((*clusters_NOFTB_pre)[j]);
            charge_data.charge_T.push_back(ChargeT0 * DTS_sampling + (*clusters_PeakTime)[vector_index] * CRO_sampling + CRO_time_offset);
            charge_data.charge_E.push_back((*clusters_E)[vector_index]);
            charge_data.charge_PeakTime.push_back((*clusters_PeakTime)[vector_index]);
            charge_data.charge_Z.push_back((*clusters_Z)[vector_index]);
            charge_data.charge_Y.push_back((*clusters_Y)[vector_index]);
            
            // assign APA number to clusters, 
            // account for incorrectly formatted NearOrFarToTheBeam branch
            if ((*clusters_Z)[vector_index] < 230) {
                charge_data.charge_apa.push_back(clusters_NOFTB[vector_index] == 1 ? 1 : 3);
            } else {
                charge_data.charge_apa.push_back(clusters_NOFTB[vector_index] == 1 ? 2 : 4);
            }
            vector_index++;
        }
        charge_data.charge_event_indices.push_back({lastLength, static_cast<int>(charge_data.charge_Z.size()) - 1});
        lastLength = clusters_NOFTB.size();
    }
    clusters_file->Close(); 
    return charge_data;
}

void processWaveforms(const std::vector<raw::OpDetWaveform>& rawops, pdsData& pds_data, const std::unordered_map<int, pdschannelmap>& PDSChannelMap) {
    // Get waveform data and store in a struct for later use
    constexpr double DTS_sampling           = 16e-3; // DTS sampling, us
    constexpr double CRO_sampling        = 0.512; // CRO sampling, us
    constexpr int selftrig_length = 1024;
    constexpr int baseline_ticks  = 50; // how many ticks at start of wvfm to avg for baseline
    for (const auto& rawop : rawops) {
        auto const& ADC_Count = rawop.Waveform();
        if (ADC_Count.size() > selftrig_length) continue; // ignore streaming waveforms
        double baseline = std::accumulate(ADC_Count.begin(), ADC_Count.begin() + baseline_ticks, 0.0) / baseline_ticks;
        std::vector<double> waveform(ADC_Count.size());
        std::transform(ADC_Count.begin(), ADC_Count.end(), waveform.begin(), [baseline](double val) { return val - baseline; });
        double amplitude = std::abs(*std::min_element(ADC_Count.begin(), ADC_Count.end()) - baseline);
        double ts = rawop.TimeStamp() * DTS_sampling;
        int ch = rawop.ChannelNumber();

        if (PDSChannelMap.find(ch) != PDSChannelMap.end()) {
            const auto& values = PDSChannelMap.at(ch);
            pds_data.pdsAPA.push_back(values.APA);
            pds_data.pdsTPC.push_back(values.TPC);
            pds_data.pdsX.push_back(values.X);
            pds_data.pdsY.push_back(values.Y);
            pds_data.pdsZ.push_back(values.Z);
            pds_data.pdsChannel.push_back(ch);
            pds_data.pdsTimestamp.push_back(ts);
            pds_data.pdsAmplitude.push_back(amplitude);
            pds_data.pdsBaseline.push_back(baseline);
            pds_data.pdsWaveform.push_back(waveform);
        } else {
            std::cout << "Could not look up channel in channelmap, skipping this waveform." << std::endl;
        }
    }
}

void groupPDSHits1D(const pdsData& pds_data, std::vector<std::vector<int>>& grouped_ophit_indices_1D, double op_t_window) {
    // group PDS hits in 1D
    int op_hit_limit = 10;
    for (size_t i = 0; i < pds_data.pdsTimestamp.size(); i++) {
        std::vector<int> ophit_indices_1D = {static_cast<int>(i)};
        for (size_t j = 0; j < pds_data.pdsTimestamp.size(); j++) {
            if (j == i) continue;
            if (std::abs(pds_data.pdsTimestamp[j] - pds_data.pdsTimestamp[i]) < op_t_window) {
                ophit_indices_1D.push_back(static_cast<int>(j));
            }
        }
        if (ophit_indices_1D.size() < op_hit_limit) {
            grouped_ophit_indices_1D.push_back(ophit_indices_1D);
        }
    }
}

void groupPDSHits2D(const pdsData& pds_data, const std::vector<std::vector<int>>& grouped_ophit_indices_1D, std::vector<std::vector<int>>& grouped_ophit_indices_2D, double op_z_window, double op_y_window) {
    // group 1D groups in 2D, look for isolated ophit groups
    std::vector<int> indices_to_skip;
    int total_nhit1 = 0, total_nhit2 = 0, total_nhit3 = 0;

    for (const auto& ophit_indices_1D : grouped_ophit_indices_1D) {
        std::vector<int> ophit_indices_2D;
        for (size_t i = 0; i < ophit_indices_1D.size(); i++) {
            if (std::find(indices_to_skip.begin(), indices_to_skip.end(), ophit_indices_1D[i]) != indices_to_skip.end()) continue;
            ophit_indices_2D.push_back(ophit_indices_1D[i]);
            for (size_t k = i+1; k < ophit_indices_1D.size(); k++) {
                if (std::abs(pds_data.pdsZ[ophit_indices_1D[k]] - pds_data.pdsZ[ophit_indices_1D[i]]) < op_z_window * 0.5 &&
                    std::abs(pds_data.pdsY[ophit_indices_1D[k]] - pds_data.pdsY[ophit_indices_1D[i]]) < op_y_window * 0.5 &&
                    std::abs(pds_data.pdsX[ophit_indices_1D[k]] - pds_data.pdsX[ophit_indices_1D[i]]) < op_z_window * 0.5) {
                    ophit_indices_2D.push_back(ophit_indices_1D[k]);
                }
            }
            if (ophit_indices_2D.size() == 1) {
                grouped_ophit_indices_2D.push_back(ophit_indices_2D);
                indices_to_skip.push_back(ophit_indices_2D[0]);
                total_nhit1++;
            } else if (ophit_indices_2D.size() == 2) {
                if (std::abs(pds_data.pdsChannel[ophit_indices_2D[0]] - pds_data.pdsChannel[ophit_indices_2D[1]]) == 10) {
                    grouped_ophit_indices_2D.push_back(ophit_indices_2D);
                    indices_to_skip.push_back(ophit_indices_2D[0]);
                    indices_to_skip.push_back(ophit_indices_2D[1]);
                    total_nhit2++;
                }
            } else if (ophit_indices_2D.size() == 3) {
                bool cond1 = std::abs(pds_data.pdsChannel[ophit_indices_2D[0]] - pds_data.pdsChannel[ophit_indices_2D[1]]) == 10;
                bool cond2 = std::abs(pds_data.pdsChannel[ophit_indices_2D[1]] - pds_data.pdsChannel[ophit_indices_2D[2]]) == 10;
                bool cond3 = std::abs(pds_data.pdsChannel[ophit_indices_2D[0]] - pds_data.pdsChannel[ophit_indices_2D[2]]) == 10;
                if ((cond1 && cond2 && !cond3) || (cond1 && !cond2 && cond3) || (!cond1 && cond2 && cond3)) {
                    grouped_ophit_indices_2D.push_back(ophit_indices_2D);
                    indices_to_skip.push_back(ophit_indices_2D[0]);
                    indices_to_skip.push_back(ophit_indices_2D[1]);
                    indices_to_skip.push_back(ophit_indices_2D[2]);
                    total_nhit3++;
                }
            }
            ophit_indices_2D.clear();
        }
    }

    std::cout << "Total number of op hits = " << pds_data.pdsTimestamp.size() << std::endl;
    std::cout << "# of isolated op hit groups w/ 1 hit = " << total_nhit1 << std::endl;
    std::cout << "# of isolated op hit groups w/ 2 hit = " << total_nhit2 << std::endl;
    std::cout << "# of isolated op hit groups w/ 3 hit = " << total_nhit3 << std::endl;
}

void matchChargeAndLight(const chargeData& charge_data, const pdsData& pds_data, const std::vector<std::vector<int>>& grouped_ophit_indices_2D, 
    double q_z_window, double q_y_window, double upper_t_window, double lower_t_window, int evt_index, 
    std::vector<std::vector<int>>& charge_indices, std::vector<std::vector<int>>& light_indices) {
    // match isolated charge clysters with isolated op hit groups
    double upper_z_window = q_z_window / 2;
    double lower_z_window = q_z_window / 2;
    double upper_y_window = q_y_window / 2;
    double lower_y_window = q_y_window / 2;

    for (const auto& ophit_indices : grouped_ophit_indices_2D) {
        std::vector<int> charge_indices_matched;
        int NearBiSource = 0;
        for (int ophit_index : ophit_indices) {
            NearBiSource += (pds_data.pdsChannel[ophit_index] == 89 || pds_data.pdsChannel[ophit_index] == 9);
        }
        if (NearBiSource) {
            upper_z_window = q_z_window;
            lower_z_window = q_z_window;
            upper_y_window = q_y_window;
            lower_y_window = q_y_window;
        } else {
            upper_z_window = q_z_window / 2;
            lower_z_window = q_z_window / 2;
            upper_y_window = q_y_window / 2;
            lower_y_window = q_y_window / 2;
        }
        for (int k = charge_data.charge_event_indices[evt_index][0]; k < charge_data.charge_event_indices[evt_index][1]; k++) {
            int is_matched = 0;
            for (int ophit_index : ophit_indices) {
                is_matched += (charge_data.charge_T[k] < pds_data.pdsTimestamp[ophit_index] + upper_t_window) &&
                              (charge_data.charge_T[k] > pds_data.pdsTimestamp[ophit_index] - lower_t_window) &&
                              (charge_data.charge_Z[k] < pds_data.pdsZ[ophit_index] + upper_z_window) &&
                              (charge_data.charge_Z[k] > pds_data.pdsZ[ophit_index] - lower_z_window) &&
                              (charge_data.charge_Y[k] < pds_data.pdsY[ophit_index] + upper_y_window) &&
                              (charge_data.charge_Y[k] > pds_data.pdsY[ophit_index] - lower_y_window) &&
                              (charge_data.charge_apa[k] == pds_data.pdsAPA[ophit_index]);
            }
            if (is_matched) {
                charge_indices_matched.push_back(k);
            }
        }
        if (charge_indices_matched.size() == 0) {
            continue;
        }
        charge_indices.push_back(charge_indices_matched);
        light_indices.push_back(ophit_indices);
    }
}

void FillChargeTree(TFile* outputFile, std::string ChargeTreeName, chargeData charge_data) {
    // Fill charge tree with data of isolated charge clusters
    TTree *charge_tree = (TTree*) outputFile->Get(ChargeTreeName.c_str());
    //chargeData *results_charge = new chargeData();

    double charge_Z;
    double charge_Y;
    double charge_T;
    int charge_apa;
    int charge_nplanes;
    double charge_E;

    if (charge_tree){
        charge_tree->SetBranchAddress("z", &charge_Z);
        charge_tree->SetBranchAddress("y", &charge_Y);
        charge_tree->SetBranchAddress("t", &charge_T);
        charge_tree->SetBranchAddress("apa", &charge_apa);
        charge_tree->SetBranchAddress("nplanes", &charge_nplanes);
        charge_tree->SetBranchAddress("E", &charge_E);
    } else {
        charge_tree = new TTree(ChargeTreeName.c_str(), "Low energy isolated clusters");
        charge_tree->Branch("z", &charge_Z, "z/D");
        charge_tree->Branch("y", &charge_Y, "y/D");
        charge_tree->Branch("t", &charge_T, "t/D");
        charge_tree->Branch("apa", &charge_apa, "apa/I");
        charge_tree->Branch("nplanes", &charge_nplanes, "nplanes/I");
        charge_tree->Branch("E", &charge_E, "E/D");
    }

    //results_charge_events.offsets.push_back(offset_start + results_charge.pds_z.size()); 
    for (int charge_index; charge_index < charge_data.charge_Z.size(); charge_index++){
        charge_Z = charge_data.charge_Z[charge_index];
        charge_Y = charge_data.charge_Y[charge_index];
        charge_T = charge_data.charge_T[charge_index];
        charge_apa = charge_data.charge_apa[charge_index];
        charge_nplanes = charge_data.charge_nplanes[charge_index];
        charge_E = charge_data.charge_E[charge_index];
        charge_tree->Fill();
    }
    
    charge_tree->Write("", TObject::kOverwrite);
}

void FillChargeTree_Matched(TFile* outputFile, std::string ChargeTreeName, std::string offsetsTreeName, chargeData charge_data, std::vector<std::vector<int>> all_group_indices) {
    // Fill charge tree with data of isolated charge clusters
    //struct results_charge_events *charge_events = new results_charge_events();
    TTree *charge_tree = (TTree*) outputFile->Get(ChargeTreeName.c_str());
    TTree *offsets_tree = (TTree*) outputFile->Get(offsetsTreeName.c_str());

    int start_index;
    int stop_index;
    int groupsize;

    if (offsets_tree){ 
        offsets_tree->SetBranchAddress("start_index", &start_index);
        offsets_tree->SetBranchAddress("stop_index", &stop_index);
        offsets_tree->SetBranchAddress("groupsize", &groupsize);
    } else {
        offsets_tree = new TTree(offsetsTreeName.c_str(), "Starting index of each group");
        offsets_tree->Branch("start_index", &start_index, "start_index/I");
        offsets_tree->Branch("stop_index", &stop_index, "stop_index/I");
        offsets_tree->Branch("groupsize", &groupsize, "groupsize/I");
    }
   
    //chargeData *results_charge = new chargeData();

    double charge_Z;;
    double charge_Y;
    double charge_T;
    int charge_apa;
    int charge_nplanes;
    double charge_E;

    if (charge_tree){
        charge_tree->SetBranchAddress("z", &charge_Z);
        charge_tree->SetBranchAddress("y", &charge_Y);
        charge_tree->SetBranchAddress("t", &charge_T);
        charge_tree->SetBranchAddress("apa", &charge_apa);
        charge_tree->SetBranchAddress("nplanes", &charge_nplanes);
        charge_tree->SetBranchAddress("E", &charge_E);
    } else {
        charge_tree = new TTree(ChargeTreeName.c_str(), "Low energy isolated clusters");
        charge_tree->Branch("z", &charge_Z, "z/D");
        charge_tree->Branch("y", &charge_Y, "y/D");
        charge_tree->Branch("t", &charge_T, "t/D");
        charge_tree->Branch("apa", &charge_apa, "apa/I");
        charge_tree->Branch("nplanes", &charge_nplanes, "nplanes/I");
        charge_tree->Branch("E", &charge_E, "E/D");
    }

    int offset = charge_tree->GetEntries();
    for (const auto& group_indices : all_group_indices){
        start_index = offset;
        stop_index = offset + group_indices.size();
        groupsize = group_indices.size();
        offsets_tree->Fill();
        offset = offset + group_indices.size();
        for (int charge_index : group_indices){
            charge_Z = charge_data.charge_Z[charge_index];
            charge_Y = charge_data.charge_Y[charge_index];
            charge_T = charge_data.charge_T[charge_index];
            charge_apa = charge_data.charge_apa[charge_index];
            charge_nplanes = charge_data.charge_nplanes[charge_index];
            charge_E = charge_data.charge_E[charge_index];
            charge_tree->Fill();
        }
    }
    charge_tree->Write("", TObject::kOverwrite);
    offsets_tree->Write("", TObject::kOverwrite);
}

void FillLightTree(TFile* outputFile, pdsData pds_data, std::string LightTreeName, std::string WaveformsTreeName, 
    std::string OffsetsTreeName, std::vector<std::vector<int>> grouped_ophit_indices) {
    TTree *pds_tree = (TTree*) outputFile->Get(LightTreeName.c_str());
    TTree *offsets_tree = (TTree*) outputFile->Get(OffsetsTreeName.c_str());

    //pdsData *results_light = new pdsData();

    double pdsZ;
    double pdsY;
    double pdsX;
    double pdsTimestamp;
    int pdsChannel;
    int pdsAPA;
    double pdsAmplitude;
    double pdsBaseline;
    std::vector<double>* pdsWaveform = nullptr;

    if (pds_tree){
        pds_tree->SetBranchAddress("z", &pdsZ);
        pds_tree->SetBranchAddress("y", &pdsY);
        pds_tree->SetBranchAddress("x", &pdsX);
        pds_tree->SetBranchAddress("t", &pdsTimestamp);
        pds_tree->SetBranchAddress("ch", &pdsChannel);
        pds_tree->SetBranchAddress("apa", &pdsAPA);
        pds_tree->SetBranchAddress("amplitude", &pdsAmplitude);
        pds_tree->SetBranchAddress("baseline", &pdsBaseline);
        pds_tree->SetBranchAddress("waveform", &pdsWaveform);
    } else {
        pds_tree = new TTree(LightTreeName.c_str(), "Low energy charge-light matching results");
        pds_tree->Branch("z", &pdsZ, "z/D");
        pds_tree->Branch("y", &pdsY, "y/D");
        pds_tree->Branch("x", &pdsX, "x/D");
        pds_tree->Branch("t", &pdsTimestamp, "t/D");
        pds_tree->Branch("ch", &pdsChannel, "ch/I");
        pds_tree->Branch("apa", &pdsAPA, "apa/I");
        pds_tree->Branch("amplitude", &pdsAmplitude, "amplitude/D");
        pds_tree->Branch("baseline", &pdsBaseline, "baseline/D");
        pds_tree->Branch("waveform", &pdsWaveform);
    }
    //struct results_light_events *light_events = new results_light_events();
    int start_index;
    int stop_index;
    int groupsize;
    if (offsets_tree){ 
        offsets_tree->SetBranchAddress("start_index", &start_index);
        offsets_tree->SetBranchAddress("stop_index", &stop_index);
        offsets_tree->SetBranchAddress("groupsize", &groupsize);
    } else {
        offsets_tree = new TTree(OffsetsTreeName.c_str(), "Indices of each group");
        offsets_tree->Branch("start_index", &start_index, "start_index/I");
        offsets_tree->Branch("stop_index", &stop_index, "stop_index/I");
        offsets_tree->Branch("groupsize", &groupsize, "groupsize/I");
    }

    //results_light->pdsZ.clear();
    //results_light->pdsY.clear();
    //results_light->pdsX.clear();
    //results_light->pdsChannel.clear();
    //results_light->pdsTimestamp.clear();
    //results_light->pdsAPA.clear();
    //results_light->pdsBaseline.clear();
    //results_light->pdsAmplitude.clear();
    //light_events->offsets.clear();
    //light_events->groupsize.clear();

    int offset = pds_tree->GetEntries();
    
    for (const auto& ophit_indices : grouped_ophit_indices){
        start_index = offset;
        stop_index = offset + ophit_indices.size();
        groupsize = ophit_indices.size();
        offsets_tree->Fill();
        offset = offset + ophit_indices.size();
        //light_events->offsets.push_back(offset + pds_data.pdsZ.size());  // Store the starting index
        for (int ophit_index : ophit_indices){
            pdsZ = pds_data.pdsZ[ophit_index];
            pdsY = pds_data.pdsY[ophit_index];
            pdsX = pds_data.pdsX[ophit_index];
            pdsChannel = pds_data.pdsChannel[ophit_index];
            pdsTimestamp = pds_data.pdsTimestamp[ophit_index];
            pdsAPA = pds_data.pdsAPA[ophit_index];
            pdsBaseline = pds_data.pdsBaseline[ophit_index];
            pdsAmplitude = pds_data.pdsAmplitude[ophit_index];
            pdsWaveform = &pds_data.pdsWaveform[ophit_index];
            pds_tree->Fill();
        }
    }
    
    pds_tree->Write("", TObject::kOverwrite);
    offsets_tree->Write("", TObject::kOverwrite);
}

void CL_lowenergy_analysis(const int RUN, const int PROCESS_ID, const char* SINGLEHIT_FILELIST, 
    const char* RECO_FILELIST, const char* OUTFILE, const char* CHANNELMAP, const char* OPTAG, const char* OPTWINDOW, 
    const char* OPZWINDOW, const char* OPYWINDOW, const char* MATCHTWINDOWLOW, const char* MATCHTWINDOWHIGH, 
    const char* QZWINDOW, const char* QYWINDOW){
    
    //// Script for matching charge and light for low energy TPC clusters
    // Inputs: RUN - Run number
    //         PROCESS_ID - Index specifying file to process in filelist
    //         SINGLEHIT_FILELIST - Text file with paths to singlehit results
    //         RECO_FILELIST - Text file with paths to reconstructed ART ROOT files (found with rucio/metacat)             
    //         OUTFILE - Name of the output file
    //         CHANNELMAP - Path to channel map file
    //         OPTAG - Product tag for optical data in ART ROOT file
    //         OPTWINDOW - Time window in usec to group self-triggered waveforms
    //         OPZWINDOW - Z window in cm to group self-triggered waveforms
    //         OPYWINDOW - Y window in cm to group self-triggered waveforms
    //         MATCHTWINDOWLOW - Time window lower bound in usec to match charge and light (can be negative)
    //         MATCHTWINDOWHIGH - Time window upper bound in usec to match charge and light
    //         QZWINDOW - Z window in cm around isolated op hit to select nearby charge clusters
    //         QYWINDOW - Y window in cm around isolated op hit to select nearby charge clusters
    //// For standalone setup: setup xrootd, kx509, dunesw, and Gallery.
    //// Run `kx509;voms-proxy-init -noregen -rfc -voms dune:/dune/Role=Analysis`, and then
    //// precede the root command with `LD_PRELOAD=$XROOTD_LIB/libXrdPosixPreload.so`.

    // example inputs:
    //const int run = 28086;
    //const int process_offset = 0;
    //const int process_id = 0 + process_offset;
    //std::string singlehit_filelist = "/exp/dune/data/users/sfogarty/DUNEFD_39Ar/singlehit_lists/run28086_r248_20000files.list";
    //std::string reco_filelist = "/exp/dune/data/users/sfogarty/DUNEFD_39Ar/runlists/28086_filepaths.txt";
    //std::string OUTFILE = "PDHD_CLMATCHING_RESULTS_28086_jobnum0.root";
    //std::string channelmap_path = "/pnfs/dune/persistent/users/sfogarty/JUSTIN/PDHD_PDS_ChannelMap.csv";
    //std::string optag = "pdhddaphne:daq";
    //InputTag rawop_tag(optag);
    //const double op_t_window = 0.1;
    //const double op_z_window = 600.0;
    //const double op_y_window = 600.0;
    //const double match_t_window_low = -100;
    //const double match_t_window_high = 2200;
    //const double q_z_window = 55.0;
    //const double q_y_window = 12.0;

    // get inputs from command line
    const int run = RUN;
    const int process_id = PROCESS_ID;
    std::string singlehit_filelist = SINGLEHIT_FILELIST;
    std::string reco_filelist = RECO_FILELIST;
    std::string outfilename = OUTFILE;
    std::string channelmap_path = CHANNELMAP;
    std::string optag = OPTAG;
    InputTag rawop_tag(optag); // product tag for op waveforms
    const double op_t_window = std::atof(OPTWINDOW);
    const double op_z_window = std::atof(OPZWINDOW);
    const double op_y_window = std::atof(OPYWINDOW);
    const double match_t_window_low = std::atof(MATCHTWINDOWLOW);
    const double match_t_window_high = std::atof(MATCHTWINDOWHIGH);
    const double q_z_window = std::atof(QZWINDOW);
    const double q_y_window = std::atof(QYWINDOW);

    auto PDSChannelMap = loadPDSChannelMap(channelmap_path);
    if (PDSChannelMap.empty()){
        handleError("Failed to load the channel map.");
    }

    // get singlehit file corresponding to process id
    std::string singlehit_file = getIthLine(singlehit_filelist, process_id);
    std::filesystem::path fs_path(singlehit_file);
    singlehit_file = fs_path.filename().string();
    std::cout << "Processing singlehit file: " << singlehit_file << std::endl;
    // load clusters data
    chargeData charge_data = LoadClustersData(singlehit_file); // struct to store charge data
    // find ART ROOT filename corresponding to singlehit file
    //std::string identifier = findRecoFileIdentifier(singlehit_file);
    //std::cout << "identifier = " << identifier << std::endl;
    std::string reco_filename = findRecoFilename(singlehit_file, reco_filelist);
    std::vector<std::string> filenames = {reco_filename};
    
    

    // event loop
    int evt_index = 0;
    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) { 

        pdsData pds_data; // struct to store PDS data
        
        auto const& rawops = *ev.getValidHandle<vector<raw::OpDetWaveform>>(rawop_tag);
        processWaveforms(rawops, pds_data, PDSChannelMap);
    
        // do grouping of PDS hits in 1D
        std::vector<std::vector<int>> grouped_ophit_indices_1D;
        groupPDSHits1D(pds_data, grouped_ophit_indices_1D, op_t_window);
        std::cout << "Total dt clusters found = " << grouped_ophit_indices_1D.size() << std::endl;

        // group 1D groups in 2D, look for isolated ophit groups
        std::vector<std::vector<int>> grouped_ophit_indices_2D;
        groupPDSHits2D(pds_data, grouped_ophit_indices_1D, grouped_ophit_indices_2D, op_z_window, op_y_window);

        // match charge and light
        std::vector<std::vector<int>> charge_indices;
        std::vector<std::vector<int>> light_indices;
        matchChargeAndLight(charge_data, pds_data, grouped_ophit_indices_2D, q_z_window, q_y_window, match_t_window_high, match_t_window_low, evt_index, charge_indices, light_indices);

        std::cout << "# of candidate events found = " << charge_indices.size() << std::endl;

        // save results to ROOT file
        TFile *outputFile = new TFile(OUTFILE, "UPDATE");
        if (!outputFile) {
            handleError("Error: TFile not found.");
        }

        // save all clustered results, without matching
        FillChargeTree(outputFile, "IsolatedChargeClusters", charge_data);
        FillLightTree(outputFile, pds_data, "IsolatedPDSGroups_All", "PDSWaveforms_All", "Offsets_Light_All", grouped_ophit_indices_2D);
        
        // save matched results
        FillChargeTree_Matched(outputFile, "IsolatedChargeClusters_Matched", "Offsets_Charge_CLMatched", charge_data, charge_indices);
        FillLightTree(outputFile, pds_data, "IsolatedPDSGroups_Matched", "PDSWaveforms_Matched", "Offsets_Light_CLMatched", light_indices);
                    
    outputFile->Close();
    delete outputFile;
    evt_index++;
    } // end gallery loop 
}

