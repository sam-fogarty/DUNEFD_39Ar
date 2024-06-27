#include <TFile.h>
#include <TTree.h>
#include <iostream>

void read_root() {
    // Open the ROOT file
    TFile *file = TFile::Open("pdhd_daphne_decodermodule.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    TDirectory* dir = file->GetDirectory("pdhddaphne");
    if (!dir) {
       std::cerr << "Error: Directory pdhddaphne not found in file " << std::endl;
       file->Close();
       return;
    }
    TTree* tree = dynamic_cast<TTree*>(dir->Get("WaveformTree"));
    if (!tree) {
        std::cerr << "Error retrieving TTree" << std::endl;
        return;
    }
    file->Print();
    // mean baseline
    Int_t Baseline;
    tree->SetBranchAddress("Baseline", &Baseline);
    
    // adc counts for a waveform
    Short_t adc_value[1024];
    tree->SetBranchAddress("adc_channel", adc_value);
    
    // timestamp
    ULong64_t TimeStamp;
    tree->SetBranchAddress("TimeStamp", &TimeStamp);
    
    // frame timestamp
    ULong64_t FrameTimestamp;
    tree->SetBranchAddress("FrameTimestamp", &FrameTimestamp);

    // daphne channel
    Int_t DaphneChannel;
    tree->SetBranchAddress("DaphneChannel", &DaphneChannel);

    // offline channel
    Int_t OfflineChannel;
    tree->SetBranchAddress("OfflineChannel", &OfflineChannel);

    // TriggerSampleValue
    Int_t TriggerSampleValue;
    tree->SetBranchAddress("TriggerSampleValue", &TriggerSampleValue);

    // Threshold
    Int_t Threshold;
    tree->SetBranchAddress("Threshold", &Threshold);

    // Crate
    Int_t Crate;
    tree->SetBranchAddress("Crate", &Crate);

    // Slot
    Int_t Slot;
    tree->SetBranchAddress("Slot", &Slot);
   
    // Run
    Int_t Run;
    tree->SetBranchAddress("Run", &Run);

    // Event
    Int_t Event;
    tree->SetBranchAddress("Event", &Event);

    // TriggerNumber
    Int_t TriggerNumber;
    tree->SetBranchAddress("TriggerNumber", &TriggerNumber);


    //TCanvas* canvas = new TCanvas("canvas", "ADC Values", 800, 600);
    //TGraph* graph = new TGraph(); 
 
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Number of entries in tree = " << nEntries << std::endl; 
    for (Long64_t i = 0; i < 3; ++i) {
       tree->GetEntry(i);
       std::cout << "Baseline = " << Baseline << std::endl;
       std::cout << "TimeStamp = " << TimeStamp << std::endl;   
       std::cout << "FrameTimestamp = " << FrameTimestamp << std::endl;
       std::cout << "DaphneChannel = " << DaphneChannel << std::endl;
       std::cout << "OfflineChannel = " << OfflineChannel << std::endl;
       std::cout << "TriggerSampleValue = " << TriggerSampleValue << std::endl;
       std::cout << "Threshold = " << Threshold << std::endl;
       std::cout << "Crate = " << Crate << std::endl;
       std::cout << "Slot = " << Slot << std::endl;
       std::cout << "Run = " << Run << std::endl;
       std::cout << "Event = " << Event << std::endl;
       std::cout << "TriggerNumber = " << TriggerNumber << std::endl;                
       std::cout << std::endl;      
    }

    //tree->Print();
    //for (Long64_t i = 1; i < 2; ++i) {
    //    tree->GetEntry(i);
    //    for (int j = 0; j < 1024; ++j) {
             //std::cout << adc_value[j] << " ";
             //hist->Fill(j, adc_value[j]);
             //graph->SetPoint(j, j, adc_value[j] - Baseline); 

    //    }
        //for (const auto& value : *adcs) {
        //     std::cout << value << " ";
        //}
        //std::cout << std::endl;

    //}
    //graph->Draw("AL");
    

    // Close the file
    file->Close();
}

