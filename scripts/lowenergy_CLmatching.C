R__ADD_INCLUDE_PATH("gallery/Event.h")

#include "tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace art;

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

void lowenergy_CLmatching(){
string filelist = "runlist.txt";
std::vector<std::string> filenames = readFilePaths(filelist);
std::string line;
//for (const auto& path : paths) {
//    std::cout << path << std::endl;
//    }
string run = "028850";
int nfiles = 1;
//string singlehit_runlist="/pnfs/dune/scratch/users/lavaut/03684/1/run28850.list";
string singlehit_runlist="/exp/dune/app/users/lavaut/workflow/run28086_r248_10000files.list";
string optag = "pdhddaphne:daq";
InputTag rawop_tag(optag);

float dt    = 16e-3; // us
float dtcrp = 0.512;
int selftrig_length = 1024;
int evt = 0;
int Nlimit=1;

double pdsTimestamp;
int pdsChannel;
size_t start;
size_t end;

//std::ifstream singlehit_files(singlehit_runlist); 
//if (!singlehit_files) {
//   cout << "Unable to open runlist file!" << endl;
//}
std::string identifier;
for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
   
   //start = filenames[evt].find("raw_");
   //end = filenames[evt].find("_datawriter");
   //if (start != std::string::npos && end != std::string::npos && start < end) {
        // extract the substring between "raw_" and "_datawriter"
   //     identifier = filenames[evt].substr(start + 4, end - (start + 4));
   //     } else {
   //           std::cout << "Error getting charge cluster file.!" << std::endl;
   //           continue;
   //     }
   //string singlehit_file = "";
   //while (std::getline(singlehit_files, line)) { 
   //     if (line.find(identifier) != std::string::npos) {
   //         singlehit_file = line;
   //         break;
   //     }
   // }
   //if (singlehit_file == ""){
   //        cout << "Cannot locate singlehit clusters file." << endl;
   //        continue;
   //} else {
   //        cout << "Using singlehit clusters file: " << singlehit_file << endl;
   //}
   //if(evt > Nlimit){continue;}
   //auto const& rawops = *ev.getValidHandle<vector<raw::OpDetWaveform>>(rawop_tag);
   //for (int i = 0; i < rawops.size(); i++){
   //   auto const& ADC_Count = rawops[i].Waveform();
   //   if (ADC_Count.size() > selftrig_length) {continue;} // ignore streaming wvfms
   //   pdsTimestamp = rawops[i].TimeStamp() * dt;
   //   pdsChannel = rawops[i].ChannelNumber();
      
   //}
   //evt++;

   //exit 1;
}
 
