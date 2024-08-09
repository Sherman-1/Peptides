#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <regex>
#include <string>
#include <iostream>
#include <cstring>


char mapSSType(const std::string& rawType) {

    /* 

     As provided by the official dssp documentation :
     https://github.com/PDB-REDO/dssp/blob/trunk/doc/mkdssp.md
    
     DSSP Code      mmCIF Code    Description
     H         HELX_RH_AL_P   Alphahelix
     B             STRN       Betabridge
     E             STRN         Strand
     G         HELX_RH_3T_P     Helix_3
     I         HELX_RH_PI_P     Helix_5
     P         HELX_LH_PP_P   Helix_PPII
     T          TURN_TY1_P       Turn
     S             BEND          Bend
     ' ' (space)      OTHER          Loop

    */

    if (rawType.find("HELX") != std::string::npos) {
        return 'H';
    } else if (rawType.find("STRN") != std::string::npos) {
        return 'S';
    } else if (rawType.find("TURN") != std::string::npos) {
        return 'T';
    } else if (rawType.find("BEND") != std::string::npos) {
        return 'B';
    }

    return 'L'; // Default to loop if no match is found

}

std::vector<std::string> splitLine(const std::string& line) {

    std::regex non_whitespaces("\\S+");
    std::sregex_iterator words_begin = std::sregex_iterator(line.begin(), line.end(), non_whitespaces);
    std::sregex_iterator words_end = std::sregex_iterator();

    std::vector<std::string> words;

    for (std::sregex_iterator i = words_begin; i != words_end; ++i) {

        words.push_back((*i).str());

    }

    return words;
}

std::unordered_map<std::string, std::string> createHashMap(const std::vector<std::string>& keys, const std::vector<std::string>& values) {
    std::unordered_map<std::string, std::string> hashmap;

    if (keys.size() != values.size()) {
        std::cerr << "Error: Size of keys and values lists must be equal" << std::endl;

        std::cerr << "Keys size : " << keys.size() << std::endl;
        std::cerr << "Values size : " << values.size() << std::endl;

        std::cerr << "Keys : " << std::endl;
        for (size_t i = 0; i < keys.size(); ++i) {
            std::cerr << keys[i] << std::endl;
        }

        std::cerr << "Values : " << std::endl;
        for (size_t i = 0; i < values.size(); ++i) {
            std::cerr << values[i] << std::endl;
        }

        return hashmap; // Return an empty hashmap if sizes don't match
    }

    for (size_t i = 0; i < keys.size(); ++i) {
        hashmap[keys[i]] = values[i];
    }
    return hashmap;
}

int parseCIF(std::ifstream& cifFile, std::ofstream& outFile) {

    std::string line;
    std::vector<std::string> keyList;

    // Read until we find the _struct_conf section
    // When the first line of the section is found
    // start to populate the keyList with the headers ( one header per line 
    // before the data actually start )

    while (std::getline(cifFile, line)) {
        if (line.find("_struct_conf.") != std::string::npos) {

            size_t pos = line.find("_struct_conf.");
            std::string extractedString = line.substr(pos + strlen("_struct_conf."));
            size_t lastChar = extractedString.find_last_not_of(" \n\r\t");
            if (lastChar != std::string::npos) {
                extractedString.erase(lastChar + 1);
            } else {
                extractedString.clear();
            }
            keyList.push_back(extractedString);
            break;
        }
    }

    // We reached the first header, now we populate the keyList
    // with the headers until we reach the data
    while (std::getline(cifFile, line)) {
        if (line.empty() ) {
            
            // Skip potential empty lines
            continue;

        } else if (line.find("_struct_conf.") != std::string::npos) {

            size_t pos = line.find("_struct_conf.");
            std::string extractedString = line.substr(pos + strlen("_struct_conf."));
            size_t lastChar = extractedString.find_last_not_of(" \n\r\t");
            if (lastChar != std::string::npos) {
                extractedString.erase(lastChar + 1);
            } else {
                extractedString.clear();
            }
            keyList.push_back(extractedString);

        } else {

            std::string buffline = line;
            break;
        }

    }

    // Now we have the keyList populated with the headers
    // But we also have the first line of the data ! 
    // Parse it before re entering the loop to parse the rest of the data

    std::vector<std::string> tokens = splitLine(line);
    std::unordered_map<std::string, std::string> hashmap = createHashMap(keyList, tokens);


    int start = std::stoi(hashmap["beg_auth_seq_id"]);
    int end = std::stoi(hashmap["end_auth_seq_id"]);
    std::string chain_ID = hashmap["beg_auth_asym_id"];

    char one_letter_SS = mapSSType(hashmap["conf_type_id"]);

    // chain_id << '\t' << res_number << '\t' << secondary_structure << '\n';
    for (int i = start; i <= end; ++i) {
        outFile << chain_ID << '\t' << i << '\t' << one_letter_SS << '\n';
    }

    // Now we parse the rest of the data
    while (std::getline(cifFile, line)) {

        if (line.empty()) {
            continue;
        }

        if (line.find('#') != std::string::npos || line.front() == '_') {
            break;
        }

        std::vector<std::string> tokens = splitLine(line);
        std::unordered_map<std::string, std::string> hashmap = createHashMap(keyList, tokens);

        int start = std::stoi(hashmap["beg_auth_seq_id"]);
        int end = std::stoi(hashmap["end_auth_seq_id"]);
        std::string chain_ID = hashmap["beg_auth_asym_id"];

        char one_letter_SS = mapSSType(hashmap["conf_type_id"]);
        // chain_id << '\t' << res_number << '\t' << secondary_structure << '\n';
        for (int i = start; i <= end; ++i) {
            outFile << chain_ID << '\t' << i << '\t' << one_letter_SS << '\n';
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {


    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        return 1;
    }

    std::ifstream cifFile(argv[1]);
    if (!cifFile) {
        std::cerr << "Error opening input file: " << argv[1] << std::endl;
        return 1;
    }

    std::ofstream outFile(argv[2]);
    if (!outFile) {
        std::cerr << "Error opening output file: " << argv[2] << std::endl;
        return 1;
    }

    parseCIF(cifFile, outFile);

    cifFile.close();
    outFile.close();
    return 0;
}

