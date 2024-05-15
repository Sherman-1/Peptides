#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>  // Include for std::cerr and std::endl

// Function to split string by whitespace
std::vector<std::string> split(const std::string& s) {
    std::istringstream stream(s);
    std::string token;
    std::vector<std::string> tokens;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file_name>" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];
    std::ifstream infile(inputFileName);
    std::ofstream outfile("secondary_coord.tsv");
    std::string line;

    if (!infile.is_open() || !outfile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    while (std::getline(infile, line)) {
        std::vector<std::string> tokens = split(line);

        // Assuming the structure of each line matches the provided Python script logic
        if (tokens.size() < 16) continue; // Safety check for index access
        std::string start = tokens[12];
        std::string end = tokens[15];
        std::string chain_id = tokens[14];
        std::string secondary_structure = tokens[1];

        if (secondary_structure[0] == 'H') {
            secondary_structure = "H";
        } else if (secondary_structure.substr(0, 2) == "ST") {
            secondary_structure = "S";
        } else {
            secondary_structure = "O";
        }

        int start_num = std::stoi(start);
        int end_num = std::stoi(end);
        for (int res_number = start_num; res_number <= end_num; res_number++) {
            std::cout << chain_id << '\t' << res_number << '\t' << secondary_structure << '\n';
        }
    }

    infile.close();
    outfile.close();
    return 0;
}
