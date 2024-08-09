#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <set>
#include <vector>

namespace fs = std::filesystem;

void checkResidueGaps(const std::string& directoryPath) {
    // Iterate over all files in the directory
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.path().extension() == ".pdb") {
            std::ifstream pdbFile(entry.path());
            std::string line;
            std::unordered_map<char, std::set<int>> chainResidues;

            while (getline(pdbFile, line)) {
                if (line.substr(0, 4) == "ATOM") {
                    int residueNumber = std::stoi(line.substr(22, 4));
                    char chainID = line[21];

                    chainResidues[chainID].insert(residueNumber);
                }
            }

            pdbFile.close();

            // Check for gaps in each chain
            for (const auto& [chainID, residues] : chainResidues) {
                int previousResidue = -1000; // Start with an impossible residue number
                for (int residue : residues) {
                    if (previousResidue != -1000 && residue - previousResidue > 3) {
                        std::cout << "Gap found in chain " << chainID << " in file " << entry.path().filename()
                                  << ": between residues " << previousResidue << " and " << residue << std::endl;
                    }
                    previousResidue = residue;
                }
            }
        }
    }
}

int main() {
    std::string directoryPath = "/Users/k-archi/Desktop/Peptides/pdbs"; // Change to your directory path
    checkResidueGaps(directoryPath);

    return 0;
}
