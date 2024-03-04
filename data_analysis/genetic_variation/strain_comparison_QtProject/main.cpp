/*
 * File: main.cpp
 * --------------
 * Blank C++ project configured to use Stanford cslib and Qt
 */

#include "console.h"
#include "simpio.h"
#include <iostream>
#include <fstream>
#include "SimpleTest.h"
#include "map.h"
#include "set.h"
#include "strlib.h"

using namespace std;

void analyzeAlignment(string inputFile, string outputFile);
void baseByBaseAlignment(string inputFile, string outputFile, string refStrain);
void fillGenesMap(Map<string, Map<string, string>>& map, string path);
void compareSequences(string seq1, string seq2, Vector<string>& vec);
string compareSequencesBaseByBase(string refSeq, string seq2);
Vector<Vector<string>> mapAlignmentAnalysis(Map<string, Map<string, string>>& map);
void testing();
void testing2();
void testing3(string inputFile);
string testing4(string teststring);
void generateSubAlignmentFiles(string inputFile);
void generateSubAlignmentFilesHandleDups(string inputFile);
Vector<string> geneList;

int main()
{
    // FOR ANALYZING COMPRESSED ALIGNMENTS:

   // analyzeAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-02_fliS_compressedAlignment_nt.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-02_fliS_compressedAlignment_nt_outout.txt");
   // analyzeAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-02_fliS_compresseddAlignment_aa.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-02_fliS_compressedAlignment_aa_output.txt");

   // analyzeAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-10_fliC_compressedAlignment_nt.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-11_fliC_seqID_nt_output.txt");
   // analyzeAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-10_fliC_compressedAlignment_aa.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-11_fliC_seqID_aa_output.txt");


    //analyzeAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-06-28_nt_compressed_flag+hk_promoter_alignments_all.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-06-28_nt_compressed_flag+hk_promoter_alignments_all_output.txt");
    //analyzeAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-09-26_nt_compressed_lac+phostrans+transporters_promoter_alignments_all.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-09-26_nt_compressed_lac+phostrans+transporters_promoter_alignments_all_output.txt");
    //generateSubAlignmentFiles("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-09-26_nt_compressed_lac+phostrans+transporters_promoter_alignments_all_output.txt");

    // 2023-11-08:
    //analyzeAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-11-08_nt_compressed_ALL_promoter_alignments.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-11-08_nt_compressed_ALL_promoter_alignments_output.txt");
    // generateSubAlignmentFiles("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-11-08_nt_compressed_ALL_promoter_alignments_output.txt");

    // FOR MAKING A SINGLE SEQUENCE NT-BY-NT ALIGNMENT:
    // baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-02_fliS_compressedAlignment_nt.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-10_fliS_compressedAlignmentBBB_nt_output.txt", "MG1655");

    //baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-10_fliC_compressedAlignment_nt.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-10_fliC_compressedAlignmentBBB_nt_QOXR01.1ref_output.txt", "QOXR01.1");


    // 2023-11-10
    Vector<string> geneList = { "lacZ", "lacY", "lacA", "lacI", "cyaA", "hns", "dksA", "crp", "dnaA", "dnaB", "dnaC", "dnaE", "dnaG", "dnaN", "dnaQ", "polA", "polB", "ssb", "gyrA", "gyrB", "parC", "parE", "ligA", "rpoA", "rpoZ", "rpoB", "rpoC", "rpoD", "rpsR", "rpsB", "rplS", "rpsO", "rplO", "rpmJ", "rpsL", "rpmF", "rplL", "rpsI", "rpmC", "rplF", "rpsF", "rplX", "rplC", "rplU", "rpsS", "rpsC", "rplK", "rpsP", "rplQ", "rpsM", "rpmG", "rplM", "rpsJ", "rpmD", "rplI", "rpsG", "rplY", "rplD", "rpmA", "rpsT", "rpsD", "rplV", "rplA", "ykgO", "rpsQ", "rplR", "rpsN", "rpmH", "rplN", "rpsK", "rpmE", "rplJ", "rpmB", "rplE", "rplP", "rpsU", "rpsE", "rplW", "rplB", "ykgM", "rpmI", "rpsA", "rplT" };
   // Vector<string> geneList = { 'string1'};

    for (string gene : geneList) {
        baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-11-08_allPromotersTest/BBBAlignments/" + gene + "_compressedAlignment.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-11-08_allPromotersTest/BBBAlignments/" + gene + "_compressedAlignmentBBB_nt_MGref_output.txt", "MG1655");
    }

//    // NT:
//     baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_fliD_compressedAlignment_nt.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_fliD_compressedAlignmentBBB_nt_MGref_output.txt", "MG1655");
//     baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_tap_compressedAlignment_nt.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_tap_compressedAlignmentBBB_nt_MGref_output.txt", "MG1655");
//     baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_trg_compressedAlignment_nt.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_trg_compressedAlignmentBBB_nt_MGref_output.txt", "MG1655");
//    // AA:
//     baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_fliD_compressedAlignment_aa.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_fliD_compressedAlignmentBBB_aa_MGref_output.txt", "MG1655");
//     baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_tap_compressedAlignment_aa.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_tap_compressedAlignmentBBB_aa_MGref_output.txt", "MG1655");
//     baseByBaseAlignment("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_trg_compressedAlignment_aa.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-26_trg_compressedAlignmentBBB_aa_MGref_output.txt", "MG1655");
//    // testing3();

    // TESTING AS I DEVELOP THE SEPARATE-FILE GENERATION:
//testing3("/Users/rlporter/Documents/gut_native_ec_motility/testalignment.txt");
  //  generateSubAlignmentFilesHandleDups("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2022-12-29_nt_alignment_output_all.txt");
  //  generateSubAlignmentFilesHandleDups("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-04_aa_compressed_alignments_output_repeat12.29.txt");

   // generateSubAlignmentFiles("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2022-12-29_nt_alignment_output_all.txt");
   // generateSubAlignmentFiles("/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2023-01-04_aa_compressed_alignments_output_repeat12.29.txt");

    //cout << output << endl;
    return 0;
}

/* This is a wrapper function around the functions fillGenesMap, mapAlignmentAnalysis, and
 * compareSequences. It accepts as input a compressed CLW file, and uses its contents to fill in a
 * genemap, analyze all alignemnts, and output a file with results from each pairwise analysis for
 * each gene. This is written to the specified text file, but the function itself does not return
 * anything.
 */
void analyzeAlignment(string inputFile, string outputFile) {
    // read in compressed file & use to construct a map:
    Map<string, Map<string, string>> genemap;
    fillGenesMap(genemap, inputFile);

    // processes alignments within genemap and output to vector:
    Vector<Vector<string>> vec = mapAlignmentAnalysis(genemap);

    // output vector to outfile:
    ofstream out;
    out.open(outputFile);
    for (Vector<string> vector : vec) {
        for (string str : vector) {
            out << str << "," ;
        }
        out << endl;
    }
    out.close();
}

/* This function accepts an input file of the compressedCLW format, which should
 * only contain a single gene. A map of strains to sequences is created and each
 * sequence is parsed through, and compared base by base to the reference stain.
 * If the sequence contains a gap, the outout shows '-' for that position.
 * If the bases match, they are marked with a '0', if they vary, they are marked
 * with a '1'.
 *
 * ******************************************************************************* *
 * NOTE: I need to figure out what to do with duplicates in this procedure, and if
 * it makes more sense to handle them here/in python steps, etc;
 * ******************************************************************************* *
 */
void baseByBaseAlignment(string inputFile, string outputFile, string refStrain) {
    string gene;
    string refSeq;
    string seq2;
    string outStr;
    Vector<string> outPairVec;
    Vector<Vector<string>> outVec;
    // read in compressed file and use to construct a map (should have only 1 gene!)
    Map<string, Map<string, string>> genemap;
    fillGenesMap(genemap, inputFile);
    if (genemap.size() != 1) {
        cout << "WARNING: input file may contain more than one gene." << endl;
        // error?
    }

    // find the single key for the main map:
    gene = genemap.firstKey();

    // compare strains & build vector with outputs:
    // check if the value map for the gene contains the reference strain, and output
    // warning message if not:
    if (!genemap[gene].containsKey(refStrain)) {
        cout << "WARNING: input file does not contain reference strain." << endl;
        // error?
    } else {
        // obtain reference sequence, and compare all sequences included in map against this. The
        // output of this is written along with the strain name (key) to the output file.
        refSeq = genemap[gene].get(refStrain);
        for (string key : genemap[gene].keys()) {
            seq2 = genemap[gene].get(key);
            outStr = compareSequencesBaseByBase(refSeq, seq2);
            outPairVec = {key, outStr};
            outVec.add(outPairVec);
        }

        // write outVec to output file:
        ofstream out;
        out.open(outputFile);
        for (Vector<string> vector : outVec) {
         //   for (string str : vector) {
         //       out << str << "," ;
         //   }

            out << vector[0] << ",";
            for (char i : vector[1]) {
                out << i << "," ;
            }
            out << endl;
        }
        out.close();
    }
}

/* This function will read in from a text file and create a map with keys
 * corresponding to strain names, and values corresponding to that strain's
 * sequence for the provided gene.
 * fillMap accepts as input an empty map by reference and a string to the
 * file to read in and process. The file is read in, iterated through character
 * by character, and the map is filled with gene/map pairings, where each map
 * for a gene contains strain/sequence pairings.
 */
void fillGenesMap(Map<string, Map<string, string>>& map, string path) {
    //initialize temp variables:
    string gene;
    string strain;
    string sequence;
    char current;

    //set up ifstream to read file:
    ifstream inputFile;
    inputFile.open(path);

    // while the file is nonempty, parse through each character to grab genes,
    // strain names, and sequences
    while ( inputFile ) {
        current = inputFile.get();

        // if exit delimiter is hit, terminate loop:
        if (current == '!') {
            break;
        }
        // get gene name & add as a key in main map:   // <-- 'find gene' function
        if (current == '@') {
            gene = "";                          // reinitialize gene to empty string
            current = inputFile.get();          // move to next character (first character in gene name)
            while ( current != '@' ) {
                gene += charToString(current);
                current = inputFile.get();
            }
            cout << endl << "Gene is : " << gene;

            // if gene is already in map:
            if (map.containsKey(gene)) {
                string geneAppended = gene + "_1";
                int i = 2;
                while (map.containsKey(geneAppended)) {
                    //gene = geneAppended.substr(0, geneAppended.length() - 1);
                    //geneAppended = gene + integerToString(i);
                    geneAppended = gene + "_" + integerToString(i);
                    i ++;
                }
                cout << "WARNING: " << gene << " already exists in map. Will be added as " << geneAppended << " instead." << endl;
                gene = geneAppended;
            }
            if (!map.containsKey(gene)) {
                map[gene];
            } else {
                cout << "WARNING: " << gene << " was not successfully added to map." << endl;
            }
            current = inputFile.get();      // this will set current to first character after '@'
        }

        // continue parsing until ':' to get strain name:
        strain = "";
        //current = inputFile.get();
        while ( current != ':' ) {
            strain += charToString(current);
            current = inputFile.get();
        } // at end of while, current = ':'

        // continue parsing until '*' to get sequence:
        sequence = "";
        current = inputFile.get();
        while ( current != '*' ) {
            sequence += charToString(current);
            current = inputFile.get();
        } // at end of while, current = '*'

        // check if strain is already in map & update strain name if so:
        if (map[gene].containsKey(strain)) {
            string strainAppended = strain + "_1";
            int i = 2;
            while (map[gene].containsKey(strainAppended)) {
                //strain = strainAppended.substr(0, strainAppended.length() - 1);
                //strainAppended = strain + integerToString(i);
                strainAppended = strain + "_" + integerToString(i);
                i ++;
            }
            cout << "WARNING: " << strain << " already exists in " << gene << " map. Added as " << strainAppended << " instead." << endl;
            strain = strainAppended;
        }
        // add strain/sequence pairing to map for gene:
        if (! map[gene].containsKey(strain)) {
            map[gene][strain] = sequence;
        } else {
            cout << "WARNING: " << strain << " was not successfully added to map for gene " << gene << " ." << endl;
        }



        // move to next character (after '*'). This will either be the start
        // of a new gene (if '@'), or the start of another strain name.
        //current = inputFile.get();
    }
    return;
}

/* This function accepts an alignment map by reference, which should contain genes as
 * keys, and maps of strains to sequences as values. For each gene in the main map,
 * all pairwise alignments between strains are made, without repeats.
 *
 */
Vector<Vector<string>> mapAlignmentAnalysis(Map<string, Map<string, string>>& map){
    // initialize variables
    string strain1;
    string strain2;
    string seq1;
    string seq2;
    Vector<string> rowEntry = {"", "", "", "", "", ""};
    Vector<Vector<string>> output = {};

    // for each gene, get a list of strains, update rowEntry template, and assign a pointer
    // to the map corresponding to the gene:
    for (string gene : map.keys() ){
        rowEntry[0] = gene;
        cout << "Gene is " << gene << endl;
        Map<string, string>* geneMap = &map[gene];
        Vector<string> strains = geneMap->keys();
        // for each unique strain combination, update RowEntry and compareSequences:
        int strainSize = strains.size();
        for (int i = 0; i < strainSize; i++) {
            for (int j = i; j < strainSize; j++) {
                strain1 = strains[i];
                strain2 = strains[j];
                rowEntry[1] = strain1;
                rowEntry[2] = strain2;
                seq1 = geneMap->get(strain1);
                seq2 = geneMap->get(strain2);
                if (seq1.length() != seq2.length()) {
                    error("WARNING: Sequences for gene " + gene + " are not all of the same length!");
                } else {
                    // add comparison output to rowEntry, and add this to main output vec:
                    compareSequences(seq1, seq2, rowEntry);
                    output.add(rowEntry);
                }
            }
        }
    }
    return output;
}
//Vector<string> analyzeAlignment(Map<string, Map<string, string>>& genemap, string outpath, ) {

//}

/* This function accepts as input two sequences and compares them on a character
 * by character basis, calculating the number of times one has a hyphen and the other
 * doesn't (gaps) and the number of times each has a non-hyphen character, but these
 * don't match (snps). The total number of times both have hyphens is also tabulated
 * to find the specific length of each pairwise alignment.
 *
*/
void compareSequences(string seq1, string seq2, Vector<string>& vec) {
    // initialize variables:
    int subLength = 0;
    int snps = 0;
    int gaps = 0;
    //string strain1 = "";              // <-- these can actually just be inputs, no new var needed
    //string strain2 = "";
    //string gene = "";

    // first, check that the two sequences have the same length:
//    if (seq1.length() != seq2.length()) error("WARNING: Sequence lengths do not match for gene.");   // <-- actually, makes sense to do this before compareSeq call!
    for (int i = 0; i < seq1.length(); i ++ ) {
        // if there's a gap in both:
        if (seq1[i] == seq2[i] && seq1[i] == '-') subLength ++;
        // if sequence characters don't match:
        if (seq1[i] != seq2[i]) {
            if (seq1[i] == '-' || seq2[i] == '-') gaps ++;          // if char aren't the same, and either has '-' --> gaps +1
            else snps ++;                                           // otherwise --> snps +1
        }
    }
    vec[3] = integerToString(gaps);
    vec[4] = integerToString(snps);
    vec[5] = integerToString(seq1.length() - subLength);
}

/* This function accepts two string sequences and outputs a single string representing similarities of the
 * two inputs. If the second sequence has a gap ('-'), the output is appended with '-'. Otherwise, the
 * characcters at this position in both strains are compared and if they match, the output is appended with
 * a '0'. If they do not match, the output is appended with a '1':
 *
 *  ref:     seq:          char:
 *  gap     gap             -
 *  gap     non-gap         1
 *  base    gap             -
 *  base    right base      0
 *  base    wrong base      1
 *
 */
string compareSequencesBaseByBase(string refSeq, string seq2) {
    if (refSeq.length() != seq2.length()) cout << "WARNING: sequence lengths do not match!" << endl;
    // error?
    else {
        string output = "";
        // THIS WAS THE ORIGINAL CODE, WHERE ALL GAPS ARE CONSIDERED THE SAME:
//        for (int i = 0; i < refSeq.length(); i ++) {
//            if (seq2[i] == '-') output += "-";
//            else if (refSeq[i] == seq2[i]) output += "0";
//            else output += "1";
//        }
        // NEW CODE WHERE CORRECT/INCORRECT GAPS ARE DIFFERENTIATED:
        for (int i = 0; i < refSeq.length(); i++) {
            if ((seq2[i] == '-') && (refSeq[i] == '-')) output += "-";         // right gap
            else if (seq2[i] == '-') output += "&";                            // wrong gap
            else if (seq2[i] == refSeq[i]) output += "0";                      // right monomer
            else output += "1";                                                // wrong monomer
            }


        if (output.length() != refSeq.length()) cout << "WARNING: sequence and output lengths do not match!" << endl;
        // error?
        return output;
    }

        // WORKING HERE!

    return "";
}

/* This function accepts as input a '.txt' file of the format generated using analyzeAlignments. From this input,
 * it parses through each line, identifying each strain that is compared in that line, and adds the line to a map
 * connecting each strain to a set of lines that include the strain. A check is added to ensure that if both
 * strains are the same (strain1 == strain2), it won't be added to the set again (but this is actually unecessary
 * based on the definition of a set...).
 * After constructing the map, each key (strain) included in the map is iterated through, and an output file
 * containing each line in the associated set is generated, using the input file name as a prefix, followed by
 * '_*AsRef.txt' where * is the strain name.
 */
void generateSubAlignmentFiles(string inputFile) {
    // do shit here
    Map<string, Set<string>> strainMap;
    string strain1;
    string strain2;
    string outputFile;
    string prefix;
    string line;
    Vector<string> entries;
    bool skipStrain2;

    // read in file...
    ifstream input(inputFile);
    while(getline(input, line)) {
        entries = stringSplit(line, ",");
        strain1 = trim(entries[1]);                         // extract strain1 and strain 2 names, removing spaces
        strain2 = trim(entries[2]);
        skipStrain2 = false;

        // check if strain1 == strain2 and handle accordingly:
        if (strain1 == strain2) {
            skipStrain2 = true;
        }

        // add line to strain1 set in map:
        if (strainMap.containsKey(strain1)) {
            strainMap[strain1].add(line);                   // append the current line to the set strainMapp[strain1]
        }
        else {
            strainMap[strain1].add(line);                   // append the current line to the set strainMapp[strain1] (these are the same...)
        }

        // add line to strain2 set in map, if skipStrain2 == false:
        if (skipStrain2 == false) {
            if (strainMap.containsKey(strain2)) {
                strainMap[strain2].add(line);               // append the current line to the set strainMapp[strain1]
            }
            else {
                strainMap[strain2].add(line);               // append the current line to the set strainMapp[strain1] (these are the same...)
            }
        }
      //  cout << line << endl;  // actually do something different here...
    }
    cout << "Done filling map!" << endl;

    // convert each set of lines into an output text file...
    prefix = inputFile.erase(inputFile.length() - 4, 4);
    for (string curStrain : strainMap.keys()) {
        ofstream out;
        outputFile = prefix + "_" + curStrain + "AsRef.txt";            // append filename with "_*AsRef.txt"
        out.open(outputFile);
//        if(outputFile.fail()){
//            // do some shit if file doesn't exist?
//        }
        for (string line : strainMap[curStrain]) {
            out << line << endl;
        }
        out.close();
    }
    cout << "Done making files!" << endl;
    return;
}

//void testmakemap(string name) {
//    Map<string, string> map;
//}

// THIS IS A COPY OF THE ABOVE THAT IS MEANT TO HANDLE DUPLICATES!!!
/* This function accepts as input a '.txt' file of the format generated using analyzeAlignments. From this input,
 * it parses through each line, identifying each strain that is compared in that line, and adds the line to a map
 * connecting each strain to a set of lines that include the strain. A check is added to ensure that if both
 * strains are the same (strain1 == strain2), it won't be added to the set again (but this is actually unecessary
 * based on the definition of a set...).
 * After constructing the map, each key (strain) included in the map is iterated through, and an output file
 * containing each line in the associated set is generated, using the input file name as a prefix, followed by
 * '_*AsRef.txt' where * is the strain name.
 */
void generateSubAlignmentFilesHandleDups(string inputFile) {
    // do shit here
    Map<string, Set<string>> strainMap;
    string strain1;
    string strain2;
    string outputFile;
    string prefix;
    string line;
    Vector<string> entries;
    Vector<string> lineVector;
    bool skipStrain2;

    // read in file...
    ifstream input(inputFile);
    while(getline(input, line)) {
        entries = stringSplit(line, ",");
        strain1 = trim(entries[1]);                         // extract strain1 and strain 2 names, removing spaces
        strain2 = trim(entries[2]);
        skipStrain2 = false;

        // This is where the duplicate handling occurs: I will automatically split all strings with the delimiter '_',
        // and take only the first resulting string BEFORE comparing strings 1 and 2...
        strain1 = stringSplit(strain1, "_")[0];
        strain2 = stringSplit(strain2, "_")[0];
        lineVector = stringSplit(line, ",");
        lineVector[1] = strain1;
        lineVector[2] = strain2;
        line = stringJoin(lineVector,",");

        // check if strain1 == strain2 and handle accordingly:
        if (strain1 == strain2) {
            skipStrain2 = true;
        }

        // add line to strain1 set in map:
        if (strainMap.containsKey(strain1)) {
            strainMap[strain1].add(line);                   // append the current line to the set strainMapp[strain1]
        }
        else {
            strainMap[strain1].add(line);                   // append the current line to the set strainMapp[strain1] (these are the same...)
        }

        // add line to strain2 set in map, if skipStrain2 == false:
        if (skipStrain2 == false) {
            if (strainMap.containsKey(strain2)) {
                strainMap[strain2].add(line);               // append the current line to the set strainMapp[strain1]
            }
            else {
                strainMap[strain2].add(line);               // append the current line to the set strainMapp[strain1] (these are the same...)
            }
        }
      //  cout << line << endl;  // actually do something different here...
    }
    cout << "Done filling map!" << endl;

    // convert each set of lines into an output text file...
    prefix = inputFile.erase(inputFile.length() - 4, 4);
    for (string curStrain : strainMap.keys()) {
        ofstream out;
        outputFile = prefix + "_" + curStrain + "AsRef_HD.txt";            // append filename with "_*AsRef.txt"
        out.open(outputFile);
//        if(outputFile.fail()){
//            // do some shit if file doesn't exist?
//        }
        for (string line : strainMap[curStrain]) {
            out << line << endl;
        }
        out.close();
    }
    cout << "Done making files!" << endl;
    return;
}

//void testmakemap(string name) {
//    Map<string, string> map;
//}




/***********************************/
/* * * * * Testing section * * * * */
/***********************************/

/* Function to test functions */
void testing() {
    Vector<string> vec = {"", "", "", "", "", ""};

    cout << "Test case 1: basic with snips, gaps, sublength != 0" ;
    compareSequences("AAATTTCCCGGG---", "AATTTTCCCGGGAA-", vec);
    cout << vec << endl;

    cout << "Test case 2: all sublen" ;
    compareSequences("---", "---", vec);
    cout << vec << endl;

    cout << "Test case 3: all gap" ;
    compareSequences("---", "AAA", vec);
    cout << vec << endl;

    cout << "Test case 4: all snp" ;
    compareSequences("AAA", "CCC", vec);
    cout << vec << endl;

    cout << "Test case 5: empty" ;
    compareSequences("", "", vec);
    cout << vec << endl;
}

/* This is a testing function that verifies the compareSequencesBaseByBase outputs are as expected:
 */
void testing2() {
    //baseByBaseAlignment( "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2022-12-31_tap_compressed_alignment.txt", "/Users/rlporter/Documents/gut_native_ec_motility/strain_seq_variation/Alignment Analysis Testing/C++ version/2022-12-31_tap_compressed_alignments.txt", "MG1655");
    //2022-12-31_tap_compressed_alignment.txt
    cout << compareSequencesBaseByBase("ABCDEFG", "ABCDEFG") << " = 0000000" << endl;
    cout << compareSequencesBaseByBase("ABC----", "ABC----") << " = 000----" << endl;
    cout << compareSequencesBaseByBase("-----", "-----") <<  " = -----" << endl;

    cout << compareSequencesBaseByBase("ABCDEFG", "AbCDEFG") << " = 0100000" << endl;
    cout << compareSequencesBaseByBase("ABC----", "AbC--B-") << " = 010--1-" << endl;
    cout << compareSequencesBaseByBase("-----", "AA---") <<  " = 11---" << endl;

}

/* This is a testing function for the BaseByBaseAlignment function:
 */
void testing3(string inputFile) {
    string line;

    // read in file...
    ifstream input(inputFile);
    while (getline(input, line)) {
        cout << line << endl;
    }
    return;
}

string testing4(string inputFile) {
//    string output = teststring.erase(teststring.length() - 4, 4);
//    return output;
    string line;

    // read in file...
    ifstream input(inputFile);
    while (getline(input, line)) {
        cout << line << endl;
    }
    return "the end";
}
