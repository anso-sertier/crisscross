/*==========================================================
 * // This program extract:
 * //   - read pairs from sample bam which have at least one of these properties:
 * //		- wrong read orientations
 * //		- mapped on different chromosome
 * //       - insert size distance is greater than 1000 kb
 * //   - soft-clip reads
 * //
 * // Many outputs are created :
 * //		- a bamfile of abnormal reads
 * //       - a bamfile of soft-clip reads
 * //	    -
 * //
 * // AUTHOR : anso
 * // DATE : 10/06/2016
 *
 *
 * // MODIFICATION :
 * //		- partie soft-clip : suppression de la partie recherche paire anormale : le bonus associé est rarissime (1 fois sur 7 paires)
 * // MODIFICATION 22/08/2017:
 * //		- partie soft-clip : regarde si la position clippée est la bonne relativement à la référence (parfois quelques nt sont égaux à la référence)
 *
 * // g++ -o bam_extract_ab_sc_check bam_extract_ab_sc_check.cpp -I /data-ddn/software/bamtools/current/include -L /data-ddn/software/bamtools/current/lib -I /data-ddn/software/libraries/boost_1_55_0/include/ -L /data-ddn/software/libraries/boost_1_55_0/lib -lbamtools -lz -lboost_regex -lboost_program_options
 * ==========================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <math.h>
#include <algorithm>

#include "boost/regex.hpp"
#include <boost/algorithm/string.hpp>
#include "boost/program_options.hpp"

#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;

namespace po = boost::program_options;


/*---------------------- Global variable definition ABNORMAL  -----------*/
/* ABNORMAL ONLY */
int MIN_MAP_QUAL     = 20;
/* SOFT-CLIP ONLY */
int MIN_PRC_IDENTITY = 90;
int MIN_SC_LENGTH    = 3;
int MAX_PRC_REF      = 10;
/* BOTH */
int MIN_BASE_QUAL    = 25;
int MIN_PRC_HIGHQUAL = 80;
int MAX_IS           = 1000;
/*---------------------- General fonctions -------------------------------*/

int string2int(string s){
	int number;
	istringstream ss( s );
	ss >> number;
	return(number);
}

unsigned long int string2ulint(string s){
	unsigned long int number;
	istringstream ss( s );
	ss >> number;
	return(number);
}

const char *convert_to_cstr(const std::string & s){
   return s.c_str();
}

bool process_command_line(int argc, char** argv, string& inputFileName,string& refFileName, string& outA, string& outB, string& outC, string& outD, string& outE){
	po::options_description config("Mandatory");
	config.add_options()
		("input,i", po::value<string>(&inputFileName)->required(),  "input bam file name")
		("ref,r",   po::value<string>(&refFileName)->required(),  "reference genome filename")
		("bamAB,a", po::value<string>(&outA)->required(),"output bam file name of abnormal pairs")
		("bamSC,b", po::value<string>(&outB)->required(),"output bam file name of soft-clip reads")
		("fasta,c", po::value<string>(&outC)->required(),  "output fasta file name of soft-clip sequences")
		("tab,d",   po::value<string>(&outD)->required(),  "output text file name of clipped positions summary")
		("tab2,e",  po::value<string>(&outE)->required(),  "output text file name of clipped positions that partially matched reference");
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help,h", "produce help message")
		("insertsize", po::value<int>(&MAX_IS)->default_value(1000),         "set min insert size for deletion detection")
		("sclength",   po::value<int>(&MIN_SC_LENGTH)->default_value(3),     "set min length of soft-clipped sequence at a given position")
	    ("basequal",   po::value<int>(&MIN_BASE_QUAL)->default_value(25),    "set min base quality threshold")
		("prcbases",   po::value<int>(&MIN_PRC_HIGHQUAL)->default_value(80), "set min percent of high quality bases per read")
		("mapqual",    po::value<int>(&MIN_MAP_QUAL)->default_value(20),     "set min mapping quality threshold")
		("identity",   po::value<int>(&MIN_PRC_IDENTITY)->default_value(90), "set min percent of identity in aligned region of soft-clipped reads");
		("prcref",     po::value<int>(&MAX_PRC_REF)->default_value(10),      "set max percent of diff between sc bases and reference (to remove sc identical to ref)");
	po::options_description cmdline_options;
	cmdline_options.add(desc).add(config);
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, cmdline_options), vm);

	try{
		if (vm.count("help")) {
			cout << config << "\n" << desc << "\n";
			return 0;
		}
		po::notify(vm);
	}catch(exception& e){
        cerr << "Error: " << e.what() << "\n" << config << "\n" << desc << "\n";
        return false;
    }catch(...){
        cerr << "Unknown error!" << "\n";
        return false;
    }
	return true;
}

bool validateBaseQuality(string qualities){
    int i,qual;
	int nbGoodQual = 0;
    int nbTot = qualities.length();
 	for(i = 0 ; i < qualities.length() ; i++){
		qual = int(qualities[i]) - 33;
		if(qual >= MIN_BASE_QUAL) nbGoodQual ++;
	}
	// CHECK THRESHOLDS
	if (nbTot == 1){
		if(nbGoodQual != 1) return false;
		else return true;
	}else if (nbTot <= 10){
        if(nbGoodQual < nbTot - 1) return false;
        else return true;
    }else if(nbTot <= 20){
        if(nbGoodQual < nbTot - 2) return false;
        else return true;
    }else{
        if(nbGoodQual*100 < nbTot*MIN_PRC_HIGHQUAL) return false;
        else return true;
    }
}


/*---------------------- ABNORMAL specific fonctions -------------------*/
int testPair(BamAlignment read,BamAlignment mate){
	int hitOptRead,hitOptMate;
	// TESTS UNICITY OF ALIGNMENT
	if (read.GetTag("X0", hitOptRead) && mate.GetTag("X0", hitOptMate)){
		if (hitOptRead > 1 || hitOptMate > 1){ // one of them is non uniq
			return 1;
		}
	}
	// TESTS MAPPING QUALITIES
	if(read.MapQuality < MIN_MAP_QUAL || mate.MapQuality < MIN_MAP_QUAL){
		return 1;
	}
	// TESTS BASE QUALITIES
	if( (! validateBaseQuality(read.Qualities)) && (! validateBaseQuality(mate.Qualities))){ // if both are bad quality then remove pair
		return 1;
	}
	return 0;
}

int analysePair(vector<BamAlignment> abnormalpair, BamWriter& writer){
	BamAlignment read = abnormalpair[0];
	BamAlignment mate = abnormalpair[1];
	int flagWrite = testPair(read, mate);
	if(flagWrite == 0){
		writer.SaveAlignment(read);
		writer.SaveAlignment(mate);
		return 2;
	}
	return 0;
}

/*---------------------- Reference genome handling fonctions -------------------*/
map<string, vector<unsigned long int> > ReadFaidx(string filename){
	ifstream fic;
	fic.open(convert_to_cstr(filename));
	string line;
	map<string, vector<unsigned long int> > idx;
	vector<string> strs;
	vector<string> chrom;
	while(getline(fic,line)){
		boost::split(strs, line, boost::is_any_of("\t"));
		boost::split(chrom, strs[0], boost::is_any_of(" "));
		idx[chrom[0]].push_back(string2ulint(strs[1])); // total length of ref seq in bases
		idx[chrom[0]].push_back(string2ulint(strs[2])); // offset in the file of this sequence first base
		idx[chrom[0]].push_back(string2ulint(strs[3])); // number of bases on each line
		idx[chrom[0]].push_back(string2ulint(strs[4])); // number of bytes in each line including new line
	}
	fic.close();
	return idx;
}

string getSubsequence(string chrom, int start, int lg, map<string, vector<unsigned long int> > idx, string filename){
	int newlines_before = start > 0 ? (start - 1) / idx[chrom][2] : 0;
    int newlines_by_end = (start + lg - 1) / idx[chrom][2];
    int newlines_inside = newlines_by_end - newlines_before;
    int seqlen = lg + newlines_inside;
    char* seq = (char*) calloc (seqlen + 1, sizeof(char));
	FILE * is;
	is = fopen((convert_to_cstr(filename)),"r");
    fseek64(is, (off_t) (idx[chrom][1] + newlines_before + start), SEEK_SET);
    fread(seq, sizeof(char), (off_t) seqlen, is);
	fclose(is);
    seq[seqlen] = '\0';
    char* pbegin = seq;
    char* pend = seq + (seqlen/sizeof(char));
    pend = remove(pbegin, pend, '\n');
    pend = remove(pbegin, pend, '\0');
    string s = seq;
    free(seq);
    s.resize((pend - pbegin)/sizeof(char));
    return s;
}

bool validateReference(string ref, string test){
	int diff = 0;
	for (unsigned i = 0; i < test.size(); i++ ) {
		if (ref[i] != test[i]) {
            diff ++;
		}
	}
	if ( diff*100 / test.size() <= MAX_PRC_REF ){
		return false; // equal to reference (max 5% mismatch)
	}else{
		return true;
	}
}

int nbEqualToReference(string ref, string test, string side){
	int nb = 0;
	if (side == "LEFT"){  // check from end to start
		for (unsigned i = ref.size()-1 ; i >= 0 ; i--){
			if (ref[i] == test[i]) nb ++;
			else break;
		}
	}else{ // check from start to end
		for (unsigned i = 0 ; i < ref.size() ; i++){
			if (ref[i] == test[i]) nb ++;
			else break;
		}
	}
	return nb;
}


/*---------------------- SOFT-CLIP specific fonctions -------------------*/
string parseMD(string s){
	vector<string> alpha, num;
	string matchPos;
	boost::sregex_token_iterator endIt;
	boost::regex re1("[0-9]+");
	boost::sregex_token_iterator it1(s.begin(), s.end(), re1, -1);
	*it1++; // first always empty
	while (it1 != endIt) {
		alpha.push_back(*it1++);
	}
	boost::regex re2("[A-Z^]+");
	boost::sregex_token_iterator it2(s.begin(), s.end(), re2, -1);
	while (it2 != endIt) {
		num.push_back(*it2++);
	}
	int i,k, value;
	for (k = 0 ; k < alpha.size(); k++){
		value = string2int(num[k]);
		for(i = 0; i < value ; i++){
			matchPos += "i"; // identity to reference
		}
		if(alpha[k][0] != '^'){
			matchPos += "s"; // substitution
		}else{
            for(i = 0; i < alpha[k].length() - 1 ; i ++)
            matchPos += "d"; // deletion
        }
	}
	if (num.size() > alpha.size()){
		value = string2int(num[k]);
		for(i = 0; i < value ; i++){
			matchPos += "i";
		}
	}
	return(matchPos);
}

string parseCigar(vector<CigarOp> cigar){
	string cigarStr;
	int i,j;
	for (i = 0 ; i < cigar.size(); i++){
		for(j = 0; j < cigar[i].Length ; j++){
			cigarStr += cigar[i].Type;
		}
	}
	return(cigarStr);
}

bool validateAlignedBases(string cigarStr,string mdStr){
	int nIdentity = 0, alignLength = 0;
	int mdCompteur = 0;   // md string = cigar string minus insertions and clipping
    int flag = 0;
    if(mdStr.length() == 0) flag = 1;
	for (int i = 0 ; i < cigarStr.length() ; i++) {
		if (cigarStr[i] == 'M'){
			alignLength ++;
            if (flag == 1) nIdentity ++;
            else if (mdStr[mdCompteur] == 'i') nIdentity ++;
			mdCompteur ++;
		}else if (cigarStr[i] == 'D') {
			alignLength ++;
            mdCompteur ++;
        }else if (cigarStr[i] == 'I') {
			alignLength ++;
        }
	}
	// CHECK THRESHOLDS
	if (alignLength <= 20){
		if(nIdentity < alignLength - 1) return false;
		else return true;
	}else{
		if (nIdentity * 100 / alignLength < MIN_PRC_IDENTITY) return false;
		else return true;
	}
}

bool validateRead(BamAlignment ali, string side, vector<CigarOp> cigar){
	if(ali.MapQuality == 0) return false;
 	string baseQualities = ali.Qualities;
	string clipQual;
    string alignedQual;
	if (side == "LEFT"){
		clipQual = baseQualities.substr(0, cigar[0].Length);
        alignedQual = baseQualities.substr(cigar[0].Length, string::npos);
        if(cigar[cigar.size()-1].Type == 'S'){
            alignedQual = alignedQual.substr(0, ali.Length - cigar[0].Length - cigar[cigar.size()-1].Length);
        }
	}else if (side == "RIGHT"){
		clipQual = baseQualities.substr(ali.Length - cigar[cigar.size()-1].Length, string::npos);
        alignedQual = baseQualities.substr(0, ali.Length - cigar[cigar.size()-1].Length);
        if(cigar[0].Type == 'S'){
            alignedQual = alignedQual.substr(cigar[0].Length, string::npos);
        }
	}
    // validate clipped base qualities
    if (! validateBaseQuality(clipQual)) return false;
	// validate aligned base qualities
    if (! validateBaseQuality(alignedQual)) return false;
	// validate alignment quality
    string MDtag;
 	string cigarStr = parseCigar(cigar);
    string mdStr = "";
    if (ali.GetTag("MD",MDtag)) mdStr = parseMD(MDtag);
    if (!validateAlignedBases(cigarStr,mdStr)) return false;
    // Soft clip PASS all filters => KEEP
	return true;
}

struct clipping {
	int nbr;
	int max_all;
	string seq;
};

/*------------------------------------------------------------------------*/
/*------------------------------ MAIN PROGRAM ----------------------------*/
/*------------------------------------------------------------------------*/
int main(int argc, char *argv[]){

	/********************************/
	/* PROCESS COMMAND LINE */
	/********************************/
	string inputFileName, refFileName, outputBamABName, outputBamSCName, outputTabName, outputSeqName, outputPartial;
	bool flag = process_command_line(argc,argv,inputFileName, refFileName, outputBamABName, outputBamSCName,
									 outputSeqName, outputTabName, outputPartial);
	if (! flag){
		return 1;
	}

	/********************************/
	/* OPEN INPUT AND OUTPUT FILES */
	/********************************/
	BamReader reader;
	BamWriter writerAB,writerSC;

	if (!reader.Open(convert_to_cstr(inputFileName))) {
		cerr << "Could not open input BAM file." << endl;
		return 1;
	}
	SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();

	if (!writerAB.Open(convert_to_cstr(outputBamABName),header,references)) {
		cerr << "Could not open output BAM file for abnormal pairs." << endl;
		return 1;
	}
	if (!writerSC.Open(convert_to_cstr(outputBamSCName),header,references)) {
		cerr << "Could not open output BAM file for soft_clip reads." << endl;
		return 1;
	}

	ofstream out,outseq,outpartial;
	out.open(convert_to_cstr(outputTabName));
	outseq.open(convert_to_cstr(outputSeqName));
	outpartial.open(convert_to_cstr(outputPartial));
	if(!out.is_open()){
		cerr << "Could not open output text file for soft-clip summary." << endl;
		return 1;
	}
	if(!outseq.is_open()){
		cerr << "Could not open output fasta file for soft-clip sequences." << endl;
		return 1;
	}
	if(!outpartial.is_open()){
		cerr << "Could not open output text file for soft-clip partial." << endl;
		return 1;
	}
	out << "chrom\tpos\tside\tnbReads\tmaxSC\n";
	outpartial << "chrom\tpos\tside\tnbidentical\tseq\n";
	
	ifstream ref, refidx;
	ref.open(convert_to_cstr(refFileName));
	refidx.open(convert_to_cstr(refFileName + ".fai"));
	if(! ref.is_open()){
		cerr << "Could not open input reference genome." << endl;
		return 1;
	}
	if(! refidx.is_open()){
		cerr << "Could not open reference genome index, please check syntax : ref.fasta.fai" << endl;
		return 1;
	}
	ref.close();
	refidx.close();
	
	map<string, vector<unsigned long int> > idx = ReadFaidx(refFileName + ".fai");

	/*********************************/
	/* GET BAM SEQUENCE INFORMATIONS */
	/*********************************/
	cout << "Starts on " << inputFileName << endl;

	SamSequenceDictionary refsequences = header.Sequences ;
	SamSequenceIterator itseq;
	SamSequence currSeq;
	vector<string> referenceNames;
	string name;
	for (itseq = refsequences.Begin(); itseq < refsequences.End() ; ++itseq){
		currSeq = *itseq;
		if (currSeq.HasName()){
			name = currSeq.Name;
			referenceNames.push_back(name);
			if(name == "Y" || name == "chrY"){
				break;
			}
		}
	}

	/***********************/
	/* READ INPUT BAM FILE */
	/***********************/
	BamAlignment ali;
	string strand1,strand2,strand;
	int pos,matePos;
	unsigned long int nbSee = 0;

	/*--- variable for AB ---*/
	map<string, vector<BamAlignment> > abnormalpairs;
	int nbABkeep = 0;

	/*--- variable for SC---*/
	int clipSeqLg, nbIdent;
	int scPos;
	int flagWritten;
	string clipSeq,refSeq;
	vector<CigarOp> cigar;
	map<int, map<int, pair<clipping,clipping> > > posCounts; // chrom : pos : (clipping_left, clipping_right)
	int nbSCkeep = 0;
	int nbIdentical = 0;
	int nbPartialIdentical = 0;
	
	while ( reader.GetNextAlignment(ali) ) {
		nbSee ++;
		/**********************************/
		/***** COMMON CHECKS        *******/
		/**********************************/
		if(ali.RefID >= referenceNames.size() || (!ali.IsMapped()) || ali.IsDuplicate() || ali.IsFailedQC() || (! ali.IsPrimaryAlignment())){
			continue;
		}
		/**********************************/
		/***** PROCESS SOFT-CLIP    *******/
		/**********************************/
		cigar = ali.CigarData;
		flagWritten = 0;
		if (cigar[0].Type == 'S' && validateRead(ali, "LEFT", cigar)){
            scPos     = ali.Position; // ali.Position is 0-based -> last soft-clipped position is 1-based
			if(scPos == 0) scPos = 1;
			clipSeqLg = cigar[0].Length;
			clipSeq   = ali.QueryBases.substr(0,cigar[0].Length);
			refSeq    = getSubsequence(referenceNames[ali.RefID], scPos - clipSeqLg , clipSeqLg, idx, refFileName);
			if (validateReference(refSeq, clipSeq)){
				nbIdent   = nbEqualToReference(refSeq, clipSeq, "LEFT");
				if (nbIdent > 0){
					nbPartialIdentical ++;
					outpartial << referenceNames[ali.RefID] << "\t" << scPos << "\t5'\t" << nbIdent << "\t" << clipSeq << endl;
				}else{
					nbSCkeep ++;
					posCounts[ali.RefID][scPos].first.nbr     ++;
					posCounts[ali.RefID][scPos].first.max_all = max(posCounts[ali.RefID][scPos].first.max_all, clipSeqLg);
					posCounts[ali.RefID][scPos].first.seq     += clipSeq + "\n";
				}
				writerSC.SaveAlignment(ali);
				flagWritten = 1;
			}else {
				nbIdentical ++ ;
			}
		}
		if (cigar[cigar.size()-1].Type == 'S' && validateRead(ali, "RIGHT",cigar)){
			scPos     = ali.GetEndPosition(); // 1-based
			clipSeqLg = cigar[cigar.size()-1].Length;
			clipSeq   = ali.QueryBases.substr(ali.Length - clipSeqLg , string::npos);
			refSeq    = getSubsequence(referenceNames[ali.RefID], scPos , clipSeqLg, idx, refFileName);
			if (validateReference(refSeq, clipSeq)){
				nbIdent   = nbEqualToReference(refSeq, clipSeq, "RIGHT");
				if (nbIdent > 0){
					nbPartialIdentical ++;
					outpartial << referenceNames[ali.RefID] << "\t" << scPos << "\t3'\t" << nbIdent << "\t" << clipSeq << endl;
				}else{
					nbSCkeep ++;
					posCounts[ali.RefID][scPos].second.nbr     ++;
					posCounts[ali.RefID][scPos].second.max_all = max(posCounts[ali.RefID][scPos].second.max_all, clipSeqLg);
					posCounts[ali.RefID][scPos].second.seq     += clipSeq + "\n";
				}
				if(flagWritten == 0) writerSC.SaveAlignment(ali); // otherwise, already saved in previous if
			}else{
				nbIdentical ++ ;
			}
		}


		/***********************************/
		/***** PROCESS ABNORMAL PAIR *******/
		/***********************************/
		
		/*** ADDITIONNAL CHECKS ************/
		if(ali.MateRefID >= referenceNames.size() || (! ali.IsMateMapped())){
			continue;
		}
		/*** INTERCHROM PAIRS TESTS ********/
		if(ali.RefID != ali.MateRefID){
			abnormalpairs[ali.Name].push_back(ali);
			if(abnormalpairs[ali.Name].size() == 2){
				nbABkeep += analysePair(abnormalpairs[ali.Name], writerAB);
				abnormalpairs.erase(ali.Name);
			}
			continue;
		}
		/*** DEFINE STRAND AND POSITION  ***/
		if (ali.IsReverseStrand()) strand1 = "R";
		else strand1 = "F";
		if (ali.IsMateReverseStrand()) strand2 = "R";
		else strand2 = "F";
		if (ali.Position < ali.MatePosition){
			strand = strand1 + strand2;
			pos = ali.Position;
			matePos = ali.MatePosition;
		}else{
			strand = strand2 + strand1;
			pos = ali.MatePosition;
			matePos = ali.Position;
		}
		
		/*** STRAND TESTS ***/
		if( (! ((strand == "FR" && pos < matePos) || (strand == "RF" && matePos < pos))) || (abs(ali.InsertSize) >= MAX_IS) ){
			abnormalpairs[ali.Name].push_back(ali);
			if(abnormalpairs[ali.Name].size() == 2){
				nbABkeep += analysePair(abnormalpairs[ali.Name], writerAB);
				abnormalpairs.erase(ali.Name);
			}
		}
	}
	reader.Close();
	writerAB.Close();
	cout << "Read processed               : " << nbSee << endl;
	cout << "Nb abnormal reads kept       : " << nbABkeep << endl;
	cout << "Nb sc identical to reference : " << nbIdentical << endl;
	cout << "Nb partially ident to ref    : " << nbPartialIdentical << endl;
	cout << "Nb soft-clip keep            : " << nbSCkeep << endl;

	/**********************************/
	/***** WRITE TEXT FILES     *******/
	/**********************************/
	map<int, map<int, pair<clipping,clipping> > >::iterator itchrom;
	map<int, pair<clipping,clipping> >::iterator itpos;
	clipping left,right;
	int chrom;
	int nbPos = 0;
	int nbPosKept = 0;
	for(itchrom = posCounts.begin() ; itchrom != posCounts.end() ; itchrom++){
		chrom = itchrom->first;
		for(itpos = itchrom->second.begin() ; itpos != itchrom->second.end() ; itpos++){
			pos = itpos->first;
			left = posCounts[chrom][pos].first;
			right = posCounts[chrom][pos].second;
			if(left.max_all >= MIN_SC_LENGTH ){
				out    << referenceNames[chrom] << "\t" << pos << "\t5'\t" << left.nbr << "\t" << left.max_all << endl;
				outseq << "> " << referenceNames[chrom] << ":" << pos << "\t5'\n" << left.seq;
				nbPosKept ++;
			}
			if(right.max_all >= MIN_SC_LENGTH ){
				out << referenceNames[chrom] << "\t" << pos << "\t3'\t" << right.nbr   << "\t" << right.max_all << endl;
				outseq << "> " << referenceNames[chrom] << ":" << pos << "\t3'\n" << right.seq;
				nbPosKept ++;
			}
			nbPos ++;
		}
	}
	out.close();
	outseq.close();
	outpartial.close();
	cout << "\tNb position found : " << nbPos << endl;
	cout << "\tNb position kept  : " << nbPosKept << endl;
	cout << "\ndone." << endl;
	return 0;
}
