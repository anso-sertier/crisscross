/*==========================================================
// This program create cluster of reads of the same type
//
// Two outputs are created :
//		- a tab file discribing the blocs
//		- a bed file of blocs
//
// AUTHOR : anso
// DATE : 03/08/2016


// g++ -o ab_clustering ab_clustering.cpp -I /data-ddn/software/bamtools/current/include -L /data-ddn/software/bamtools/current/lib -I /data-ddn/software/libraries/boost_1_55_0/include/ -L /data-ddn/software/libraries/boost_1_55_0/lib -lbamtools -lz -lboost_program_options
==========================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "boost/program_options.hpp"

#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/SamHeader.h"
#include "api/BamWriter.h"


using namespace std;
using namespace BamTools;
namespace po = boost::program_options;

int OFFSET=300;
int MAX_IS=1000;

string STRANDS[] = {"FR","RF","FF","RR"};

/* PROCESS COMMAND LINE*/
const char *convert_to_cstr(const std::string & s){
   return s.c_str();
}

bool process_command_line(int argc, char** argv, string& input, string& outA, string& outB){
	po::options_description config("Mandatory");
	config.add_options()
		("input,i", po::value<string>(&input)->required(), "input abnormal bam file name")
		("tab,a",   po::value<string>(&outA)->required(),"output tab file discribing blocs")
		("bed,b",   po::value<string>(&outB)->required(),"output bed file discribing blocs");
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help,h", "produce help message")
		("insertsize", po::value<int>(&MAX_IS)->default_value(1000), "set min insert size for deletion detection")
		("delta",      po::value<int>(&OFFSET)->default_value(300),  "set distance around clusters to include new read");
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


/* STRUCTURES DEFINITION */
struct bloc{
	int chr1;
	int chr2; // same as chr1 in case of intra chrom pair
	int start1;
	int start2;
	int end1;
	int end2;
	int insert; // sum of pair insert size, 0 in case of interchr
	int nbp;    // abnormal pairs
	int nbo;    // overlapping pairs
	string strands;
};

struct readinfo{
	int refchr;
	int start;
	int end;
};

map<string,int> initialize_counters(){
	map<string,int> m;
	m["interFR"] = 1;
	m["interRF"] = 1;
	m["interFF"] = 1;
	m["interRR"] = 1;
	m["intraFR"] = 1;
	m["intraRF"] = 1;
	m["intraFF"] = 1;
	m["intraRR"] = 1;
return m;
}

string int2string(int i){
        string s;
        stringstream ss;
        ss << i;
        s = ss.str();
        return s;
}

/* FUNCTIONS */

readinfo fillReadinfo(BamAlignment ali){
   readinfo newread;
   newread.refchr = ali.RefID;
   newread.start  = ali.Position;
   newread.end    = ali.GetEndPosition();
   return(newread);
}

bloc fillBloc(int chr1, int chr2, int start1, int start2, int end1, int end2, int insert, int nbp, int nbo,string strand){
	bloc newbloc;
	newbloc.chr1   = chr1;
	newbloc.chr2   = chr2;
	newbloc.start1 = start1;
	newbloc.start2 = start2;
	newbloc.end1   = end1;
	newbloc.end2   = end2;
	newbloc.insert = insert;
	newbloc.nbp    = nbp;
	newbloc.nbo    = nbo;
	newbloc.strands = strand;
	return newbloc;
}

bloc NUL = fillBloc(0,0,0,0,0,0,0,0,0,"");

bloc mergeBloc(bloc current, bloc added){
	bloc merge;
	merge.chr1    = current.chr1;
	merge.chr2    = current.chr2;
	merge.start1  = min(current.start1, added.start1);
	merge.start2  = min(current.start2, added.start2);
	merge.end1    = max(current.end1,   added.end1);
	merge.end2    = max(current.end2,   added.end2);
	merge.insert  = max(current.insert, added.insert);
	merge.nbp     = current.nbp + added.nbp;
	merge.nbo     = current.nbo + added.nbo;
	merge.strands = current.strands;
	return merge;
}

void compareBloc(vector<bloc>& blocs, bloc newbloc){
	int flag = 0;
	bloc current;
	vector<bloc>::iterator it = blocs.end();
	it --;
	while(flag == 0 && it >= blocs.begin()){
		current = *it;
		if(newbloc.start1 <= current.end1 + OFFSET && newbloc.start2 <= current.end2 + OFFSET && newbloc.end2 >= current.start2 - OFFSET){
			*it = mergeBloc(current,newbloc);
			flag = 1;
		}else if(newbloc.start1 <= current.end1 + 1000){
			it--;
		}else{
			break;
		}
	}
	if(flag == 0){ // add new bloc
		blocs.push_back(newbloc);
	}
}

int compareAndMergeBlocs(vector<bloc>& blocs){
   vector<bloc>::iterator it = blocs.end();
   vector<bloc>::iterator itto ;
  	it --;
   bloc current,tocompare;
   int nbmerge = 0;
   while(it >= blocs.begin()){
      current = *it;
      int flag = 0;
      itto = it - 1;
      while(flag == 0 && itto >= blocs.begin()){
         tocompare = *itto;
         if( min(current.end1,tocompare.end1) - max(current.start1,tocompare.start1) > 0 && min(current.end2,tocompare.end2) - max(current.start2,tocompare.start2) > 0){
            *itto = mergeBloc(current,tocompare);
            blocs.erase(it);
            flag = 1;
            nbmerge ++;
         }else if(current.start1 <= tocompare.end1 + 3000){
            itto--;
         }else{
            break;
         }
      }
      it--;
   }
   return nbmerge;
}

void writeOutput(vector<bloc> list, string name, map<string,int>& counters, vector<string> referenceNames, ofstream& out, ofstream& outbed){
   int num = counters[name];
	for (int i = 0; i < list.size(); i++){
		if (list[i].nbp + list[i].nbo != 0){
			out << name << "\t" << num << "\t" \
			    << referenceNames[list[i].chr1] << "\t" << list[i].start1 << "\t" << list[i].end1 << "\t" \
			    << referenceNames[list[i].chr2] << "\t" << list[i].start2 << "\t" << list[i].end2 << "\t" \
				<< list[i].nbp << "\t" << list[i].nbo << "\t"<< list[i].insert << "\n"; // 0-based
			outbed << referenceNames[list[i].chr1] << "\t" << list[i].start1 << "\t" << list[i].end1 + 1 << "\t" << name << num <<"a\t" << list[i].nbp + list[i].nbo << "\n"; // 0-based
			outbed << referenceNames[list[i].chr2] << "\t" << list[i].start2 << "\t" << list[i].end2 + 1 << "\t" << name << num <<"b\t" << list[i].nbp + list[i].nbo << "\n"; // 0-based
         num ++;
		}
	}
   counters[name] = num;
}

void defineBothStrand(BamAlignment read, string &strand){
	if (read.IsReverseStrand()) strand = "R";
	else strand = "F";
	if (read.IsMateReverseStrand()) strand += "R";
	else strand += "F";
}

void treatIntraPairs(map<string, map<int, vector<string> > >& data, int chrom, map<string, vector<readinfo> >& abnormalpairs,
					 vector<string> referenceNames, ofstream& out, ofstream& outbed, map<string,int>& counters, int& nbmerge){
	int insert,overlap;
	vector<string> intra;
	readinfo read,mate;
	string name,strand;
	bloc newbloc;
	vector<bloc> blocs;
	blocs.push_back(NUL);
	string blocname;
	for (int s = 0; s <= 3; s++ ){
		strand = STRANDS[s];
		intra = data[strand][chrom];
		blocname = "intra" + strand;
		for(int i = 0 ; i < intra.size() ; i++){
			name = intra[i];
			if(abnormalpairs[name].size() != 2) continue;
			read = abnormalpairs[name][0];
			mate = abnormalpairs[name][1];
			insert   = max(read.end,mate.end) - min(read.start,mate.start) ;
			overlap  = read.end - mate.start;
			if(overlap >= 0) newbloc = fillBloc(read.refchr, mate.refchr, read.start, mate.start, read.end, mate.end, insert,0,1,strand);
			else             newbloc = fillBloc(read.refchr, mate.refchr, read.start, mate.start, read.end, mate.end, insert,1,0,strand);
			compareBloc(blocs, newbloc);
			abnormalpairs.erase(name);
		}
      nbmerge += compareAndMergeBlocs(blocs);
		writeOutput(blocs, blocname, counters, referenceNames, out, outbed);
		data[strand].erase(chrom);
		blocs.clear();
		blocs.push_back(NUL);
	}
}

void treatInterPairs(map<string, map<int, map<int, vector<string> > > >& data , int chr1, int chr2, map<string, vector<readinfo> >& abnormalpairs,
					 vector<string> referenceNames, ofstream& out, ofstream& outbed, map<string,int>& counters, int& nbmerge){
	vector<string> inter;
	readinfo read,mate;
	string name,strand;
	bloc newbloc;
	vector<bloc> blocs;
	blocs.push_back(NUL);
	string blocname;
	for (int s = 0; s <= 3; s++ ){
		strand   = STRANDS[s];
		inter    = data[strand][chr1][chr2];
		blocname = "inter" + strand;
      cout << blocname << " : " << inter.size() << endl;
		for(int i = 0 ; i < inter.size() ; i++){
			name = inter[i];
			read = abnormalpairs[name][0];
			mate = abnormalpairs[name][1];
			newbloc = fillBloc(read.refchr, mate.refchr, read.start, mate.start, read.end, mate.end, 0,1,0,strand);
			compareBloc(blocs, newbloc);
			abnormalpairs.erase(name);
		}
      nbmerge += compareAndMergeBlocs(blocs);
		writeOutput(blocs, blocname, counters, referenceNames, out, outbed);
		data[strand][chr1].erase(chr2);
		blocs.clear();
		blocs.push_back(NUL);
	}
}


int main(int argc, char* argv[]){

	/********************************/
	/* PROCESS COMMAND LINE */
	/********************************/
	string inputFileName,outputTabName, outputBedName;
	bool flag = process_command_line(argc,argv,inputFileName, outputTabName, outputBedName);
	if (! flag){
		return 1;
	}
	cout << "Start clustering abnormal pairs : " << inputFileName << endl;

	/*******************************/
	/* OPEN INPUT AND OUTPUT FILES */
	/*******************************/
	BamReader reader;
	if (!reader.Open(convert_to_cstr(inputFileName))) {
		cerr << "Could not open input BAM file." << endl;
		return 1;
	}
	SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();

	SamSequenceDictionary refsequences = header.Sequences ;
	SamSequenceIterator itseq;
	SamSequence currSeq;
	vector<string> referenceNames;
	for (itseq = refsequences.Begin(); itseq < refsequences.End() ; ++itseq){
		currSeq = *itseq;
		if (currSeq.HasName()){
			referenceNames.push_back(currSeq.Name);
		}
	}
	ofstream out, outbed;
	out.open(convert_to_cstr(outputTabName));
	outbed.open(convert_to_cstr(outputBedName));
	if(!out.is_open()){
		cerr << "Could not open output text file for clusters summary." << endl;
		return 1;
	}
	if(!outbed.is_open()){
		cerr << "Could not open output bed file of clusters positions." << endl;
		return 1;
	}
	out << "type\tnum\tchrom1\tstart1\tend1\tchrom2\tstart2\tend2\tnbpab\tnbpov\tinsert\n";

	/***********************************/
	/***** TREAT INPUT BAM FILE ******* /
	/***********************************/
	BamAlignment ali;
	string strand;
	int nbSee = 0;
	int current_chrom = 0;
   int nbmerge = 0;
 
	map<string, vector<readinfo> >                     abnormalpairs; // ali.Name : [read, mate]
	map<string, map<int, vector<string> > >            intrapairs;	  // ori : chr1 : [ali.Name]
	map<string, map<int, map<int, vector<string> > > > interpairs ;   // ori : chr1 : chr2 : [ali.Name]
	map<string,int>                                    counters = initialize_counters() ;

	while ( reader.GetNextAlignment(ali) ) {
		nbSee ++;
		if (ali.RefID > current_chrom){
			treatIntraPairs(intrapairs, current_chrom, abnormalpairs,referenceNames, out, outbed, counters, nbmerge);
			for (int i = (current_chrom - 1) ; i >= 0 ; i--){
				treatInterPairs(interpairs,i , current_chrom, abnormalpairs, referenceNames, out, outbed, counters, nbmerge);
			}
			current_chrom = ali.RefID;
		}
		abnormalpairs[ali.Name].push_back(fillReadinfo(ali));
		if(abnormalpairs[ali.Name].size() == 1){ // first read in pair seen
			defineBothStrand(ali,strand);
			if(ali.RefID < ali.MateRefID){
				interpairs[strand][ali.RefID][ali.MateRefID].push_back(ali.Name);
			}else{
				intrapairs[strand][ali.RefID].push_back(ali.Name);
			}
		}
	}
	current_chrom = ali.RefID;
	treatIntraPairs(intrapairs, current_chrom, abnormalpairs,referenceNames, out, outbed, counters, nbmerge);
	for (int i = (current_chrom - 1) ; i >= 0 ; i--){
		treatInterPairs(interpairs,i , current_chrom, abnormalpairs, referenceNames, out, outbed, counters, nbmerge);
	}
	reader.Close();
	out.close();
	outbed.close();

	cout << "\tNb reads seen           : " << nbSee << endl;
   cout << "\tMerged blocs (2nd time) : " << nbmerge << endl;
	cout << "\tStay in abnormalpairs   : " << abnormalpairs.size() << endl;
	cout << "done." << endl;
	return 0;
}
