// Parse a PDB file, identify contacts calculate P(contact|residue_type_i,residue_type_j)
// (c) Tim Nugent 2013

#include <stdio.h>
#include <math.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

#define VERSION 1.0

/* Convert AA letter to numeric code (0-21) */
int aanum(int ch, int alphabet){

    const static int aacvs[] =
    {
	999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2,
	21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
    };

    if(!alphabet){
    	return (isalpha(ch) ? aacvs[ch & 31] : 20);
    }else if (alphabet == 1){

/*

Accuracy of Sequence Alignment and Fold Assessment
Using Reduced Amino Acid Alphabets
Francisco Melo and Marc A. Marti-Renom

As it should have been expected before hand, the JO20
potential exhibited the best performance in model assess-
ment, irrespective of the sensitivity/specificity balance at
any given classification threshold. The potential based on
the novel MM5 alphabet was the second best classifier for
fold assessment. Then, the SR5 and WW5 potentials
exhibited an overall similar performance, although that
for high specificities the SR5 potential was more accurate.
Finally, the MU4 potential showed the worst performance
in fold assessment among all the reduced potentials tested.

JO20	A-C-D-E-F-G-H-I-K-L-M-N-P-Q-R-S-T-V-W-Y 	20	*
WW5		AHT-CFILMVWY-DE-GP-KNQRS 					5	3rd
SR5		AEHKQRST-CFILMVWY-DN-G-P 					5	2nd
MU4		AGPST-CILMV-DEHKNQR-FYW 					4	4th
MM5		AG-C-DEKNPQRST-FILMVWY-H 					5	1st
RD5		100 randomly reduced alphabets 				5

A reduced amino acid alphabet for understanding and designing
protein adaptation to mutation
C. Etchebest Æ C. Benros Æ A. Bornot Æ
A.-C. Camproux Æ A. G. de Brevern

IV-WFY-ALM-EQRK-ND-CHST-P-G

*/
    	// MM5	AG-C-DEKNPQRST-FILMVWY-H	5

    	static const char* g1 = "AG";
    	static const char* g2 = "C";
    	static const char* g3 = "DEKNPQRST";
    	static const char* g4 = "FILMVWY";
    	static const char* g5 = "H";

    	if (strchr(g1, ch) != NULL){
    		return 0;
    	}else if (strchr(g2, ch) != NULL){
    		return 1;
    	}else if (strchr(g3, ch) != NULL){
    		return 2;
    	}else if (strchr(g4, ch) != NULL){
    		return 3;
    	}else if (strchr(g5, ch) != NULL){
    		return 4;
    	}else{
    		return 5;
    	} 
    }else if (alphabet == 2){

    	// SR5	AEHKQRST-CFILMVWY-DN-G-P	5

    	static const char* g1 = "AEHKQRST";
    	static const char* g2 = "CFILMVWY";
    	static const char* g3 = "DN";
    	static const char* g4 = "G";
    	static const char* g5 = "P";

    	if (strchr(g1, ch) != NULL){
    		return 0;
    	}else if (strchr(g2, ch) != NULL){
    		return 1;
    	}else if (strchr(g3, ch) != NULL){
    		return 2;
    	}else if (strchr(g4, ch) != NULL){
    		return 3;
    	}else if (strchr(g5, ch) != NULL){
    		return 4;
    	}else{
    		return 5;
    	} 

    }else if (alphabet == 3){

    	// WW5	AHT-CFILMVWY-DE-GP-KNQRS	5

    	static const char* g1 = "AHT";
    	static const char* g2 = "CFILMVWY";
    	static const char* g3 = "DE";
    	static const char* g4 = "GP";
    	static const char* g5 = "KNQRS";

    	if (strchr(g1, ch) != NULL){
    		return 0;
    	}else if (strchr(g2, ch) != NULL){
    		return 1;
    	}else if (strchr(g3, ch) != NULL){
    		return 2;
    	}else if (strchr(g4, ch) != NULL){
    		return 3;
    	}else if (strchr(g5, ch) != NULL){
    		return 4;
    	}else{
    		return 5;
    	} 
    }else{
    	cout << "Unknown amino acid alphabet!" << endl;
    	exit(-1);
    }
    return -1;	
}

struct atom{	

	int atomnum_;
	string atom_;
	int resnum_;
	string res_;
	int res_alph_;
	double x_, y_, z_;
	atom(int an, string a, int rn, string r, int ra, double x, double y, double z) : atomnum_(an), atom_(a), resnum_(rn), res_(r), res_alph_(ra), x_(x), y_(y), z_(z) {}     

};  

class PDB{

public:

	PDB() : contact_count(0), alphabet_sz(20) {}
	~PDB();
	void parse_pdb(const char*, int);
	void count_contacts(bool print, int min_sep = 5, double threshold1 = 8.0, double threshold2 = 6.0);
	void print_atoms();
	void get_stats();
	bool is_contact(atom&,atom&,double,double);
	int get_length(){return residues.size();}
	int get_confreq(int i, int j){return freqs_con[i][j];}
	int get_nconfreq(int i, int j){return freqs_ncon[i][j];}
private:

	map<int, vector<atom> > residues;
	int** freqs_con;
	int** freqs_ncon;
	int contact_count, alphabet, alphabet_sz;
	
};

PDB::~PDB(){

    for(int i = 0; i < alphabet_sz; i++){
    	delete [] freqs_con[i];
    }
	delete [] freqs_con;    
    for(int i = 0; i < alphabet_sz; i++){
    	delete [] freqs_ncon[i];
    }
	delete [] freqs_ncon;    

}

void PDB::parse_pdb(const char* pdb, int a){

	alphabet = a;
	if(alphabet) alphabet_sz = 5;
	ifstream input(pdb);

	string res3[20] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
	char res1[20] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
    map <string, char> aa3_1; 
    for(int i = 0; i < 20; i++){
            aa3_1[res3[i]] = res1[i];
    }
    freqs_con = new int*[alphabet_sz];
    for(int i = 0; i < alphabet_sz; i++){
    	freqs_con[i] = new int[alphabet_sz];
    	for(int j = 0; j < alphabet_sz; j++){
    		freqs_con[i][j] = 0;
    	}	
    }
    freqs_ncon = new int*[alphabet_sz];
    for(int i = 0; i < alphabet_sz; i++){
    	freqs_ncon[i] = new int[alphabet_sz];
    	for(int j = 0; j < alphabet_sz; j++){
    		freqs_ncon[i][j] = 0;
    	}	
    }

	// Parse PDB file
	for(std::string line; getline(input, line);){
		if (strncmp(line.substr(0,4).c_str(),"ATOM",4) == 0){

			string at = line.substr(13,4);
			string::iterator end_pos = remove(at.begin(), at.end(), ' ');
			at.erase(end_pos, at.end());

			if(at.compare(0,1,"H") != 0){

				string res = line.substr(17,3);
				double x = atof(line.substr(30,8).c_str());
				double y = atof(line.substr(38,8).c_str());
				double z = atof(line.substr(46,8).c_str());		
				string atomnum = line.substr(6,5);		
				string resnum = line.substr(22,4);
				int an = stoi(atomnum);
				int rn = stoi(resnum);
				
				// Map the 3 letter aa code to a single letter aa code, then to an int depending on which alphabet is selected 
				int res_alph =  aanum(aa3_1[res],alphabet);

				atom p(an,at,rn,res,res_alph,x,y,z);

				if(residues.find(rn) == residues.end()){
					vector<atom> tmp;
					tmp.push_back(p);
					residues[rn] = tmp;	
				}else{	
					residues[rn].push_back(p);	
				}
			}
		}	
	}	
	input.close();

}

void PDB::print_atoms(){

	for(auto it = residues.begin(); it != residues.end(); it++){

		for(size_t i = 0; i < it->second.size(); i++){
			cout << it->second[i].atomnum_ << "\t"  << it->second[i].atom_ << "\t" << it->second[i].resnum_ << "\t" <<  it->second[i].res_ << "\t" << it->second[i].x_ << ","  << it->second[i].y_ << ","  << it->second[i].z_ << endl;
		}
	}

}

bool PDB::is_contact(atom& a, atom& b, double threshold1, double threshold2){

	double dist = sqrt((a.x_-b.x_)*(a.x_-b.x_) + (a.y_-b.y_)*(a.y_-b.y_) + (a.z_-b.z_)*(a.z_-b.z_));
	if( (a.res_ == "GLY" && a.atom_ == "CA" && b.res_ != "GLY" && b.atom_ == "CB") ||
		(a.res_ == "GLY" && a.atom_ == "CA" && b.res_ == "GLY" && b.atom_ == "CA") ||
		(a.res_ != "GLY" && a.atom_ == "CB" && b.res_ == "GLY" && b.atom_ == "CA") ||
		(a.res_ != "GLY" && a.atom_ == "CB" && b.res_ != "GLY" && b.atom_ == "CB") ){

		if(dist < threshold1){
			// pymol check
			//cout << "distance mydist, " << a.resnum_ << "/" << a.atom_ << "," <<  b.resnum_ << "/" << b.atom_ << "\t# " << dist << endl;
			return true;
		}else{
			return false;
		}

	}else{

		if(dist < threshold2){
			// pymol check
			//cout << "distance mydist, " << a.resnum_ << "/" << a.atom_ << "," <<  b.resnum_ << "/" << b.atom_ << "\t# " << dist << endl;
			return true;
		}else{
			return false;
		}

	}
}

void PDB::count_contacts(bool print, int min_sep, double threshold1, double threshold2){

	int r1 = 0;
	for(auto ita = residues.begin(); ita != residues.end(); ita++){
		auto itb = ita;
		itb++;
		int r2 = 0;
		for(; itb != residues.end(); itb++){
			bool contact = false;
			for(size_t i = 0; i < ita->second.size(); i++){
				for(size_t j = 0; j < itb->second.size(); j++){
					if(abs(ita->second[0].resnum_ - itb->second[0].resnum_) > min_sep){
						if(is_contact(ita->second[i],itb->second[j],threshold1,threshold2)){
							contact = true;
							break;
						}
					}else{
						break;
					}	
				}
				if(contact){
					break;
				}
			}
			if(contact){
				contact_count++;				

				freqs_con[ita->second[0].res_alph_][itb->second[0].res_alph_]++;
				freqs_con[itb->second[0].res_alph_][ita->second[0].res_alph_]++;
				//cout << tmp[0] << tmp[1] << "\t" << abs(ita->second[0].resnum_ - itb->second[0].resnum_) << endl;
				//if(print) cout << ita->first << "\t" << itb->first << "\t" << abs(ita->second[0].resnum_ - itb->second[0].resnum_) << "\t" << tmp[0] << tmp[1] << endl;				
					
			}else{

				if(abs(ita->second[0].resnum_ - itb->second[0].resnum_) > min_sep){
					// Frequencies for non-contacting pairs
					freqs_ncon[ita->second[0].res_alph_][itb->second[0].res_alph_]++;
					freqs_ncon[itb->second[0].res_alph_][ita->second[0].res_alph_]++;
				}
			}
			r2++;
		}	
		r1++;
	}
	
}

void PDB::get_stats(){

	cout << "# Length:\t\t" << residues.size() << endl;
	cout << "# Contacts (sep > 5):\t" << contact_count << endl;	
	cout << "# Density:\t\t" << setprecision(3) << (double)contact_count/((double)(residues.size()*(residues.size()-1))/2) << endl;
	cout << "# Alphabet:\t\t" << alphabet << " (size = " << alphabet_sz << ")" << endl;

	cout << "Contacting pairs:" << endl;
    for(int i = 0; i < alphabet_sz; i++){
    	for(int j = 0; j < alphabet_sz; j++){
    		cout << freqs_con[i][j] << "\t";
    	}
    	cout << endl;	
    }
    cout << "Non-contacting pairs:" << endl;
    for(int i = 0; i < alphabet_sz; i++){
    	for(int j = 0; j < alphabet_sz; j++){
    		cout << freqs_ncon[i][j] << "\t";
    	}
    	cout << endl;	
    }
}	

void usage(const char* progname){

	cout << "Usage: " << endl << progname << " PDB " << endl; 

}

void scan_directory(const char* path, int alphabet){

	int alphabet_sz = 20;
	if(alphabet) alphabet_sz = 5;

    double** priors = new double*[alphabet_sz];
    for(int i = 0; i < alphabet_sz; i++){
    	priors[i] = new double[alphabet_sz];
    	for(int j = 0; j < alphabet_sz; j++){
    		priors[i][j] = 0.0;
    	}	
    }

    int** freqs_con = new int*[alphabet_sz];
    for(int i = 0; i < alphabet_sz; i++){
    	freqs_con[i] = new int[alphabet_sz];
    	for(int j = 0; j < alphabet_sz; j++){
    		freqs_con[i][j] = 0;
    	}	
    }
    int** freqs_ncon = new int*[alphabet_sz];
    for(int i = 0; i < alphabet_sz; i++){
    	freqs_ncon[i] = new int[alphabet_sz];
    	for(int j = 0; j < alphabet_sz; j++){
    		freqs_ncon[i][j] = 0;
    	}	
    }

	size_t i = 0;
	DIR *dir;
	struct dirent *ent;
	struct stat st;
	string pdb_path = path;
	if(pdb_path.compare(pdb_path.length()-1,1,"/") != 0){
		pdb_path.append("/");
	}
	dir = opendir(pdb_path.c_str());
	if (dir != NULL) {
	  	while ((ent = readdir (dir)) != NULL) {

			string pdb = pdb_path;
			pdb.append(ent->d_name);
			lstat(pdb.c_str(), &st);

			// PDB file			
			if(S_ISREG(st.st_mode) && (pdb.compare(pdb.length()-3,3,"pdb") == 0 || pdb.compare(pdb.length()-3,3,"PDB") == 0)){
							
				cout << pdb << endl; 
				PDB* protein = new PDB;
				protein->parse_pdb(pdb.c_str(),alphabet);

				protein->count_contacts(false,5);

			    for(int i = 0; i < alphabet_sz; i++){
			    	for(int j = i; j < alphabet_sz; j++){
			    		freqs_con[i][j] += protein->get_confreq(i,j);
			    		freqs_con[j][i] = freqs_con[i][j];
			    		freqs_ncon[i][j] += protein->get_nconfreq(i,j);
			    		freqs_ncon[j][i] = freqs_ncon[i][j];
			    	}
			    }
				delete protein;	
			}
	  	}
	  	closedir (dir);
	}else{
		cout << "Couldn't open PDB path " << pdb_path << endl;
	}

	int total_contacts  = 0;
	int total_ncontacts = 0;

    for(int i = 0; i < alphabet_sz; i++){
    	for(int j = i; j < alphabet_sz; j++){
    		total_contacts += freqs_con[i][j];
    		total_ncontacts += freqs_ncon[i][j];
    	}	
    }	

    int total = total_contacts+total_ncontacts;
    double p_contact = (double)total_contacts/total;
    double p_ncontact = 1-p_contact;
    cout << "Average density: " << p_contact << endl;


	cout << "Contacting pairs:" << endl;
    for(int i = 0; i < alphabet_sz; i++){
    	for(int j = 0; j < alphabet_sz; j++){
    		cout << freqs_con[i][j] << "\t";
    	}
    	cout << endl;	
    }
    cout << "Non-contacting pairs:" << endl;
    for(int i = 0; i < alphabet_sz; i++){
    	for(int j = 0; j < alphabet_sz; j++){
    		cout << freqs_ncon[i][j] << "\t";
    	}
    	cout << endl;	
    }

    // Calculate priors here
    for(int i = 0; i < alphabet_sz; i++){
    	for(int j = i; j < alphabet_sz; j++){

    		// p(contact|aa=i,aa=j)
    		// = p(contact) * freq_con(i,j)/all_contacts
    		//   ---------------------------------------------------------
    		//   p(contact) * freq(i,j)/all_contacts   +  p(ncontact) * freq_ncon(i,j)/all_ncontact

    		double top = p_contact * ((double)freqs_con[i][j]/total_contacts);
    		double bottom = p_ncontact * ((double)freqs_ncon[i][j]/total_ncontacts);
    		priors[i][j] = top/(top+bottom);
    		priors[j][i] = priors[i][j];

    	}	
    }	

    cout << "Priors:" << endl;
    cout << "double priors["<< alphabet_sz << "][" << alphabet_sz << "] = {" << endl;
    for(int i = 0; i < alphabet_sz; i++){
    	cout << "{";	
    	for(int j = 0; j < alphabet_sz; j++){
    		cout << priors[i][j];
    		if(j != alphabet_sz-1) cout << ",";
    	}
    	cout << "}";
    	if(i != alphabet_sz-1) cout << ",";	
    	cout << endl;
    }
    cout << "};" << endl;

    for(int i = 0; i < alphabet_sz; i++){
    	delete [] priors[i];
    }
	delete [] priors;   
    for(int i = 0; i < alphabet_sz; i++){
    	delete [] freqs_con[i];
    }
	delete [] freqs_con;    
    for(int i = 0; i < alphabet_sz; i++){
    	delete [] freqs_ncon[i];
    }
	delete [] freqs_ncon;   

}

int main(int argc, const char* argv[]){
	
	if (argc < 2){
		usage(argv[0]);
		return(1);
	}

	int i = 0, alphabet = 1;
    while(i < argc){
        if(argv[i][0] == '-'){
        i++;
        switch(argv[i-1][1]){
        	case 'a' : {alphabet = atoi(argv[i]);break;}
            case 'd' : {scan_directory(argv[i],alphabet); return(1);}
            case 'h' : {usage(argv[0]); return(1);}
            default  : {usage(argv[0]); return(1);}
            }   
       	}
       	i++;
   	}


	// Parse PDB file
	PDB* protein = new PDB;
	protein->parse_pdb(argv[argc-1],alphabet);
	protein->count_contacts(true);
	protein->get_stats();
	delete protein;
	return(0);
}
