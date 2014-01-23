#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <math.h> 
#include <RInside.h>

using namespace std;

#define MAXSEQLEN 5000

void fail(string fmt){

	cout << fmt << endl;
    exit(-1);

}

struct sc_entry{
	sc_entry() : sc(0.0) {}
    double sc;
    int i, j;
};

/* Sort descending */
int cmpfn(const void *a, const void *b){

    if (((struct sc_entry *)a)->sc == ((struct sc_entry *)b)->sc)
        return 0;

    if (((struct sc_entry *)a)->sc < ((struct sc_entry *)b)->sc)
        return 1;

    return -1;
}

/*
int aanum(int ch){
    
    const static int aacvs[] =
    {
        999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2,
        21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
    };
	return (isalpha(ch) ? aacvs[ch & 31] : 20);

}
*/

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
    	fail("Unknown amino acid alphabet!");
    }
    return -1;	
}

/* Allocate matrix */
double** allocmat(int rows, int columns){
    double** p = new double*[rows];
    if(p == NULL){
        fail("allocmat failed!");
    }
    for(int i = 0; i < rows; i++){
       p[i] = new double[columns];
       if(p[i] == NULL){
            fail("allocmat failed!");
       }
       for(int j = 0; j < columns; j++){
       		p[i][j] = 0.0;
       }	
    }
    return p;
}

/* Free matrix */
void freemat(double** p, int rows, int columns){
    for(int i = 0; i < rows; i++){
       delete [] p[i];
    }   
    delete [] p;
}

int index4D(int i, int j, int a, int b, int w, int ww){
        int row = (i*w)+a;
        int column = (j*w)+b;
        return((row*(w*ww))+column);
}

int test_cholesky(Rcpp::NumericMatrix& a, const int n){

    //cout << "#In cholesky..." << endl;
    int i, j, k;
    double sum;
    double *diag = new double[n];

    for (i=0; i<n; i++){

        for (j=i; j<n; j++){

            sum = a(i,j);

            for (k=i-1; k >= 0; k--)
                sum -= a(i,k)*a(j,k);

            if (i == j){

                if (sum <= 0.0){
                    delete [] diag;
                    return TRUE;
                }

                diag[i] = sqrt(sum);
            }
            else
                a(j,i) = sum / diag[i];
        }
    }
    delete [] diag;
    return FALSE;
}

int main(int argc, char **argv){

	double wtsum, pcmean, pc, idthresh = -1.0, pseudoc = 1.0, sumfnzero = 0.0, percentile = 0.05;
	unsigned int nseqs, opt, seqlen, i, j, k, a, b, npair, nnzero, filt = 0, shrinkflg = 1, minseqsep = 5, apcflg = 1, all_theta = 1, wts = 1, alphabet = 1, alphabet_sz = 21;
	int o = 1;

	double rho = 1e-06;
    double l1  = 1e-08;
    double l2  = 1e-08;

    /*
	rho = 7e-04;
	l1  = 4e-06;
	l2  = 3e-02;	
    */

	char seq[MAXSEQLEN];
	char *native = NULL;
	string jgl_weight = "\"equal\"";
	string penalty = "group";

    // Bayesian priors
    double priors_a0[20][20] = {
        {0.0384709,0.0315522,0.0252431,0.0267246,0.0553167,0.0280882,0.0234792,0.0309125,0.0348211,0.0575771,0.0544788,0.0252531,0.0506351,0.0580766,0.0278139,0.0277613,0.03306,0.0535248,0.0569143,0.0548293},
        {0.0315522,0.0377166,0.0395014,0.0448287,0.0462737,0.0351509,0.0437371,0.0327015,0.0489691,0.0430665,0.0381715,0.0291482,0.0372603,0.0521655,0.0323253,0.0352403,0.0354181,0.0633929,0.0600957,0.0396119},
        {0.0252431,0.0395014,0.030496,0.0298982,0.0365854,0.0313936,0.025177,0.0273452,0.0345646,0.0338721,0.0308303,0.0264725,0.0329986,0.0432402,0.0295818,0.0336339,0.0309915,0.0546282,0.0525324,0.0364171},
        {0.0267246,0.0448287,0.0298982,0.0256499,0.0350737,0.0242983,0.0168017,0.0255038,0.0466577,0.0332681,0.0283254,0.0325473,0.0298229,0.039159,0.0274418,0.0288016,0.0298283,0.040271,0.0447852,0.0313199},
        {0.0553167,0.0462737,0.0365854,0.0350737,0.191235,0.0479339,0.0242558,0.0432369,0.0639594,0.0714286,0.0772233,0.0302804,0.0968137,0.115055,0.0458971,0.0424671,0.0475578,0.0919689,0.0986005,0.0721717},
        {0.0280882,0.0351509,0.0313936,0.0242983,0.0479339,0.0308642,0.0225366,0.0282611,0.040107,0.0384517,0.0365945,0.0274413,0.0398421,0.0444999,0.0320218,0.028232,0.0327219,0.0632699,0.0503027,0.0380285},
        {0.0234792,0.0437371,0.025177,0.0168017,0.0242558,0.0225366,0.0168077,0.0196585,0.0346681,0.0308996,0.0284074,0.0358452,0.0265916,0.0360999,0.0236483,0.0263646,0.0267739,0.045389,0.0471886,0.0305756},
        {0.0309125,0.0327015,0.0273452,0.0255038,0.0432369,0.0282611,0.0196585,0.0309288,0.0381158,0.0348741,0.0326037,0.0259678,0.0399194,0.0455832,0.0300998,0.0280512,0.0309509,0.047529,0.0472133,0.0343498},
        {0.0348211,0.0489691,0.0345646,0.0466577,0.0639594,0.040107,0.0346681,0.0381158,0.0644567,0.0473282,0.0461363,0.0257868,0.0596919,0.0649699,0.042328,0.0426136,0.0498899,0.0743728,0.0664414,0.0452912},
        {0.0575771,0.0430665,0.0338721,0.0332681,0.0714286,0.0384517,0.0308996,0.0348741,0.0473282,0.101794,0.0999256,0.0326071,0.0779581,0.102943,0.0371438,0.0361649,0.0495798,0.0930421,0.0847148,0.097934},
        {0.0544788,0.0381715,0.0308303,0.0283254,0.0772233,0.0365945,0.0284074,0.0326037,0.0461363,0.0999256,0.095633,0.0297359,0.0796293,0.107474,0.0334634,0.0350263,0.0482116,0.0995159,0.078081,0.0861285},
        {0.0252531,0.0291482,0.0264725,0.0325473,0.0302804,0.0274413,0.0358452,0.0259678,0.0257868,0.0326071,0.0297359,0.0247232,0.0294663,0.0382533,0.0258621,0.0268877,0.0273437,0.0463936,0.0499595,0.0349487},
        {0.0506351,0.0372603,0.0329986,0.0298229,0.0968137,0.0398421,0.0265916,0.0399194,0.0596919,0.0779581,0.0796293,0.0294663,0.0797342,0.098869,0.0350605,0.0379338,0.0435572,0.0910165,0.0748803,0.0736103},
        {0.0580766,0.0521655,0.0432402,0.039159,0.115055,0.0444999,0.0360999,0.0455832,0.0649699,0.102943,0.107474,0.0382533,0.098869,0.143694,0.0463234,0.0462518,0.0550044,0.114045,0.0935484,0.101713},
        {0.0278139,0.0323253,0.0295818,0.0274418,0.0458971,0.0320218,0.0236483,0.0300998,0.042328,0.0371438,0.0334634,0.0258621,0.0350605,0.0463234,0.037079,0.0246798,0.034871,0.0523203,0.0572042,0.0346482},
        {0.0277613,0.0352403,0.0336339,0.0288016,0.0424671,0.028232,0.0263646,0.0280512,0.0426136,0.0361649,0.0350263,0.0268877,0.0379338,0.0462518,0.0246798,0.0257632,0.0323108,0.0415786,0.045959,0.0389282},
        {0.03306,0.0354181,0.0309915,0.0298283,0.0475578,0.0327219,0.0267739,0.0309509,0.0498899,0.0495798,0.0482116,0.0273437,0.0435572,0.0550044,0.034871,0.0323108,0.0408865,0.0578888,0.0545251,0.0485699},
        {0.0535248,0.0633929,0.0546282,0.040271,0.0919689,0.0632699,0.045389,0.047529,0.0743728,0.0930421,0.0995159,0.0463936,0.0910165,0.114045,0.0523203,0.0415786,0.0578888,0.107769,0.0944882,0.0823899},
        {0.0569143,0.0600957,0.0525324,0.0447852,0.0986005,0.0503027,0.0471886,0.0472133,0.0664414,0.0847148,0.078081,0.0499595,0.0748803,0.0935484,0.0572042,0.045959,0.0545251,0.0944882,0.0743175,0.0763453},
        {0.0548293,0.0396119,0.0364171,0.0313199,0.0721717,0.0380285,0.0305756,0.0343498,0.0452912,0.097934,0.0861285,0.0349487,0.0736103,0.101713,0.0346482,0.0389282,0.0485699,0.0823899,0.0763453,0.0915821}
    };

    double priors_a1[5][5] = {
        {0.0329891,0.0493345,0.0273959,0.0467397,0.0364271},
        {0.0493345,0.191235,0.038926,0.0838025,0.0639594},
        {0.0273959,0.038926,0.029892,0.0387009,0.0402696},
        {0.0467397,0.0838025,0.0387009,0.0927607,0.0525634},
        {0.0364271,0.0639594,0.0402696,0.0525634,0.0644567}
    };

    double priors_a2[5][5] = {
        {0.0309809,0.0432743,0.0297187,0.0284433,0.0290587},
        {0.0432743,0.0923695,0.0352178,0.0377782,0.0395318},
        {0.0297187,0.0352178,0.0285302,0.0262594,0.0282917},
        {0.0284433,0.0377782,0.0262594,0.0309288,0.0300998},
        {0.0290587,0.0395318,0.0282917,0.0300998,0.037079}
    };

    double priors_a3[5][5] = {
        {0.0380645,0.0532596,0.0282864,0.032105,0.0304052},
        {0.0532596,0.0923695,0.0328692,0.0384289,0.0391278},
        {0.0282864,0.0328692,0.0186447,0.0234738,0.0322004},
        {0.032105,0.0384289,0.0234738,0.0313865,0.0284151},
        {0.0304052,0.0391278,0.0322004,0.0284151,0.0303451}
    };

	while(o < argc){
        if(argv[o][0] == '-'){
            o++;
            switch(argv[o-1][1]){
            		case 'p' : {penalty = "fused"; o--; break;} // use fused rather than group penalty (slow/might crash)
            		case 'c' : {native = argv[o]; break;} // file containing native contacts
            		case 'r' : {rho = atof(argv[o]); break;} // jgl rho param
                    case 'm' : {l1 = atof(argv[o]); break;} // jgl lambda1 param
                    case 'n' : {l2 = atof(argv[o]); break;} // jgl lambda2 param
                    case 'z' : {apcflg = atoi(argv[o]); break;} // apply apc correction
                    case 'x' : {all_theta = atoi(argv[o]); break;} // get all thetas
                    case 'a' : {alphabet = atoi(argv[o]); break;} // use reduced alphabet (1-3)
                    case 'w' : {wts = atoi(argv[o]); break;} // use weights, 1 = condition, 2 = bayesian priors
                    case 's' : {shrinkflg = atoi(argv[o]); break;} // shrinkage type, 1 = JS shrinkage, 2 = diagonal, 0 = none
                    case 'f' : {filt = atoi(argv[o]); break;} // filter low weights at whatever percentile
                    case 'g' : {percentile = atof(argv[o]); break;} // filter matrices with weights in this percentile
                    default  : {fail("Options:\t-m float lamda1 penalty\n-n float lambda2 penalty\n-w use weights\n-c file native contacts\n");}
                }   
           }
       o++;
    }

    if(alphabet) alphabet_sz = 6;
    unsigned int classes = (alphabet_sz-1)*(alphabet_sz-1);

    cout << "#penalty   = " << penalty << endl;
    cout << "#alphabet  = " << alphabet << endl;
    cout << "#weight    = " << wts << endl;
    cout << "#shrink    = " << shrinkflg << endl;
    cout << "#rho       = " << rho<< endl;
    cout << "#lambda1   = " << l1 << endl;
    cout << "#lambda2   = " << l2 << endl;
    cout << "#all_theta = " << all_theta << endl;
    cout << "#apc flag  = " << apcflg << endl;
    cout << "#aln       = " << argv[o-1] << endl;

    // Parse native contacts is provided
    vector<string> native_contacts;
    if(native != NULL){
	    FILE *ifp = fopen(native, "r");
	    if (ifp == NULL){
	    	cout << "Failed to parse native contact file!" << endl;
	    	return(1);	
	    }
	    int i, j ,t1, t2;
	    double score;	    
	    char buf[100];
	    while(fgets(buf,100,ifp) != NULL){
			if (sscanf(buf, "%d %d %d %d %lf",&i, &j, &t1, &t2, &score) == 5){
				string con = to_string((long long int)i);
				con.append("-");
				con.append(to_string((long long int)j));
				native_contacts.push_back(con);
			}
	    }
	    fclose(ifp);	
	    if(native_contacts.size()) cout << "#Finished parsing native contacts" << endl;
	}

	FILE* ifp = fopen(argv[o-1], "r");
    if (!ifp)
        fail("Unable to open alignment file!");

    for (nseqs = 0;; nseqs++)
        if (!fgets(seq, MAXSEQLEN, ifp))
            break;

    //aln = (char**)allocvec(nseqs, sizeof(char *));    
    //weight = (double*)allocvec(nseqs, sizeof(double));
    //wtcount = (unsigned int*)allocvec(nseqs, sizeof(unsigned int));
 
    char** aln = new char*[nseqs];
    double* weight = new double[nseqs];
    unsigned int* wtcount = new unsigned int[nseqs];
    rewind(ifp);
    
    if (!fgets(seq, MAXSEQLEN, ifp))
        fail("Bad alignment file!");
    
    seqlen = strlen(seq)-1;

    if (!(aln[0] = new char[seqlen]))
        fail("Out of memory!");

    for (j = 0; j < seqlen; j++){
    	aln[0][j] = aanum(seq[j],alphabet);
    }            

    for (i = 1; i < nseqs; i++){

        if (!fgets(seq, MAXSEQLEN, ifp))
            break;
        
        if (seqlen != strlen(seq)-1)
            fail("Length mismatch in alignment file!");
        
        if (!(aln[i] = new char[seqlen]))
            fail("Out of memory!");
        
        for (j = 0; j < seqlen; j++)
            aln[i][j] = aanum(seq[j],alphabet);

    }

    cout << "#Finished parsing alignment file" << endl;

	double **pa = allocmat(seqlen,alphabet_sz);
    double *pab = new double[seqlen*seqlen*alphabet_sz*alphabet_sz];

    /* Calculate sequence weights */
    //if(0){
    if (idthresh < 0.0){
        double meanfracid = 0.0;
        for (i=0; i<nseqs; i++){
            for (j=i+1; j<nseqs; j++){
                
                int nids;
                double fracid;
                
                for (nids=k=0; k<seqlen; k++)
                    if (aln[i][k] == aln[j][k])
                        nids++;
                
                fracid = (double)nids / seqlen;
                meanfracid += fracid;
            }
        }
        meanfracid /= 0.5 * nseqs * (nseqs - 1.0);
        idthresh = 0.38 * 0.32 / meanfracid;
    }

    for (i=0; i<nseqs; i++){
        for (j=i+1; j<nseqs; j++){
            int nthresh = seqlen * (int)idthresh;
            for (k =0 ; nthresh > 0 && k < seqlen; k++){
                if (aln[i][k] != aln[j][k])
                    nthresh--;
            }
            if (nthresh > 0){
                wtcount[i]++;
                wtcount[j]++;
            }
        }
    }
    for (wtsum=i=0; i<nseqs; i++){
        wtsum += (weight[i] = 1.0 / (1 + wtcount[i]));    
    }
	//}
    cout << "#Finished calculating sequence weights" << endl;

	/* Calculate singlet frequencies with pseudocount */
    for (i=0; i<seqlen; i++){
		for (a = 0; a < alphabet_sz; a++)
	    	pa[i][a] = pseudoc;
		for (k = 0; k < nseqs; k++){
	    	a = aln[k][i];
	    	if (a < alphabet_sz)
				pa[i][a] += weight[k];
		}	
		for (a = 0; a < alphabet_sz; a++)
	    	pa[i][a] /= pseudoc * alphabet_sz + wtsum;
    }

    /* Calculate pair frequencies with pseudocount */
    for (i=0; i<seqlen; i++){
		for (j=i+1; j<seqlen; j++){

		    for (a=0; a < alphabet_sz; a++){
				for (b=0; b < alphabet_sz; b++){
			    	pab[index4D(i,j,a,b,alphabet_sz,seqlen)] = pseudoc / (double)alphabet_sz;
				}
		    }
		    for (k=0; k<nseqs; k++){
				a = aln[k][i];
				b = aln[k][j];
				if (a < alphabet_sz && b < alphabet_sz){
			    	pab[index4D(i,j,a,b,alphabet_sz,seqlen)] += weight[k];
				}
		    }		    
		    for (a=0; a < alphabet_sz; a++){
				for (b=0; b< alphabet_sz; b++){
			    	pab[index4D(i,j,a,b,alphabet_sz,seqlen)] /= pseudoc * alphabet_sz + wtsum;
			    	pab[index4D(j,i,a,b,alphabet_sz,seqlen)] = pab[index4D(i,j,a,b,alphabet_sz,seqlen)];
		    	}
			}
    	}
    }

    for (i=0; i<seqlen; i++){
		for (a=0; a < alphabet_sz; a++){
	    	for (b=0; b < alphabet_sz; b++){
				pab[index4D(i,i,a,b,alphabet_sz,seqlen)] = (a == b) ? pa[i][a] : 0.0;
	    	}
	    }		
	}    
	cout << "#Finished calculating frequencies" << endl;



	/* Form the covariance matrices */
    RInside R(argc, argv);               
	R.parseEval("suppressMessages(library(JGLx))");
	R.parseEval("library(corpcor)");
	string r_code = "paste(1:" + to_string((long long unsigned int)seqlen) + ",sep=\"\")";
	Rcpp::CharacterVector labels = R.parseEval(r_code);
	Rcpp::List dimnms = Rcpp::List::create(labels,labels);
	Rcpp::List mylist(classes);
	Rcpp::NumericVector weights(classes);
    Rcpp::NumericVector prior_weights(classes);
	double max_cond = -100000.0;
    unsigned int ndim = seqlen*alphabet_sz-1;


    for (a = k = 0; a < alphabet_sz-1; a++){
        for (b = 0; b < alphabet_sz-1; b++){
    
        	Rcpp::NumericMatrix rMat(seqlen,seqlen);
        	//cout << "a = " << a << " --- b = " << b << endl;

    		for (i = 0; i < seqlen; i++){
        		for (j = 0; j < seqlen; j++){
                    //if (i != j){
                    	rMat(i,j) = pab[index4D(i,j,a,b,alphabet_sz,seqlen)] - pa[i][a] * pa[j][b];
	                //}else{
                    //    rMat(i,j) = pa[i][a];
                    //}       
    			}
    		}
    		rMat.attr("dimnames") = dimnms;			   
		    
		    // Shrinkage
		    if(shrinkflg == 1){
		    	
		    	R["X"] = rMat;
			    
			    
			    //int posdef = R.parseEval("is.positive.definite(X)");
			    //cout << "(pre)Positive definite: " << posdef << endl;
			    Rcpp::List condition = R.parseEval("rank.condition(X)");
			    int rank1 = Rcpp::as<int>(condition[0]);
			    //double cond = (condition[1]);
			    //cout << "(pre)Rank: " << rank << endl;
			    //cout << "(pre)Condition: " << cond << endl;
			    //conut << "---shrinking---" << endl;
				//time_t start,end;
			    //time (&start);
			    
			    
			    R.parseEval("s2 = cov.shrink(X,verbose=FALSE)");

			    int posdef = R.parseEval("is.positive.definite(s2)");
			    //cout << "Positive definite: " << posdef << endl;
			    //Rcpp::List condition_shrunk = R.parseEval("rank.condition(s2)");
			    //int rank2 = Rcpp::as<int>(condition_shrunk[0]);
			    //double cond = Rcpp::as<double>(condition_shrunk[1]);
    		    //cout << "Rank: " << rank << endl;
			    //cout << "Condition: " << cond << endl;
			    Rcpp::NumericMatrix M = R.parseEval("s2");
                //if(rank2 < rank1){
		    	mylist[k] = M;
                //}    

	    	}else if(shrinkflg == 2){

                //cout << "#Calculing mean of diag..." << endl;
                double smean = 0.0;
                for(smean=i=0; i<seqlen; i++){
                    //cout << i << " --- " << rMat(i,i) << endl;
                    smean += rMat(i,i);
                }
                        
                smean /= (double)seqlen;
                double lambda = 0.25;
                cout << "#Starting diag shrinkage..." << endl;
                for (;;){

                    //cout << "#Made a copy..." << endl;
                    Rcpp::NumericMatrix tempmat(rMat);
                    //memcpy(tempmat, cmat, ndim*ndim*sizeof(double));
                    //cout << "#Done..." << endl;

                    /* Test if positive definite using Cholesky decomposition */
                    //R["X"] = tempmat;
                    //int posdef = R.parseEval("is.positive.definite(X)");
                    //if (!posdef)
                    //    break;

                    if (!test_cholesky(tempmat, seqlen)){
                        break;
                    }
                    
                    for (i=0; i<seqlen; i++){
                        for (j=0; j<seqlen; j++){
                            if (i != j){
                                rMat(i,j) *= 1.0 - lambda;
                            }else{
                                rMat(i,j) = smean * lambda + (1.0 - lambda) * rMat(i,j);
                            }
                        }
                    }
                }

                //cout << "#Finished diag shrinkage..." << endl;
                mylist[k] = rMat;

            }else{
                // No shrinkage
                mylist[k] = rMat;
            }    

            // Set up weights
            R["X"] = mylist[k];
            Rcpp::List condition_shrunk = R.parseEval("rank.condition(X)");
            weights[k] = Rcpp::as<double>(condition_shrunk[1]);
            if(Rcpp::as<double>(condition_shrunk[1]) > max_cond){
                max_cond = Rcpp::as<double>(condition_shrunk[1]);
            }

            if(!alphabet){
                prior_weights[k++] = 10*priors_a0[a][b];
            }else if(alphabet == 1){    
                prior_weights[k++] = 10*priors_a1[a][b];
            }else if(alphabet == 2){    
                prior_weights[k++] = 10*priors_a2[a][b];
            }else if(alphabet == 2){    
                prior_weights[k++] = 10*priors_a3[a][b];
            }    

            //cout << "#Prior for (" << a << "," << b << ") = " << 10*priors[a][b] << endl;
    	}	
    }	
    cout << "#Finished shrinking covariance matrices" << endl;

    if(wts) jgl_weight = "weights";
    vector<double> wts_sorted;
	for(k = 0; k < classes; k++){    	
		//cout << "# " << k << " - " << weights[k] << " ---> ";    		
		if(weights[k] == max_cond){
			// Can't be zero
			weights[k] = 1e-5;
		}else{
			weights[k] = 1.0-(weights[k]/max_cond);

		}
		wts_sorted.push_back(weights[k]);
		//cout << weights[k] << " (max = " << max_cond << ")" << endl;
	}
    cout << "#Finished calculating weights" << endl;

    // Filter
    Rcpp::List mylist_filt(classes-(int)(classes*percentile));
    Rcpp::NumericVector weights_filt(classes-(int)(classes*percentile));    
    Rcpp::NumericVector prior_weights_filt(classes-(int)(classes*percentile));    
    if(filt){
    	double filt_threshold = 0.0;
	    sort(wts_sorted.begin(),wts_sorted.end());
	    for(i = 0; i < classes; i++){
	    	cout << i << " --- " << wts_sorted[i];
	    	if(i+1 == (unsigned int)(classes*percentile)){
	    		cout << " **** ";
	    		filt_threshold = wts_sorted[i];
	    	}
	    	cout << endl;	
	    }    	
	    for(k = i = 0; k < classes; k++){ 
	    	if(weights[k] > filt_threshold){
	    		weights_filt[i] = weights[k];
                prior_weights_filt[i] = prior_weights[k];
	    		mylist_filt[i++] = mylist[k];
	    	}else{
	    		cout << "Skipping weight " << weights[k] << " - below/equal threhsold " << filt_threshold << endl;
	    	}	
	    }
	    cout << "#Finished filtering matrices by weight" << endl;
	    R["data"] = mylist_filt;
        if(wts > 1){
            cout << "#Using Bayesian priors" << endl;
            R["weights"] = prior_weights_filt;
        }else{    
            R["weights"] = weights_filt;
       }

	}else{
		R["data"] = mylist;
        if(wts > 1){
            cout << "#Using Bayesian priors" << endl;
            R["weights"] = prior_weights;
        }else{    
		    R["weights"] = weights;
	   }
    }

    //sc_entry* sclist = new sc_entry[seqlen * (seqlen - 1) / 2]; 
    time_t start,end;    

	cout << "#Running JGL..." << endl;

    if(native_contacts.size()){
        cout << "# l1\tl2\trho\tncon\tdens\t1\t2\t5\t10\t20\t50\t100" << endl;
    }

    for(l1 = 1e-10; l1 <= 1e-01; l1 *= 10){
        for(l2 = 1e-10; l2 <= 1e-01; l2 *= 10){

    time (&start);

    sc_entry* sclist = new sc_entry[seqlen * (seqlen - 1) / 2]; 
    unsigned int ncon = 0;
    double pcmean = 0.0;

    if(all_theta){
        //cout << "#Returning all thetas..." << endl;

        r_code = "ggl.results = JGLx(Y=data,penalty=\"" + penalty + "\",lambda1=" + to_string((long double)l1) + ",lambda2=" + to_string((long double)l2) + ",rho=" + to_string((long double)rho) + ",weight=" + jgl_weight + ",penalize.diagonal=FALSE,maxiter=500,tol=1e-5,warm=NULL,return.whole.theta=TRUE,screening=\"fast\",truncate = 1e-5)";
        Rcpp::List ret = R.parseEval(r_code);            
        Rcpp::List theta = ret[0];

        double** pcmat = allocmat(seqlen,seqlen);
        double* pcsum = new double[seqlen];
        for (i = 0; i < seqlen; i++){
            pcsum[i] = 0.0;
        }    
        pcmean = 0.0;

        for(k = 0; k < classes; k++){
            
            SEXP ll = theta[k];
            Rcpp::NumericMatrix y(ll);    
            //cout << "#Theta : " << k << endl;
            for (i = 0; i < seqlen; i++){
                for (j = i+1; j < seqlen; j++){
                    //cout << i << "-" << j << ":" << y(i,j) << "- ";
                    pcmat[i][j] += fabs(y(i,j));
                    pcmat[j][i] = pcmat[i][j];
                }
                //cout << endl;
            }
            //cout << endl;        
        }

        for (i = 0; i < seqlen; i++){
            for (j = i+1; j < seqlen; j++){
                pcsum[i] += pcmat[i][j];
                pcsum[j] += pcmat[i][j];
                pcmean += pcmat[i][j];
            }
        }

        pcmean /= seqlen * (seqlen - 1) * 0.5;

        for (ncon=i=0; i<seqlen; i++){
            for (j=i+minseqsep; j<seqlen; j++){
                if (pcmat[i][j] > 0.0){
                    /* Calculate APC score */
                    if (apcflg){
                        sclist[ncon].sc = pcmat[i][j] - pcsum[i] * pcsum[j] / sqrt(seqlen - 1.0) / pcmean;
                    }else{
                        sclist[ncon].sc = pcmat[i][j];
                    }
                    sclist[ncon].i = i+1;
                    sclist[ncon++].j = j+1;
                }
            }    
        }
        freemat(pcmat,seqlen,seqlen);
        delete [] pcsum;

    }else{
	
    	r_code = "ggl.results = JGLx(Y=data,penalty=\"" + penalty + "\",lambda1=" + to_string((long double)l1) + ",lambda2=" + to_string((long double)l2) + ",rho=" + to_string((long double)rho) + ",weight=" + jgl_weight + ",penalize.diagonal=FALSE,maxiter=500,tol=1e-5,warm=NULL,return.whole.theta=FALSE,screening=\"fast\",truncate = 1e-5)";
    	R.parseEvalQ(r_code);    	
    	//cout << "#Finished group lasso (" << (int)difftime (end,start) << " seconds)" << endl;
        
        r_code = "dim(ggl.results[[1]][[1]])[1]";
        unsigned int dim = Rcpp::as<int>(R.parseEval(r_code));
        //cout << "Dim = " << dim << endl;

        r_code = "dimnames(ggl.results[[1]][[1]])[[1]]";
        Rcpp::List v = (R.parseEval(r_code));  
        //cout << "dimnames:" << endl;
        //for(int i = 0; i < dim; i++){
        //    cout << Rcpp::as<string>(v[i]) << endl;     
        //}   

        r_code = "ggl.results[[1]][[1]]";
        SEXP ll = R.parseEval(r_code);
        Rcpp::NumericMatrix y(ll);  
       
        for(i = 0; i < dim; i++){
        	for (j = i+minseqsep; j< dim; j++){
                //cout << y(i,j) << " ";
                if(i != j && y(i,j) != 0){
    				sclist[ncon].sc = fabs(y(i,j));
    	            sclist[ncon].i = atoi(Rcpp::as<string>(v[i]).c_str());
    	            sclist[ncon++].j = atoi(Rcpp::as<string>(v[j]).c_str());
               }
            }
            //cout << endl;    
        }

    }
    time (&end);    

    int tp100 = 0;  
    int tp50 = 0;
    int tp20 = 0;
    int tp10 = 0;
    int tp5 = 0;
    int tp2 = 0;
    int tp1 = 0;   

    sumfnzero = (double)ncon/(seqlen * (seqlen - 1) * 0.5);
    //cout << "#Contacts = " << ncon << endl;
    //cout << "#Density  = " << setprecision(4) << sumfnzero << endl;

	qsort(sclist, ncon, sizeof(struct sc_entry), cmpfn);
    for (i = 0; i < ncon; i++){
    	if(native_contacts.size()){
			string con = to_string((long long int)sclist[i].i);
			con.append("-");
			con.append(to_string((long long int)sclist[i].j));
			if (find(native_contacts.begin(),native_contacts.end(),con) != native_contacts.end()){
				if(i < 100)tp100++;	
				if(i < 50)tp50++;
				if(i < 20)tp20++;
				if(i < 10)tp10++;
				if(i < 5)tp5++;
				if(i < 2)tp2++;
				if(i < 1)tp1++;
			}
		}else{
			printf("%d %d 0 8 %f\n", sclist[i].i, sclist[i].j, sclist[i].sc);
		}        	
    }    
    delete [] sclist;

    if(native_contacts.size()){

		double p100 = 0.0;	
		double p50 = 0.0;
		double p20 = 0.0;
		double p10 = 0.0;
		double p5 = 0.0;
		double p2 = 0.0;
		double p1 = 0.0;

	    if(tp100) p100 = (double)tp100/100.0;
	    if(tp50) p50 = (double)tp50/50.0;
	    if(tp20) p20 = (double)tp20/20.0;
	    if(tp10) p10 = (double)tp10/10.0;
	    if(tp5) p5 = (double)tp5/5.0;
	    if(tp2) p2 = (double)tp2/2.0;
	    if(tp1) p1 = (double)tp1/1.0;	    

	    //cout << "# ncon = " << ncon << endl;
    	//cout << "# l1\tl2\trho\tncon\ndens.\t1\t2\t5\t10\t20\t50\t100" << endl;
		cout << "# " << setprecision(4) << l1 << "\t" << l2 << "\t" << rho << "\t" << ncon << "\t" << sumfnzero << "\t" << p1 << "\t" << p2 << "\t" << p5 << "\t" << p10 << "\t" << p20 << "\t" << p50 << "\t" << p100;
        if(all_theta){
            cout << "\t" << pcmean;
        }
        cout << endl;


    }	
    
        }
    }
	

    // Clean up
    freemat(pa,seqlen,alphabet_sz); 
    for(i = 0; i < nseqs; i++){
    	delete [] aln[i];
    }
    delete [] aln;
    delete [] pab;
    return(0);

}
