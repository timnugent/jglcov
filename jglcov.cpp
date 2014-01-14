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

int main(int argc, char **argv){

	double wtsum, pcmean, pc, idthresh = -1.0, pseudoc = 1.0, sumfnzero = 0.0;
	unsigned int nseqs, ncon = 0, opt, seqlen, i, j, k, a, b, npair, nnzero, shrinkflg = 1, minseqsep = 5, apcflg = 0, wts = 1, alphabet = 1, alphabet_sz = 21;
	int o = 1;
	double rho = 1e-06;
	double l1  = 1e-08;
	double l2  = 1e-02;
	char seq[MAXSEQLEN];
	char *native = NULL;
	string jgl_weight = "\"equal\"";
	string penalty = "group";

	while(o < argc){
        if(argv[o][0] == '-'){
            o++;
            switch(argv[o-1][1]){
            		case 'p' : {penalty = "fused"; o--; cout << "#Using fused penalty" << endl; break;}
            		case 'c' : {native = argv[o]; break;}
            		case 'r' : {rho = atof(argv[o]); break;}
                    case 'm' : {l1 = atof(argv[o]); break;}
                    case 'n' : {l2 = atof(argv[o]); break;}
                    case 'a' : {alphabet = atoi(argv[o]); cout << "#Using alphabet " << alphabet << endl; break;}
                    case 'w' : {wts = 1; o--; cout << "#Using weights" << endl; break;}
                    default  : {fail("Options:\t-m float lamda1 penalty\n-n float lambda2 penalty\n-w use weights\n-c file native contacts\n");}
                }   
           }
       o++;
    }

    if(alphabet) alphabet_sz = 6;
    unsigned int classes = (alphabet_sz-1)*(alphabet_sz-1);

    cout << "#lambda1 = " << l1 << endl;
    cout << "#lambda2 = " << l2 << endl;
    cout << "#aln     = " << argv[o-1] << endl;

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
				pab[index4D(i,j,a,b,alphabet_sz,seqlen)] = (a == b) ? pa[i][a] : 0.0;
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
	double max_cond = -100000.0;

    for (a = k = 0; a < alphabet_sz-1; a++){
        for (b = 0; b < alphabet_sz-1; b++){
    
        	Rcpp::NumericMatrix rMat(seqlen,seqlen);
        	//cout << "a = " << a << " --- b = " << b << endl;

    		for (i = 0; i < seqlen; i++){
        		for (j = 0; j < seqlen; j++){
                    if (i != j){
                    	rMat(i,j) = pab[index4D(i,j,a,b,alphabet_sz,seqlen)] - pa[i][a] * pa[j][b];
	                }        
    			}
    		}
    		rMat.attr("dimnames") = dimnms;		    
		    mylist[k] = rMat;
		    
		    // Shrinkage
		    if(shrinkflg){
		    	
		    	R["X"] = rMat;
			    
			    /*
			    int posdef = R.parseEval("is.positive.definite(X)");
			    cout << "(pre)Positive definite: " << posdef << endl;
			    Rcpp::List condition = R.parseEval("rank.condition(X)");
			    int rank = Rcpp::as<int>(condition[0]);
			    double cond = (condition[1]);
			    cout << "(pre)Rank: " << rank << endl;
			    cout << "(pre)Condition: " << cond << endl;
			    conut << "---shrinking---" << endl;
				time_t start,end;
			    time (&start);
			    */
			    
			    R.parseEval("s2 = cov.shrink(X,verbose=FALSE)");

			    int posdef = R.parseEval("is.positive.definite(s2)");
			    //cout << "Positive definite: " << posdef << endl;
			    Rcpp::List condition_shrunk = R.parseEval("rank.condition(s2)");
			    int rank = Rcpp::as<int>(condition_shrunk[0]);
			    //double cond = Rcpp::as<double>(condition_shrunk[1]);
			    weights[k] = Rcpp::as<double>(condition_shrunk[1]);
			    if(Rcpp::as<double>(condition_shrunk[1]) > max_cond){
			    	max_cond = Rcpp::as<double>(condition_shrunk[1]);
			    }
			    //cout << "Rank: " << rank << endl;
			    //cout << "Condition: " << cond << endl;
			    Rcpp::NumericMatrix M = R.parseEval("s2");
		    	mylist[k] = M;
	    	}

			/*
			cout << "Matrix " << k+1 << ":" << endl;
    		for (i = 0; i < seqlen; i++){
        		for (j = 0; j < seqlen; j++){
                    if (i != j){
                    	cout << M(i,j) << " ";
	                }        
    			}
    			cout << endl;
    		}
			*/

		    k++;
		    //time (&end);
		    //double dif = difftime (end,start);
		    //cout << "Shrinkage runtime:\t" << dif << " seconds" << endl;


    	}	
    }	
    cout << "#Finished shrinking covariance matrices" << endl;

    if(wts) jgl_weight = "weights";
	for(k = 0; k < classes; k++){    	
		//cout << "# " << k << " - " << weights[k] << " ---> ";    		
		if(weights[k] == max_cond){
			// Can't be zero
			weights[k] = 1e-5;
		}else{
			weights[k] = 1.0-(weights[k]/max_cond);
		}
		//cout << weights[k] << " (max = " << max_cond << ")" << endl;
	}
    cout << "#Finished calculating weights" << endl;

    sc_entry* sclist = new sc_entry[seqlen * (seqlen - 1) / 2]; 
    time_t start,end;
    Rcpp::List theta;
	R["data"] = mylist;
	R["weights"] = weights;


	cout << "#Running JGL..." << endl;


	/*
	time (&start);
	r_code = "JGLx(Y=data,penalty=\"" + penalty + "\",lambda1=" + to_string((long double)l1) + ",lambda2=" + to_string((long double)l2) + ",rho=" + to_string((long double)rho) + ",weight=" + jgl_weight + ",penalize.diagonal=FALSE,maxiter=500,tol=1e-5,warm=NULL,return.whole.theta=TRUE,screening=\"fast\",truncate = 1e-5)";
	Rcpp::List ret = R.parseEval(r_code);	
	time (&end);	
	//cout << "Finished group lasso (" << (int)difftime (end,start) << " seconds)" << endl;
	
	theta = ret[0];
	for(k = 0; k < classes; k++){
		//http://stackoverflow.com/questions/12719334/how-to-handle-list-in-r-to-rcpp
		SEXP ll = theta[k];
        Rcpp::NumericMatrix y(ll);				
	    double** pcmat = allocmat(seqlen,seqlen);
		double* pcsum = new double[seqlen];
		pcmean = 0.0;


		//cout << "Matrix " << k+1 << endl;
 		//for (i = 0; i < seqlen; i++){
        //    for (j = 0; j < seqlen; j++){
        //    	cout << y(i,j) << " ";
        //    }
        //    cout << endl;
        //}

 		for (npair = nnzero = i = 0; i < seqlen; i++){
            for (j = i+1; j < seqlen; j++, npair++){
                if (y(i,j) > 0.0){
                    nnzero++;
                }
                pcmat[i][j] = pcmat[j][i] = fabs(y(i,j));
                //pcmat[i][j] = pcmat[j][i] = -y(i,j)/sqrt(y(i,i)*y(j,j));
                pcsum[i] += pcmat[i][j];
    			pcsum[j] += pcmat[i][j];
				pcmean += pcmat[i][j];
            }
 		}

 			
		//cout << "PCMatrix " << k+1 << endl;
 		//for (i = 0; i < seqlen; i++){
        //    for (j = 0; j < seqlen; j++){
        //    	cout << pcmat(i,j) << " ";
        //    }
        //    cout << endl;
        //}
        //cout << endl;
		
 		pcmean /= seqlen * (seqlen - 1) * 0.5;
        double fnzero = (double) nnzero / npair;
        sumfnzero += fnzero;

        sc_entry* sclist_k = new sc_entry[seqlen * (seqlen - 1) / 2]; 

        // Calculate APC score
	    for (ncon = i = 0; i < seqlen; i++){
	        for (j=i+minseqsep; j<seqlen; j++){
	            //if (pcmat(i,j) > 0.0){			

	            	if(apcflg){
	                	sclist[ncon].sc += pcmat[i][j] - pcsum[i] * pcsum[j] / sqrt(seqlen-1.0) / pcmean;
	                	sclist_k[ncon].sc += pcmat[i][j] - pcsum[i] * pcsum[j] / sqrt(seqlen-1.0) / pcmean;
	                }else{
	                	sclist[ncon].sc += pcmat[i][j];
	                	sclist_k[ncon].sc += pcmat[i][j];
	                }
	                sclist[ncon].i = i;
	                sclist[ncon].j = j;
	                sclist_k[ncon].i = i;
	                sclist_k[ncon++].j = j;

	            //}
	        }    
	    } 
        delete [] pcsum;

        if(native_contacts.size()){

			int tp100 = 0;	
			int tp50 = 0;
			int tp20 = 0;
			int tp10 = 0;
			int tp5 = 0;
			int tp2 = 0;
			int tp1 = 0;	
			double p100 = 0.0;	
			double p50 = 0.0;
			double p20 = 0.0;
			double p10 = 0.0;
			double p5 = 0.0;
			double p2 = 0.0;
			double p1 = 0.0;

			qsort(sclist_k, ncon, sizeof(struct sc_entry), cmpfn);
		    for (i = 0; i < ncon; i++){
		        if(sclist_k[i].sc > 0.0){
		        	//printf("%d %d 0 8 %f\n", sclist_k[i].i+1, sclist_k[i].j+1, sclist_k[i].sc);
					string con = to_string(sclist_k[i].i+1);
					con.append("-");
					con.append(to_string(sclist_k[i].j+1));
					if (sclist_k[i].sc > 0.0 && find(native_contacts.begin(),native_contacts.end(),con) != native_contacts.end()){
						if(i < 100)tp100++;	
						if(i < 50)tp50++;
						if(i < 20)tp20++;
						if(i < 10)tp10++;
						if(i < 5)tp5++;
						if(i < 2)tp2++;
						if(i < 1)tp1++;
					}
		        }	
		    }

		    if(tp100) p100 = (double)tp100/100.0;
		    if(tp50) p50 = (double)tp50/50.0;
		    if(tp20) p20 = (double)tp20/20.0;
		    if(tp10) p10 = (double)tp10/10.0;
		    if(tp5) p5 = (double)tp5/5.0;
		    if(tp2) p2 = (double)tp2/2.0;
		    if(tp1) p1 = (double)tp1/1.0;

		    if(!k) cout << "#k\twt\tdens.\t1\t2\t5\t10\t20\t50\t100" << endl;
		    cout << "#" << k << "\t" << setprecision(4) << weights[k] << "\t" << fnzero << "\t" << p1 << "\t" << p2 << "\t" << p5 << "\t" << p10 << "\t" << p20 << "\t" << p50 << "\t" << p100 << endl;

		}

        delete [] sclist_k;
        freemat(pcmat,seqlen,seqlen);    
    }
    sumfnzero /= (double)classes;
    cout << "#" << l1 << " " << l2 << " - mean fnzero = " << sumfnzero << "\t(" << (int)difftime (end,start) << " seconds)" << endl;

	qsort(sclist, ncon, sizeof(struct sc_entry), cmpfn);
    for (i = 0; i < ncon; i++){
        if(sclist[i].sc > 0.0){
        	printf("%d %d 0 8 %f\n", sclist[i].i+1, sclist[i].j+1, sclist[i].sc);
        	if(native_contacts.size()){
				string con = to_string(sclist[i].i+1);
				con.append("-");
				con.append(to_string(sclist[i].j+1));
				if (find(native_contacts.begin(),native_contacts.end(),con) != native_contacts.end()){
					if(i < 100)tp100++;	
					if(i < 50)tp50++;
					if(i < 20)tp20++;
					if(i < 10)tp10++;
					if(i < 5)tp5++;
					if(i < 2)tp2++;
					if(i < 1)tp1++;
				}
			}        	
        }	
    }

    if(native_contacts.size()){
	    if(tp100) p100 = (double)tp100/100.0;
	    if(tp50) p50 = (double)tp50/50.0;
	    if(tp20) p20 = (double)tp20/20.0;
	    if(tp10) p10 = (double)tp10/10.0;
	    if(tp5) p5 = (double)tp5/5.0;
	    if(tp2) p2 = (double)tp2/2.0;
	    if(tp1) p1 = (double)tp1/1.0;
    	    cout << "# dens.\t1\t2\t5\t10\t20\t50\t100" << endl;
	    cout << "# " << setprecision(4) << sumfnzero << "\t" << p1 << "\t" << p2 << "\t" << p5 << "\t" << p10 << "\t" << p20 << "\t" << p50 << "\t" << p100 << endl;
    }	

	*/

	time (&start);
	r_code = "ggl.results = JGLx(Y=data,penalty=\"" + penalty + "\",lambda1=" + to_string((long double)l1) + ",lambda2=" + to_string((long double)l2) + ",rho=" + to_string((long double)rho) + ",weight=" + jgl_weight + ",penalize.diagonal=FALSE,maxiter=500,tol=1e-5,warm=NULL,return.whole.theta=FALSE,screening=\"fast\",truncate = 1e-5)";
	R.parseEvalQ(r_code);
	time (&end);	
	cout << "#Finished group lasso (" << (int)difftime (end,start) << " seconds)" << endl;
    
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

	int tp100 = 0;	
	int tp50 = 0;
	int tp20 = 0;
	int tp10 = 0;
	int tp5 = 0;
	int tp2 = 0;
	int tp1 = 0;	

    for(i = 0; i < dim; i++){
    	for (j = i+1; j< dim; j++){
            //cout << y(i,j) << " ";
            if(i != j && y(i,j) != 0){
		    sclist[ncon].sc += fabs(y(i,j));
	            sclist[ncon].i = atoi(Rcpp::as<string>(v[i]).c_str());
	            sclist[ncon++].j = atoi(Rcpp::as<string>(v[j]).c_str());
           }
        }
        //cout << endl;    
    }

    qsort(sclist, ncon, sizeof(struct sc_entry), cmpfn);
    for (i = 0; i < ncon; i++){
    	printf("%d %d 0 8 %f\n", sclist[i].i, sclist[i].j, sclist[i].sc);
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
		}        	
    }


    sumfnzero = (double)ncon/(seqlen * (seqlen - 1) * 0.5);

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
	    cout << "# ncon = " << ncon << endl;
    	    cout << "# dens.\t1\t2\t5\t10\t20\t50\t100" << endl;
	    cout << "# " << setprecision(4) << sumfnzero << "\t" << p1 << "\t" << p2 << "\t" << p5 << "\t" << p10 << "\t" << p20 << "\t" << p50 << "\t" << p100 << endl;
    }	

    // Clean up
    freemat(pa,seqlen,alphabet_sz); 
    for(i = 0; i < nseqs; i++){
    	delete [] aln[i];
    }
    delete [] aln;
    delete [] pab;
    delete [] sclist;
    return(0);

}
