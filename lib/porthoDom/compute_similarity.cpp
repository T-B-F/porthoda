/*
 * Compute the distance between two domain arrangement using cosine distance
 * and a similarity matrix based on hmm domain model similarity
 * 
 * Copyright (C) 2013 IEB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* informations :
 * char _author[] = "Tristan Bitard-Feildel, Carten Kemena";
 * char _email[] = "t.bitard.feildel@uni-muenster.de";
 * char _institute[] = "Insitute for Evolution and Biodiversity";
 * char _lab[] = "Evolutionary Bioinformatics";
 * char _vesion[] = "0.01";
 */

/* beginning of the program : */ // useless comment
#include <cstdio>
#include <cstdlib>
#include <cmath>
/* structures */
#include <utility>      // std::pair, std::make_pair
#include <algorithm>    // std::sort, std::min, std:;max
#include <tuple>
#include <vector>
#include <unordered_map>
#include <list>
#include <valarray>

/* string part */
#include <sstream>
#include <iostream>
#include <fstream>
/* old c part */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits>
/* boost */
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>


/* use namespace for lazyness*/
using namespace boost;

/* Some constant */
#define MINIMAL_SCORE 0.0f
#define EPSILON 0.0001f

/* CRS_MAT + function, code by carsten */

FILE *
my_fopen (const char *filename, const char *mode )
{
    FILE *tmp = fopen(filename, mode);
    if (tmp != NULL)
        return tmp;
    else
    {
        printf ("Error opening file %s: %s\n",filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
}

typedef struct
{
    char name[11];
    int n_domains;
    int n_vals;
    int *ids;
    int *row_ids;
    int *col_ids;
    short *vals;
} CRS_MAT;

int compareints (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

CRS_MAT*
read_crs_mat(char *mat_f)
{
    CRS_MAT *mat = (CRS_MAT*)malloc(sizeof(CRS_MAT));
    FILE *mat_F = my_fopen (mat_f , "rb");
    int error=0;
    error += fread(&mat->name, sizeof(char), 10, mat_F);
    error += fread(&mat->n_domains, sizeof(int), 1, mat_F);
    error += fread(&mat->n_vals, sizeof(int), 1, mat_F);
    mat->ids=(int*)malloc(sizeof(int)*mat->n_domains);
    mat->row_ids=(int*)malloc(sizeof(int)*mat->n_domains);
    mat->col_ids=(int*)malloc(sizeof(int)*mat->n_vals);
    mat->vals=(short*)malloc(sizeof(short)*mat->n_vals);
    error += fread(mat->ids, sizeof(int), mat->n_domains, mat_F);
    error += fread(mat->row_ids, sizeof(int), mat->n_domains, mat_F);
    error += fread(mat->col_ids, sizeof(int), mat->n_vals, mat_F);
    error += fread(mat->vals, sizeof(short), mat->n_vals, mat_F);
    if (error != (12+2*mat->n_domains+2*mat->n_vals))
    {
        printf("ERROR when reading matrix file!\n");
        exit(1);
    }
    fclose(mat_F);
    return mat;
}

void
free_crs_mat(CRS_MAT *crs_mat)
{
    free(crs_mat->ids);
    free(crs_mat->row_ids);
    free(crs_mat->col_ids);
    free(crs_mat->vals);
    free(crs_mat);
}

short
get_crs_val(int i, int j, CRS_MAT* mat)
{
    int *id1p, *id2p;
    if (i<j)
    {
        id1p = (int*) bsearch (&i, mat->ids, mat->n_domains, sizeof (int), compareints);
        id2p = (int*) bsearch (&j, mat->ids, mat->n_domains, sizeof (int), compareints);
    }
    else
    {
        id1p = (int*) bsearch (&j, mat->ids, mat->n_domains, sizeof (int), compareints);
        id2p = (int*) bsearch (&i, mat->ids, mat->n_domains, sizeof (int), compareints);
    }
    int id1,id2;
    if ((id1p == NULL) || (id2p == NULL))
        return -1;
    else
    {
        id1=id1p-mat->ids;
        id2=id2p-mat->ids;
    }
    int id_col=mat->row_ids[id1];
    int end=(id1!=(mat->n_domains-1)) ? mat->row_ids[id1+1] : mat->n_vals;
    for (int i =id_col; i<end; ++i)
    {
        if (mat->col_ids[i]==id2)
            return mat->vals[i];
    }
    return -1;
}

/* END of Carsten code */
/* need some comments, no ? */

/* similarity part */
float cosineSimilarity( std::vector<float> vec1, std::vector<float> vec2 ) {
    /*
    Cosine index of two vector
    \param vec1 a vector of float, ie scores
    \param vec2 a vector of float, ie scores
    \return 1 - dot(vec1,vec2)/(norm(vec1)*norm(vec2))
    */
    float  result = 0.0 ;
    assert ( vec1.size() == vec2.size() );
    float numerator = 0;
    float denoma = 0;
    float denomb = 0;
    float ai, bi ;
    for ( unsigned int i = 0 ; i < vec1.size(); i++ ) {
        ai = vec1[i] ;
        bi = vec2[i] ;
        numerator += ai*bi ;
        denoma += ai*ai ;
        denomb += bi*bi ;
    }
    result = 1.0 - ( numerator / (sqrt(denoma)*sqrt(denomb)) );
    return result ;
}

float euclideanSimilarity( std::vector<float> vec1, std::vector<float> vec2 ) {
    /*
    Cosine index of two vector
    \param vec1 a vector of float, ie scores
    \param vec2 a vector of float, ie scores
    \return 1 - dot(vec1,vec2)/(norm(vec1)*norm(vec2))
    */
    float  result = 0.0 ;
    assert ( vec1.size() == vec2.size() );
    float score = 0;
    for ( unsigned int i = 0 ; i < vec1.size(); i++ ) {
        score += ( (vec1[i]-vec2[i]) * (vec1[i]-vec2[i]) );
    }
    
    result = sqrt( score ) ;
    return result ;
}

std::vector< std::tuple< std::string,std::string > > domPairSet( std::vector<std::string> protein, int order ) {
    /*
    Create a set of domain pair of a specified order from a protein
    \param protein1
    \param order
    \return pairsOfDom
    */
    std::vector< std::tuple< std::string,std::string > > pairsOfDom ;    
    int cnt;    
    for (unsigned int i = 0 ; i < protein.size( ) - order ; i ++ ) {
        std::tuple<std::string,std::string> tmp( protein[i], protein[i+order]) ;
        cnt = std::count( pairsOfDom.begin(), pairsOfDom.end(), tmp ) ; 
        if ( cnt == 0 )
            pairsOfDom.push_back( tmp );
    }
    return pairsOfDom ;
}

std::vector< std::string> domainSet( std::vector<std::string> protein, int order ) {
    /*
    Create a set of domain pair of a specified order from a protein
    \param protein1
    \param order
    \return vectorOfDom
    */
    std::vector< std::string > vectorOfDom ;    
    int cnt;
    /*std :: cout << "@" << order << " "  ;
    for (unsigned int i = 0 ; i < protein.size( ); ++i ) {
        std::cout << protein[i] << " " << i << " "   ;
    }
    std :: cout << std::endl;
    std :: cout << "]" ;*/
    for (unsigned int i = 0 ; i < protein.size( ) - order + 1; ++i ) {
        std::string tmp = protein[i];
        //std::cout << protein[i] << " " << i << " "   ;

        for (unsigned int j = i+1; j < i+order  ; ++j ) {
            //std::cout << i << "  " << j << " " << i+order << std::endl;
            tmp += "_"+protein[j] ;
        }
        //std::cout << "> " << tmp << std::endl;
        cnt = std::count( vectorOfDom.begin(), vectorOfDom.end(), tmp ) ; 
        if ( cnt == 0 )
            vectorOfDom.push_back( tmp );
    }
    //std :: cout << std::endl;

    return vectorOfDom ;
}


float findMaxWeight( std::tuple<std::string,std::string> domA, 
                     std::vector<std::tuple<std::string,std::string>> domSet, 
                     CRS_MAT *mat ){
    /* 
    \param domA the list of domain in a protein A 
    \param domSet the list of all domain
    \param mat the distance matrix between domains 
    \return current_score
    */
    int idom1, idom2, jdom1, jdom2; //, best_jdom1 = -1  ;
    float w1, w2, w; 
    float current_score = 0.0 ;
    
    std::vector<std::tuple<std::string,std::string>>::iterator it;
    
    idom1 = atoi( &(std::get<0>(domA).c_str( ))[2] ) ;
    idom2 = atoi( &(std::get<1>(domA).c_str( ))[2] ) ;
    
    for (it = domSet.begin(); it != domSet.end(); ++it ) { 
        jdom1 = atoi( &( std::get<0>(*it).c_str( ))[2] );
        jdom2 = atoi( &( std::get<1>(*it).c_str( ))[2] );
        /* Warning : score in the matrix sould be \in [0,1] */
        w1 = get_crs_val( idom1, jdom1, mat ) ;
        if ( w1 == -1 ) {
            w1 = MINIMAL_SCORE ;
        } else  {
            w1 /= 100. ;
        }
        w2 = get_crs_val( idom2, jdom2, mat ) ;
        if (w2 == -1 ) {
            w2 = MINIMAL_SCORE ;
        } else {
            w2 /= 100. ;
        }
        w = (w1+w2)/2.0 ;
        if (w > current_score ) {
            current_score = w ;
            //best_jdom1 = jdom1 ;
        }
    }
    //std::cout << idom1 << " " << best_jdom1 << " " << current_score << std::endl;
    return current_score;
}

inline float normalize_0_1( float a, float x ) {
    return (pow(a,x) - 1.0) / (a - 1.0 ) ;
}

std::unordered_map< std::string, float > computeSetMatrix( 
                      std::vector< std::string> domain_set, 
                      int order, CRS_MAT *mat, bool with_weight ){
    /* Compute all pairwise score
     * \param domain_set
     * \param order
     * \param mat
     * \param with_weight
     * \return matrix
     */
    
    std::unordered_map< std::string, float > matrix;
    char_separator<char> char_sep(" \n\t_");
    float score_rev = 0.0, score = 0.0, score_nor = 0.0, w, w_rev;
    
    /* for each domain in the domain set used */
    /* a "domain" in a domain set is a set of domain, examples :
     *     dom1
     *     dom1_dom2
     *     dom1_dom3_dom1 
     *     etc ...
     */
    for(unsigned int i = 0; i< domain_set.size( ); i++) { 
        tokenizer< char_separator<char> > char_tokens(domain_set[i], char_sep);
        std::vector<int> dai ;
        /* cut the domain set into individual domain and fill dai vector */
        BOOST_FOREACH (const std::string& tmp, char_tokens) {
            dai.push_back( atoi( &tmp.c_str()[2] ) );
        }
        /* try to match the reverse word */
        std::vector<int> dai_rev( dai.begin(), dai.end()) ;
        std::reverse( dai_rev.begin(), dai_rev.end() );
        /* compare to all other words in the domain set */
        for(unsigned int j = i+1; j < domain_set.size() ; j++) {
            tokenizer< char_separator<char> > char_tokens(domain_set[j], char_sep);
            std::vector<int> daj ;       
            /* same  as dai, cut the domain and fill daj */
            BOOST_FOREACH (const std::string& tmp, char_tokens) {
                daj.push_back( atoi( &tmp.c_str()[2] ) );
            }
            /* set score for both reverse and normal */
            score = 0.0;
            score_nor = 0.0;
            score_rev = 0.0;
            /* check pairwise matching */
            for( int k=0; k < order; k++ ){ // order start from 0 
                //std::cout << dai[k] << " " << daj[k] << std::endl ;
                if (dai[k] == daj[k] ) {
                    w = 1.0;
                } else {
                    w = get_crs_val( dai[k], daj[k], mat );
                    if (w == -1 ) { 
                        w = MINIMAL_SCORE;
                    } else { 
                        w /= 100.;
                    }
                }
                score_nor += w ;
                if (dai_rev[k] == daj[k] ) {
                    w_rev = 1.0;
                } else {
                    w_rev = get_crs_val( dai_rev[k], daj[k], mat );
                    if (w_rev == -1 ) { 
                        w_rev = MINIMAL_SCORE;
                    } else { 
                        w_rev /= 100.;
                    }
                }
                score_rev += w_rev ;
                /*std::cout << k << " " << domain_set[i] << " " << domain_set[j] << " ";
                std::cout << dai[k] << " " << " " << dai_rev[k] << " " << daj[k] << " ";
                std::cout << w << " " << w_rev << std::endl;*/
            }
            score = std::max( score_nor, score_rev);
            score /= float(order+1); // order start from 0 
            //std::cout << order+1.0+EPSILON << " "<< score << " ";
            /* 
             * score is weighted according to the order, 
             * when the order increase (more domain in the word),
             * the stringence for domain pair matching increase also
             * 
             *     a = pow( 10, order + 1 )
             * 
             * score = (( a ^ score ) - 1 ) / ( a - 1 )
             * 
             */
            if (with_weight)
                score = normalize_0_1( pow(10,order+1), score ); // (order+1.0+EPSILON), score );
            //std::cout << score << std::endl;
            //std::cout << domain_set[i] << " "<< domain_set[j] << " " << score << std::endl;
            matrix.insert( std::make_pair( domain_set[i]+"|"+domain_set[j], score ) );
            matrix.insert( std::make_pair( domain_set[j]+"|"+domain_set[i], score ) );
        }
    }
    return matrix;    
}

float  weightedCosine( std::vector<std::string> protein1, 
                       std::vector<std::string> protein2, 
                       CRS_MAT * mat, unsigned int order,
                       bool with_weight ) { 
    /*
    Compute a score depending of the best pair similarity between
    domain pairs 
    \param protein1 list of domains in protein 1
    \param protein2 list of domains in protein 2 
    \param mat the distance matrix between domains 
    \param order the pairing order (0 mean single, 1 mean direct pair, 2 mean pair at +2 position etc... )
    \param with_weight use or not the weighting scheme
    \return score a floating score
    */
    std::vector< std::string > setOfP1 ;
    std::vector< std::string > setOfP2 ;
    std::vector< std::string >::iterator its;
    std::vector<float> vec1;
    std::vector<float> vec2;
    std::unordered_map< std::string, float > matrix_set;
    float score = 1.0;
    //std::cout << "score : " <<score << std::endl;
    setOfP1 = domainSet( protein1, order );
    setOfP2 = domainSet( protein2, order );   
    
    
    std::vector< std::string > domains_set2;
    std::vector< std::string >::iterator it2 ;
    std::vector< std::string >::iterator it2_ ;
    
    for(its = setOfP1.begin( ); its!= setOfP1.end(); ++its ) { 
        
        if ( std::count( domains_set2.begin(), domains_set2.end(), *its ) == 0 ) {
            domains_set2.push_back( *its );
            //std :: cout << *its << " ";
        }
    }
    //std::cout << std::endl;
    for(its = setOfP2.begin( ); its!= setOfP2.end(); ++its ) {
        //std::cout << *its << std::endl ;;
        if (std::count( domains_set2.begin(), domains_set2.end(), *its ) == 0 ) {
            domains_set2.push_back( *its );
        }
    }
    
    /* TODO compute matrix sim of domains set and store them into an hash */
    matrix_set = computeSetMatrix( domains_set2, order, mat, with_weight );
    float max_score = 1.0;
    vec1.clear( );
    vec2.clear( );
    for ( it2 = domains_set2.begin( ); it2!= domains_set2.end() ; ++it2) {
            
        if ( std::count( setOfP1.begin(), setOfP1.end(), *it2 ) > 0 ) {
            vec1.push_back( 1.0 ) ;
        } else {
            max_score = 0.0;
            for( it2_ = setOfP1.begin(); it2_ != setOfP1.end(); ++it2_ ){
                if ( matrix_set[ (*it2)+"|"+(*it2_) ] > max_score )
                    //std::cout<< "1 "+(*it2)+"|"+(*it2_) << std::endl;
                    max_score = matrix_set[ (*it2)+"|"+(*it2_) ] ;
            }
            vec1.push_back(  max_score  );
        }
        //std::cout << ">" << *it2 << " "<< std::count( setOfP2.begin(), setOfP2.end(), *it2 ) << std::endl ;
        if ( std::count( setOfP2.begin(), setOfP2.end(), *it2 ) > 0 ) {
            vec2.push_back( 1.0 ) ;
        } else {
            max_score = 0.0;
            for( it2_ = setOfP2.begin(); it2_ != setOfP2.end(); ++it2_ ){
                if ( matrix_set[ (*it2)+"|"+(*it2_) ] > max_score )
                    //std::cout<< "2 "+(*it2)+"|"+(*it2_) << " " << max_score ; 
                    max_score = matrix_set[ (*it2)+"|"+(*it2_) ] ;
                    //std::cout << " " << max_score << std::endl;
            }
            vec2.push_back(  max_score  );
        }
    }
    /*
    for( unsigned int l=0 ; l<vec1.size(); l++) 
        std::cout<< vec1[l] << " " ;
    std::cout<<std::endl;
    for( unsigned int l=0 ; l<vec2.size(); l++) 
        std::cout<< vec2[l] << " " ;
    std::cout<<std::endl;
    */
    score = cosineSimilarity( vec1, vec2 ) ;
    //std::cout << "score : " <<score << std::endl;
    return 1.0 - score ;
}
    
void usage( ) {
    std::cout <<"compute_similarity [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "   -i string  path to a file containing in two columns format the pairwise DA to compare" << std::endl;
    std::cout << "   -m string  path to a hmm model matrix stored in a Compressed Row Storage format" << std::endl;
    std::cout << "   -c float   only display similarity superior to this cutoff in the output" << std::endl;
    std::cout << "   -o order   order of magnitude use for domain pair comparison (default:0, no order, self domain comparison) " << std::endl;
    std::cout << "              for a protein with domain ABCD, order 0 is A B C D, order 1 : AB BC CD, order 2 : ABC BCD, etc ..." << std::endl;
    std::cout << "              order cannot bit larget than the smallest length of protein in domain" << std::endl;
    std::cout << "   -v verbose choose to display or not the messages [T/F] default = F " << std::endl;
    std::cout << "   -w weight  choose to add weight or not to the final score [T/F] default = F " << std::endl;
    std::cout << "   -h         display this help" << std::endl ;
}

int main ( int argc, char **argv ) {
    /* argument parser */
    int c ;
    unsigned int order=0, ori_order;
    char *pathlist=NULL, *pathmatrix=NULL;
    float cutoff=-1.0;
    /* main */
    CRS_MAT *mat;
    std::ifstream listfile;
    std::string buffer;
    char_separator<char> char_sep(" \n\t");
    char_separator<char> dom_sep(";");
    float sim;
    bool weight=false;
    bool verbose=false;
    
    /* parse arguments */
    while ( ( c = getopt (argc, argv, "i:m:c:o:v:w:h:" ) ) != -1 ) {
        switch (c) {
        case 'c':
            cutoff = atof(optarg) ;
            break ;
        case 'i' :
            pathlist = optarg ;
            break ;
        case 'm':
            pathmatrix = optarg ;
            break ;
        case 'o':
            order = atoi(optarg) ;
            break ; 
        case 'v':
            if (*optarg =='T')
                verbose = true;
            break ;
        case 'w' :
            if(*optarg=='T')
                weight=true;
            break;
        case 'h' :
            usage();
            return EXIT_SUCCESS ;
        default:
            usage( ) ;
            return EXIT_FAILURE ;
        }
    }
    
    /* some check for the existence of the files */
    if ( pathlist==NULL || pathmatrix==NULL || cutoff==-1.0 ) {
        std::cerr << "Parsing option problem, make sure that the paths to the list, to the matrix and the cutoff paramterers are initialized" << std::endl;
        usage( ) ;
        return EXIT_FAILURE;
    }
    
    /* memorize the order */
    ori_order = order ;
    
    /* read the matrix file */
    mat = read_crs_mat( pathmatrix ) ;
    
    /* read the pairwise file */
    listfile.open ( pathlist ) ; 
    if ( !listfile.is_open( ) ) { 
        std::cerr << "Error : Unable to open domain architecture "<<pathlist<< std::endl ;
        exit(1) ;
    }
    /* for each pairwise compute the similarity */
    
    while ( getline( listfile, buffer ) ) {
        tokenizer< char_separator<char> > char_tokens(buffer, char_sep);
        std::vector<std::string> das;
        BOOST_FOREACH (const std::string& tmp, char_tokens) {
            das.push_back( tmp );
        }
        std::vector<std::string> dai ;
        std::vector<std::string> daj ;
        tokenizer< char_separator<char> > domi_tokens(das[0], dom_sep);
        BOOST_FOREACH (const std::string& tmp, domi_tokens) {
            dai.push_back( tmp );
        }
        tokenizer< char_separator<char> > domj_tokens(das[1], dom_sep);
        BOOST_FOREACH (const std::string& tmp, domj_tokens) {
            daj.push_back( tmp );
        }

        /* check if order is appropriate for protein */
        if( order > dai.size( ) || order > daj.size( ) ) { 
            if ( verbose ) { 
                std::cerr << "Warning, the order pair is superior to the proteins size analysed" << std::endl ;
                std::cerr << "         using minimal size order instead" << std::endl ;
            }
            order = std::min( dai.size(), daj.size() ) ;
        }
        /* compute similarity */
        sim = weightedCosine( dai, daj, mat, order, weight);
        sim *= (order+1.0) / (ori_order+1.0);
        /* only display if sim > cutoff */
        if ( sim > cutoff )
            std::cout << das[0] << " " << das[1] << " " << sim << " " << order << std::endl ;
        /* back to original order */
        order = ori_order;
    }

    listfile.close( ) ;
    /* free matrix */
    free_crs_mat( mat );

    return EXIT_SUCCESS;  
} 

