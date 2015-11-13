#include<cstdlib>
#include<cmath>
#include<fstream>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;

/**
 * Author: Carsten Kemena
 * 
 * This program was developed to turn a pfam domain vs domain similarity matrix into a sparse version. Only values above the threshold are kept and only the 
 * upper triangle matrix. Values are turned into integer values. The binary format of the ouput file looks as following: 
 * 
 * The file contains 6 fields:
 * Field_id	type	#elements	name		description
 * 1		char	10			name		The name of this matrix (and or version, etc)
 * 2 		int 	1 			n_domains 	Contains the number of domains in the matrix (=size of the matrix)
 * 3 		size_t	1			n_vals		The number of values stored
 * 4		int		n_domains	ids			The PFAM id as integer (eg: PF00001 = 1)
 * 5		int		n_domains	row_ids		The position from wich the row can be found in col_ids
 * 6		int		n_vals		col_ids		The columns position in the original matrix
 * 7		short	n_vals		values		The actual matrix entry
 * 
 */

int
main(int argc, char *argv[])
{
	bool show_help=false;
	for (int i=1; i<argc; ++i)
		if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help"))
			show_help=true;
	if (show_help || (argc!=5))
	{
		printf("USAGE: mat2bincrs.cpp <matrix file> <output file> <filter threshold> <name/identifier>\n\n");
		printf("This program converts the input matrix into a binary sparse matrix in CRS (Compressed Row Storage) format.\n\n Output file description:\n Field_id name type number\n 1 name char 10\n 2 n_rows int 1\n 3 n_vals int 1\n 4 names int n_rows\n 5 row_ids int n_rows\n 6 col_ids int n_vals\n 7 vals values short n_vals\n");
		if (show_help)
			exit(EXIT_SUCCESS);
		else
			exit(EXIT_FAILURE);
	}
	char *mat_f = argv[1];
	char *out_f = argv[2];
	int threshold = atoi(argv[3]);
	
	char name[11]="\0\0\0\0\0\0\0\0\0\0";
	strncpy(&name[0], argv[4], 11);
	ifstream mat_F(mat_f);
	string line;

	// read the pfam ids
	getline(mat_F, line);
	char_separator<char> sep(" \n\t");
	tokenizer< char_separator<char> > tokens(line, sep);
	vector<int> ids;
	BOOST_FOREACH (const string& t, tokens)
	{
		int acc=atoi(&t.c_str()[2]);
		ids.push_back(acc);
	}
	
	//read the values of the matrix row by row
	vector<int> row_ids, col_ids;
	vector<short> vals;
	row_ids.reserve(ids.size());
	int row_num=0;
	int col_num;
	float val;
	while (getline(mat_F, line))
	{
		col_num=0;
		row_ids.push_back(col_ids.size());
		tokenizer< char_separator<char> > tokens(line, sep);
		BOOST_FOREACH (const string& t, tokens)
		{
			if (col_num>=row_num)
			{
 				val=atof(t.c_str());
				// keep only matrix entries larger or equal to the threshold 
				if (val>=threshold)
				{
 					col_ids.push_back(col_num);
					vals.push_back(static_cast<short>(floor(val + 0.5f)));
				}
			}
			++col_num;
		}
		++row_num;
	}
	
	mat_F.close();

	// write data to binary file
	ofstream out_F;
    out_F.open (out_f, ios::out | ios::binary);
	int n_domains=ids.size();
	int n_vals=col_ids.size();

	out_F.write((char*)&name[0], sizeof(char)*10);
	out_F.write((char*)&n_domains, sizeof(int));
	out_F.write((char*)&n_vals,sizeof(int));
	out_F.write((char*)&ids[0],sizeof(int)*n_domains);
	out_F.write((char*)&row_ids[0],sizeof(int)*n_domains);
	out_F.write((char*)&col_ids[0],sizeof(int)*col_ids.size());
	out_F.write((char*)&vals[0],sizeof(short)*vals.size());
	out_F.close();
	
	return(EXIT_SUCCESS);
}
