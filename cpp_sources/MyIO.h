#ifndef MYIO_H
#define MYIO_H
#include <vector>
#include <string> // "" for string, '' for char
#include <iostream> // cout/cin, ofstream: Stream class to write on files, ifstream : Stream class to read from files, istringstream is for input, ostringstream for output
#include <fstream> // fstream : Stream class to both read and write from / to files

#include <H5Cpp.h>
#include <hdf5_hl.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

using namespace std;

const char delim = ',';  /// Bad practice: this is hard-coded here!

void write_vector_HDF5(Group & group, const vector<bool> & v, const string & v_name);
void write_vector_HDF5(Group & group, const vector<unsigned int> & v, const string & v_name);
void write_vector_HDF5(Group & group, const vector<int> & v, const string & v_name);
void write_vector_HDF5(Group & group, const vector<double> & v, const string & v_name);
void append_vector_to_matrix_HDF5(DataSet & dataset_tmp, const vector<double> & v, const int colNum);
void append_vector_to_matrix_HDF5(DataSet & dataset_tmp, const vector<int> & v, const int colNum);
void write_matrix_HDF5(Group & group, const vector< vector<double> > & m, const string & m_name);
void write_matrix_HDF5(Group & group, const vector< vector<int> > & m, const string & m_name);
void write_matrix_HDF5(Group & group, const vector< vector<bool> > & m, const string & m_name);
	
void write_string_HDF5(Group & group, const string & s, const string & s_name);
void write_scalar_HDF5(Group & group, unsigned int s, const string & v_name);
void write_scalar_HDF5(Group & group, int s, const string & v_name);
void write_scalar_HDF5(Group & group, double s, const string & v_name);
	
void read_string_HDF5(const H5File & file, const string &  s_name,  string & s);
void read_matrix_HDF5(const H5File & file, const string & name, vector< vector<int> > & m_tmp);
void read_matrix_HDF5(const H5File & file, const string & name, vector< vector<bool> > & m_tmp);
void read_matrix_HDF5(const H5File & file, const string & name, vector< vector<double> > & m_tmp);
void read_vector_HDF5(const H5File & file, const string & name, vector<unsigned int> & v_tmp);
void read_vector_HDF5(const H5File & file, const string & name, vector<int> & v_tmp);
void read_vector_HDF5(const H5File & file, const string & name, vector<bool> & v_tmp);
void read_vector_HDF5(const H5File & file, const string & name, vector<double> & v_tmp);
bool dataset_exist_HDF5(const H5File & file, const string & name);
bool group_exist_HDF5(const H5File & file, const string & name);
bool group_exist_HDF5(const string & filename, const string & name);
template < typename Type > Type read_scalar_HDF5(const H5File & file, const string & name);

void write2file(ofstream& output_file, const vector< vector<int> >& v); /// write integer matrix to output file
void write2file(ofstream& output_file, const vector< vector<double> >& v); /// write double matrix to output file
void write2file(ofstream& output_file, const vector<int>& v); /// write integer vector to output file
void write2file(ofstream& output_file, const vector<double>& v); /// write double vector to output file

#endif