#include "MyIO.h"

using namespace std;

// Use function templates when you want to perform the same action on types that can be different.
// Use function overloading when you want to apply different operations depending on the type.
// In this case, just save yourself the trouble and use overloading.
void write2file(ofstream& output_file,  const vector< vector<int> >& v){
	if (!v.empty()){
		for (unsigned int i = 0; i < v.size(); ++i){
			//for (double f : v[i]){ output_file << f << delim; } // range-based "for" in C++11
			for (unsigned int j = 0; j < v[i].size(); ++j){
				output_file << v[i][j] << delim;
			}
			output_file << endl;
		}
	}
	else {output_file << " " << endl;}
}

void write2file(ofstream& output_file, const vector< vector<double> >& v){
	if (!v.empty()){
		for (unsigned int i = 0; i < v.size(); ++i){
			//for (double f : v[i]){ output_file << f << delim; } // range-based "for" in C++11
			for (unsigned int j = 0; j < v[i].size(); ++j){
				output_file << v[i][j] << delim;
			}
			output_file << endl;
		}
	}
	else {output_file << " " << endl;}
}


void write2file(ofstream& output_file, const vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}


void write2file(ofstream& output_file, const vector<double>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}






/*--------------------------------------- HDF5 --------------------------------------------*/
#ifdef HDf5
void write_scalar_HDF5(Group & group, int s, const string & v_name){
	vector<int> v_tmp;
	v_tmp.push_back(s);
	write_vector_HDF5(group, v_tmp, v_name);
}

void write_scalar_HDF5(Group & group, double s, const string & v_name){
	vector<double> v_tmp;
	v_tmp.push_back(s);
	write_vector_HDF5(group, v_tmp, v_name);
}


void write_vector_HDF5(Group & group, const vector<int> & v, const string & v_name){
	hsize_t dims[1]; 
	dims[0] = v.size();
	DataSpace fspace(1, dims); 
	DataSet v_dataset = group.createDataSet(v_name, PredType::NATIVE_INT32, fspace);
	v_dataset.write( v.data(), PredType::NATIVE_INT32, fspace, fspace );	
}

void write_vector_HDF5(Group & group, const vector<double> & v, const string &   v_name){
	hsize_t dims[1]; 
	dims[0] = v.size();
	DataSpace fspace(1, dims); 
	DataSet v_dataset = group.createDataSet(v_name, PredType::NATIVE_DOUBLE, fspace);
	v_dataset.write( v.data(), PredType::NATIVE_INT32, fspace, fspace );	
}

void write_string_HDF5(Group & group, const string & s, const string &  s_name){
   // HDF5 only understands vector of char* :-(
   vector<const char*> arr_c_str;
   arr_c_str.push_back(s.c_str());

   hsize_t str_dimsf[1] {arr_c_str.size()};
   DataSpace dataspace(1, str_dimsf);

   // Variable length string
   StrType datatype(PredType::C_S1, H5T_VARIABLE); 
   DataSet str_dataset = group.createDataSet(s_name, datatype, dataspace);

   str_dataset.write(arr_c_str.data(), datatype);
}

void write_matrix_HDF5(Group & group, const vector< vector<double> > & m, const string &  m_name){
	hsize_t dims[2]; 
	dims[0] = m.size();
	dims[1] = m[0].size();
	DataSpace fspace(2, dims); 
	DataSet m_dataset = group.createDataSet(m_name, PredType::NATIVE_DOUBLE, fspace);
	for (int i = 0; i < int(m.size()); i++){
		append_vector_to_matrix_HDF5(m_dataset, m[i],  i);
	}
}


void append_vector_to_matrix_HDF5(DataSet & dataset_tmp, const vector<double> & v, const int colNum){
   	hsize_t offset[2];
    offset[1] = 0;
    offset[0] = colNum;
	
    hsize_t fdims[2];            // new data dimensions 
	fdims[0] = 1;
	fdims[1] = v.size();
	DataSpace mspace( 2, fdims );

	DataSpace fspace = dataset_tmp.getSpace();
    fspace.selectHyperslab( H5S_SELECT_SET, fdims, offset );
    dataset_tmp.write(  v.data(), PredType::NATIVE_DOUBLE, mspace, fspace );	
}


void read_vector_HDF5(const H5File & file, const string & name, vector<int> & v_tmp){
	const H5std_string dataset_name( name );
	DataSet dataset = file.openDataSet( dataset_name );
	
	DataSpace dataspace = dataset.getSpace();
	hsize_t dims_out[1];
	dataspace.getSimpleExtentDims( dims_out, NULL);

	v_tmp.resize((unsigned long)(dims_out[0]));
	dataset.read( v_tmp.data(), PredType::NATIVE_INT, dataspace, dataspace );
	
	/*
	cout << "The length of the vector (int) is: " << (unsigned long)(dims_out[0]) << endl;
	cout << "The vector is: ";
	for (unsigned i = 0; i < v_tmp.size(); i++){
		cout << v_tmp[i] << ",";
	}
	cout << endl;
	*/
	
}

	
void read_vector_HDF5(const H5File & file, const string & name, vector<bool> & v_tmp){
	
	vector<int> v_tmp_int;
	read_vector_HDF5(file, name, v_tmp_int);
	
	v_tmp.resize(v_tmp_int.size());
	for (unsigned i = 0; i < v_tmp_int.size(); i++){
		v_tmp[i] = bool(v_tmp_int[i]);
	}
	
}

void read_vector_HDF5(const H5File & file, const string & name, vector<double> & v_tmp){
	const H5std_string dataset_name( name );
	DataSet dataset = file.openDataSet( dataset_name );
	
	DataSpace dataspace = dataset.getSpace();
	hsize_t dims_out[1];
	dataspace.getSimpleExtentDims( dims_out, NULL);

	v_tmp.resize((unsigned long)(dims_out[0]));
	dataset.read( v_tmp.data(), PredType::NATIVE_DOUBLE, dataspace, dataspace );
	
	/*
	cout << "The length of the vector (double) is: " << (unsigned long)(dims_out[0]) << endl;
	cout << "The vector is: ";
	for (unsigned i = 0; i < v_tmp.size(); i++){
		cout << v_tmp[i] << ",";
	}
	cout << endl;
	*/
}

void write_string_HDF5(Group & group, const string & s, const string &  s_name){
   // HDF5 only understands vector of char* :-(
   vector<const char*> arr_c_str;
   arr_c_str.push_back(s.c_str());

   hsize_t str_dimsf[1] {arr_c_str.size()};
   DataSpace dataspace(1, str_dimsf);

   // Variable length string
   StrType datatype(PredType::C_S1, H5T_VARIABLE); 
   DataSet str_dataset = group.createDataSet(s_name, datatype, dataspace);

   str_dataset.write(arr_c_str.data(), datatype);
}

bool group_exist_HDF5(const H5File & file, const string & name){
	// how to suppress the error message??
	bool r = true;
	try {  // to determine if the group exists in the file
		file.openGroup( name );
	}
	catch( FileIException& not_found_error ) {
	  r = false;
	}
	return r;
}

bool group_exist_HDF5(const string & filename, const string & name){
	hid_t file_id = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	//hbool_t check_object_valid = true;
	bool r = H5LTpath_valid(file_id, name.c_str(), true);
	
	H5Fclose(file_id);
	return r;
}

template < typename Type > Type read_scalar_HDF5(const H5File & file, const string & name){
	vector<Type> v_tmp;
	read_vector_HDF5(file, name, v_tmp);
	return v_tmp[0];
}

#endif