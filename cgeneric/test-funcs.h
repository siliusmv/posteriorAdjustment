
#pragma once 

#include "cgeneric.h"

// Given a directory dir, locate all the files in the directory, read the data from
// the files (who should all contain a matrix or a sparse matrix on a prespecified format),
// and save the data into an inla_cgeneric_data_tp object
inla_cgeneric_data_tp read_cgeneric_data_from_dir(char const * dir);

// This is a struct used for keeping track of a list of n filenames
typedef struct {
  char **files;
  int n;
} file_list;

// Locate all filenames in a directory, and return them as a file_list
file_list files_in_dir(char const * dirname);

// Read an inla_cgeneric_mat_tp from the file filename, and save it into the object res
void read_mat_from_file(char const * filename, inla_cgeneric_mat_tp * res);

// Read an inla_cgeneric_smat_tp from the file filename, and save it into the object res
void read_smat_from_file(char const * filename, inla_cgeneric_smat_tp * res);

// Create an inla_cgeneric_mat_tp that is a copy of an inla_cgeneric_smat_tp
inla_cgeneric_mat_tp smat_to_mat(inla_cgeneric_smat_tp const A);

// Print an inla_cgeneric_smat_tp
void print_smat(inla_cgeneric_smat_tp const A);

// Print an inla_cgeneric_mat_tp
void print_mat(inla_cgeneric_mat_tp const A);

// Print an inla_cgeneric_smat_tp by first converting it into an inla_cgeneric_mat_tp
void print_smat_as_mat(inla_cgeneric_smat_tp const A);
