#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <dirent.h>

#include "cgeneric.h"
#include "test-funcs.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))

inla_cgeneric_data_tp read_cgeneric_data_from_dir(const char *dir) {

  // Preallocate the result
  inla_cgeneric_data_tp res;
  res.n_mat = 0;
  res.n_smat = 0;
  res.n_ints = 0;
  res.n_doubles = 0;
  res.n_chars = 0;

  // Find all the files in the directory
  file_list filenames = files_in_dir(dir);

  // Open all files and read the first line, which states whether the file contains
  // an inla_cgeneric_mat_tp or an inla_cgeneric_smat_tp
  FILE * file;
  char line[1024];
  for (int i = 0; i < filenames.n; ++i) {
    file = fopen(filenames.files[i], "r");
    assert(fgets(line, 1024, file) != NULL);
    fclose(file);
    if (!strcmp(line, "mat\n")) {
      res.n_mat += 1;
    } else if (!strcmp(line, "smat\n")) {
      res.n_smat += 1;
    }
  }

  // Allocate the correct number of mats/smats in res
  res.mats = Calloc(res.n_mat, inla_cgeneric_mat_tp *);
  res.smats = Calloc(res.n_smat, inla_cgeneric_smat_tp *);
  for (int i = 0; i < res.n_mat; ++i) {
    res.mats[i] = Calloc(1, inla_cgeneric_mat_tp);
    res.smats[i] = Calloc(1, inla_cgeneric_smat_tp);
  }

  // Read all the mats/smats from the files and into the res object
  int mat_count = 0, smat_count = 0;
  for (int i = 0; i < filenames.n; ++i) {
    file = fopen(filenames.files[i], "r");
    assert(fgets(line, 1024, file) != NULL);
    fclose(file);
    if (!strcmp(line, "mat\n")) {
      read_mat_from_file(filenames.files[i], res.mats[mat_count]);
      mat_count += 1;
    } else if (!strcmp(line, "smat\n")) {
      read_smat_from_file(filenames.files[i], res.smats[smat_count]);
      smat_count += 1;
    }
  }

  return res;
}

file_list files_in_dir(const char * dirname) {

  // Preallocate the result
  file_list res;
  int i = 0, max_size = 5;
  res.files = Calloc(max_size, char *);

  // Open the directory and loop through all available files
  DIR *dir = opendir(dirname);
  struct dirent *file;
  while ((file = readdir(dir)) != NULL) {
    // Ignore the "." and ".." files
    if (strcmp(file->d_name, ".") && strcmp(file->d_name, "..")) {
      // Possibly reallocate res.files, if there is no more available space
      if (i == max_size) {
	max_size *= 2;
	res.files = realloc(res.files, max_size * sizeof(char *));
      }
      // Copy the filename into res.files[i]
      res.files[i] = Calloc(1024, char);
      strcpy(res.files[i], dirname);
      strcat(res.files[i], "/");
      strcat(res.files[i], file->d_name);
      i += 1;
    }
  }
  closedir(dir);
  res.n = i;

  // Reallocate res.files to free unused memory
  res.files = realloc(res.files, res.n * sizeof(char *));

  return res;
}

void read_mat_from_file(char const * filename, inla_cgeneric_mat_tp * res) {
  FILE *file = fopen(filename, "r");
  char line[1024], *tmp;
  assert(fgets(line, 1024, file) != NULL);
  assert(!strcmp(line, "mat\n"));
  assert(fgets(line, 1024, file) != NULL);
  res->nrow = atoi(line);
  assert(fgets(line, 1024, file) != NULL);
  res->ncol = atoi(line);
  res->x = Calloc(res->nrow * res->ncol, double);
  for (int i = 0; i < res->nrow * res->ncol; i++) {
    assert(fgets(line, 1024, file) != NULL);
    res->x[i] = strtod(line, &tmp);
  }
  fclose(file);
}

void read_smat_from_file(char const * filename, inla_cgeneric_smat_tp * res) {
  FILE *file = fopen(filename, "r");
  char line[1024];
  char *token, *tmp;
  assert(fgets(line, 1024, file) != NULL);
  assert(!strcmp(line, "smat\n"));
  assert(fgets(line, 1024, file) != NULL);
  res->n = atoi(line);
  assert(fgets(line, 1024, file) != NULL);
  res->nrow = atoi(line);
  assert(fgets(line, 1024, file) != NULL);
  res->ncol = atoi(line);
  res->x = Calloc(res->n, double);
  res->i = Calloc(res->n, int);
  res->j = Calloc(res->n, int);
  for (int i = 0; i < res->n; ++i) {
    assert(fgets(line, 1024, file) != NULL);
    token = strtok(line, ";");
    res->i[i] = atoi(token);
    token = strtok(NULL, ";");
    res->j[i] = atoi(token);
    token = strtok(NULL, ";");
    res->x[i] = strtod(token, &tmp);
  }
  fclose(file);
}

inla_cgeneric_mat_tp smat_to_mat(inla_cgeneric_smat_tp const A) {
  inla_cgeneric_mat_tp res;
  res.nrow = A.nrow;
  res.ncol = A.ncol;
  res.name = A.name;
  res.x = Calloc(res.nrow * res.ncol, double);
  for (int i = 0; i < A.n; i++) {
    res.x[A.j[i] + A.i[i] * res.ncol] = A.x[i];
  }
  return res;
}

void print_smat(inla_cgeneric_smat_tp const A) {
  printf("%s: (n: %d, nrow: %d, ncol: %d):\n", A.name, A.n, A.nrow, A.ncol);
  printf("i j x\n");
  for (int i = 0; i < A.n; i++) {
    printf("%d %d %f\n", A.i[i], A.j[i], A.x[i]);
  }
  printf("\n");
}

void print_mat(inla_cgeneric_mat_tp const A) {
  printf("%s: (nrow: %d, ncol: %d):\n", A.name, A.nrow, A.ncol);
  for (int i = 0; i < A.nrow; i++) {
    for (int j = 0; j < A.ncol; j++) {
      printf("%f ", A.x[j + i * A.ncol]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_smat_as_mat(inla_cgeneric_smat_tp const A) {
  inla_cgeneric_mat_tp tmp = smat_to_mat(A);
  print_mat(tmp);
}
