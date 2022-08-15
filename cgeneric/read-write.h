
#pragma once 

#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <dirent.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

typedef struct {
  char **files;
  int n;
} file_list;

file_list files_in_dir(const char *dirname) {
  file_list res;
  int i = 0, max_size = 5;
  res.files = Calloc(max_size, char *);
  DIR *dir = opendir(dirname);
  struct dirent *file;
  while ((file = readdir(dir)) != NULL) {
    if (strcmp(file->d_name, ".") && strcmp(file->d_name, "..")) {
      if (i == max_size) {
	max_size *= 2;
	res.files = realloc(res.files, max_size * sizeof(char *));
      }
      res.files[i] = Calloc(1024, char);
      strcpy(res.files[i], dirname);
      strcat(res.files[i], "/");
      strcat(res.files[i], file->d_name);
      i += 1;
    }
  }
  closedir(dir);
  res.n = i;
  return res;
}

void read_mat_from_file(const char *filename, inla_cgeneric_mat_tp *res) {
  FILE *file = fopen(filename, "r");
  char line[1024], *tmp;
  fgets(line, 1024, file);
  assert(!strcmp(line, "mat\n"));
  fgets(line, 1024, file);
  res->nrow = atoi(line);
  fgets(line, 1024, file);
  res->ncol = atoi(line);
  res->x = Calloc(res->nrow * res->ncol, double);
  for (int i = 0; i < res->nrow * res->ncol; i++) {
    fgets(line, 1024, file);
    res->x[i] = strtod(line, &tmp);
  }
  fclose(file);
}

void read_smat_from_file(const char *filename, inla_cgeneric_smat_tp *res) {
  FILE *file = fopen(filename, "r");
  char line[1024];
  char *token, *tmp;
  fgets(line, 1024, file);
  assert(!strcmp(line, "smat\n"));
  fgets(line, 1024, file);
  res->n = atoi(line);
  fgets(line, 1024, file);
  res->nrow = atoi(line);
  fgets(line, 1024, file);
  res->ncol = atoi(line);
  res->x = Calloc(res->n, double);
  res->i = Calloc(res->n, int);
  res->j = Calloc(res->n, int);
  for (int i = 0; i < res->n; i++) {
    fgets(line, 1024, file);
    token = strtok(line, ";");
    res->i[i] = atoi(token);
    token = strtok(NULL, ";");
    res->j[i] = atoi(token);
    token = strtok(NULL, ";");
    res->x[i] = strtod(token, &tmp);
  }
  fclose(file);
}

inla_cgeneric_data_tp read_cgeneric_data_from_dir(const char *dir) {
  file_list filenames = files_in_dir(dir);
  inla_cgeneric_data_tp res;
  res.n_mat = 0;
  res.n_smat = 0;
  res.n_ints = 0;
  res.n_doubles = 0;
  res.n_chars = 0;
  FILE *file;
  char line[1024];
  for (int i = 0; i < filenames.n; i++) {
    file = fopen(filenames.files[i], "r");
    fgets(line, 1024, file);
    fclose(file);
    if (!strcmp(line, "mat\n")) {
      res.n_mat += 1;
    } else if (!strcmp(line, "smat\n")) {
      res.n_smat += 1;
    }
  }
  res.mats = Calloc(res.n_mat, inla_cgeneric_mat_tp *);
  res.smats = Calloc(res.n_smat, inla_cgeneric_smat_tp *);
  for (int i = 0; i < res.n_mat; i++) {
    res.mats[i] = Calloc(1, inla_cgeneric_mat_tp);
    res.smats[i] = Calloc(1, inla_cgeneric_smat_tp);
  }
  int mat_count = 0, smat_count = 0;
  for (int i = 0; i < filenames.n; i++) {
    file = fopen(filenames.files[i], "r");
    fgets(line, 1024, file);
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
