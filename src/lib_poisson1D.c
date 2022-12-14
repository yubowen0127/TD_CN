/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
	
	for(int i = 0; i < *la; i++) {
		
		for(int j = 0; j < *(lab); j++) {
			
			if (j < *kv) {
				AB[i * *lab + j] = 0.0;
			}
			else if(j == *kv + 1) {
				AB[i * *lab + j] = 2.0;
			}
			else {
				B[i * *lab + j] = -1,0;
			}
		}
	}
	
	AB[*kv] = 0.0;
	AB[(*lab) * (*la) -1] = 0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
	
	for(int i = 0; i < *la; i++) {
		
		for(int j = 0; j < *lab; j++) {
			
			if(j == 1) {
				AB[i * *lab + j] = 1.0;
			}
			else {
				AB[i * *lab + j] = 0.0;
			}
		}
	}
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
	
	RHS[0] = *BC0;
	for(int i = 1; i < *la - 1; i++) {
		RHS[i] = 0;
	}
	RHS[*la - 1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
	
	for(int i = 0; i < *la; i++) {
		
		EX_SOL[i] = *BC0 + X[i] * (*BC1 - *BC0);
	}
}  

void set_grid_points_1D(double* x, int* la){
	
	const double h = 1.0/ (double) (*la + 1); 
	
	for(int i = 1; i < *la + 1; i++) {
		
		X[i-1] = i*h;
	}
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
	
	if(*la < 0) {
		*info = -1;}
	else if(*n < 0) {
		*info = -2;}
	else if(*kl != 1) {
		*info = -3;}
	else if(*ku != 1) {
		*info = -4;}
	else if(AB == NULL) {
		*info = -5;}
	else if(*lab < 3) {
		*info = -6;}
	else if(ipiv == NULL) {
		*info = -7;}
	else{
		AB[*lab - 1] /= AB[*lab-2];
		for(int i = 1; i < *n; i++) {
			AB[(i+1) * (*lab) - 2] += AB[i *     (*lab) - 1];
			if(AB[(i+1) * (*lab) -2] == 0) {
				*info = i;
				break;
			}
			AB[(i+1) * (*lab) - 1] /= AB[(i+1) * (*lab) - 2];
			if(i == *n -1) *info = 0;
		}
	}
	return *info;

}
