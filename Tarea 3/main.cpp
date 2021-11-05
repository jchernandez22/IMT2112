#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <math.h>
using namespace std;

double alpha(double i, double j, double h_x, double h_y){
    // Esta función calcula el valor alpha dados i y j
    double x=i*h_x;
    double y=j*h_y;
    return x*(x-1)*y*(y-1) + 1;
}

double calculateN(double i, double j, double h_x, double h_y, int firstIndex){
    // Esta función calcula el valor de N dados i y j
    return -alpha(i, j + firstIndex + 0.5, h_x, h_y)/(h_y*h_y);
}

double calculateS(double i, double j, double h_x, double h_y, int firstIndex){
    // Esta función calcula el valor de S dados i y j
    return -alpha(i, j + firstIndex - 0.5, h_x, h_y)/(h_y*h_y);
}
double calculateE(double i, double j, double h_x, double h_y, int firstIndex){
    // Esta función calcula el valor de E dados i y j
    return -alpha(i + 0.5, j + firstIndex, h_x, h_y)/(h_x*h_x);
}
double calculateW(double i, double j, double h_x, double h_y, int firstIndex){
    // Esta función calcula el valor de W dados i y j
    return -alpha(i - 0.5, j + firstIndex, h_x, h_y)/(h_x*h_x);
}

double calculateC(double i, double j, double h_x, double h_y, int firstIndex){
    // Esta función calcula el valor de C dados i y j
    return -calculateW(i,j,h_x,h_y,firstIndex) - calculateE(i,j,h_x,h_y,firstIndex) - calculateS(i,j,h_x,h_y,firstIndex) - calculateN(i,j,h_x,h_y,firstIndex) + 1;
}

void print_matrix(double** matrix, int n, int m) {
  printf("\n");
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
         printf("%f ", matrix[i][j]); 
      }
    printf("\n");
  }
}

int main(){
    // Obtenemos el nro de procesadores y rank del proceso
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Definimos los valores n_x n_y h_x y h_x
    int n_x = 10;
    int n_y = 7;
    
    // Usamos pasos equidistantes
    double h_x = (double)(1/(n_x-1));
    double h_y = (double)(1/(n_y-1));

    // Ahora definimos el número de filas x procesador
    // El dominio se divide horizontalmente así que usamos n_y
    int firstIndex, localRowSize, err;
    localRowSize = n_y/world_size;
    firstIndex = world_rank*localRowSize;
    
    // Las cargas restantes se agregan al último procesador
    if(world_rank==world_size-1){
        localRowSize += (n_y%world_size);
    }

    // Ahora, creamos los 5 arrays de 2 dimensiones correspondientes
    // a C_ij E_ij W_ij N_ij y S_ij
    
    // Creamos el array 2D: C
    double localC[n_x][localRowSize] = {0}; // Creamos un array de 0's

    for(int i=1; i<n_x-1; i++){
        for(int j=1; j<localRowSize-1; j++){
            localC[i][j] = calculateC(i, j, h_x, h_y, firstIndex);
        }
    }

    // Creamos el array 2D: N
    double localN[n_x][localRowSize] = {0}; // Creamos un array de 0's

    for(int i=1; i<n_x-1; i++){
        for(int j=1; j<localRowSize-1; j++){
            localN[i][j] = calculateN(i, j, h_x, h_y, firstIndex);
        }
    }

    // Creamos el array 2D: S
    double localS[n_x][localRowSize] = {0}; // Creamos un array de 0's

    for(int i=1; i<n_x-1; i++){
        for(int j=1; j<localRowSize-1; j++){
            localS[i][j] = calculateS(i, j, h_x, h_y, firstIndex);
        }
    }

    // Creamos el array 2D: E
    double localE[n_x][localRowSize] = {0}; // Creamos un array de 0's

    for(int i=1; i<n_x-1; i++){
        for(int j=1; j<localRowSize-1; j++){
            localE[i][j] = calculateE(i, j, h_x, h_y, firstIndex);
        }
    }

    // Creamos el array 2D: W
    double localW[n_x][localRowSize] = {0}; // Creamos un array de 0's

    for(int i=1; i<n_x-1; i++){
        for(int j=1; j<localRowSize-1; j++){
            localW[i][j] = calculateW(i, j, h_x, h_y, firstIndex);
        }
    }

    print_matrix(localC, n_x, localRowSize);

    MPI_Finalize();
}