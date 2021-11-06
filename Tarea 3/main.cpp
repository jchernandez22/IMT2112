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

void print_2DArray(double *array, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << *((array+i*n) + j) << " ";
        }
        cout << endl;
  }
}

int main(){
    MPI_Init(NULL,NULL);
    // Obtenemos el nro de procesadores y rank del proceso
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Obtenemos el nombre del procesador y el largo del nombre
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Definimos los valores n_x n_y h_x y h_x
    int n_x = 10;
    int n_y = 7;
    
    // Usamos pasos equidistantes
    double h_x = (double)1/(n_x-1);
    double h_y = (double)1/(n_y-1);

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
    double localC[n_x][localRowSize];

    for(int i=0; i<n_x; i++){
        for(int j=0; j<localRowSize; j++){
            // Calculamos el valor expresado en el stencil
            localC[i][j] = calculateC(i, j, h_x, h_y, firstIndex);

            // Si nos paramos en los bordes de los lados:
            if(i==0 || i==n_x-1){
                localC[i][j]=0;
            }

            // Si somos el primer proceso y estamos en la primera fila:
            if(world_rank==0 && j==0){
                localC[i][j]=0;
            }

            // Si somos el último proceso y estamos en la última fila:
            if(world_rank==world_size-1 && j==localRowSize-1){
                localC[i][j]=0;
            }
        }
    }

    // Creamos el array 2D: N
    double localN[n_x][localRowSize];

    for(int i=0; i<n_x; i++){
        for(int j=0; j<localRowSize; j++){
            localN[i][j] = calculateN(i, j, h_x, h_y, firstIndex);
        }
    }

    // Creamos el array 2D: S
    double localS[n_x][localRowSize];

    for(int i=0; i<n_x; i++){
        for(int j=0; j<localRowSize; j++){
            localS[i][j] = calculateS(i, j, h_x, h_y, firstIndex);
        }
    }

    // Creamos el array 2D: E
    double localE[n_x][localRowSize];

    for(int i=0; i<n_x; i++){
        for(int j=0; j<localRowSize; j++){
            localE[i][j] = calculateE(i, j, h_x, h_y, firstIndex);
        }
    }

    // Creamos el array 2D: W
    double localW[n_x][localRowSize]; 

    for(int i=0; i<n_x; i++){
        for(int j=0; j<localRowSize; j++){
            localW[i][j] = calculateW(i, j, h_x, h_y, firstIndex);
        }
    }

    cout << "Name:" << processor_name << " Rank:" << world_rank << " Array C: " << endl;
    print_2DArray((double *)localC, n_x, localRowSize);

    MPI_Finalize();
}