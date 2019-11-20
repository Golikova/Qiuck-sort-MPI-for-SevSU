#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>
#include <string.h>

using namespace std;

void swap(int* a, int*b){
	int t = *a;
	*a = *b;
	*b = t;
}

int partition(int arr[], int low, int high){

	int pivot = arr[high];
	int i = (low-1);
	int j=0;
	for(j=low;j<=high-1;j++){
		if(arr[j]<=pivot){
			i++;
			swap(&arr[i],&arr[j]);
		}
	}
	swap(&arr[i+1], &arr[high]);
	return(i+1);
}

void quickSort(int arr[], int low, int high){
	if(low<high){
		int pi = partition(arr,low,high);
		quickSort(arr,low,pi-1);
		quickSort(arr,pi+1,high);
	} 
}

void print(int * array, int size) {
	for (int i = 0; i < size; ++i)
	{
		cout << array[i] << " " ;
	}
	cout << endl;
}


int main(int argc, char* argv[]){


	int rank, size;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	int globalSize=16;
	int* globalArray ;

	globalArray = new int[globalSize];

	for (int i = 0; i < globalSize; i++) {
			globalArray[i] = (rand() % 100);
	}


	if (rank == 0) {
		cout << "Создан массив: ";
		print(globalArray,globalSize);
	}


	int elementsPerProcess,pivot,localPivot,breakPoint,flag=0;

	// Initialize local array and compute the elementsPerProcess
	pivot = ceil(globalSize/2);
	pivot = globalArray[pivot];
	elementsPerProcess = globalSize / size ;

	int localArray[elementsPerProcess],i=0,j=0;
	int SampleArr[size];
	int FullSampleArr[size*size];


	/* ------------------------------- BUILD LOCAL ARRAYS ------------------------------------------------------*/

	MPI_Scatter(globalArray,elementsPerProcess,MPI_INT,localArray,elementsPerProcess,MPI_INT,0,MPI_COMM_WORLD);

	cout << "Процесс №" << rank << " работает с фрагментом " ;
	print(localArray, elementsPerProcess);

	// Sort the individual arrays

	quickSort(localArray,0,elementsPerProcess-1);

	MPI_Barrier(MPI_COMM_WORLD);

	/* -------------------------- BUILD SAMPLE SET ---------------------------------------------------------------*/

	// Collect Sample set and sort it
	for(i=0;i<size;i++)
		SampleArr[i]=localArray[i*globalSize/(size*size)];


	/*------------------------- AGGREGATE SAMPLE SET AND COLLECT INTO A SINGLE ARRAY ------------------------------*/
	MPI_Gather(SampleArr,size,MPI_INT,FullSampleArr,size,MPI_INT,0,MPI_COMM_WORLD);
	quickSort(FullSampleArr,0,size*size-1);

	MPI_Bcast(FullSampleArr,size*size,MPI_INT,0,MPI_COMM_WORLD);

	localPivot = FullSampleArr[(int)ceil(sizeof(FullSampleArr)/sizeof(FullSampleArr[0]))/2];


	/*---------------------------- BREAK POINT AND SUBARRAY GENERATION -------------------------------------------*/

	// Calculate the break point of each local array
	for(breakPoint=0;breakPoint<elementsPerProcess;breakPoint++){ 
		if(localPivot<=localArray[breakPoint]){
			if(localPivot==localArray[breakPoint])
				flag=1;
			break;
		}
	}


	int leftSub[elementsPerProcess];
	int rightSub[elementsPerProcess];


	// SET VALUES FOR THE LEFT SUBARRAY
	for(i=0;i<breakPoint;i++)
		leftSub[i]=localArray[i];	// Set the values upto breakpoint for the initial values

	for(i=breakPoint ; i<elementsPerProcess ; i++)
		leftSub[i]=-99999;	// Set remaining values to -99999


	// SET VALUES FOR THE RIGHT SUBARRAY

	if(flag==0){
		for(i=breakPoint;i<elementsPerProcess;i++)
			rightSub[i-breakPoint]=localArray[i];
		for(i=elementsPerProcess-breakPoint;i<elementsPerProcess;i++)
			rightSub[i]=-99999;
	}

	else{
		breakPoint++; // Skip the broadcasted pivot
		for(i=breakPoint ; i<elementsPerProcess ; i++){
			rightSub[i-breakPoint]=localArray[i];
		}
		for(i=elementsPerProcess-breakPoint;i<breakPoint+1;i++)
			rightSub[i]=-99999;
	}


	MPI_Barrier(MPI_COMM_WORLD);


	int leftSubSize = sizeof(leftSub)/sizeof(int);
	int rightSubSize = sizeof(rightSub)/sizeof(int);


	MPI_Barrier(MPI_COMM_WORLD);

	cout << endl << "Левая часть массива (ранг "<< rank << "): ";
	for(i=0;i<leftSubSize;i++) {
		if (leftSub[i]!= -99999)
		printf("%d ",leftSub[i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	rightSubSize = sizeof(rightSub)/sizeof(int);

	cout << endl << "Правая часть массива (ранг "<< rank << "): ";
	for(i=0;i<rightSubSize;i++) {
		if (rightSub[i]!= -99999)
		printf("%d ",rightSub[i]);
	}
	cout << endl;

/*------------------------------ ALL to ALL SECTION ---------------------*/


	int Sorted1[globalSize];
	int Sorted2[globalSize];
	int Final[globalSize];

	int counter1=0;

	MPI_Gather(leftSub,elementsPerProcess,MPI_INT,Sorted1,elementsPerProcess,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(rightSub,elementsPerProcess,MPI_INT,Sorted2,elementsPerProcess,MPI_INT,0,MPI_COMM_WORLD);

	quickSort(Sorted1,0,globalSize-1);	
	quickSort(Sorted2,0,globalSize-1);


	if(rank ==0){
		for(i=0;i<globalSize;i++){
			if(Sorted1[i]!=-99999){
				Final[counter1]=Sorted1[i];	
				counter1++;
			}
		}
		Final[counter1]=localPivot;
		counter1++;
		for(j=0;j<globalSize;j++){
			if(Sorted2[j]!=-99999){
				Final[counter1]=Sorted2[j];
				counter1++;
			}
		}
	}


	if(rank==0){
		print(Final, globalSize);
	}



	MPI_Finalize();
} 