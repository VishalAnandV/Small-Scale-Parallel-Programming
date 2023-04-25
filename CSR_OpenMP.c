#include <stdio.h>
#include <stdlib.h>
#include<sys/time.h>
#include<string.h>
#include<omp.h>


// Structure to represent a triplet
typedef struct 
{
    int row;
    int col;
    double val;
} Triplet;

// Structure to represent a CSR matrix
typedef struct 
{
    int *rowptr;
    int *colidx;
    double *val;
} CSRMatrix;

//Matrix information are initialised
int nrow,ncol,nnz;

//Structure declaration for matrix data in triplet and CSR formats
Triplet *triplets;
CSRMatrix csr;

/*
 * Function: compare
 * -----------------
 *   Compares row indices of triplet data
 *
 *   returns: 0 if row indices of triplets x and y match;
 *            1 if row index of triplet x greater than row index of y
 *           -1 if row index of triplet x lesser than row index of y
 */
int compare(const void *x, const void *y)
{
    Triplet *tripletsx=(Triplet*)x;
	Triplet *tripletsy=(Triplet*)y;
	
	if(tripletsx->row > tripletsy->row)
	{
        	return 1;
    	}
	else if(tripletsx->row < tripletsy->row) 
	{
        	return -1;
   	}
	else
	{
		return 0;
	}
}

/*
 * Function: read_matrix
 * ---------------------
 *   Reads the matrix market file and processes the type of matrix
 *   and stores the data in a triplet structure
 *
 *   filename: Name of matrix market file to be read
 *
 *   returns: None
 */
void read_matrix(char* filename)
{
    char line[512];
    FILE* file=fopen(filename,"r");

    // Read the first line to identify the type of matrix
    fgets(line,sizeof(line),file);

	char *token;
	const char delimiter[2]=" ";
	int word=0;

	//get the fourth and fifth word
	token = strtok(line,delimiter);
	
	int k=0;	
	int mtype=0; // if type of matrix is real -> 0, if pattern -> 1
	int mformat=0; // if general -> 0, if symmetric ->1
	while(token != NULL)
	{
        k++;
        if(k == 4) 
		{
			if(strncmp(token, "pattern", 7) == 0)
			{
			    mtype=1;
          printf("\nThe type of the Matrix Market is Pattern\n");
        	}
		} 
		else if(k == 5) 
		{
			if(strncmp(token, "symmetric", 9) == 0)
			{                   
                mformat=1;
                printf("\nThe type of the Matrix Market is Symmetric\n");
            }

        }
        token = strtok(NULL, delimiter);
    }
	
	//Read comments
    do{
        fgets(line,sizeof(line),file);
    }while(line[0]=='%');

    // Read basic matrix information
    sscanf(line,"%d %d %d",&nrow,&ncol,&nnz);
	
	//memory allocation for triplet storage
    if(mformat==1)
    {
        triplets=(Triplet*)malloc(nnz*2*sizeof(Triplet));
    }
    else
	{
        triplets=(Triplet*)malloc(nnz*sizeof(Triplet));
    }
	
    //Checks whether the type of matrix is pattern and reads the data accordingly	
	if(mtype==1)
	{
	    for(int i=0;i<nnz;i++)
        {
            fgets(line,sizeof(line),file);
            sscanf(line,"%d %d",&triplets[i].row,&triplets[i].col);
			triplets[i].val=1.0;
            triplets[i].row--;
            triplets[i].col--;
        }
	}
	else
	{
		for(int i=0;i<nnz;i++)
        {
            fgets(line,sizeof(line),file);
            sscanf(line,"%d %d %lf",&triplets[i].row,&triplets[i].col,&triplets[i].val);
            triplets[i].row--;
            triplets[i].col--;
        }
	}
	
	int index=nnz;
	printf("\nNo. of non-zero elements before reconstruction of the matrix: %d\n",nnz);
	
	//Checks whether the matrix format is symmetric to construct the full symmetric matrix
	if(mformat==1)
	{
		for(int i=0;i<nnz;i++)
		{
			if(triplets[i].row!=triplets[i].col)
			{
				triplets[index].row=triplets[i].col;
				triplets[index].col=triplets[i].row;
				triplets[index].val=triplets[i].val;
				index+=1;
			}
		}
	}

	nnz=index;
	printf("\nNo. of non-zero elements after reconstruction of the matrix: %d\n",nnz);
    qsort(triplets,nnz,sizeof(Triplet),compare);
	
	fclose(file);
}

/*
 * Function: convertToCSR
 * --------------------------
 *   Converts the triplet matrix structure to CSR format
 *
 *   returns: None
 *
 */
void convertToCSR()
{
    // Allocate memory for csr format matrix
    csr.rowptr = (int *)calloc((nrow + 1),sizeof(int));
    csr.colidx = (int *)malloc(nnz * sizeof(int));
    csr.val = (double *)malloc(nnz * sizeof(double));

    // Counts the number of non-zero elements in each row
    for (int i = 0; i < nnz; i++)
	{
        csr.rowptr[triplets[i].row]++;
    }

	// Compute the sum of the row pointer array
    int sum=0;
    for (int i = 0; i < nrow; i++)
    {
        int temp = csr.rowptr[i];
        csr.rowptr[i]=sum;
        sum+=temp;
    }
	csr.rowptr[nrow]=nnz;

    // Copy the column indices and values into the CSR matrix
    for (int i = 0; i < nnz; i++) 
    {
        int row = triplets[i].row;
        int idx = csr.rowptr[row]++;
        csr.colidx[idx] = triplets[i].col;
        csr.val[idx] = triplets[i].val;
    }
	//  The row pointer array is shifted to the right by one to set the first element to 0
    for (int i = nrow-1; i >= 0; i--) 
    {
        csr.rowptr[i+1] = csr.rowptr[i];
    }
	csr.rowptr[0] = 0;

}

/*
 * Function: serial_CSR
 * ------------------------
 *   Multiplies the csr matrix with a vector in serial
 *
 *   *x contains the vector data; *y stores the results obtained during serial multiplication
 *   returns: None
 *
 */
void serial_CSR(int *x,double *y)
{
	for(int i=0;i<nrow;i++)
    {
        double t=0.0;
        for(int j=csr.rowptr[i];j<csr.rowptr[i+1];j++)
        {
            t+=csr.val[j]*x[csr.colidx[j]];
        }
        y[i]=t;
    }
}


/*
 * Function: parallel_CSR
 * --------------------------
 *   Multiplies the csr matrix with a vector by parallelization
 *   
 *   *x contains the vector data; *z stores the results obtained during parallel multiplication
 *   
 *   returns: None
 *
 */
void parallel_CSR(int *x,double *z)
{
	int i,j;
        double t2;


        #pragma omp parallel for private(i,j) shared(csr,x,z)
        for (i = 0; i < nrow; i++)
        {
                t2=0.0;
                for (j = csr.rowptr[i]; j < csr.rowptr[i+1]; j++)
                {
                        t2+= csr.val[j] * x[csr.colidx[j]];
                }
                z[i]=t2;
        }
}


/*
 * Function: validate
 * ------------------
 *   Validation of parallel multiplication result vector by comparing it with serial
 *   multiplication result vector
 *
 *   *y contains the serial results; *z contains the parallel results
 * 
 *   returns: None
 *
 */
void validate(double *y,double *z)
{
    int valid=0; // Stores the number of matched values
    for(int i=0;i<nrow;i++)
    {
        if(y[i]==z[i])
        {
            valid++;
        }
    }
    printf("\nNumber of elements in the result vector: %d\n",nrow);
    printf("Number of matching values: %d\n",valid);
    printf("Value of the last element in the result vector in serial: %f  in parallel: %f\n",y[valid-1],z[valid-1]);
    if(valid==nrow)
    {
        printf("\nValidation Successful\n");
    }
    else
    {
        printf("\nValidation Failed\n");
    }
}

/*
 * Function: csr_matrix_multiply
 * ----------------------------------------
 *   Multiplies the csr matrix with a vector in both serial and parallel methods
 *
 *   returns: None
 *   
 */
void csr_matrix_multiply()
{
    // Declare vector to be multiplied with the ellpack matrix
    int *x=(int *)malloc(nrow*sizeof(int));
    double *y=(double *)malloc(nrow*sizeof(double));//Serial result vector
    double *z=(double *)malloc(nrow*sizeof(double));//Parallel result vector

	//Input data into x vector
	for(int i=0; i<nrow; i++)
    {
        x[i] = i;
    }

	struct timeval start_time, end_time;
  	double comp_time;

	gettimeofday(&start_time, NULL);

	// serial sparse matrix multiplication
	serial_CSR(x,y);	

	gettimeofday(&end_time, NULL);

	/*
	printf("\nSerial results:\n");
	for(int i=0;i<nrow;i++)
    {
        printf("\n%f\n",y[i]);
    }
    */

    //Computation time for performing serial matrix vector multiplication is computed
    comp_time = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec)/1e6;
    printf("\nComputation time for CSR serial implementation: %f seconds\n",comp_time);
    double flops = (2*nnz)/comp_time;
    printf("FLOPS for CSR serial implementation: %.2f\n",flops);

	double start_ptime = omp_get_wtime();
	
	// parallel sparse matrix multiplication
    parallel_CSR(x,z);
	
	double end_ptime = omp_get_wtime();

	/*
	printf("\nParallel results:\n");
	for(int i=0;i<nrow;i++)
    {
        printf("\n%f\n",z[i]);
    }
    */
	
	//Computation time for performing parallel matrix vector multiplication is computed
	double comp_ptime = end_ptime - start_ptime;
    printf("\nComputation time for CSR parallel implementation: %f seconds\n",comp_ptime);
    double pflops = (2*nnz)/comp_ptime;
    printf("FLOPS for CSR parallel implementation: %.2f\n",pflops);

	//Validation of results
	validate(y,z);

}

int main()
{
    read_matrix("cage4.mtx");
    convertToCSR();
    // Print the CSR matrix
        
	/*
	printf("rowptr:");
    for (int i = 0; i <= nrow; i++)
    {
        printf(" %d", csr.rowptr[i]);
    }
	
	printf("\ncolidx:");
    for (int i = 0; i < nnz; i++)
    {
        printf(" %d", csr.colidx[i]);
    }
	printf("\nval:\n");
    for (int i = 0; i < nnz; i++)
    {
        printf("\n %f\n", csr.val[i]);
    }
	printf("\n");
	*/

    csr_matrix_multiply();
    printf("\nCSR Successful\n\n");
    
	/*
	printf("Sorted matrix");
	for(int i=0;i<nnz;i++)
	{
		printf("\n%d  %d  %f\n",triplets[i].row,triplets[i].col,triplets[i].val);
	}
	*/
    
    // Free memory
    return 0;
}