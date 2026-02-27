#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int decomp1d(int n,int p,int myid,int *s,int *e);

int main(int argc,char**argv)
{
    MPI_Init(&argc,&argv);

    int myid,p;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&p);

    int n;
    double *vec=NULL;

    if(myid==0){
        FILE *f=fopen(argv[1],"r");
        fscanf(f,"%d",&n);
        vec=malloc(n*sizeof(double));
        for(int i=0;i<n;i++) fscanf(f,"%lf",&vec[i]);
        fclose(f);
    }

    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

    int s,e;
    decomp1d(n,p,myid,&s,&e);
    int local_n=e-s+1;

    double *local=malloc(local_n*sizeof(double));

    MPI_Scatter(vec,local_n,MPI_DOUBLE,local,local_n,MPI_DOUBLE,0,MPI_COMM_WORLD);

    double local_sum=0;
    for(int i=0;i<local_n;i++)
        local_sum+=local[i]*local[i];

    double global_sum;
    MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(myid==0)
        printf("Norm^2 = %lf\n",global_sum);

    MPI_Finalize();
    return 0;
}

