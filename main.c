
#include <stdio.h>
#include <math.h>
#include <mpi.h>

void MPI_FlattreeColectiva(void * buff, void * recvbuff, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
void MPI_BinomialColectiva(void* buff, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

int main(int argc, char *argv[]){
    int i, done = 0, n;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x;

    int rank, numprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    while (!done) {
        if (rank == 0) {
            printf("Enter the number of intervals: (0 quits) \n");
            scanf("%d", &n);
        }
        MPI_BinomialColectiva(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (n == 0) break;

        h = 1.0 / (double) n;
        sum = 0.0;
        for (i = rank + 1; i <= n; i += numprocs) {
            x = h * ((double) i - 0.5);
            sum += 4.0 / (1.0 + x * x);
        }
        double resultado;
        MPI_FlattreeColectiva(&sum, &resultado, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            pi = h * resultado;
            printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }
    MPI_Finalize();
}

void MPI_FlattreeColectiva(void * buff, void * recvbuff, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    int size, myrank, i;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &myrank);

    if(myrank != root){
        MPI_Send(buff, count, datatype, root, 0, comm);
    }
    else{
        double aux;
        *(double*)recvbuff = *(double*)buff;
        for (i = 1; i < size; i++) {
            MPI_Recv(&aux, count, datatype, i, 0, comm, MPI_STATUS_IGNORE);
            *(double*)recvbuff += aux;
        }
    }
}

void MPI_BinomialColectiva(void* buff, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int size, myrank, i;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &myrank);
    if(myrank!=root){
        MPI_Recv(buff, 1, datatype, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
    }

    int offset= pow(2,0);

    for(i=0;(myrank+offset) < size;i++) {
        offset = pow(2, i);
        if (myrank < offset) {
            int destino = myrank + offset;
            if (destino < size) {
                MPI_Send(buff, 1, datatype, destino, 0, comm);
            }
        }
    }
}