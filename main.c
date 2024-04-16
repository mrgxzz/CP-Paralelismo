#include <stdio.h>
#include <math.h>
#include <mpi.h>

// Definición de MPI_FlattreeColectiva para la recolección de sum
void MPI_FlattreeColectiva(double *mypi, double *pi, int numprocs, int myid) {
    MPI_Status status;
    int mask = 1;
    while (mask < numprocs) {
        if ((myid & mask) == 0) {
            int partner = myid | mask;
            if (partner < numprocs) {
                double recv_pi;
                MPI_Recv(&recv_pi, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &status);
                *pi += recv_pi;
            }
        } else {
            int partner = myid & (~mask);
            MPI_Send(mypi, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
        }
        mask <<= 1;
    }
}

// Definición de MPI_BinomialColectiva para la distribución de n
void MPI_BinomialColectiva(int *n, int numprocs, int myid) {
    int mask = 1;
    while (mask < numprocs) {
        int partner = myid ^ mask;
        if (partner < numprocs) {
            if (myid & mask) {
                MPI_Recv(n, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Send(n, 1, MPI_INT, partner, 0, MPI_COMM_WORLD);
            }
        }
        mask <<= 1;
    }
}

int main(int argc, char *argv[]) {
    int i, n, numprocs, myid;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x, mypi;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    while (1) {
        if (myid == 0) {
            printf("Enter the number of intervals: (0 quits) \n");
            scanf("%d", &n);
        }
        MPI_BinomialColectiva(&n, numprocs, myid);

        if (n == 0) break;

        h = 1.0 / (double)n;
        sum = 0.0;
        for (i = myid + 1; i <= n; i += numprocs) {
            x = h * ((double)i - 0.5);
            sum += 4.0 / (1.0 + x * x);
        }
        mypi = h * sum;

        // Utilizar MPI_FlattreeColectiva para recolectar la estimación de PI de cada proceso
        pi = 0.0;
        MPI_FlattreeColectiva(&mypi, &pi, numprocs, myid);

        if (myid == 0) {
            printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize();
    return 0;
}
