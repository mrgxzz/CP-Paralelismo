#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int i, n, numprocs, myid;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x, mypi;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid == 0) {
        while (1) {
            printf("Enter the number of intervals: (0 quits) \n");
            scanf("%d", &n);

            for (int dest = 1; dest < numprocs; dest++) {
                MPI_Send(&n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
            }

            if (n == 0) break;

            h = 1.0 / (double)n;
            sum = 0.0;
            for (i = myid + 1; i <= n; i += numprocs) {
                x = h * ((double)i - 0.5);
                sum += 4.0 / (1.0 + x * x);
            }
            pi = h * sum;

            for (int source = 1; source < numprocs; source++) {
                MPI_Recv(&mypi, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pi += mypi;
            }

            printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    } else {
        while (1) {
            MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (n == 0) break;

            h = 1.0 / (double)n;
            sum = 0.0;
            for (i = myid + 1; i <= n; i += numprocs) {
                x = h * ((double)i - 0.5);
                sum += 4.0 / (1.0 + x * x);
            }
            mypi = h * sum;

            MPI_Send(&mypi, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}
