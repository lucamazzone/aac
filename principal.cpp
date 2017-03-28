#include <mpi.h>

using namespace std;

int main() {
    int numProcs, rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int* data;
    int num;

    data = new int[5];
    data[0] = 0;
    data[1] = 1;
    data[2] = 2;
    data[3] = 3;
    data[4] = 4;
//    MPI_Scatter(data, 5, MPI_INT, &num, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(data, 1, MPI_INT, &num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cout << rank << " recieved " << num << endl; 

    MPI_Finalize();
    return 0;
}
