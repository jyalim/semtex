#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define ROOTONLY if (ip == 0)

const int NB = 25600;		/* -- Code fails for NB = 16000! */


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Test block transposition across 2 processors.  Run: mpirun -np 2 a.out
 * ------------------------------------------------------------------------- */
{
  int         i, ip, np;
  double*     data;
  MPI_Status  status;
  MPI_Request send_req;
  MPI_Request recv_req;

  MPI_Init (&argc, &argv);
 
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);

  data = (double*) malloc (np*NB*sizeof (double));
  for (i = 0; i < np*NB; i++) data[i] = 0.0;

  MPI_Barrier (MPI_COMM_WORLD);
  ROOTONLY printf ("-- Initially:\n");
  MPI_Barrier (MPI_COMM_WORLD);

  for (i = 0; i < np; i++)
    if (i != ip) {
      MPI_Isend (data + i*NB, NB, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
      MPI_Irecv (data + i*NB, NB, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_req);
      MPI_Wait  (&send_req, &status);
      MPI_Wait  (&recv_req, &status);
    }

  MPI_Barrier (MPI_COMM_WORLD);
  ROOTONLY printf ("-- Finally:\n");
  MPI_Barrier (MPI_COMM_WORLD);


  MPI_Finalize();

  return EXIT_SUCCESS;
}
