#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "proto.h"

#define MAXBLOCKS 80000

static void **Table;

static size_t *BlockSize;

static int Nblocks = 0;

static int calls_of_free = 0;

static int calls_of_malloc = 0;


void hdf5_memory_cleanup(void)
{
  int i, count;

  size_t n;

  for(i = 0, n = 0, count = 0; i < Nblocks; i++)
    {
      free(Table[i]);
      n += BlockSize[i];
      count++;
    }

  if(n)
    {
      printf("\nFreed %d bytes on Task=%d in %d blocks that were not freed by HDF5.\n\n", (int) n, ThisTask,
	     count);
    }

  if(ThisTask == 0)
    {
      printf("Task=%d HDF5 calls_of_malloc = %d\n", ThisTask, calls_of_malloc);
      printf("Task=%d HDF5 calls_of_free   = %d\n", ThisTask, calls_of_free);
    }

  calls_of_malloc = 0;
  calls_of_free = 0;
  Nblocks = 0;
}


void hdf5_mymalloc_init(void)
{
  BlockSize = (size_t *) malloc(MAXBLOCKS * sizeof(size_t));
  Table = (void **) malloc(MAXBLOCKS * sizeof(void *));
  Nblocks = 0;
}


void hdf5_mymalloc_release(void)
{
  free(Table);
  free(BlockSize);
  Nblocks = 0;
}


void *hdf5_mymalloc(size_t n)
{
  calls_of_malloc++;

  if(Nblocks >= MAXBLOCKS)
    {
      printf("Task=%d: No blocks left in hdf5_mymalloc().\n", ThisTask);
      endrun(813);
    }
  Table[Nblocks] = malloc(n);

  if(n == 0)
    endrun(1231);

  if(Table[Nblocks])
    {
      BlockSize[Nblocks] = n;
      Nblocks++;
      return Table[Nblocks - 1];
    }

  return NULL;
}


void hdf5_myfree(void *p)
{
  int i;

  calls_of_free++;

  if(p)
    {
      for(i = 0; i < Nblocks; i++)
	if(Table[i] == p)
	  {
	    free(Table[i]);

	    Table[i] = Table[Nblocks - 1];
	    BlockSize[i] = BlockSize[Nblocks - 1];

	    Nblocks--;
	    return;
	  }

      printf("Free failed on Task=%d\n", ThisTask);
      endrun(3123);

      free(p);
    }
}


void *hdf5_myrealloc(void *p, size_t n)
{
  int i;

  if(p == NULL)
    return hdf5_mymalloc(n);

  if(n == 0)
    {
      hdf5_myfree(p);
      return NULL;
    }

  for(i = 0; i < Nblocks; i++)
    if(Table[i] == p)
      {
	Table[i] = realloc(p, n);
	BlockSize[i] = n;
	return Table[i];
      }
  printf("Realloc failed on Task=%d\n", ThisTask);
  endrun(3124);
  return NULL;
}


void *hdf5_mycalloc(size_t n, size_t m)
{
  void *p;

  char *c;

  size_t i, j;

  p = hdf5_mymalloc(n * m);

  c = (char *) p;

  for(i = 0; i < n; i++)
    for(j = 0; j < m; j++)
      *c++ = 0;

  return p;
}


char *hdf5_mystrdup(const char *s)
{
  char *p;

  p = (char *) hdf5_mymalloc(strlen(s) + 1);
  strcpy(p, s);

  return p;
}
