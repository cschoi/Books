#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "somap.h"
#include "randomNumber.h"

static double **allocInput(int ninputs, int depth)
{
  int i, j;
  double **input;

  input = (double **) malloc(ninputs * sizeof(double **));
  if (!input)
    return NULL;

  for (i = 0; i < ninputs; i++)
  {
    input[i] = (double *) malloc(depth * sizeof(double *));
    if (!input[i])
      return NULL;

    for (j = 0; j < depth; j++)
      input[i][j] = randomNumber();
  }

  return input;
}

static double **allocSpread(int width, int height)
{
  int i;
  double **spread;

  spread = (double **) malloc(width * sizeof(double *));
  if (!spread)
    return NULL;

  for (i = 0; i < width; i++)
  {
    spread[i] = (double *) malloc(height * sizeof(double));
    if (!spread[i])
      return NULL;
  }

  return spread;
}

int main(int argc, char **argv)
{
  //int nrows = 100, ncols = 100, depth = 3;
  int nrows = 50, ncols = 50, depth = 3;
  int width = 5, height = 5;
  int nsteps = 10, ninputs = nrows*ncols;
  double **input = allocInput(ninputs, depth);
  double **spread = allocSpread(width, height);
  double ***somap;
  time_t t0, t1;

  if (!input || !spread)
  {
    printf("Memory allocation problem\n");
    exit(1);
  }

  spread[0][0] = spread[0][4] = spread[4][0] = spread[4][4] = 0.0;
  spread[0][1] = spread[0][3] = spread[1][0] = spread[1][4] = 0.1;
  spread[3][0] = spread[3][4] = spread[4][1] = spread[4][3] = 0.1;
  spread[0][2] = spread[2][0] = spread[2][4] = spread[4][2] = 0.2;
  spread[1][1] = spread[1][3] = spread[3][1] = spread[3][3] = 0.35;
  spread[1][2] = spread[2][1] = spread[2][3] = spread[3][2] = 0.5;
  spread[2][2] = 1.0;

  t0 = time(NULL);

  somap = selfOrganisingMap(input, ninputs, depth, nrows, ncols, spread, width, height, nsteps);

  t1 = time(NULL);

  printf("Time taken = %ld\n", t1-t0);
}

