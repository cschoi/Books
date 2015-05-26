#include "somap.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "randomNumber.h"

static double ***randomMap(int nrows, int ncols, int depth)
{
  int i, j, k;
  double ***x;

  x = (double ***) malloc(nrows * sizeof(double **));
  if (!x)
    return NULL;

  for (i = 0; i < nrows; i++)
  {
    x[i] = (double **) malloc(ncols * sizeof(double *));
    if (!x[i])
      return NULL;

    for (j = 0; j < ncols; j++)
    {
      x[i][j] = (double *) malloc(depth * sizeof(double));
      if (!x[i][j])
        return NULL;

      for (k = 0; k < depth; k++)
        x[i][j][k] = randomNumber();
    }
  }

  return x;
}

static void updateSomap(double ***somap, double *input, int depth, int nrows, int ncols,
                        double **spread, int width, int height, double decay)
{
  int halfWidth, halfHeight, i, j, k, l, m, imin, jmin;
  double diff, diff2, diff2min, lambda;

  imin = jmin = -1;
  diff2min = 0; // will change in first pass
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      diff2 = 0.0;
      for (k = 0; k < depth; k++)
      {
        diff = somap[i][j][k] - input[k];
        diff2 += diff * diff;
      }

      if ((imin == -1) || (diff2 < diff2min))
      {
        imin = i;
        jmin = j;
        diff2min = diff2;
      }
    }
  }

  halfWidth = (width-1) / 2;
  halfHeight = (height-1) / 2;

  for (k = 0; k < width; k++)
  {
    i = (imin + k - halfWidth + nrows) % nrows;
    for (l = 0; l < height; l++)
    {
      j = (jmin + l - halfHeight + ncols) % ncols;
      lambda = decay * spread[k][l];
      for (m = 0; m < depth; m++)
        somap[i][j][m] = (1.0-lambda) * somap[i][j][m] + lambda * input[m];
    }
  }
}

double ***selfOrganisingMap(double **inputs, int ninputs, int depth,
	int nrows, int ncols, double **spread, int width, int height, int nsteps)
{
  int step, i;
  double decay;
  double ***somap;

/*
printf("%d %d %d %d %d %d %d\n", ninputs, depth, nrows, ncols, width, height, nsteps);
*/
  somap = randomMap(nrows, ncols, depth);
  if (!somap)
    return NULL;

  for (step = 0; step < nsteps; step++)
  {
    decay = exp(-step / (double) nsteps);
    for (i = 0; i < ninputs; i++)
      updateSomap(somap, inputs[i], depth, nrows, ncols, spread, width, height, decay);
  }

  return somap;
}

