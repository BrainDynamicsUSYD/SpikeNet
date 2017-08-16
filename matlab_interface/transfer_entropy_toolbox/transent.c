/*=============================================================================
Copyright (c) 2011, The Trustees of Indiana University
All rights reserved.

% Authors: Michael Hansen (mihansen@indiana.edu), Shinya Ito (itos@indiana.edu)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

  3. Neither the name of Indiana University nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=============================================================================*/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <memory.h>

typedef int TimeType;

void transent_1
(const mxArray *all_series, const mwSize series_count,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result);

void transent_ho
(const mxArray *all_series, const mwSize series_count,
 const unsigned int x_order, const unsigned int y_order,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

  mwSize series_count;
  unsigned int x_order, y_order;
  TimeType y_delay, duration;
  double *data_ptr;
  mxArray *array_ptr;

	if (nlhs > 1) {
		mexErrMsgTxt("Expected only one output argument");
	}

  /* Arguments: all_series, y_delay */
  if (nrhs < 2) {
		mexErrMsgTxt("Expected at least two input arguments");
  }

  if (!mxIsCell(prhs[0])) {
		mexErrMsgTxt("First argument must be a cell array");
  }

  /* Extract arguments */
  series_count = mxGetNumberOfElements(prhs[0]) - 2; /* Last two cells hold info */

  data_ptr = mxGetPr(prhs[1]);
  y_delay = (TimeType)data_ptr[0];

  array_ptr = mxGetCell(prhs[0], series_count + 1);
  data_ptr = mxGetPr(array_ptr);
  duration = (TimeType)data_ptr[1];

  if (nrhs == 4) {
    data_ptr = mxGetPr(prhs[2]);
    x_order = (unsigned int)data_ptr[0];

    data_ptr = mxGetPr(prhs[3]);
    y_order = (unsigned int)data_ptr[0];
  }
  else {
    x_order = 1;
    y_order = 1;
  }

  /* Create result matrix */
	plhs[0] = mxCreateDoubleMatrix(series_count, series_count, mxREAL);

  /* Do calculation */
  data_ptr = mxGetPr(plhs[0]);

  if ((x_order == 1) && (y_order == 1)) {
    transent_1(prhs[0], series_count,
               y_delay, duration,
               data_ptr);
  } else {
    transent_ho(prhs[0], series_count,
                x_order, y_order,
                y_delay, duration,
                data_ptr);
  }
}

double Log2(double n) {
	return log(n) / log(2.0);
}

/* Computes the first-order transfer entropy matrix for all pairs. */
void transent_1
(const mxArray *all_series, const mwSize series_count,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result) {

  /* Constants */
  const unsigned int x_order = 1, y_order = 1,                
               num_series = 3,
               num_counts = 8,
               num_x = 4,
               num_y = 2;

  /* Locals */
  TimeType counts[8];
  unsigned long code;
  long k, l, idx, c1, c2;
  double te_final, prob_1, prob_2, prob_3;

  double *ord_iter[3];
  double *ord_end[3];

  TimeType ord_times[3];
  TimeType ord_shift[3];

  const unsigned int window = y_order + y_delay;
  const TimeType end_time = duration - window + 1;
  TimeType cur_time, next_time;

  /* Calculate TE */
  mxArray *array_ptr;
  double *i_series, *j_series;
  mwSize i_size, j_size;
  mwIndex i, j;

  /* MATLAB is column major */
  for (j = 0; j < series_count; ++j) {
    for (i = 0; i < series_count; ++i) {

      /* Extract series */
      array_ptr = mxGetCell(all_series, i);
      i_size = mxGetNumberOfElements(array_ptr);
      i_series = mxGetPr(array_ptr);

      array_ptr = mxGetCell(all_series, j);
      j_size = mxGetNumberOfElements(array_ptr);
      j_series = mxGetPr(array_ptr);

      if ((i_size == 0) || (j_size == 0)) {
        te_result[(i * series_count) + j] = 0;
		continue;
      }

      /* Order is x^(k+1), y^(l) */
      idx = 0;

      /* x^(k+1) */
      for (k = 0; k < (x_order + 1); ++k) {
        ord_iter[idx] = i_series;
        ord_end[idx] = i_series + i_size;
        ord_shift[idx] = (window - 1) - k;

        while ((TimeType)*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* y^(l) */
      for (k = 0; k < y_order; ++k) {
        ord_iter[idx] = j_series;
        ord_end[idx] = j_series + j_size;
        ord_shift[idx] = -k;
        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* Count spikes */
      memset(counts, 0, sizeof(TimeType) * num_counts);

      /* Get minimum next time bin */
      cur_time = ord_times[0];
      for (k = 1; k < num_series; ++k) {
        if (ord_times[k] < cur_time) {
          cur_time = ord_times[k];
        }
      }

      while (cur_time <= end_time) {

        code = 0;
        next_time = end_time + 1;

        /* Calculate hash code for this time bin */
        for (k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {      
            code |= 1 << k;

            /* Next spike for this neuron */
            ++(ord_iter[k]);

            if (ord_iter[k] == ord_end[k]) {
              ord_times[k] = end_time + 1;
            }
            else {
              ord_times[k] = (TimeType)*(ord_iter[k]) - ord_shift[k];
            }
          }

          /* Find minimum next time bin */
          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code]);
        cur_time = next_time;

      } /* while spikes left */

      /* Fill in zero count */
      counts[0] = end_time;
      for (k = 1; k < num_counts; ++k) {
        counts[0] -= counts[k];
      }

      /* ===================================================================== */

      /* Use counts to calculate TE */
      te_final = 0;

      /* Order is x^(k), y^(l), x(n+1) */
      for (k = 0; k < num_counts; ++k) {
        prob_1 = (double)counts[k] / (double)end_time;

        if (prob_1 == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        c1 = 0;
        c2 = 0;

        for (l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += (prob_1 * Log2(prob_2 / prob_3));
      }

      /* MATLAB is column major, but flipped for compatibility */
      te_result[(i * series_count) + j] = te_final;

    } /* for i */

  } /* for j */
 
} /* transent_1 */

/* Computes the higher-order transfer entropy matrix for all pairs. */
void transent_ho
(const mxArray *all_series, const mwSize series_count,
 const unsigned int x_order, const unsigned int y_order,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result) {

  /* Constants */
  const unsigned int num_series = 1 + y_order + x_order,
               num_counts = (unsigned int)pow(2, num_series),
               num_x = (unsigned int)pow(2, x_order + 1),
               num_y = (unsigned int)pow(2, y_order);

  /* Locals */
  TimeType *counts = (TimeType*)malloc(sizeof(TimeType) * num_counts);
  unsigned long code;
  long k, l, idx, c1, c2;
  double te_final, prob_1, prob_2, prob_3;

  double **ord_iter = (double**)malloc(sizeof(double*) * num_series);
  double **ord_end = (double**)malloc(sizeof(double*) * num_series);

  TimeType *ord_times = malloc(sizeof(TimeType) * num_series);
  TimeType *ord_shift = malloc(sizeof(TimeType) * num_series);

  const unsigned int window = (y_order + y_delay) > (x_order + 1) ? (y_order + y_delay) : (x_order + 1);
  const TimeType end_time = duration - window + 1;
  TimeType cur_time, next_time;

  /* Calculate TE */
  mxArray *array_ptr;
  double *i_series, *j_series;
  mwSize i_size, j_size;
  mwIndex i, j;

  /* MATLAB is column major */
  for (j = 0; j < series_count; ++j) {
    for (i = 0; i < series_count; ++i) {

      /* Extract series */
      array_ptr = mxGetCell(all_series, i);
      i_size = mxGetNumberOfElements(array_ptr);
      i_series = mxGetPr(array_ptr);

      array_ptr = mxGetCell(all_series, j);
      j_size = mxGetNumberOfElements(array_ptr);
      j_series = mxGetPr(array_ptr);

      if ((i_size == 0) || (j_size == 0)) {
        te_result[(i * series_count) + j] = 0;
		continue;
      }

      /* Order is x^(k+1), y^(l) */
      idx = 0;

      /* x^(k+1) */
      for (k = 0; k < (x_order + 1); ++k) {
        ord_iter[idx] = i_series;
        ord_end[idx] = i_series + i_size;
        ord_shift[idx] = (window - 1) - k;

        while ((TimeType)*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* y^(l) */
      for (k = 0; k < y_order; ++k) {
        ord_iter[idx] = j_series;
        ord_end[idx] = j_series + j_size;
        ord_shift[idx] = (window - 1) - y_delay - k;

        while ((TimeType)*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* Count spikes */
      memset(counts, 0, sizeof(TimeType) * num_counts);

      /* Get minimum next time bin */
      cur_time = ord_times[0];
      for (k = 1; k < num_series; ++k) {
        if (ord_times[k] < cur_time) {
          cur_time = ord_times[k];
        }
      }

      while (cur_time <= end_time) {

        code = 0;
        next_time = end_time + 1;

        /* Calculate hash code for this time bin */
        for (k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {       
            code |= 1 << k;

            /* Next spike for this neuron */
            ++(ord_iter[k]);

            if (ord_iter[k] == ord_end[k]) {
              ord_times[k] = end_time + 1;
            }
            else {
              ord_times[k] = (TimeType)*(ord_iter[k]) - ord_shift[k];
            }
          }

          /* Find minimum next time bin */
          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code]);
        cur_time = next_time;

      } /* while spikes left */

      /* Fill in zero count */
      counts[0] = end_time;
      for (k = 1; k < num_counts; ++k) {
        counts[0] -= counts[k];
      }

      /* ===================================================================== */

      /* Use counts to calculate TE */
      te_final = 0;

      /* Order is x^(k), y^(l), x(n+1) */
      for (k = 0; k < num_counts; ++k) {
        prob_1 = (double)counts[k] / (double)end_time;

        if (prob_1 == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        c1 = 0;
        c2 = 0;

        for (l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += (prob_1 * Log2(prob_2 / prob_3));
      }

      /* MATLAB is column major, but flipped for compatibility */
      te_result[(i * series_count) + j] = te_final;

    } /* for i */

  } /* for j */

  /* Clean up */
  free(counts);
  free(ord_iter);
  free(ord_end);
  free(ord_times);
  free(ord_shift);
 
} /* transent_ho */

