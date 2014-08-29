#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define bi_spike_diffs_future_t_out plhs[0]
    
    #define num_pairs_in prhs[0]
    #define run_time_lengths_ruc_in prhs[1]
    #define num_trains_in prhs[2]
    #define folli_in prhs[3]
    #define run_time_range_in prhs[4]
    #define foll_indy_in prhs[5]
    #define udists_in prhs[6]
    
    int *num_pairs, *run_time_lengths_ruc, *num_trains, *foll_indy, *run_time_range, pac = 0, sac, trac1, trac2, M;
    double *bi_spike_diffs_future_t, *folli, *udists, *udists2, fnormy;
    const mxArray *udistsPr, *udists2Pr;
    
    num_pairs = mxGetPr(num_pairs_in);
    run_time_lengths_ruc = mxGetPr(run_time_lengths_ruc_in);
    num_trains = mxGetPr(num_trains_in);
    folli = mxGetPr(folli_in);
    run_time_range = mxGetPr(run_time_range_in);
    foll_indy = mxGetPr(foll_indy_in);
    
    bi_spike_diffs_future_t_out = mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(bi_spike_diffs_future_t_out, *num_pairs);
    mxSetN(bi_spike_diffs_future_t_out, *run_time_lengths_ruc);
    mxSetData(bi_spike_diffs_future_t_out, mxMalloc(sizeof(double) * *num_pairs * *run_time_lengths_ruc));
    bi_spike_diffs_future_t = mxGetPr(bi_spike_diffs_future_t_out);
    
    M  = mxGetM(udists_in);
    
    for(trac1 = 0; trac1 < *num_trains-1; ++trac1)
        for(trac2 = trac1 + 1; trac2 < *num_trains;  ++trac2) {
            pac++;
            
            udistsPr = mxGetCell(udists_in, trac2 * M + trac1);
            udists2Pr = mxGetCell(udists_in, trac1 * M + trac2);
            
            udists = mxGetPr(udistsPr);
            udists2 = mxGetPr(udists2Pr);
            
            for(sac = 0; sac < *run_time_lengths_ruc; ++sac) {
                fnormy = folli[trac1 + *num_trains * sac] + folli[trac2 + *num_trains * sac];
                if (folli[trac1 + *num_trains * sac ] < folli[trac2 + *num_trains * sac])
                    bi_spike_diffs_future_t[(pac-1) + *num_pairs * sac] = (folli[trac2 + *num_trains * sac] - folli[trac1 + *num_trains * sac]
                    + udists2[foll_indy[trac2 + *num_trains * sac] - 1])/(2*fnormy + !fnormy);
                else
                    bi_spike_diffs_future_t[(pac-1) + *num_pairs * sac] = (folli[trac1 + *num_trains * sac] - folli[trac2 + *num_trains * sac]
                    + udists[foll_indy[trac1 + *num_trains * sac] - 1])/(2*fnormy + !fnormy);
            }
        }
    return;
}