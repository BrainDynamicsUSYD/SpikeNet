#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define bi_spike_diffs_t_out plhs[0]
    
    #define num_pairs_in prhs[0]
    #define run_time_lengths_ruc_in prhs[1]
    #define num_trains_in prhs[2]
    #define folli_in prhs[3]
    #define previ_in prhs[4]
    #define isi_indy_in prhs[5]
    #define prev_indy_in prhs[6]
    #define run_isi_starts_ruc_in prhs[7]
    #define ints_in prhs[8]
    #define run_isi_range_in prhs[9]
    #define udists_in prhs[10]
    
    int *num_pairs, *run_time_lengths_ruc, *run_isi_starts_ruc, *num_trains, *prev_indy, *isi_indy, *run_isi_range, pac = 0, sac, trac1, trac2, M;
    double *bi_spike_diffs_t, *folli, *previ, *ints, *udists, *udists2;
    const mxArray *udistsPr, *udists2Pr;
    
    num_pairs = mxGetPr(num_pairs_in);
    run_time_lengths_ruc = mxGetPr(run_time_lengths_ruc_in);
    num_trains = mxGetPr(num_trains_in);
    folli = mxGetPr(folli_in);
    previ = mxGetPr(previ_in);
    isi_indy = mxGetPr(isi_indy_in);
    prev_indy = mxGetPr(prev_indy_in);
    run_isi_starts_ruc = mxGetPr(run_isi_starts_ruc_in);
    ints = mxGetPr(ints_in);
    run_isi_range = mxGetPr(run_isi_range_in);
    
    bi_spike_diffs_t_out = mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(bi_spike_diffs_t_out, *num_pairs);
    mxSetN(bi_spike_diffs_t_out, *run_time_lengths_ruc);
    mxSetData(bi_spike_diffs_t_out, mxMalloc(sizeof(double) * *num_pairs * *run_time_lengths_ruc));
    bi_spike_diffs_t = mxGetPr(bi_spike_diffs_t_out);
    
    M  = mxGetM(udists_in);
    
    for(trac1 = 0; trac1 < *num_trains-1; ++trac1)
        for(trac2 = trac1 + 1; trac2 < *num_trains;  ++trac2) {
            pac++;
            
            udistsPr = mxGetCell(udists_in, trac2 * M + trac1);
            udists2Pr = mxGetCell(udists_in, trac1 * M + trac2);
            
            udists = mxGetPr(udistsPr);
            udists2 = mxGetPr(udists2Pr);
            
            for(sac = 0; sac < *run_time_lengths_ruc; ++sac)
                bi_spike_diffs_t[(pac-1) + *num_pairs * sac] = ((udists[prev_indy[trac1 + *num_trains * sac] - 1]*folli[trac1 + *num_trains * sac]
                + udists[prev_indy[trac1 + *num_trains*sac]]*previ[trac1 + *num_trains*sac])
                /ints[trac1 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)]*ints[trac2 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)]
                + (udists2[prev_indy[trac2 + *num_trains*sac] - 1]*folli[trac2 + *num_trains*sac]
                +  udists2[prev_indy[trac2 + *num_trains*sac]]*previ[trac2 + *num_trains*sac])
                /ints[trac2 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)]*ints[trac1 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)])
                /((ints[trac1 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)] + ints[trac2 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)])
                *(ints[trac1 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)] + ints[trac2 + *num_trains*(isi_indy[sac] - *run_isi_starts_ruc)])/2);
        }
    return;
}