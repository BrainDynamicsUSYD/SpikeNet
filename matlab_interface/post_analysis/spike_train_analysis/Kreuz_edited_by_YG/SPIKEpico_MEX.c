#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define bi_spike_diffs_out plhs[0]
    
    #define num_pairs_in prhs[0]
    #define run_isi_lengths_ruc_in prhs[1]
    #define num_trains_in prhs[2]
    #define foll_spikes_in prhs[3]
    #define isi_pos_in prhs[4]
    #define prev_spikes_indy_in prhs[5]
    #define run_isi_starts_ruc_in prhs[6]
    #define ints_in prhs[7]
    #define run_isi_range_in prhs[8]
    #define foll_spikes_indy_in prhs[9]
    #define prev_spikes_in prhs[10]
    #define udists_in prhs[11]
    
    int *num_pairs, *run_isi_lengths_ruc, *run_isi_starts_ruc, *num_trains, *prev_spikes_indy, *foll_spikes_indy, *run_isi_range, pac = 0, sac, trac1, trac2, M;
    double *bi_spike_diffs, *isi_pos, *foll_spikes, *prev_spikes, *ints, *udists, *udists2;
    const mxArray *udistsPr, *udists2Pr;
    
    num_pairs = mxGetPr(num_pairs_in);
    run_isi_lengths_ruc = mxGetPr(run_isi_lengths_ruc_in);
    num_trains = mxGetPr(num_trains_in);
    foll_spikes = mxGetPr(foll_spikes_in);
    prev_spikes_indy = mxGetPr(prev_spikes_indy_in);
    isi_pos = mxGetPr(isi_pos_in);
    run_isi_starts_ruc = mxGetPr(run_isi_starts_ruc_in);
    ints = mxGetPr(ints_in);
    run_isi_range = mxGetPr(run_isi_range_in);
    foll_spikes_indy = mxGetPr(foll_spikes_indy_in);
    prev_spikes = mxGetPr(prev_spikes_in);
        
    bi_spike_diffs_out = mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(bi_spike_diffs_out, *num_pairs);
    mxSetN(bi_spike_diffs_out, *run_isi_lengths_ruc);
    mxSetData(bi_spike_diffs_out, mxMalloc(sizeof(double) * *num_pairs * *run_isi_lengths_ruc));
    bi_spike_diffs = mxGetPr(bi_spike_diffs_out);
    
    M  = mxGetM(udists_in);
    for(trac1 = 0; trac1 < *num_trains-1; ++trac1)
        for(trac2 = trac1 + 1; trac2 < *num_trains;  ++trac2) {
            pac++;
            
            udistsPr = mxGetCell(udists_in, trac2 * M + trac1);
            udists2Pr = mxGetCell(udists_in, trac1 * M + trac2);
            
            udists = mxGetPr(udistsPr);
            udists2 = mxGetPr(udists2Pr);
            
            for(sac = 0; sac < *run_isi_lengths_ruc; ++sac)
                bi_spike_diffs[(pac-1) + *num_pairs * sac] =
                ((udists[prev_spikes_indy[trac1 + *num_trains * sac] - 1]*(foll_spikes[trac1 + *num_trains * sac] - isi_pos[sac])
                + udists[foll_spikes_indy[trac1 + *num_trains * sac] - 1]*(isi_pos[sac] - prev_spikes[trac1 + *num_trains*sac]))
                /ints[trac1 + *num_trains*sac]*ints[trac2 + *num_trains*sac]
                + (udists2[prev_spikes_indy[trac2 + *num_trains * sac] - 1]*(foll_spikes[trac2 + *num_trains*sac] - isi_pos[sac])
                +  udists2[foll_spikes_indy[trac2 + *num_trains * sac] - 1]*(isi_pos[sac] - prev_spikes[trac2 + *num_trains*sac]))
                /ints[trac2 + *num_trains*sac]*ints[trac1 + *num_trains*sac])
                /((ints[trac1 + *num_trains*sac] + ints[trac2 + *num_trains*sac])*(ints[trac1 + *num_trains*sac] + ints[trac2 + *num_trains*sac])/2);
        }
    return;
}