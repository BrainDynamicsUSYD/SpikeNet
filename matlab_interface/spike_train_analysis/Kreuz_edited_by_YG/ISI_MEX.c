#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define bi_isi_ratio_out plhs[0]
    
    #define num_pairs_in prhs[0]
    #define run_time_lengths_ruc_in prhs[1]
    #define num_trains_in prhs[2]
    #define run_isi_starts_ruc_in prhs[3]
    #define ints_in prhs[4]
    
    int *num_pairs, *run_time_lengths_ruc, *num_trains,  pac = 0, sac, trac1, trac2;
    double *bi_isi_ratio, *ints;
    
    num_pairs = mxGetPr(num_pairs_in);
    run_time_lengths_ruc = mxGetPr(run_time_lengths_ruc_in);
    num_trains = mxGetPr(num_trains_in);
    ints = mxGetPr(ints_in);
    
    bi_isi_ratio_out = mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(bi_isi_ratio_out, *num_pairs);
    mxSetN(bi_isi_ratio_out, *run_time_lengths_ruc);
    mxSetData( bi_isi_ratio_out, mxMalloc(sizeof(double) * *num_pairs * *run_time_lengths_ruc));
    bi_isi_ratio = mxGetPr(bi_isi_ratio_out);
    
    for(trac1 = 0; trac1 < *num_trains-1; ++trac1)
        for(trac2 = trac1 + 1; trac2 < *num_trains;  ++trac2) {
            pac++;
            
            for(sac = 0; sac < *run_time_lengths_ruc; ++sac) {
                if (ints[trac1 + *num_trains * sac] < ints[trac2 + *num_trains * sac])
                    bi_isi_ratio[(pac-1) + *num_pairs * sac] = ints[trac1 + *num_trains * sac]/ints[trac2 + *num_trains * sac] - 1;
                else
                    if (ints[trac1 + *num_trains * sac] == 0)
                        bi_isi_ratio[(pac-1) + *num_pairs * sac] = 0;
                    else
                        bi_isi_ratio[(pac-1) + *num_pairs * sac] = 1 - ints[trac2 + *num_trains * sac]/ints[trac1 + *num_trains * sac];
            }
            
        }
    return;
}