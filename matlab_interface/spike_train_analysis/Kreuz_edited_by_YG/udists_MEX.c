#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define udists_out plhs[0]
    
    #define num_trains_in prhs[0]
    #define num_uspikes_in prhs[1]
    #define uspikes_in prhs[2]
    
    int *num_trains, sac, trac1, trac2, *num_uspikes, *max_num_ispikes, indx1, usac, changeSign;
    double *uspikes, *udists1, *temp, normy;
    mxArray *udists, *udists1Pr, *tempPr;
    
    num_trains = mxGetPr(num_trains_in);
    num_uspikes = mxGetPr(num_uspikes_in);
    uspikes = mxGetPr(uspikes_in);
    
    udists_out = mxCreateCellMatrix(*num_trains, *num_trains);
    
    for(trac1 = 0; trac1 < *num_trains; ++trac1)
        for(trac2 = 0; trac2 < *num_trains;  ++trac2)
            if (trac1 != trac2) {
                
                udists1Pr = mxCreateDoubleMatrix(1, num_uspikes[trac1], mxREAL);
                tempPr =  mxCreateDoubleMatrix(1, num_uspikes[trac2], mxREAL);
                
                udists1 = mxGetPr(udists1Pr);
                temp =  mxGetPr(tempPr);
                
                for(sac = 0; sac < num_uspikes[trac2]; ++sac)
                    temp[sac] = uspikes[trac1] - uspikes[trac2 + *num_trains*sac];
                
                udists1[0] = temp[0]; indx1 = 0;
                if (udists1[0] < 0)
                    udists1[0] = -udists1[0];
                
                for(usac = 1; usac < num_uspikes[trac2]; ++usac) {
                    if (temp[usac] < 0)
                        changeSign = -1;
                    else
                        changeSign = 1;
                    if (udists1[0] > changeSign*temp[usac]) {
                        udists1[0] = changeSign*temp[usac];
                        indx1 = usac;
                    }
                    else
                        break;
                }
                
                for(sac = 1; sac < num_uspikes[trac1]; ++sac) {
                    for(usac = 0; usac < num_uspikes[trac2]; ++usac)
                        temp[usac] = temp[usac] + uspikes[trac1 + *num_trains*sac] - uspikes[trac1 + *num_trains*(sac - 1)];
                    
                    udists1[sac] = temp[indx1];
                    if (udists1[sac] < 0)
                        udists1[sac] = -udists1[sac];
                    
                    for(usac = indx1+1; usac < num_uspikes[trac2]; ++usac) {
                        if (temp[usac] < 0)
                            changeSign = -1;
                        else
                            changeSign = 1;
                        if (udists1[sac] > changeSign*temp[usac]) {
                            udists1[sac] = changeSign*temp[usac];
                            indx1 = usac;
                        }
                        else
                            break;
                    }
                }
                
                mxSetCell(udists_out, trac2 * *num_trains + trac1, udists1Pr);
                
            }
    return;
}