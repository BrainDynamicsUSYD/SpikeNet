Package for Transfer Entropy Calculation (Ver. 0.5)
Authors:
Michael Hansen <mihansen@indiana.edu> (Indiana University)
Shinya Ito <itos@indiana.edu> (Indiana University)
Date: June 11, 2013

How to use TEpackage
Index
1. Preparation
2. Another Spiking Data Format (ASDF)
3. Functions
4. About the order and delay of TE
5. Tested environment


1. Preparation
==============
Compile the mex file (C program) with gcc or lcc (Matlab default C compiler).
>> mex transent.c % if you use mex for the first time, you might be prompted to choose a compiler.
Put all the programs into your Matlab path (or current directory).
That's it!
For demonstration of the code, please download Izhik_100_0.mat from our
project web page, open and run 'simple_demo.m'.


2. Another Spiking Data Format (ASDF)
=====================================
Another Spike Data Format is basically cell array of spike timing of each neuron.
In order to calculate TE correctly, I recommend to use only integer for the timing.
To ensure that, you can do
>> asdf = ChangeBinning(asdf, 1);

Last two cells contains special information of the data.

asdf{end-1}: Binning size in unit of ms (milisecond). (e.g. 1.2 -> 1.2ms/bin, 10 -> 10ms/bin etc...)
asdf{end}: Array of number of neurons, number of total bins of the data
  (e.g. [10 300000] -> 10 neurons, 300000 time bins)


3. Functions (Please refer to help of each function for details)
================================================================
* TE Calculation
ASDFTE.m: (requires transentmex) delayed higher order TE calculator for ASDF.
ASDFTE_parallel.m: (requires transentmex, and Parallel Computing Toolbox)
Parallel version of ASDFTE.m. Call 'matlabpool(N)' before running.

* Changing Data Format
SparseToASDF.m: Convert matrix form of raster to ASDF.
ASDFToSparse.m: Convert ASDF to sparse matrix of size (nNeu, duration)

* ASDF utilities
ASDFSubsample.m: Subsample specified neurons from ASDF.
ASDFChooseTime.m: Crop a time segment from larger ASDF.
ASDFGetfrate.m: Get firing rate (per bin) of all the neurons.
ASDFChangeBinning.m: Change binning size of ASDF.

* Supporting functions (not to be excuted directly)
transent.c: Mex file for rapid calculation of TE.


4. About the order and delay of TE
==================================
Order is 'length of time you include in the past'.
If you use 3rd order, you include 3 bins before the time you want to predict.

* - the time bin you want to predict
o - the time bin you use for prediction
1st order
train i:   ....o*....
train j:   ....o.....

3rd order for both i and j
train i:   ..ooo*....
train j:   ..ooo.....

Delay is 'distance of time from what you want to predict to what you use for predict'.

delay=1 (1st order)
train i:   ....o*....
train j:   ....o.....

delay=3 (1st order)
train i:   ....o*....
train j:   ..o.......

You can also combine these two.

delay=3 (2nd order for j, 3rd order for i)
train i:   ..ooo*....
train j:   .oo.......


5. Tested environment
=====================
These programs are tested on these environments
Environment 1
Spec: Core 2 Duo E8200 (2.66GHz) + 4GB RAM
OS: Ubuntu Linux 11.04 64-bit
Compiler: GCC 4.4.3
Matlab: R2010b
Octave: 3.2.4

Environment 2
Spec: Core 2 Duo E8200 (2.66GHz) + 4GB RAM
OS: Windows XP 32-bit
Compiler: LCC
Matlab: R2006a

Environment 3
Spec: Mac Pro 5,1
OS: Mac OS X Mountain Lion
Compiler: GCC 4.2.1 (LLVM)
Matlab: R2011a
