% This is a demonstrative script of how to use TEpackage
% This program should work if you have everything set up correctly.
%
% Author: Shinya Ito, Indiana University
% Last Updated: 5/5/2011

% Make sure you have...
% * transent.c compiled.
% * all the programs in your Matlab path. (or current directory)


%% 1. Prepare some spike trains
% I'm makin in a way that train2 tends to fire 1 bin after train1.
% Also, train1 tends to fires 2 bins after train2.
disp('Here''re sample trains')
train1 = [0 1 0 0 1 1 0 1 0 0 0 0 1 0 0 1]
train2 = [0 0 1 0 0 0 1 1 0 0 1 0 0 1 0 0]


%% 2. Combine them together to make matrix representation of spike trains
trains = [train1 ; train2];


%% 3. Make proper data format called ASDF (Another Spiking Data Format)
asdf = SparseToASDF(trains, 1); % the last 1 means to use 1ms bin (doesn't matter in this case)

disp('You get TE out of these');
[peakTE, CI, TEdelays] = ASDFTE(asdf, 1:3) % using all the pairs in asdf, with delay of 1 to 3 bins
% peakTE is a peak value over delay. CI is coincidence index.
% TEdelays is 3-d matrix; each dimension means (sendind train, receiving train, delay), respectively.
% You should see TEdelays(1,2,1) and TEdelays(2,1,2) strong as I made them correlated.

%% 4. Some Higher Order TE calculations
disp('You can also calculate Higher Order TE with arbitrary delays');
TE133 = ASDFTE(asdf, 1, 3, 3) % Higher order TE with delay 1, order k=l=3
TE333 = ASDFTE(asdf, 3, 3, 3) % Higher order TE with delay 3, order k=l=3
TE1to312 = ASDFTE(asdf, 1, 1, 2) % Higher order TE with delay 1 to 3, order k=1, l=2


%% 5. Let's do the same thing with large spike train.
% One of the simulations from the paper. 100 subsampled neurons with 1.8 million bins. (30 minutes with 1ms bins)
load Izhik_100_0

disp('Calculating TE with 30 different delays for the 30 minute data... It will take a few minutes.')
tic
[peakTE, CI, TEdelays_big] = ASDFTE(asdf, 1:30); % Now it has delay of 1ms to 30ms
toc
% Using one core of my Core 2 duo 2.66GHz, it takes about 135 seconds.


%% 6. Show scatter plots
% removing self connections
peakTE = peakTE - diag(diag(peakTE));
CI = CI - diag(diag(CI));
% plot
from0to10 = find(conmat>0 & conmat<10);
subplot(1,2,1); semilogy(conmat(from0to10), peakTE(from0to10),'r.'); title 'weight vs TEPk'
subplot(1,2,2); plot(conmat(from0to10), CI(from0to10),'r.'); title 'weight vs TECI'


%% 7. Check if the values are consistent with my values.
cond = 0; % no error
if TEdelays(1,1,1) ~= 0
	cond = 1;
elseif abs(TEdelays(1,2,1) - 0.1953) > 0.001
	cond = 1;
elseif abs(TE333(1,2) - 0.5009) > 0.001
	cond = 1;
end

if cond == 1
	error('Some of your values are not correct. Please check settings!');
end
