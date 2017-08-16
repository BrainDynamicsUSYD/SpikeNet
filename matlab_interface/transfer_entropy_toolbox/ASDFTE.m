% [te_result, ci_result, all_te] = ASDFTE(asdf, j_delay, i_order, j_order, windowsize)
% Parameters:
%   asdf        - Time series in Another Spike Data Format (ASDF)
%   j_delay     - Number of bins to lag sender (j series) or a vector [default 1]
%   i_order     - Order of receiver [default 1]
%   j_order     - Order of sender [default 1]
%   windowsize  - window size used for Coincidence Index calculation (odd number only)
%
% Returns:
%   te_result - (nNeu, nNeu) NxN matrix where N(i, j) is the transfer entropy from i->j
%               If multiple j_delay is given, this is a peak value of TE over delay.
%   ci_result - (nNeu, nNeu) NxN matrix where N(i, j) is the Coincidence Index from i->j
%               Multiple delays are necessary to calculate it.
%   all_te    - (nNeu, nNeu, delays) NxNxd matrix where N(i, j, k) is the transfer entropy
%               from i->j at delay of jdelay(k). For those who need all the delays.
%
% Examples:
% >> te = ASDFTE(asdf); % gives you delay 1 TE.
% >> [maxte, ci] = ASDFTE(asdf, 1:30); % gives you TEPk and TECI based on delay of 1 to 30.
% >> [maxte, ci] = ASDFTE(asdf, 1:30, 2, 3); % gives you HOTEPk and HOTECI (k=2, l=3) based on delay of 1 to 30.
% >> [maxte, ci, te] = ASDFTE(asdf, 1:30, 1, 1, 3); % gives you TEPk and TECI with CI window of 3, also gives you TE at all delays.

%==============================================================================
% Copyright (c) 2011, The Trustees of Indiana University
% All rights reserved.
%
% Authors: Michael Hansen (mihansen@indiana.edu), Shinya Ito (itos@indiana.edu)
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%   1. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%   3. Neither the name of Indiana University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================

function [te_result, ci_result, all_te] = ASDFTE(asdf, j_delay, i_order, j_order, windowsize)
% Set defaults
if nargin < 2
	j_delay = 1;
end

if nargin < 3
	i_order = 1;
end

if nargin < 4
	j_order = 1;
end

if nargin < 5
	windowsize = 5;
end

if length(j_delay) == 1
	te_result = transent(asdf, j_delay, i_order, j_order); % Single delay
	return;
else
	% Multiple delays
	num_delays = length(j_delay);
	info = asdf{end};
	num_neurons = info(1);

	% Allocate space for all matrices
	all_te = zeros(num_neurons, num_neurons, num_delays);

	% Compute TE for delay times
	for d = 1:num_delays % Change this for to parfor for parallelization.
		all_te(:, :, d) = transent(asdf, j_delay(d), i_order, j_order);
	end

	% Reduce to final matrix
	te_result = max(all_te, [], 3); % reduction in 3rd dimension

	ci_result = zeros(num_neurons);
	if nargout > 1
		for i = 1:num_neurons
			for j = 1:num_neurons
				ci_result(i, j) = CIReduce(all_te(i, j, :), windowsize);
			end
		end
	end
end % if multiple delays
