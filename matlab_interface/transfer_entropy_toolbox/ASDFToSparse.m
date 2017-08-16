% [raster, binunit] = ASDFToSparse(asdf)
%
%    asdf - {n_neu + 2, 1} ASDF data.
%
% Returns :
%    raster - (n_neu, duration) (Sparse) time raster expressed as a sparse matrix.
%    binunit - (string) the unit of time in the data. (length of a bin in real scale)
%
% Description :
%    This function converts the time raster of ASDF to sparse matrix.

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

function [raster, binunit] = ASDFToSparse(asdf)
info = asdf{end};
n_neu = info(1);
duration = info(2);
binunit = asdf{end - 1};


% very simple check of validity of ASDF
if n_neu ~= size(asdf,1) - 2
	error('Invalid n_neu information is contained in this ASDF');
	return
end


total_firings = 0;
for neu = 1:n_neu
	asdf{neu} = ceil(asdf{neu}); % I need to ceil them since an index can't be a float number
	total_firings = total_firings + length(asdf{neu});
end

neuron_data = zeros(total_firings,1);
time_data = zeros(total_firings,1);
spikes = ones(total_firings,1);

current = 1;

for neu = 1:n_neu
	leng = length(asdf{neu});
	if leng ~= 0
		range = current:current+leng-1;
		neuron_data(range) = neu;
		time_data(range) = asdf{neu};
		current = current + leng;
	end
end

raster = sparse(neuron_data, time_data, spikes, n_neu, duration); 
