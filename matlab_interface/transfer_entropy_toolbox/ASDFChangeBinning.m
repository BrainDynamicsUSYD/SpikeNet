% new_asdf = ASDFChangeBinning(asdf, factor)
%    asdf - {n_neu,n_bin} the old asdf binned with different bin size
%    factor - (1,1) the number of bins that you want to consider as one bin
%                   in the new asdf. must be an integer.
%
% Returns:
%    new_asdf - {n_neu,n_bin_new} the new asdf binned with different bin size
%
% Description :
%    Change binning size. New asdf will have the size n_bin_new = floor(n_bin / factor)
%    Since it "bins" the spiking, the resulting timing becomes an integer (and unique).

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

function new_asdf = ASDFChangeBinning(asdf, factor)
n_neu = asdf{end}(1);
original_duration = asdf{end}(2);
new_duration = ceil(original_duration / factor);

new_asdf = cell(size(asdf));
for i = 1:n_neu
	new_asdf{i} = unique(ceil(asdf{i} / factor));
end

% Change the binning size
if isnumeric(asdf{end-1})
	new_asdf{end-1} = asdf{end-1} * factor;
else
	new_asdf{end-1} = 'void';
end
new_asdf{end} = [n_neu, new_duration];
