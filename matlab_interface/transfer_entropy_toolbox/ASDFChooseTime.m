% cutasdf = ASDFchoosetime(asdf, startTime, endTime)
%
%    asdf - {nNeu+2,1} ASDF to subsample from
%    startTime - (scalar) Starting time of new asdf
%                  if only two arguments are given, this is the time that you wanna choose.
%    endTime - (scalar) Ending time of new asdf
%    
%
% Returns:
%    cutasdf - {nNeu+2,1} ASDF with new time and duration

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

function cutasdf = ASDFChoosetime(asdf, startTime, endTime)
nNeu = asdf{end}(1);
duration = asdf{end}(2);
if nargin > 2
	newDuration = endTime - startTime + 1; % edges are included
	cutasdf = cell(nNeu+2,1);

	for i = 1:nNeu
		% throw away the unnecessary tails
		nontails = find(asdf{i} <= endTime);
		newTrain = asdf{i}(nontails) - startTime + 1;

		% throw away the heads.
		positives = find(newTrain > 0);
		cutasdf{i} = newTrain(positives);
	end

else
	newDuration = length(startTime);
	cutasdf = cell(nNeu+2,1);
	for i = 1:nNeu
		inds = find(ismember(asdf{i}, startTime));
		cutasdf{i} = asdf{i}(inds);
	end
end

cutasdf{end-1} = asdf{end-1}; % binning size doesn't change
cutasdf{end} = [nNeu, newDuration]; % new [nNeu duration]
