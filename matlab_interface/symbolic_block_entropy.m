function [ H, word_count ] = symbolic_block_entropy(seq, k)
% k: finite alphabet size
% order: order of block entropy (length of word, currently only 2nd order
% is supported
%
% seq: symbolic sequence, can be integer vector, ideally from 1 to k
% ref: Entropy estimation of very short symbolic sequences, Phys Rev E, 79,
% 046208 (2009)
%
%

if isempty(seq) || length(seq) == 1
    H = 0;
    word_count = 0;
else
    
    order = 2;
    
    seq = seq(:); % column vector
    
    n = length(seq);
    words = zeros(n,order); % Pre Allocate Memory
    
    for i = 1:order;
        words(:,i) = circshift(seq,1-i);
    end
    words = words(1:end-order+1,:); % Truncate Last Part
    
    if min(seq) >= 1 && max(seq) <= k
        word_count = zeros(k,k);
        for i = 1:(n-order+1)
            word_count( words(i,1), words(i,2) ) = word_count( words(i,1), words(i,2) ) + 1;
        end
        word_freq = word_count/sum(sum(word_count));
        p = word_freq(word_freq > 0);
        H = -sum(p.*log2(p));
    else
        error('Need to reformat the symbolic sequence with alphabet from 1 to k')
    end
    
end

end

