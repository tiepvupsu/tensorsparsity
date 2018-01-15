function Mi = get_block_col(M, C, col_range)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 2:39:32 AM
	% Last modified	: 1/27/2016 2:39:36 AM
	% Description	: Get column blocks of a big matrix M, blocks indexed by C 
	% 	INPUT
	%		M        : the big matrix. M = [M1 , M2 , ... , MC]
	%		i        : block index 
	%		col_range: a vector store the last index of each block. col_range(1) = 0.
	%					i-th block is indexed by col_range(i)+1: col_range(i+1).
	% 	OUTPUT 
	%		Mi: output block matrix  
	%
	%% ================== end File info ==========================
	Mi = [];
	for i = 1: numel(C)
		c = C(i);
		range = col_range(c) + 1: col_range(C+1);
        tmp =  size(size(M));
        if tmp(2) == 3
            Mi  = [Mi, M(:, range, :)];
        else
            Mi = [Mi, M(:, range)];
        end
	end 
end