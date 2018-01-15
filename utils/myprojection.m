function X = myprojection(U, opts)
% function X = myprojection(U, opts)
% Description: Solve: 
% if opts.regul == 'l1':
% 	xi = \arg\min 0.5*||xi - ui||_F^2 + opts.lambda ||xi||_1
% if opts.regul == 'l2':
% xi = \arg\min 0.5*||xi - ui||_F^2 + opts.lambda ||xi||_2
% where xi and ui are the i-th columns of X and U 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 6/8/2016 3:36:06 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0 
		d = 1000;
		n = 1000;
		U = rand(d, n);
		opts.lambda = .15;
		opts.pos = true;
	end 
	lambda = opts.lambda;
	[d, n, T] = size(U);
	%%
	% [~, n] = size(U);
	if ~isfield(opts, 'pos')
		opts.pos = false;
	end 

	if ~isfield(opts, 'regul')
		opts.regul = 'l1';
	end 
	%%
	if opts.lambda == 0 
		X = U;
	else 
		switch opts.regul
			%% ========================== l1 =============================
			case 'l1'
				if numel(opts.lambda) > 1 && size(opts.lambda, 2) == 1 % column vector 
					lambda = repmat(lambda, 1, size(U, 2));
				end 
				%%
				if opts.pos 
					X = max(0, U - lambda);
				else 
					X = max(0, U - lambda) + min(0, U + lambda);
				end
			%% ---------------------- end of l1 --------------------------
			%% ========================== l2 =============================
			case 'l2'
				norm2_cols = sqrt(sum(U.^2, 1));
				k = max(1 - opts.lambda./norm2_cols, 0);
				if opts.pos 
					X = repmat(k, [d, 1, 1]).*max(0, U);
				else 
					X = repmat(k, [d, 1, 1]).*U;
				end 
			%% ---------------------- end of l2 --------------------------
			%% ========================== l2-tensor =============================
			case 'l2_tensor'
				norm2_vertical_slice = sqrt(sum(sum(U.^2, 3), 1));
				k = max(1 - opts.lambda./norm2_vertical_slice, 0);
				if opts.pos 
					X = repmat(k, [d, 1, T]).*max(0, U);
				else 
					X = repmat(k, [d, 1, T]).*U;
				end
			%% ---------------------- end of l2-tensor --------------------------
			%% ========================== group =============================
			case 'group'
				U_range = opts.range;
				X = zeros(size(U));
				n_ranges = numel(U_range) - 1;
				for g = 1: n_ranges
					Ug = get_block_row(U, g, U_range);
					opts_tmp = opts;
					opts_tmp.regul = 'l2';
					X(U_range(g)+1: U_range(g+1), :, :) = myprojection(Ug, opts_tmp);
				end  
			%% ---------------------- end of group --------------------------
			%% ========================== group tensor =============================
			case 'group_tensor'
				U_range = opts.range;
				X = zeros(size(U));
				n_ranges = numel(U_range) - 1;
				for g = 1: n_ranges
					Ug = get_block_row(U, g, U_range);
					opts_tmp = opts;
					opts_tmp.regul = 'l2_tensor';
					X(U_range(g)+1: U_range(g+1), :, :) = myprojection(Ug, opts_tmp);
				end  
			%% ---------------------- end of group tensor --------------------------
			%% ========================== row sparsity =============================
			case 'l12'
				norm2_rows = sqrt(sum(U.^2, 2));
				k = max(1 - opts.lambda./norm2_rows, 0);
				if opts.pos 
					X = repmat(k, [1, n, 1]).*max(0, U);
				else 
					X = repmat(k, [1, n, 1]).*U;
				end 
			 
			%% ---------------------- end of row sparsity --------------------------
			%% ========================== row sparsity =============================
			case 'tube'
				norm2_tubs = sqrt(sum(U.^2, 3));
		        lambda = opts.lambda;
		        % if numel(lambda) > 1 && size(lambda, 2) == 1
		        %     lambda = repmat(lambda, 1, size(U, 2));
		        % end
				k = max(1 - lambda./norm2_tubs, 0);
				if opts.pos 
					X = repmat(k, [1, 1, T]).*max(0, U);
				else 
					X = repmat(k, [1, 1, T]).*U;
				end 
			%% ---------------------- end of row sparsity --------------------------
			otherwise
				X = U;
		end
    end 
    if nargin == 0 
        X = [];
    end 
end
	