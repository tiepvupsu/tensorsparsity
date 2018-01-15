function [X, res] = general_sparse_coding(Y, D, Xinit, opts)
% function [X, res] = SHIRC_SC_l12_multi_2(Y, D, Xinit, opts)
% Solving SHIRC problem 
% `Y`: (d * n * T) - T representations of a signal 
% `D`: (d * k * T) - T dictionaries 
% `opts`:
%		opts.eps: tolerance 
%		opts.L: number of nonzero columns of X 
%       opts.regul: regularization function:
%           - 'l1'
%           - 'l2'
%           - 'l2_tensor'
% output: 
%		`X`: (k*C) and the corresponding residual
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 6/6/2016 2:14:08 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    addpath('utils/');
	if nargin == 0
        close all;
		d = 400;
		k = 20;
		T = 4;
		n = 1;
		Y = normc(rand(d, n, T));
		D = zeros(d, k, T);
		for c = 1: T 
			D(:, :, c) = normc(rand(d, k));
        end 
        Xinit = zeros(size(D, 2), size(Y, 2), size(Y, 3));
		opts.eps = .000001;
		opts.L = 10;
        opts.p = 1.5;
        opts.lambda = 0.5;
        opts.verbal = 1;
        opts.pos = false;
        opts.regul = 'tube';
        opts.range = [0, 5, 10, 15, 20];
    end 

    % D = normc(D);
    % Y = normc(Y);

    %%
    [d, k, T] = size(D);
    n = size(Y, 2);
    if numel(Xinit) == 0 
        Xinit = zeros(size(D, 2), n, size(Y, 3));
    end
    %% 
    if ~isfield(opts, 'p')
    	opts.p = 2;
    end 
    if ~isfield(opts, 'eps')
    	opts.eps = 0.0001;
    end
    if ~isfield(opts, 'weighted')
        opts.weighted = false;
    end 
    if ~isfield(opts, 'W')
        opts.W = ones(size(Xinit));
    end 
    %% ========================== l1 (Cumulative Residual) =============================
    if strcmp(opts.regul, 'l1')
        D = normc(D);
        Y = normc(Y);
        X = zeros(size(Xinit));
        for t = 1: T 
            X(:, :, t) = lasso_fista(Y(:, :, t), D(:, :, t), Xinit(:, :, t), opts);
        end 
        % X
        if nargin == 0
            imagesc(reshape(X, size(X, 1), size(X, 2)*size(X, 3)));
            X = [];
        end 
        return;
    end 
    %% ---------------------- end of l1 (Cumulative residual) --------------------------
    %% ========================== SRC-CC concatenation =============================
    if strcmp(opts.regul, 'concat')
        % Yconcat = normc(concat_tensor_channels(Y));
        % Dconcat = normc(concat_tensor_channels(D));
        Yconcat = Y; % concatenation was executed in SRC_local inside general_SRC_w_roc.m
        Dconcat = D; 
        X1 = full(lasso_fista(Yconcat, Dconcat, [], opts));
        % X = zeros(size(X1, 1), size(X1, 2), size(Y, 3));
        X = repmat(X1, [1, 1, T]);
        return
    end
    % ---------------------- end of SRC-CC concatenation --------------------------
    %%
    lambda = opts.lambda;
    % if numel(lambda) > 1 && size(lambda, 2) == 1 
    %     lambda = repmat(lambda, 1, size(Xinit, 2));
    % end 
    %%
    % function cost = calc_cost(X)
    %     cost = 0.5*norm(vec(Y - tensor_mult(D,X)))^2;
    %     norm2X = sum(vec(sqrt(sum(X.^2, 3))));
    %     cost = cost + lambda*norm2X;
    % end 

    %% 
    function cost = calc_f(X)
    	cost =  0.5*norm(vec(Y - tensor_mult(D,X)))^2; 
    end
    %% 
    function cost = calc_F(X) 
        calcfx = calc_f(X);
        if opts.weighted 
            X = opts.W.*X;
        end 
        switch opts.regul
            case 'l2' % column-sparsity
                cost = calcfx + lambda*sum(vec(sqrt(sum(X.^2, 1))));
            case 'l12' % row-sparsity
                cost = calcfx + lambda*sum(vec(sqrt(sum(X.^2, 2))));
            case 'tube' % tube-sparsity
                cost = calcfx + lambda*sum(vec(sqrt(sum(X.^2, 3))));
            case 'l2_tensor' % vertical slice sparsity 
                cost = calcfx + lambda*sum(vec(sqrt(sum(sum(X.^2, 3), 1))));
            case 'group'
                cost_group = 0;
                n_groups = numel(opts.range) - 1;
                for g = 1: n_groups
                    Xg = get_block_row(X, g, opts.range);
                    cost_group = cost_group + sum(vec(sqrt(sum(Xg.^2, 1))));
                end 
                cost = calcfx + lambda*cost_group;
            case 'group_tensor'
                cost_group = 0;
                n_groups = numel(opts.range) - 1;
                for g = 1: n_groups
                    Xg = get_block_row(X, g, opts.range);
                    cost_group = cost_group + sum(sqrt(sum(sum(Xg.^2, 3), 1)));
                end 
                cost = calcfx + lambda*cost_group;
        end
    end 
    %% gradient 
    D = normc(D);
    Y = normc(Y);
    Dt = permute(D, [2, 1, 3]);
    DtD = tensor_mult(Dt, D);
    DtY = tensor_mult(Dt, Y);
    function res = grad(X) 
    	res = tensor_mult(DtD, X) - DtY;
    end 
    %%
	max_eig = zeros(1, T);
	for i = 1: T 
		max_eig(i) = max(eig(DtD(:, :, i)));
	end 
	L = max(max_eig);
	[X, ~] = general_fista(@grad, Xinit, L, opts, @calc_F);
    % end 

	if nargin == 0
        switch opts.regul
            case {'l1', 'l2', 'group', 'l12'}
                imagesc(reshape(X, size(X, 1), size(X, 2)*size(X, 3)));
            % case 'l12'
            case {'tube', 'l2_tensor', 'group_tensor'}
                X = permute(X, [1, 3, 2]);
                imagesc(reshape(X, size(X, 1), size(X, 2)*size(X, 3)));

        end
		X = [];
	end
end 
