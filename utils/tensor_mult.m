function C = tensor_mult(A, B)
	n_channels = size(A, 3); % = size(B, 3);
    if n_channels == 1 
        C = A*B;
    else 
    	C = zeros(size(A, 1), size(B, 2), n_channels);
    	for i = 1: n_channels 
    		C(:, :, i) = A(:, :, i) * B(:, :, i);
    	end
    end 
end 
	