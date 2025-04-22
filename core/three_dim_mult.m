function result = three_dim_mult(A, B, Adim, Bdim)
% Get the size of the matrices A and B
[xA, yA, zA] = size(A);
[xB, yB, zB] = size(B);

% Check for valid dimensions for multiplication
if Adim == Bdim
    error('Cannot multiply the same dimension from both matrices');
end

% Initialize the result matrix
% The size of the result depends on the other dimensions of A and B
if Adim == 2 && Bdim == 1
    if yA == xB
        result = zeros(xA,yB,zA*zB);
        for j=1:
        for i=1:yA
          result()
        end
        result = A * reshape(B, [xB, yB * zB]);
    else
        error('Dimensions of A and B do not match for multiplication along Adim=2 and Bdim=1');
    end
end

