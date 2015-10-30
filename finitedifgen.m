%% 

% finite difference approximation generator
prompt = 'ENTER THE ORDER OF DERIVATIVE YOU WANT TO APPROXIMATE:';
m = input(prompt);
prompt ='ENTER p:';
p = input(prompt);
prompt = 'ENTER q:'; % the derivative  will be approximated by a linear combination of u_i-p, u_i-p+1, ... u_i+q
q = input(prompt);
if ((p+q)<m)
    error('IMPOSSIBLE TO COMPUTE. WRONG INPUT FOR P AND Q.');
end

%%

% initialize matrix for taylor table
% taylor tables looks as follows:
%           u_i     h*u_i(1) ...    h^k*u_i(k) ... h^(p+q)*u_i(p+q) ... h^(p+q+2)*u_i(p+q+2) (2 extra columns for errors) 
% 
% ap*u_i-p  1       (-p)^1/1!       (-p)^k/k!  ...     
% ...       ...      ...           ...
% an*u_i+n  ...      ...            (n)^k/k!               
% ...       ...      ...
% a0*u_i    ...      ...                           
% ...       ...      ...
% aq*u_i+q  1        ... 


taylor = zeros(p+q+1,p+q+3);

b = zeros(p+q+1,1);
b(m+1,1) = 1;

for n = (-p):1:q
    for k = 0:1:p+q+2 % indices are messy, see comment above for reference
        taylor(n+p+1,k+1)=((n)^(k))/(factorial(k));
    end
end
%%

taylor =(taylor)';
solve = taylor(1:1:p+q+1,:); % extract submatrix without error terms unnecessary for solving
x = solve\b;
error = taylor*x;
i=1;
display(x','COEFFICIENT VECTOR:');
while(round(m+1+i)~=round(p+q+4))
    if (abs(error(m+1+i))>10^(-9))
        display((i),'ORDER OF ERROR');
        break;
    end
    i=i+1;
end

    
%%
