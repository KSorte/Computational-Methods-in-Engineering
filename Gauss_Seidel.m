%% 
% Name    : KAUSHAL SORTE
% Roll no : BT18MEC044
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   GAUSS-SIEDEL ALGORITHM DEMONSTRATION   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NECESSARY DATA
format long
%A = [20 6 7; 1 9 7; 7 9 20];    % Example of convergent matrix 
%A = [5 6 7; 6 3 9; 7 9 10]     % Example of divergent matrix
%B = [18; 18; 26];
%A = [9 1 2 3; 1 7 2 3; 1 5 10 2; 1 2 3 9] 
%B = [1;2;3;4]
%Tolerance = 0.001;
A = input('Enter coefficient matrix');
B = input('Enter right hand side vector');
Tolerance = input('Enter tolerance for the first norm of unknown vector');


%% FUNCTION CALL
[x,Data, Norm_Values] = GaussSiedel(A,B,Tolerance);
% This will give 3 outputs.

%% POST-PROCESSING   
plot(Norm_Values,'-o')
ylabel('First norm values')
xlabel('Iterations')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAIN PROGRAM ENDS HERE. GuassSiedel FUNCTION WILL BE DEFINED BELOW.  
%% FUNCTION DEFINITION 
function [x,Data,Norm_Values] = GaussSiedel(A,B,Tolerance)
n = size(B);
x = zeros([n,1]);      %Guess vector default set as zero for ease of operation.
Norm_Values = [];      % Storing first norm values.
K = 25;                % Limit on number of iterations.    
Data = zeros([n,K]);   % Storing output of all iterations for analysis.
%% FORMING EQUATIONS
X = sym('X',[1 n]);
digits(6)      % To limit the number of significant digits of symbolic operations to 6.
sym sigma;
for i = 1:n
sigma =0;
  for j = 1:n
      sigma  = vpa(sigma + X(j)*A(i,j));
  end
  eqns(i) = sigma == B(i);       %Equation formed. 
end
fprintf('The input matrix and vector pertain to a system of linear equations as given below:\n')
for i = 1:n
    disp(eqns(i))
end
%% CONVERGENCE CRITERIA
 
% -We will be comparing diagonal element value with sum of non-diagonal
% ones for a particular row. 
flag = 0;
for i = 1:n
    if norm(A(i,[1:i-1 i+1:n]),1)>abs(A(i,i))
    flag = 1;          % Setting flag 1 if the summation of off-diagonal elements is larger than diagonal element
    end
end
% -norm(A(i,[1:i-1 i+1:n]),1) sums up the absolute values of non-diagonal
% row elements. 
if flag == 1
    disp('WARNING : Matrix is not diagonally dominant. Iteration may not converge.')
end
% -If the convergence criteria is not satisfied the matrix may or may not
% converge. It is a sufficient condition not necessary
%% CONDITION NUMBER
ConditionNumber = norm(A)*norm(inv(A))
if (ConditionNumber>=100)
    disp('WARNING : Condition number is too large.');
end
%% ITERATIONS
redflag = 0;           % This will become 1 if the iteration diverges. Convergence indicator.
for k = 1:K
    T = x;             % Storing previous x vector
    for i = 1:n
        x(i) = (B(i) - A(i,[1:i-1 i+1:n])*x([1:i-1 i+1:n]))/A(i,i);
    end
                       % x would automatically update as Gauss-Siedel algorithm entails.
Norm_Values(k) = norm(x,1);
Data(:,k) = x;     % Storing vector  
fprintf('Iteration %d \n',k)
x'                     % Printing outcome of each iteration
if k>=3
    if( norm(Data(:,k)-Data(:,k-1),1) >  norm(Data(:,k-1)-Data(:,k-2),1))        % Checking for divergence.
        redflag = 1; 
        break
    end
end
if (norm((x-T),1) <  Tolerance)                                                  % Checking for convergence.
    fprintf('GaussSeidel has converged in %d iterations.\n\n',k)
    break
end
end
if (redflag==1)
    fprintf('ERROR : The solution of the given system of linear equation by Guass-Siedel algorithm diverges !\n ')
    fprintf('This may be mitigated by rearrangement of the matrix. It is advised that the matrix be diagonally dominant for convergence.\n')
    c = input('Do you want to re-enter a rearranged matrix?[Enter y to continue, otherwise function will exit.] ', 's');
    if c == 'y'
        A = input('Enter new matrix')                                              %Recursion of GaussSiedel. 
        [x,Data,Norm_Values] = GaussSiedel(A,B,Tolerance);
    else
      disp('GaussSiedel exiting')  
    
    end
else
    fprintf('Solution vector by the Gauss-Siedel algorithm is:\n');
    for i = 1:n
        fprintf('X%d = %f\n',i,x(i))
    end
    fprintf('The exact solution by other methods is given below for comparision:\n')
    Exact = A\B    % As calculated by Matlab
end
end



