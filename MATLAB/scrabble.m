function dataRecovered = scrabble(data,parameter,nIter,nIter_inner,error_inner_threshold,error_outer_threshold)

%% DESCRIPTION
% SCRABBLE: scrabble is used to impute the scRNAseq combining with
% bulk RNAseq of thesame tissue of cell population

% data: data is a structure data type. It contains two types of data,
% scRNA-Seq and bulk RNA-Seq data, which are consistently normalized data.
% The data matrix of scRNAseq data has the gene as the row and the cell as 
% the column. The bulk RNA-seq data is the column vector of the average 
% expression of genes across all the bulk samples. The number of rows of 
% the scRNAseq and bulk RNAseq is the same. To reduce the computation cost,
% we recommend to reduce the number of genes in the data using the experesion
% levels of the genes as the standard. In the real data, we could get 
% around 10,000 genes after setting the threshold as 10 UMIs. 
% 
% parameter: parameter is a vector containing two elements. The first element
% is the alpha parameter, which is the weight for matrix rank. 
% The second element is the beta parameter, which is the weight for 
% the consistency between the aggregated scRNA-Seq data and 
% the bulk RNA-Seq data. The recommendation range of the \alpha is from 1 to 500.
% The recommendation range of \beta is proportional to \alpha and the size of the 
% input data matrix. 
% 
% nIter: nIter is the maximal number of iterations. The recommmdation times of
% iteration are 20.

% AUTHOR: Tao Peng
% Email: pengt@email.chop.edu

%% Main Codes
% Check the inputs of SCRABBLE
if nargin < 2
   disp('The inputs of SCRABBLE are not enough. Please check your input.')
   return;
end

if nargin == 2
   nIter = 20;
   nIter_inner = 20;  % the maximum iteration times in the inner optimization
   error_inner_threshold = 1e-4; % set up the threshold for the inner loop
   error_outer_threshold = 1e-4; % set up the threshold for the outer loop
end

if nargin == 3
   nIter_inner = 20;  % the maximum iteration times in the inner optimization
   error_inner_threshold = 1e-4; % set up the threshold for the inner loop
   error_outer_threshold = 1e-4; % set up the threshold for the outer loop
end


if nargin == 4
   error_inner_threshold = 1e-4; % set up the threshold for the inner loop
   error_outer_threshold = 1e-4; % set up the threshold for the outer loop
end

if nargin == 5
   error_outer_threshold = 1e-4; % set up the threshold for the outer loop
end

% Determine if the data is consistent or not.
if ~isempty(data.data_sc) && ~isempty(data.data_bulk)
    disp('The scRNAseq data and bulk RNAseq data is loaded.')
    if size(data.data_sc,1) ~= length(data.data_bulk)
        disp('The scRNAseq data and bulk RNAseq data is not consistent.')
        return;
    end
end

% Define the parameter 
% generate the observation data
Y = data.data_sc';                   % prepare the scRNAseq for the model
zones = (Y > 0);                     % calculate the projection operator
alpha = parameter(1);                % define the weight of rank
beta = parameter(2);                 % define the weight of consistency of bulk RNAseq
gamma = parameter(3);           % define the regularization coefficient

% generate the bulk RNAseq.
% Here we consider two cases. One is the input including the bulk RNAseq
% The other is the input without the bulk RNAseq
if isempty(data.data_bulk)       
    % the case without bulk RNAseq
    beta = 0*beta;
    Z = ones(1,size(Y,2));
    disp('The inputs are without bulk RNAseq data.')
else
    % the case with bulk RNAseq
    Z = data.data_bulk'*size(Y,2);
end

% Prepare the matrices for the following iteration
n = size(Y,1);
D = ones(1,n);   % Construct the matrix to calculate the aggregated bulk RNAseq
A = beta*(D'*D) + gamma*eye(n); % Construct the iteration matrix A  
B = beta*D'*Z + Y;              % Construct the iteration matrix B

% initailize the Y, X and Lambda for iterations
X = Y;                % set up the initial value for the solution
newX = X;             % set up the iterative value for the solution
newY = Y;             % set up the iterative value Y in the model
Lambda = zeros(size(Y)); % set up the intial value for Lambda in the model

% set up the thresholds for iteration times and the errors
k = 0;         % initialize the iteration times
error = 1;     % initialize the error
% Since we save the time of optimization, we found that the 2 can get good 
% performance. If we could not get good performance, we could increase the
% iteration times here.

gamma = double(gamma);        % make the scalar double 
Y = double(Y);                % make the matrix double
B = double(B);                % make the matrix double
A = double(A);                % make the matrix double
zones = double(zones);        % make the matrix double
Lambda = double(Lambda);      % make the matrix double
newX = double(newX);          % make the matrix double
[m1,n1] = size(X);            % get the size of the data matrix
disp('SCRABBLE begins to run. Please wait ...');
while k < nIter && error > error_outer_threshold 
    % update the X
    X = newX;
    Y = newY;
    l = 1;                        % initialize the inner iteration times
    error_inner = 1;              % initialize the inner iteration errors
    X1 = newX;
    while error_inner > error_inner_threshold && l < nIter_inner           
        newX = cDescent(gamma,Y,B,Lambda,A,zones,newX);
        l = l + 1;                % count the times of the inner iterations
        error_inner = norm(log10(X1 + 1) - log10(newX + 1),'fro')/(n1*m1);  % compute the inner error
        X1 = newX;
        fprintf('The %d-th INNNER iteration and the error is %1.4e\n',l,error_inner)
    end
    % SVT
    S = newX + Lambda/gamma;             % calculate the matrix for SVT
    tau = alpha/gamma;                   % calculate the threshold for SVT
    tic
    [u,s,v] = svt(S,'lambda',tau);
    toc
    newY = u*max(s-tau,0)*v';            % SVT updating
    error = norm(log10(X+1) - log10(newX+1),'fro')/(m1*n1);% calculate the iteration error
    % To modify the error after the first iteration
    if k == 0
        error = 1;
    end
    k = k + 1;
    % update Lambda
    Lambda = Lambda+gamma*(newX - newY);
    fprintf('The %d-th iteration and the error is %1.4e\n',k,error)
end
disp('SCRABBLE finishs the imputation!');

% Return the imputed data matrix
dataRecovered = newX';