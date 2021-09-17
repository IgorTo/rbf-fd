clearvars;
close all;

N = 20;
M = 10*N;

X.nds = linspace(0,1,N)';
Y.nds = linspace(0,1,M)';

%
% Parameters.
%
bf = Basis.basisF('PHS', '1d'); % Choose basis functions.
p = 3; % PHS power (r^p).
polydeg = 2; % Augmented polynomial degree.
n = 2*nchoosek(polydeg+1,1); % Stencil size.


%
% Compute evaluation/differentiation matrices.
%
E = Matrices.generate_1d(X.nds,Y.nds,p,n,polydeg,bf,[]);

% Find a cardinal function closest to the middle of the domain.
idx = knnsearch(X.nds,0.5,'k',1);

figure;
scatter(X.nds, zeros(size(X.nds)), 'o', 'filled')
hold on;
plot(Y.nds, E(:,idx), '-b');