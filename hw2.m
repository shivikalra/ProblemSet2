% Shivi Kalra Homework 2 Macro 696G
clear all
close all
tic
global v0 beta delta a kmat k sigma alpha A
%Set up Parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A = [0.977 0.023; 0.074 0.926];
a = [1.1,0.678];
% Specify tolerance for when to stop searching over value fns
tolerance = 1e-06; % maximum tolerance
dis = 1;
num_k = 1000;

% Solve for steady state
kstar = (alpha/(1/beta - (1-delta))^(1/(1-alpha)));
%cstar = a*kstar^(alpha) - delta*kstar;
%istar = delta*kstar;
%ystar = a*kstar^(alpha);
%Set the minumum and maximum value of k between 25% and 200%
kmin = 0.25*kstar; 
kmax = 2.00*kstar;
k = linspace(kmin, kmax, num_k);
kmat = repmat(k', [1 num_k]);
 
%guess a value fn
v0 =zeros(2,num_k);

% Set up consumption and return function
cons(:,:,1) = a(1)*kmat .^ alpha + (1 - delta) * kmat - kmat';
cons(:,:,2) = a(2)*kmat .^ alpha + (1 - delta) * kmat - kmat';
ret = cons .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret(cons < 0) = -Inf;

while dis > tolerance
             
    % compute the utility value for all possible combinations of k and k':
    value_mat(:,:,1) = ret(:,:,1) + beta * (....
        A(1,1)*repmat(v0(1,:),[num_k 1])+...
        A(1,2)*repmat (v0(2,:),[num_k 1]));
    value_mat(:,:,2) = ret(:,:,2) + beta * (....
        A(2,1)*repmat(v0(1,:),[num_k 1])+...
        A(2,2)*repmat (v0(2,:),[num_k 1]));
       
    % find the optimal k' for every k:
    [vfn1, pol_indx1] = max(value_mat(:,:,1), [], 2);
    [vfn2, pol_indx2] = max(value_mat(:,:,2), [], 2);
    vfn1 = vfn1';
    vfn2 = vfn2';
    
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn1 - v0(1,:)));
    dis = max(abs(vfn2 - v0(2,:)));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v0(1,:) = vfn1;
    v0(2,:) = vfn2;
end

g1 = k(pol_indx1);
g2 = k(pol_indx2); % policy function
subplot(2,2,1)
plot(k,vfn1)
title('Ahigh')
subplot(2,2,2)
plot(k,vfn2)
title('A low')

figure
subplot(2,2,1)
plot(k,g1)
title('policy function for Ahigh')
subplot(2,2,2)
plot(k,g2)
title('policy function for Alow')


