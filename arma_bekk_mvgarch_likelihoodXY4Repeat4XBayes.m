function outputLLF = arma_bekk_mvgarch_likelihoodXY4Repeat4XBayes(parameters, data, p, q, kx, ky, Nr, Nl,priors, sig,indexX, indexY)
% input: paramter vector
%          data      -- time series
%          p         -- ar_order
%         q         -- arch_order
%         Nr, Nl  -- number of repeat experiments and the index of the measurements for each
%         experiment
% output: LLF  -- log likelihood function
%            likelihoods --- log likelihood function for each time point
%            Ht -- conditional variances
%          err -- residuals
[T,k] = size(data);
m = max(p,q);
errors = zeros(T,k); % the first m errors are zeros, since we are calculating the conditional QML
hx = zeros(kx,kx); 
hy = zeros(ky,ky);
LLF_x = 0;
LLF_y = 0;
likelihoods_x = zeros(T,1);
likelihoods_y = zeros(T,1);


para_armabekk = reshapeparasXY(parameters, p, q, k, kx, ky);
A = para_armabekk.A;
BX = para_armabekk.BX;
BY = para_armabekk.BY;
constx = para_armabekk.constx;
consty = para_armabekk.consty;


for i = 1 : Nr
    for j = Nl(i,1)+m : Nl(i,2)
        errors(j,:) = data(j,:) - data(j-1,:) * A(:,:,1);
        for l = 2 : p
            errors(j,:) = errors(j,:)- data(j-l,:) * A(:,:,l);
        end
    end
end
uncond = cov(errors);
for i = 1 : Nr
    for j = Nl(i,1) : Nl(i,2)
        if j < Nl(i,1)+m
            hx = uncond(indexX,indexX);
            hy = uncond(indexY,indexY);
         else   
             hx = constx;
             hy = consty;
             for l = 1 : q
                 Z = (data(j-l,:))'*(data(j-l,:));
                 hx = hx + BX(:,:,l)* Z * BX(:,:,l)';
                 hy = hy + BY(:,:,l)* Z * BY(:,:,l)';
             end
             likelihoods_x(j) = kx*log(2*pi)+(log(det(hx)) + errors(j,indexX)*hx^(-1)*errors(j,indexX)');
             LLF_x = LLF_x + likelihoods_x(j);
             likelihoods_y(j) = ky*log(2*pi)+(log(det(hy)) + errors(j,indexY)*hy^(-1)*errors(j,indexY)');
             LLF_y = LLF_y + likelihoods_y(j);
        end
    end
end

sig = sig .* eye(size(parameters,1));% with Bayesian priors (priors, sig)\  here we take diagonal sig matrix
sig_armabekk = reshapeparasXY(diag(sig), p, q, k, kx, ky);
sigA = sig_armabekk.A;
sigBX = sig_armabekk.BX;
sigBY = sig_armabekk.BY;
sigconstx = sig_armabekk.constx^0.5;
sigconsty = sig_armabekk.consty^0.5;

pri_armabekk = reshapeparasXY(priors, p, q, k, kx, ky);
priA = pri_armabekk.A;
priBX = pri_armabekk.BX;
priBY = pri_armabekk.BY;
priconstx = pri_armabekk.constx^0.5;
priconsty = pri_armabekk.consty^0.5;

para_x = [reshape(A(:,indexX,:), k*kx*p, 1); reshape(BX, k*kx*q, 1); reshape(constx^0.5, kx*kx, 1)];
sig_x = [reshape(sigA(:,indexX,:), k*kx*p, 1); reshape(sigBX, k*kx*q, 1); reshape(sigconstx, kx*kx, 1)];
pri_x = [reshape(priA(:,indexX,:), k*kx*p, 1); reshape(priBX, k*kx*q, 1); reshape(priconstx, kx*kx, 1)];
sig_x = diag(sig_x);

para_y = [reshape(A(:,indexY,:), k*ky*p, 1); reshape(BY, k*kx*q, 1); reshape(consty^0.5, ky*ky, 1)];
sig_y = [reshape(sigA(:,indexY,:), k*ky*p, 1); reshape(sigBY, k*ky*q, 1); reshape(sigconsty, ky*ky, 1)];
pri_y = [reshape(priA(:,indexY,:), k*ky*p, 1); reshape(priBY, k*ky*q, 1); reshape(priconsty, ky*ky, 1)];
sig_y = diag(sig_y);

LLF_x = 0.5*(LLF_x);
likelihoods_x = 0.5*likelihoods_x;
LLF_x = LLF_x + 0.5 * log(det(sig_x)) + 0.5 * (para_x - pri_x)' * sig_x^(-1) *(para_x - pri_x);
if isnan(LLF_x)
    LLF_x = 1e6;
end
Bayes_x = 0.5 * log(det(sig_x)) + 0.5 * (para_x - pri_x)' * sig_x^(-1) *(para_x - pri_x);

LLF_y = 0.5*(LLF_y);
likelihoods_y = 0.5*likelihoods_y;
LLF_y = LLF_y + 0.5 * log(det(sig_y)) + 0.5 * (para_y - pri_y)' * sig_y^(-1) * (para_y - pri_y);
if isnan(LLF_y)
    LLF_y = 1e6;
end
Bayes_y = 0.5 * log(det(sig_y)) + 0.5 * (para_y - pri_y)' * sig_y^(-1) *(para_y - pri_y);

outputLLF.LLF_x = LLF_x;
outputLLF.LLF_y = LLF_y;
outputLLF.likelihoods_x = likelihoods_x;
outputLLF.likelihoods_y = likelihoods_y;
outputLLF.Bayes_x = Bayes_x;
outputLLF.Bayes_y = Bayes_y;





