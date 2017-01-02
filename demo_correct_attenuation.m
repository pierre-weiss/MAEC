rng(1);
addpath(genpath('../utils/'));

%% Parameters
dynamics=100;
modStat=0;
n_pix=256;
max_alpha=0.03;
sigma=1;

beta_name='images/lena_gray_256.tif';
alpha_name='images/cameraman_256.png';

%% Preparing data
beta=double((imread(beta_name)));
beta=imresize(double(beta(:,:)),[n_pix,n_pix]);
beta=beta/max(beta(:))*dynamics;
beta(beta<0)=0;

alpha=double(imread(alpha_name));
alpha=imresize(alpha(:,:,1),[n_pix,n_pix]);
alpha(alpha<0)=0;
alpha=(alpha/max(alpha(:)))*max_alpha;

figure(1);imagesc(beta); colormap(gray);colorbar;axis equal;title('beta')
figure(2);imagesc(alpha); colormap(gray);colorbar;axis equal;title('alpha')

%% Image formation model
c1=cumsum(alpha,1);
c2=cumsum(alpha(end:-1:1,:),1);c2=c2(end:-1:1,:);

u1=poissrnd(beta.*exp(-c1));
u2=poissrnd(beta.*exp(-c2));

figure(3);imagesc(u1); colormap(gray);colorbar;axis equal;title('u_1')
figure(4);imagesc(u2); colormap(gray);colorbar;axis equal;title('u_2')

%% TV Minimization
lambda=0.07;
nit=100;
u0=u1+u2;
s0=exp(-c1)+exp(-c2);
uu=u0./s0;
c1=0.1;c2=0.1;
tic;[beta_best_TV,CF]=Min_beta_SDMM(alpha,uu,u1,u2,lambda,nit,c1,c2);toc;

figure(5); colormap gray;imagesc(beta_best_TV);
title(sprintf('TV -- SNR:%1.2f',SNR(beta,beta_best_TV)));colorbar;axis equal

%% Direct estimate
beta_best_directestimate=uu;

%This is for the MLE estimate
%uu1=u1./exp(-c1);Sigma1=sigma./exp(-c1);
%uu2=u2./exp(-c2);Sigma2=sigma./exp(-c2);
%beta_best_directestimate=(uu1./Sigma1+uu2./Sigma2)./(1./Sigma1+1./Sigma2);
figure(6);colormap gray;imagesc(beta_best_directestimate);axis equal;title(sprintf('Direct estimate -- SNR:%1.2f',SNR(beta,beta_best_directestimate)));colorbar;axis equal

%% Saves result
imwrite(rescaleUINT8(beta),'XP3_beta.png')
imwrite(rescaleUINT8(alpha),'XP3_alpha.png')
imwrite(rescaleUINT8(u1),'XP3_u1.png')
imwrite(rescaleUINT8(u2),'XP3_u2.png')
imwrite(rescaleUINT8(beta_best_directestimate),'XP3_MLE.png')
imwrite(rescaleUINT8(beta_best_TV),'XP3_betaTV_.png')
