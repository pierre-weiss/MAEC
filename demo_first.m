rng(1);
%% Parameters
dynamics=1000;
modStat=0;
n_pix=256;
max_alpha=0.1;

beta_name='images/insect.png';
alpha_name='images/att.png';

%% Preparing data
beta=double((imread(beta_name)));
beta=imresize(beta,[n_pix,n_pix]);
beta(beta<0)=0;
beta=beta/max(beta(:))*dynamics;

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

%% Restoration
% First, we find a warm start estimate of the attenuation
lambdaa=8;
gamma=5;
c3=100; c4=50;
nita=1000;
a0=zeros(size(u1));
tic;[beta_direct,alpha_TV,~]  = Warm_Start_SDMM(u1,u2,dynamics,nita,lambdaa,gamma,c3,c4,a0);toc;
figure(6); colormap gray;imagesc(alpha_TV);title(sprintf('alpha TV warm start -- SNR: %1.2f',SNR(alpha,alpha_TV)));colorbar;

% Then denoise the density
lambdab=0.005;
nitb=100;
u0=u1+u2;
s0=exp(-c1)+exp(-c2);
uu=u0./s0;
c1=0.1;c2=0.1;
tic;[beta_TV,CF]=Min_beta_SDMM(alpha_TV,uu,u1,u2,lambdab,nitb,c1,c2);toc;

figure(7); colormap gray;imagesc(beta_TV);
title(sprintf('TV -- SNR:%1.2f',SNR(beta,beta_TV)));colorbar;axis equal


%% Saving the results
imwrite(rescaleUINT8(beta),'XP1_beta.png')
imwrite(rescaleUINT8(alpha),'XP1_alpha.png')
imwrite(rescaleUINT8(u1),'XP1_u1.png')
imwrite(rescaleUINT8(u2),'XP1_u2.png')
imwrite(rescaleUINT8(alpha_TV),'XP1_alpha_TV.png')
imwrite(rescaleUINT8(beta_TV),'XP1_beta_TV.png')