rng(1);
%% Parameters
dynamics=100;
modStat=0;
n_pix=256;
max_alpha=0.03;

beta_name='images/lena_gray_256.tif';
alpha_name='images/cameraman_256.png';

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

%% Direct estimate
uu1=u1;uu2=u2;
uu1(u1<1e-16)=1e-16;
uu2(u2<1e-16)=1e-16;
v=log(uu2./uu1);
alpha_direct=1/2*drond1(v);
figure(5); colormap gray;imagesc(alpha_direct);title(sprintf('alpha direct estimate -- SNR: %1.2f',SNR(alpha,alpha_direct)));colorbar;

%% Restoration
lambda=10;
gamma=5;
c3=100; c4=50;
nit=1000;
a0=zeros(size(u1));
tic;[beta_TV,alpha_TV,~]  = Warm_Start_SDMM(u1,u2,dynamics,nit,lambda,gamma,c3,c4,a0);toc;
figure(6); colormap gray;imagesc(alpha_TV);title(sprintf('alpha TV estimate -- SNR: %1.2f',SNR(alpha,alpha_TV)));colorbar;

%% Saving the results
imwrite(rescaleUINT8(beta),'XP2_beta.png')
imwrite(rescaleUINT8(alpha),'XP2_alpha.png')
imwrite(rescaleUINT8(u1),'XP2_u1.png')
imwrite(rescaleUINT8(u2),'XP2_u2.png')
imwrite(rescaleUINT8(alpha_direct),'XP2_alpha_direct.png')
imwrite(rescaleUINT8(alpha_TV),'XP2_alpha_TV.png')