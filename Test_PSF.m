clc
clear all
%X =   im2double( imread('cameraman.tif'));
%X = im2double(imread('4.1.05.tiff'));
X = im2double(imread('peppers.png'));
%X = im2double(imread('papav256.pbm'));
[N,~]=size(X);
A=blur(N,6,4);
AA(:,:,1)=0.80*A; AA(:,:,2)=0.10*A; AA(:,:,3)=0.10*A;
BB(:,:,1)=A; BB(:,:,2)=zeros(N); BB(:,:,3)=zeros(N);
AA=fft(AA,[],3);
BB=fft(BB,[],3);
% AA=ndsparse(AA);
% BB=ndsparse(BB);
%-------------------------------------------------------------------------
B_exact=operator(AA,BB,X); 
E = randn(size(B_exact));
E = E/tnorm(E); 
nu=0.01*tnorm(B_exact)*E;
B = B_exact+nu;
%-------------------------------------------------------------------------
s=tnorm(nu);
eta=1.1;
mu0=10;
sigma=10^-3; % 10^-1; % % Noise level
lambda=1;
%x0=zeros(n,1); % Initial guess
ni=4;
ne=4;
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% restrt=5;
% max_it=10;
% X0=zeros(size(double(X)));
%-------------------------------------------------------------------------
tic
%[Xmu,mu]=Gl_GMRES(AA,BB,B,X0,restrt,max_it);
%[xapprox,Errel,Res,volume,TotIt,stopReSt] = ReStart_GMRES(AA,BB,vec(B),sigma,eta,7,ne,ni,vec(X));
%[xapprox,Errel,Discr,Lambda,volume,TotIt,stopReSt] = ReStart(AA,BB,vec(B),sigma,eta,ne,ni,vec(X));
[xapprox,Lambda] = RestGmres(AA,BB,vec(B),ne,ni);
toc
Xmu=reshape(xapprox(:,ne),[N,N,3]);
%toc
% Xmu(:,:,1)=Ymu(:,:,1)/A;
% Xmu(:,:,2)=Ymu(:,:,2)/A;
% Xmu(:,:,3)=Ymu(:,:,3)/A;
  RE=tnorm(X-Xmu)/tnorm(X)
  fprintf("RE=\n");disp(RE)
  SNR=snr(double(X), double(Xmu));
  fprintf("SNR=\n");disp(SNR)
%  fprintf("Tensor GGKB steps=\n");%disp(k)
%  fprintf("Regularization parameter=\n");disp(Lambda)
% % subplot(131); imshow(double(X),[]); title('Original Image');
% subplot(132);imshow(double(B),[]);title('Blurred and Noisy Image');
% subplot(133);imshow(double(Xmu),[]);title('Restored Image');
% h1=figure;
% imshow(double(X),[]);
% set(h1,'PaperSize',[6.5 6]); %set the paper size to what you want  
% print(h1,'papavo','-dpdf') % then print it
% h2=figure;
% imshow(double(B),[]);
% set(h2,'PaperSize',[6.5 6]); %set the paper size to what you want  
% print(h2,'papavb','-dpdf') % then print it
% h3=figure;
% imshow(double(Xmu),[]);
% set(h3,'PaperSize',[6.5 6]); %set the paper size to what you want  
% print(h3,'papavr','-dpdf') % then print it