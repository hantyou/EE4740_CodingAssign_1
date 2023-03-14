clc,clear,close all;
rng(990611)
%% image read
I=imread("MESSI_image.jpeg");
I=rgb2gray(I);
I=imresize(I,0.5);
I=double(I)/255;
[H,W]=size(I);
figure,
imshow(I);
sigma_n=sqrt(1e-3);
h=16;w=16; % patch height and width
n=h*w;
N_K=3;
Ks=round(2.^linspace(8,10,N_K));
Ks=round(h*w.*[0.2, 0.4,0.6,0.8]);
% Ks=round(h*w.*[0.2, 0.25]);
N_K=length(Ks);
h_group=H/h; w_group=W/w;
I_Patches=cell(h_group,w_group);
for i=1:h_group
    for j=1:w_group
        I_Patches{i,j}=I(((i-1)*w+1):i*w,((j-1)*h+1):j*h);
    end
end
%% parpool startup
p = gcp('nocreate');
if isempty(p)
    poolNum=min(w_group*h_group,16);
    parpool(poolNum);
end
for lambda=[10]
    %% OMP sampling and recovery exp
    D=dctmtx(h);
    method="OMP";
    [I_Recovered_OMP,I_Recovered_OMP_complete]=doExp(I,D,h,w,Ks,method);
    %% ISTA sampling and recovery exp
    D=dctmtx(h);
    method="ISTA";
    [I_Recovered_ISTA,I_Recovered_ISTA_complete]=doExp(I,D,h,w,Ks,method,lambda);
    
    %% Calculate Err OMP
    Err_OMP=zeros(N_K,1);
    for k=1:N_K
        err=0;
        for i=1:h_group
            for j=1:w_group
                err=err+norm(I_Recovered_OMP{i,j,k}-I_Patches{i,j},2);
            end
        end
        Err_OMP(k)=err/h_group/w_group;
    end
    Err_OMP_db=mag2db(Err_OMP);
    
    %% Calculate Err ISTA
    Err_ISTA=zeros(N_K,1);
    for k=1:N_K
        err=0;
        for i=1:h_group
            for j=1:w_group
                err=err+norm(I_Recovered_ISTA{i,j,k}-I_Patches{i,j},2);
            end
        end
        Err_ISTA(k)=err/h_group/w_group;
    end
    Err_ISTA_db=mag2db(Err_ISTA);
    %% Plot results
    gcf=figure;
    for k=1:N_K
        subplot(ceil(N_K/2),2,k);
        imshow(I_Recovered_OMP_complete{k});
        title(sprintf("S=%d, Err=%.4f dB",Ks(k),Err_OMP_db(k)));
    end
    fname="results_noGit\OMP recover results";
    fname=strcat(fname," lambda ",num2str(lambda),".png");
    sgtitle("OMP recover results");
    
    fig_export=hgexport('factorystyle');
    fig_export.Resolution=300;
    fig_export.Format='png';
    fig_export.Width=10;
    fig_export.Height=5;
    hgexport(gcf,fname,fig_export);
%     savefig(gcf,fname);
    close gcf;
    
    gcf=figure;
    for k=1:N_K
        subplot(ceil(N_K/2),2,k);
        imshow(I_Recovered_ISTA_complete{k});
        title(sprintf("S=%d, Err=%.4f dB",Ks(k),Err_ISTA_db(k)));
    end
    fname="results_noGit\ISTA recover results";
    fname=strcat(fname," lambda ",num2str(lambda),".png");
    sgtitle("ISTA recover results");
    hgexport(gcf,fname,fig_export);
    % savefig(gcf,fname);
    close gcf;
    %% Plot Err comparison
    gcf=figure;
    plot(Ks,Err_OMP_db,'-*');
    hold on ;
    plot(Ks,Err_ISTA_db,'-^');
    hold off
    legend("OMP","ISTA");
    fname="results_noGit\Error comparison of OMP and ISTA";
    fname=strcat(fname," lambda ",num2str(lambda),".png");
    title("Error comparison of OMP and ISTA");
    ylabel("MSE (dB)")
    xlabel("Sparsity S")
    fig_export.Width=5;
    fig_export.Height=3;
    hgexport(gcf,fname,fig_export);
    % savefig(gcf,fname);
    close gcf
    %% Calculate PSNR OMP
    PSNR_OMP=zeros(N_K,1);
    for k=1:N_K
        err=0;
        for i=1:h_group
            for j=1:w_group
                err=err+psnr(I_Recovered_OMP{i,j,k},I_Patches{i,j});
            end
        end
        PSNR_OMP(k)=err/h_group/w_group;
    end
    
    %% Calculate PSNR ISTA
    PSNR_ISTA=zeros(N_K,1);
    for k=1:N_K
        err=0;
        for i=1:h_group
            for j=1:w_group
                err=err+psnr(I_Recovered_ISTA{i,j,k},I_Patches{i,j});
            end
        end
        PSNR_ISTA(k)=err/h_group/w_group;
    end
    %% Plot PSNR results
    gcf=figure;
    for k=1:N_K
        subplot(ceil(N_K/2),2,k);
        imshow(I_Recovered_OMP_complete{k});
        title(sprintf("S=%d, PSNR=%.4f dB",Ks(k),PSNR_OMP(k)));
    end
    fname="results_noGit\OMP recover results PSNR";
    fname=strcat(fname," lambda ",num2str(lambda),".png");
    sgtitle("OMP recover results");
    
    fig_export=hgexport('factorystyle');
    fig_export.Resolution=300;
    fig_export.Format='png';
    fig_export.Width=10;
    fig_export.Height=5;
    hgexport(gcf,fname,fig_export);
    % savefig(gcf,fname);
    close gcf;
    
    gcf=figure;
    for k=1:N_K
        subplot(ceil(N_K/2),2,k);
        imshow(I_Recovered_ISTA_complete{k});
        title(sprintf("S=%d, PSNR=%.4f dB",Ks(k),PSNR_ISTA(k)));
    end
    fname="results_noGit\ISTA recover results PSNR";
    fname=strcat(fname," lambda ",num2str(lambda),".png");
    sgtitle("ISTA recover results");
    hgexport(gcf,fname,fig_export);
    % savefig(gcf,fname);
    close gcf;
    %% Plot PSNR comparison
    gcf=figure;
    plot(Ks,PSNR_OMP,'-*');
    hold on ;
    plot(Ks,PSNR_ISTA,'-^');
    hold off
    legend("OMP","ISTA",'Location','southeast');
    fname="results_noGit\PSNR comparison of OMP and ISTA";
    fname=strcat(fname," lambda ",num2str(lambda),".png");
    title("PSNR comparison of OMP and ISTA");
    ylabel("PSNR (dB)")
    xlabel("Sparsity S")
    fig_export.Width=5;
    fig_export.Height=3;
    hgexport(gcf,fname,fig_export);
    % savefig(gcf,fname);
    close gcf
    
    
end
save("data4plot.mat");
%%
function [s_hat,H_hat,idices] = OMP(y,H,sigma_n)
[M,N]=size(H);
maxIter=N;
iterflag=1;
epsilon=1e-4;
r_i=y;
y_i=zeros(size(y));
idices=zeros(1,N);
iterNum=0;
H_i=zeros(size(H));
while iterflag
    iterNum=iterNum+1;
    %%% OMP goes here
    r=abs(r_i'*H); % calculate correlation
    idx=find(r==max(r)); % find maximum idx
    idices(idx)=1;
    H_i(:,idx)=H(:,idx); % update H_i
    P_i=H_i*pinv(H_i); % calculate projection onto current signal subspace
    r_i=(eye(M)-P_i)*r_i; % calculate residual
    %%% OMP ends here
    if norm(r_i,2)<epsilon+sigma_n||iterNum>=maxIter
        iterflag=0;
    end
end
H_hat=H_i;
s_hat=pinv(H_hat)*y;
end

function [s_hat,H_hat,idices] = ISTA(y,H,lambda,sigma_n)
[M,N]=size(H);
maxIter=N;
iterflag=1;
epsilon=1e-4;

idices=zeros(1,N);
iterNum=0;

J=2*H'*H;
[~,J_D]=eig(J);
J_D=abs(diag(J_D));
beta=J_D(1);
ita=1/2/beta;


s_hat=zeros(N,1);
p=zeros(N,1);

while iterflag
    iterNum=iterNum+1;
    %%% ISTA goes here
    p=s_hat-ita*2*H'*(H*s_hat-y);
    c=lambda*ita;
    s_hat=softThreshold(p,c);
    y_est=H*s_hat;
    %%% ISTA ends here
    if norm(y-y_est,2)<epsilon+sigma_n||iterNum>=maxIter
        iterflag=0;
    end
end
H_hat=[];
idices=[];
end

function u=softThreshold(x,c)
[a,b]=size(x);
u=zeros(size(x));
for i=1:a
    for j=1:b
        if x(i,j)>c
            u(i,j)=x(i,j)-c;
        elseif x(i,j)<-c
            u(i,j)=x(i,j)+c;
        else
            u(i,j)=0;
        end
    end
end
end

function [I_Recovered,I_Recovered_complete]=doExp(I,D,h,w,Ks,method,lambda)
if nargin<=6
    lambda=10;
end
[H,W]=size(I);
sigma_n=sqrt(1e-3);
N_K=length(Ks);
h_group=H/h; w_group=W/w;
I_Patches=cell(h_group,w_group);
J_Patches=cell(h_group,w_group);
J_Patches_vec=cell(h_group,w_group);
for i=1:h_group
    for j=1:w_group
        I_Patches{i,j}=I(((i-1)*w+1):i*w,((j-1)*h+1):j*h);
        J_Patches{i,j}=D'*I_Patches{i,j}*D;
        J_Patches_vec{i,j}=vec(J_Patches{i,j});
    end
end
% Take Sparse Measurement and recover
U_kron=kron(D,D);
I_Recovered=cell(h_group,w_group,N_K);
I_Recovered_complete=cell(N_K,1);
for k=1:N_K
    K=Ks(k);
    D_m=[];
    for p=1:K
        d_m=vec(randi([0 1],h,w));
        d_m(d_m==0)=-1;
        D_m=[D_m;d_m'];% Take several measurement
    end
    H=D_m*U_kron; % concatenated everything into a measurement matrix
    
    tic
    % Prepare the cells that are in proper structure for parallel processing
    J_Patches_vec_par=reshape(J_Patches_vec,[h_group*w_group,1]);
    J_hat_vec_par=cell(h_group*w_group,1);
    I_Recovered_par_k=cell(h_group*w_group,1);
    % Start Parallel Processing
    parfor iter=1:h_group*w_group
        y=H*J_Patches_vec_par{iter}+sigma_n*randn(K,1); % take measurement
        if method=="OMP"
            [J_hat_vec_par{iter},~,~] = OMP(y,H,sigma_n);
        else
            [J_hat_vec_par{iter},~,~] = ISTA(y,H,lambda,sigma_n)
        end
        J_hat_vec_par{iter}=reshape(J_hat_vec_par{iter},[h,w]);
        I_Recovered_par_k{iter}=D*J_hat_vec_par{iter}*D';
    end
    % Dum the results back into matrix structured cells
    J_hat_vec=reshape(J_hat_vec_par,[h_group,w_group]);
    I_Recovered_k=reshape(I_Recovered_par_k,[h_group,w_group]);
    for i=1:h_group
        for j =1:w_group
            I_Recovered{i,j,k}=I_Recovered_k{i,j};
        end
    end
    toc
    I_Recovered_complete{k}=[];
    for i=1:h_group
        I_Recovered_complete_row=[];
        for j=1:w_group
            I_Recovered_complete_row=[I_Recovered_complete_row, I_Recovered{i,j,k}];
        end
        I_Recovered_complete{k}=[I_Recovered_complete{k};I_Recovered_complete_row];
    end
end






end