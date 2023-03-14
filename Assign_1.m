clc,clear,close all;

%% image read
I=imread("MESSI_image.jpeg");
I=rgb2gray(I);
I=double(I)/255;
[H,W]=size(I);
figure,
imshow(I);
sigma_n=sqrt(1e-3);
h=32;w=32; % patch height and width
n=h*w;
N_K=16;
Ks=round(2.^linspace(1,10,N_K));

U=dctmtx(h);
[I_Compressed_DCT,I_Compressed_DCT_complete,Err_DCT]=myCompression(I,U,h,w,H,W,N_K);
U=eye(h,w);
[I_Compressed_IDT,I_Compressed_IDT_complete,Err_IDT]=myCompression(I,U,h,w,H,W,N_K);


U=dctmtx(h);
U(U>0)=1;
U(U<0)=-1;

%% Plot results
for k=1:N_K
    subplot(ceil(N_K/4),4,k);
    imshow(I_Compressed_DCT_complete{k});
    title(sprintf("K=%d, Err=%.4f",Ks(k),Err_DCT(k)));
end
sgtitle("DCT maxk compression results");

%% Plot results
for k=1:N_K
    subplot(ceil(N_K/4),4,k);
    imshow(I_Compressed_IDT_complete{k});
    title(sprintf("K=%d, Err=%.4f",Ks(k),Err_IDT(k)));
end
sgtitle("IDT maxk compression results");
%% Plot err comparison
figure,
plot(Ks,Err_DCT);
hold on 
plot(Ks,Err_IDT);
hold off
legend("DCT compression error","I compression error");
set(gca,'XScale','log')





%% Functions 

function [I_Compressed,I_Compressed_complete,Err]=myCompression(I,U,h,w,H,W,N_K)
h_group=H/h; w_group=W/w;

Ks=round(2.^linspace(1,10,N_K));
I_Patches=cell(h_group,w_group);
J_Patches=cell(h_group,w_group);

% processing
I_Compressed=cell(h_group,w_group,N_K);
I_Compressed_complete=cell(N_K,1);
figure;
D=U;
for k=1:N_K
    K=Ks(k);
    for i=1:h_group
        for j=1:w_group
            I_Patches{i,j}=I(((i-1)*w+1):i*w,((j-1)*h+1):j*h);
            I_ij=I_Patches{i,j};
            
            J_ij=D'*I_ij*D;


            J_Patches{i,j}=J_ij;
            vec_J_ij=reshape(J_ij,[w*h,1]);
            [~,idx]=maxk(abs(vec_J_ij),K);
            vec_J_ij_maxk=zeros(size(vec_J_ij));
            vec_J_ij_maxk(idx)=vec_J_ij(idx);
            J_ij_maxk=reshape(vec_J_ij_maxk,[h,w]);
            I_maxk=D*J_ij_maxk*D';
            I_Compressed{i,j,k}=I_maxk;
        end
    end
    I_Compressed_complete{k}=[];
    for i=1:h_group
        I_Compressed_complete_row=[];
        for j=1:w_group
            I_Compressed_complete_row=[I_Compressed_complete_row, I_Compressed{i,j,k}];
        end
        I_Compressed_complete{k}=[I_Compressed_complete{k};I_Compressed_complete_row];
    end
end
% Calculate Err
Err=zeros(N_K,1);
for k=1:N_K
    err=0;
    for i=1:h_group
        for j=1:w_group
            err=err+norm(I_Compressed{i,j,k}-I_Patches{i,j},2);
        end
    end
    Err(k)=err/h_group/w_group;
end
end
