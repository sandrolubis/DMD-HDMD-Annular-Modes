clc
clear all
close all

M = 4;           %% sampling freq 1:6h; 4:24h
P = 2;           %% \tau of DMD or HDMD (day)
q = 5;           %% delay
r = 500;         %% DMD truncation value 

NumEns = 14;     %% # of datasets used for developing ROM

for ens=1:NumEns
    tic
    ens
    load(['Data/CTLrun' num2str(ens) '/ZonalWind.mat'], 'u4xDaily', 'y', 'pf')
    load(['Data/CTLrun' num2str(ens) '/Temp.mat'], 'T4xDaily')
    np=size(u4xDaily,1)/2*size(u4xDaily,2);
    nd=size(u4xDaily,3); %number of samples in one set
    
    N = floor(nd/M);
    if(ens==1)
        XX=zeros(2*np*q,(N-q+1-P)*NumEns*2,'single');
        YY=zeros(2*np*q,(N-q+1-P)*NumEns*2,'single');
    end        
    
    uN=zeros(48,39,floor(nd/M),'single');  
    uS=uN;    
    TS=uN;    %% Northern hemisphere zonal wind
    TN=uN;
    for n=M:M:nd
        uN(:,:,n/M)= u4xDaily(49:end,:,n);    %% Northern hemisphere zonal wind
        TN(:,:,n/M)= T4xDaily(49:end,:,n);    %% Northern hemisphere temperature
        uS(:,:,n/M)= u4xDaily(48:-1:1,:,n);   %% Southern hemisphere zonal wind
        TS(:,:,n/M)= T4xDaily(48:-1:1,:,n);   %% Southern hemisphere temperature
    end
    clear u4xDaily T4xDaily
    days=size(uN,3);
    
    disp('Remove mean ...')
    uNave=mean(uN,3);
    uSave=mean(uS,3);   
    TNave=mean(TN,3);
    TSave=mean(TS,3);
    uNa=uN-uNave;
    uSa=uS-uSave;
    TNa=TN-TNave;
    TSa=TS-TSave;
    clear uN uS TN TS
    
    disp('std ...')
    sduN=zeros(48,39,'single');
    sduS=sduN;
    sdTS=sduN;
    sdTN=sduN;
    for j=1:48
        for k=1:39
            sduN(j,k)= std(squeeze(uNa(j,k,:)));     
            sdTN(j,k)= std(squeeze(TNa(j,k,:)));   %% calculating zonal-wind and temperature
            sduS(j,k)= std(squeeze(uSa(j,k,:)));   %% standard deviations at different grid points
            sdTS(j,k)= std(squeeze(TSa(j,k,:)));   
        end
    end
        
    disp('Weighting ...')
    for k=1:39    %% calculating area-averaged standard deviation at different pressure levels
        sdU(k,ens)=mean(0.5*squeeze(sduN(:,k)+sduS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
        sdT(k,ens)=mean(0.5*squeeze(sdTN(:,k)+sdTS(:,k)).*cos(y(49:end)*pi/180.0)/cos(45.0*pi/180.0));
    end
    for n=1:days 
        for j=1:48       %% normalization of u and T     
            uNa(j,:,n)=squeeze(uNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
            uSa(j,:,n)=squeeze(uSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdU(:,ens))';
            TNa(j,:,n)=squeeze(TNa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
            TSa(j,:,n)=squeeze(TSa(j,:,n))*sqrt(cos(y(48+j)*pi/180.0))./squeeze(sdT(:,ens))';
        end
    end
    
    disp('Vectorizing ...')
    XNa=zeros(2*np,days);
    XSa=XNa;
    for n=1:days    %% Building state vectors 
        XNa(:,n)=[(reshape(uNa(:,:,n),np,1));(reshape(TNa(:,:,n),np,1))];
        XSa(:,n)=[(reshape(uSa(:,:,n),np,1));(reshape(TSa(:,:,n),np,1))];
    end
    clear uNa uSa TNa TSa
    
    HNa=zeros(2*np*q,N-q+1,'single');
    HSa=HNa;
    for j=1:N-q+1
        for i=1:q    %% Arranging data into Hankel matrices
            HNa((i-1)*2*np+1:i*2*np,j)=XNa(:,i+(j-1));
            HSa((i-1)*2*np+1:i*2*np,j)=XSa(:,i+(j-1));
        end
    end
    clear XNa XSa
    
    XN0=HNa(:,1:end-P);
    XNp=HNa(:,1+P:end);

    XS0=HSa(:,1:end-P);
    %size(XS0)
    XSp=HSa(:,1+P:end);
    clear HSa HNa
    
    XX(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XN0 XS0];
    YY(:,(ens-1)*2*(N-q+1-P)+1:(ens)*2*(N-q+1-P)) = [XNp XSp];

    clear XN0 XS0 XNp XSp
    toc
end

whos

disp('SVD ...')
tic
[U, S, V] = svd(XX, 'econ');
disp('SVD done')
clear XX
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);
AHDMD = U'*YY*V/S;     %% Hankel-DMD execuation
[EVec,EVal] = eig(AHDMD);
disp('EIG done')
G = YY*V/S;
clear YY V
for n=1:size(EVal,1)
    phi(:,n)=G*EVec(:,n)/EVal(n,n);   %% Hankel-DMD modes in delay coordinates
end   
[EVal, I] = sort(abs(diag(EVal)), 'descend');
phi = phi(:, I);          %% sorting Hankel-DMD modes

phi_n = zeros(size(phi, 1)/q, size(phi, 2));
for(i = 1:size(phi, 2))  %% Hankel-DMD modes in physical space (= their last block in delay coordinates)
    phi_n(:, i) = phi(end-size(phi, 1)/q+1:end, i);   
end
phi = phi_n;           

toc

disp('Saving ...')
save(['GCM_Modes_3M/Modes_Lagday' num2str(P) '_q' num2str(q) '_r' num2str(r) '.mat'], ... 
    'AHDMD', 'r', 'q', 'M', 'y', 'pf', 'U', 'S', 'phi', 'sdU', 'sdT', '-v7.3')
