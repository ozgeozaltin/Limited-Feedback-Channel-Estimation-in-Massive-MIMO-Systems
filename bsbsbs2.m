clear all 
clc

Nof=[64 128 256 512 ]; % Different number of transmitted antennas
NofQ_bit=[2 4 5 8];
%%---------------azimuth angle chosen in uniform PT and PR dictionaries


for iyte=1:4
Mt=Nof(iyte);  
for ip=1:4
Q_bit=NofQ_bit(ip);

Mr=2;
L=randi([5 10],1);
fc=2*10^9; % Carrier frequency=2GHz
c=3*10^8; % Speed of light
lambda=c/fc; % 0.15 Carrier wavelenght(m) 
d=lambda/2;
Gdb=8 ; %dB max directional gain
phi3db =55; % degree 3dB beamwidth
Am=30 ;%db front to back ratio
Gt=140;
Gr=16;
phi= randi([0 360],16,140);
a=-180;
b=180;
Ntr=64;
for r=1:Gr
PR(r)=a+r*(b-a)/(Gr+1);
end
PR_trans=transpose(PR);

 for z=1:Gt
PT(z)=a+z*(b-a)/(Gt+1);
 end


Ct=1;

Cr=1;

%%

%%---------------*********** Azimuth Angle *************-----------------------

%---***/*/*------------BS azimuth Angle---------

for yth=1:Mt

          for j=1:1
              for zth=1:Gt
              u1(j,zth) = (d/lambda)*sin(PT(j,zth).*pi/180);
              aT(yth,zth) = exp(-1i*2*pi*u1(j,zth)*(yth-1)')/sqrt(Mt); % Azimuth Scanning Steering Vector.
              end
    end
end
At=Ct.*aT;
At_h=ctranspose(At);
%------------------------------------------------
%--------------UE azimuth Angle---------
Mr=2;
for yth2=1:Mr

          for j2=1:1
              for zth2=1:Gr
              u2(j2,zth2) = (d/lambda)*sin(PR(j2,zth2).*pi/180);
              aR(yth2,zth2) = exp(-1i*2*pi*u2(j2,zth2)*(yth2-1)')/sqrt(Mr); % Azimuth Scanning Steering VeCror.
              end
    end
end
Ar=Cr.*aR;

%%
%-------------------*******************Channel*********-------------------------


k1=randi([0 40],16,140);
Kx=(sqrt(k1./(k1+1)));
Ky=(1./(k1+1));
%------------------------

for j=1:16
    
    for k=1:140
     alpha(j,k)= Kx(j,k) + 1i.*(Ky(j,k));
    end
end
alpha_f=transpose(alpha.*exp(1i.*phi));

H_compact=Ar*transpose(alpha_f)*At_h;

SS=randi([0 1],Mt,64); % training signal
TY=transpose(SS)*conj(At); 
QQ=kron(TY,Ar);    

YY1=H_compact*SS;
Received=YY1(:);


epsilon=2.25;
L_prime=log2(2240)+2*Q_bit;


%% OMp 
[N,G] = size(QQ);
g_hat= zeros(G,1);
residual=Received;
S_g=[];

for  t= 1:L_prime   
    p= QQ'*residual;  %% the complex transpose and residual 
    
    for ii=1:G
        fx(ii)=abs(p(ii));
    end
    
    if fx<epsilon
        break; 
    end
    
    
    n_star = find(fx == max(fx));
    S_g = union(S_g,n_star);
    g_hat(S_g)=pinv(QQ(:,S_g))* Received;  %%coeffecient update
    residual=Received-QQ*g_hat;   %%residu update
    
  
end


%%
%-LLOYDD algorthm------
% Quantization data: thresholds and reconstruction levels

%%%%%%%%%%%%%%%%%%%Quantization Process of REAL PART
n1=1024;
L=2^Q_bit;
 vmax=max(real(g_hat));
 vmin=min(real(g_hat));
 del=(vmax-vmin)/L;
 [part,code] = lloyds(real(g_hat),L);
 
 [ind_real,q_Real]=quantiz(real(g_hat),part,code);                     % Quantization process
                                                                      % ind contain index number and q contain quantized  values
 l1=length(ind_real);
 l2=length(q_Real);
  
 for i=1:l1
    if(ind_real(i)~=0)                                            % To make index as binary decimal so started from 0 to N
       ind_real(i)=ind_real(i)-1;
    end 
    i=i+1;
 end   
  for i=1:l2
     if(q_Real(i)==vmin-(del/2))                          % To make quantize value inbetween the levels
         q_Real(i)=vmin+(del/2);
     end
  end

%---------- Encoding Process
 code=de2bi(ind_real,'left-msb');             % Cnvert the decimal to binary
 k=1;
for i=1:l1
    for j=1:Q_bit
        coded(k)=code(i,j);                  % convert code matrix to a coded row vector
        j=j+1;
        k=k+1;
    end
    i=i+1;
end


 %  -------- Demodulation Of PCM signal
 qunt=reshape(coded,Q_bit,length(coded)/Q_bit);
 index=bi2de(qunt','left-msb');                    % Getback the index in decimal form
 q_Real=del*index+vmin+(del/2);                       % getback Quantized values

 
 %%%%%%%%%%%%%%% QUANTIZATION IMAGINNARY PART 슬슬슬슬슬�
 vmax_imag=max(imag(g_hat));
 vmin_imag=min(imag(g_hat));
 del_imag=(vmax_imag-vmin_imag)/L;
 [part_imag,code_imag] = lloyds(real(g_hat),L);
 [ind_imag,q_imag]=quantiz(imag(g_hat),part_imag,code_imag);                     % Quantization process
                                                                      % ind contain index number and q contain quantized  values
 l1_imag=length(ind_imag);
 l2_imag=length(q_imag);
  
 for i=1:l1_imag
    if(ind_imag(i)~=0)                                            % To make index as binary decimal so started from 0 to N
       ind_imag(i)=ind_imag(i)-1;
    end 
    i=i+1;
 end   
  for i=1:l2_imag
     if(q_imag(i)==vmin_imag-(del_imag/2))                          % To make quantize value inbetween the levels
         q_imag(i)=vmin_imag+(del_imag/2);
     end
  end

  
 % --------- Encoding Process
 code_imag=de2bi(ind_imag,'left-msb');             % Cnvert the decimal to binary
 k=1;
for i=1:l1_imag
    for j=1:Q_bit
        coded_imag(k)=code_imag(i,j);                  % convert code matrix to a coded row vector
        j=j+1;
        k=k+1;
    end
    i=i+1;
end

 %   Demodulation Of PCM signal
 
 qunt_imag=reshape(coded_imag,Q_bit,length(coded_imag)/Q_bit);
 index_imag=bi2de(qunt_imag','left-msb');                    % Getback the index in decimal form
 q_imag=del_imag*index_imag+vmin_imag+(del_imag/2);                       % getback Quantized values

 
 gg_hat_BS=q_Real +1i*q_imag;
 
 for i=1:2240
 if  gg_hat_BS(i) == gg_hat_BS(1)
 gg_hat_BS_arrived(i)=0;
 end 
 if  gg_hat_BS(i) ~=gg_hat_BS(1)
  gg_hat_BS_arrived(i)=gg_hat_BS(i);    
 end
 end
 

GG_arrived_BS=reshape(gg_hat_BS_arrived,16,140);
 
 
H_estimated=Ar*GG_arrived_BS*At_h;

SNR=10;

B=(H_compact-H_estimated).*conj(H_compact-H_estimated);
NMSE=mean(B./((norm(H_estimated,2))*(norm(H_compact,2))));
NMSEfin(iyte,ip)=mean(NMSE);


CAP(iyte,ip)=log2(det(eye(Mr)+(SNR/Mt)*H_estimated*H_estimated'));
CAPwithPerfect(iyte,ip)=log2(det(eye(Mr)+((10.^(SNR/10))/Mt)*(H_compact)*(H_compact')));


end
end


figure 
for als=1:4
    
    plot( Nof, NMSEfin(als,:));
    hold on
xlabel('Number of transmit antennas')
ylabel('Normalized Mean Square Error')
legend('\qbit=2', 'qbit=4', '\qbit=5', '\qbit=8', 'location', 'northeast')
 
end

