clc;clear;close all;
%% 
% ����ɢ����ʴ������

T_available = GenerateTM(1200);
%% 
% ������λͼ
% 
% $$E_m^{\mathrm{out}} =\sum nk_{\mathrm{mn}} E_n^{\mathrm{in}}$$
%%
N=32;%ͨ����һ������
Nmode = N^2;%SLM��CCD��ģʽ
N1=1024; %CCD��SLM�ķֱ���

Phase=zeros(1200,1200);
I0 = Calculation_Fresnel_diffraction_C(T_available,Phase);
ControlPhase=Phase(88:1112-1,88:1112-1);%����1024
gama=length(ControlPhase(:))/N1^2;

%������任
H = hadamard(N);%����ControlPhase��С�ĵ�Hadamard�����
H1 =hadamard(Nmode);
H2 = inv(H);
Phase_N=mat2cell(ControlPhase,ones(N,1).*N1/N,ones(N,1).*N1/N);%�������ֿ�
Phase_input=Phase;

Kobs=complex(zeros(32*32,Nmode));
h=waitbar(0,'�������');
Phase_i=Phase_N;
Phaseinit=zeros(N1,N1,4);
Phaseinit(:,:,2)=ones(N1)*pi/2;
Phaseinit(:,:,3)=ones(N1)*pi;
Phaseinit(:,:,4)=ones(N1)*pi/2*3;

%% 
% $$I_m^{\alpha } =2\mathrm{Re}\left(e^{i\alpha } \overset{-}{s_m } \sum_n 
% k_{\mathrm{mn}} E_n^{\mathrm{in}} \right)$$
% 
% $\frac{\left(I_m^0 -I_m^{\pi } \right)}{4}+i$$\frac{\left(I_m^{3\frac{\pi 
% }{2}} -I_m^{\frac{\pi }{2}} \right)}{4}$=$s_m k_{\mathrm{mn}}$

parfor i=1:Nmode
     Phase_i=Phase_N;
     
%%�Ȼ�ñ��ι��������ķֲ�
      smB=zeros(32,32);
      If=zeros(1200,1200,4);%4������Ľ��
%       buffer=zeros(N,N);
%       buffer(i)=1;
%       Phase_buffer=H*buffer*H;  
%       Phase_buffer=(Phase_buffer+1)/2;
      buffer= H1(:,i);
      buffer= reshape(buffer,N,N);
      buffer = (buffer+1)/2
      for k=1:Nmode
        % Phase_i(k)=mat2cell(ones(N1/N)*Phase_buffer(k),N1/N,N1/N);
         Phase_i(k)=mat2cell(ones(N1/N)*buffer(k),N1/N,N1/N);
      end
      Phase_HB = cell2mat(Phase_i);%��ñ�����λ 
      Phase_input=zeros(1200);
    for j=0:3%4������     
      Phase_HB1 = Phase_HB*j*pi/2;
      Phase_input(88:1112-1,88:1112-1)=Phase_HB1;%���������λ
      If(:,:,j+1)=Calculation_Fresnel_diffraction_C(T_available,Phase_input);
    end
    smK=(If(:,:,1)-If(:,:,3))/4+1j*(If(:,:,4)-If(:,:,2))/4;
    %ȡ�м�128*128һ��
    smB=smK(536:568-1,536:568-1); 

    Kobs(:,i)=reshape(smB,1,32*32);
  
  %  waitbar(i/Nmode,h,sprintf('%.3f%%',i/Nmode*100));
end
%figure,imagesc(Phase_input);
Kobs=Kobs/(H1);
close(h)
figure,imagesc(abs(Kobs));
%% 
% ��λ����
%%

KobsT=Kobs';
Knorm=Kobs./(abs(Kobs)+eps);
Ofoc_norm=Kobs*Knorm';
figure,imagesc(abs(Ofoc_norm));
%%
clear Ofoc_norm

Etarget=zeros(32,32);
Etarget=complex(Etarget);
Etarget(15,5)=1;
Etarget(15,16)=1;
Etarget(15,27)=1;

figure,imagesc(abs(Etarget))
Ein=KobsT*reshape(Etarget,length(Etarget(:)),1)./(abs(KobsT*reshape(Etarget,length(Etarget(:)),1))+eps);
Ein1=angle(reshape(Ein,N,N));
%Ein1=H*Ein1*H;
%% 
% �������λ�����
% 
% $$E^{\mathrm{in}} =K_{\mathrm{obs}}^{\dagger } E^{\mathrm{target}} /\left|K_{\mathrm{obs}}^{\dagger 
% } E^{\mathrm{target}} \right|$$

out_Phase = Phase_N;
for i = 1:Nmode
    out_Phase(i)=mat2cell(ones(N1/N).*Ein1(i),N1/N,N1/N);
end
out_Phase = cell2mat(out_Phase);
%ע������ӳ��Ӧ���ǰ�С��0�ò��ּ�pi������0�ò��ֲ���
out_Phase=mat2gray(out_Phase)*2*pi;
Phase(88:1112-1,88:1112-1)=out_Phase;
figure,imagesc(Phase);
%% 
% 


I=Calculation_Fresnel_diffraction_C(T_available,Phase);
Iout=I(536:568-1,536:568-1);
figure,subplot(1,2,1),imagesc(Iout),axis equal;
subplot(1,2,2),imagesc(I0(536:568-1,536:568-1)),axis equal;