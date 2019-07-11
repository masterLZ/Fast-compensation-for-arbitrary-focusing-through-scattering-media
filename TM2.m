clear;
N=16;
Nmode = N^2;%SLM和CCD的模式
NL=1080; NH=1920;%SLM的分辨率
NCL=1024;NCH=1024;%SLM受控部分
Phase=zeros(NL,NH);
Phase_input=Phase;
ControlPhase=Phase(28:1052-1,448:1472-1);
H1 =(hadamard(Nmode)+1)/2*pi; 
H = hadamard(N);%产生N*N大小的的Hadamard矩阵基

Phase_N=mat2cell(ControlPhase,ones(N,1).*NCL/N,ones(N,1).*NCH/N);%将向量分块

Phaseinit=zeros(NCL,NCH,4);
Phaseinit(:,:,2)=ones(NCL,NCH)*pi/2;
Phaseinit(:,:,3)=ones(NCL,NCH)*pi;
Phaseinit(:,:,4)=ones(NCL,NCH)*pi/2*3;
Phase_i=Phase_N;
h=waitbar(0,'移相计算');
for i=1:Nmode
%     buffer=zeros(N,N);
%     buffer(i)=1;
%     Phase_buffer=H*buffer*H;
%     Phase_buffer=(Phase_buffer+1)/2;
  for j=0:3%4步移相 
        buffer= H1(:,i);
        buffer(buffer>0)=angle(exp(1i*(pi+j*pi/2)))+pi;
        buffer= reshape(buffer,N,N);
        for k=1:Nmode
             Phase_i(k)=mat2cell(ones(NCL/N,NCH/N)*buffer(k),NCL/N,NCH/N);
        end
         Phase_HB = cell2mat(Phase_i);%获得本次相位  .         
         Phase_input(28:1052-1,448:1472-1)=Phase_HB/(2*pi);%获得输入相位
         str=['C:\Users\hp\Desktop\代码\TM\image256\',num2str(i),'_',num2str(j),'.bmp'];
         imwrite(Phase_input,str);
   end
    waitbar(i/Nmode,h,sprintf('%.3f%%',i/Nmode*100));
end
close(h)