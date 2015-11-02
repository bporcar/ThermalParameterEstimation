# ThermalParameterEstimation


clear a(I,J); clear P; clear A; clear B; clear dA; clear dB; clear dAB; clear Q, clear dUA1, clear dUA2, clear dgA;
clear dUA; clear QUA; clear UA1UA2; clear R;

a='';b='';da='';db='';A='';B='';dA='';dB='';dAB='';P='';Q='';QA='';dUA='';dUA1='';dUA2='';dgA='';UA1UA2='';

for I=1:no,for J=1:no,if I==J a(I,J)=1;else a(I,J)=0;end;for K=1:no*no*na+no*ni*nb,da(I,J,K)=0;end;end;end;
for I=1:no,for J=1:ni,b(I,J)=0;for K=1:no*no*na+no*ni*nb,db(I,J,K)=0;end;end;end;


for w=1:(no*no*na+no*ni*nb);
   
   clear dAB(:,:,w);
   end

for I=1:no,for J=1:no,for K=0:na-1,a(I,J)=a(I,J)-M(3,(I-1)*no*na+(I-1)*ni*nb+J+no*K);da(I,J,(I-1)*no*na+(I-1)*ni*nb+J+no*K)=-1;end;end;end
for I=1:no,for J=1:ni,for K=0:nb-1,b(I,J)=b(I,J)+M(3,I*no*na+(I-1)*ni*nb+J+ni*K);db(I,J,I*no*na+(I-1)*ni*nb+J+ni*K)=1;end;end;end

P=M(4:no*no*na+no*ni*nb+3,1:no*no*na+no*ni*nb);

A=sum(a,1);
B=sum(b,1);

dA=sum(da,1);
dB=sum(db,1);

for w=1:(no*no*na+no*ni*nb);
   
   dAB(:,:,w)=[dA(:,:,w) dB(:,:,w)];
   end

dAB=reshape(dAB,(no+ni),(no*no*na+no*ni*nb));
Q=dAB*P*transpose(dAB);

UA1=-A(1)/A(3)
UA2=A(2)/A(3)
gA=-B(1)/A(3)

for hhh=1:(ni+no) dUA1(hhh)=0;   dUA2(hhh)=0; dgA(hhh)=0; end

dUA1(1)=-1/A(3); dUA1(3)=-UA1/A(3);
dUA2(2)=1/A(3); dUA2(3)=-UA2/A(3); 
dgA(3)=-gA/A(3); dg(no+1)=-1/A(3);


dUA=[dUA1 dUA2];

dUA=transpose(reshape(dUA,ni+no,2));

QUA=dUA*Q*transpose(dUA);

m=[1 -1];
UA1UA2=transpose([UA1 UA2])-QUA*transpose(m)*inv(m*QUA*transpose(m))*(m*transpose([UA1 UA2]));
R=QUA-QUA*transpose(m)*inv(m*QUA*transpose(m))*(m*QUA);

UA1;
error_UA1=sqrt(dUA1*Q*transpose(dUA1));
UA2;
error_UA2=sqrt(dUA2*Q*transpose(dUA2));

UA=UA1UA2(1);
error_UA=sqrt(R(1));

gA
error_gA=sqrt(dgA*Q*transpose(dgA));
