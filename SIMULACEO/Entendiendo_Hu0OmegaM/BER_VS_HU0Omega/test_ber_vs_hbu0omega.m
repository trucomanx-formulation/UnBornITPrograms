addpath('functions-m');

addpath(genpath('~/lib/octmat/pdsit'));

Ps=[0:40]/80;
L=length(Ps);



M=[1 2 4 8 16];
MSTEPS = length(M);
BER    = zeros(MSTEPS,L);
H      = zeros(MSTEPS,L);
BERofH = zeros(MSTEPS,L);

for II=1:MSTEPS
	II
	BER(II,:)        = pds_bersbceo(Ps,M(II));
	H(II,:)          = pds_shu0omega(Ps,M(II));
	BERofH(II,:)     = pds_hbinv(H(II,:) );
end 

hf1=figure(1);
hp1=semilogy(	Ps,BER(4,:)-BERofH(4,:),'-+',
				Ps,BER(5,:)-BERofH(5,:),'-+'
			); grid on
ha1=gca();
axis(ha1,[0 0.5 1.0e-6 1]);
hx1=xlabel('Ps');
hy1=ylabel('BER- BER of H(u0|Omega)');
hl1=legend(	['M=',num2str(M(4))], ...
			['M=',num2str(M(5))], ...
			4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_BER  = zeros(MSTEPS,L);
DELTAH = zeros(MSTEPS,L);
MX=2;

for II=1:MSTEPS
	H_BER(II,:)  = pds_hb(BER(II,:));
	DELTAH(II,:) = H_BER(II,:)-H(II,:);
end
hf2=figure(2);
hp2=plot(	Ps,DELTAH(1,:),'-o', ...
			Ps,DELTAH(2,:),'-p', ...
			Ps,DELTAH(3,:),'-s', ...
			Ps,DELTAH(4,:),'-v', ...
			Ps,DELTAH(5,:),'-<', ...
			Ps,funcion(Ps,MX),'->' ...
		); grid on
ha2=gca();
%axis(ha2,[0 0.5 1.0e-6 1]);
hx2=xlabel('Ps');
hy2=ylabel('H(BER) - H(u0|Omega)');
hl2=legend(	['M=',num2str(M(1))], ...
			['M=',num2str(M(2))], ...
			['M=',num2str(M(3))], ...
			['M=',num2str(M(4))], ...
			['M=',num2str(M(5))], ...
			['M=',num2str(MX)], ...
			3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FONTSIZE=20;
MZ=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(hp1,'markersize',MZ);
set(ha1,'fontsize',FONTSIZE,'GridLineStyle','--'); % sets font of numbers on axes
set(hx1,'fontsize',FONTSIZE);
set(hy1,'fontsize',FONTSIZE);
set(hl1,'fontsize',FONTSIZE);

print(hf1,'test_ber_vs_hbu0omega.png','-dpng',['-F:',num2str(FONTSIZE)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(hp2,'markersize',MZ);
set(ha2,'fontsize',FONTSIZE,'GridLineStyle','--'); % sets font of numbers on axes
set(hx2,'fontsize',FONTSIZE);
set(hy2,'fontsize',FONTSIZE);
set(hl2,'fontsize',FONTSIZE);

print(hf2,'test_ber_vs_hbu0omega2.png','-dpng',['-F:',num2str(FONTSIZE)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

