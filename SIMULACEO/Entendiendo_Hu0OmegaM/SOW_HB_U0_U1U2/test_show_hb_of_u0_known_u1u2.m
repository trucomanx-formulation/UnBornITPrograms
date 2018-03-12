addpath('functions-m');
addpath(genpath('~/lib/octmat/pdsit'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% x axes : p1
%%				H(U1|U0)=Hb(p1)
%% y axes : p2
%%				H(U2|U0)=Hb(p2)
%%
%% show 3D graph of:	f(p1,p2)=Hb(U0|U1U2)
%%
%%						Hb(U0|U1U2)=Hb(p1) + Hb(p2) - Hb(p1||p2)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=40;
p1 = [1:(M-1)]/(M);
p2 = [1:(M-1)]/(M);
[P1,P2] = meshgrid(p1,p2);

C=pds_hb(P1)+pds_hb(P2)-pds_hb(P1+P2-2*P1.*P2);

figure(1)
contour(P1,P2,C),grid;
title('H(u0|u1u2)=pds_hb(p1)+pds_hb(p1)-pds_hb(p1||p2)');
print('test_show_hb_of_u0_known_u1u2_1.png','-dpng')

figure(2)
surf(P1,P2,C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% p1=alpha1+0.5
% p2=alpha2+0.5
% h(alpha+0.5)=func(alpha)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=40;
alpha1 = [-(M-1):(M-1)]/(2*M);
alpha2 = [-(M-1):(M-1)]/(2*M);
[ALPHA1,ALPHA2] = meshgrid(alpha1,alpha2);

C=func(ALPHA1)+func(ALPHA2)-func(-2*ALPHA1.*ALPHA2);

%figure(3)
%contour(ALPHA1,ALPHA2,C)

%figure(4)
%surf(ALPHA1,ALPHA2,C)

FAC=sqrt(2);
ALPHA = [-(M-1):(M-1)]/(2*M);
BETA  = [-(M-1):(M-1)]*FAC/(2*M);
D=2*func(BETA/FAC)-func(-2*(BETA/FAC).^2);

figure(5)
plot(BETA,D,'-s',ALPHA,func(ALPHA),'-o')
legend('otro','hb(alpha+0.5)')

print('test_show_hb_of_u0_known_u1u2_2.png','-dpng')
