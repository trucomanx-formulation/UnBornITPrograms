addpath(genpath('/home/fernando/lib/octmat/pdsit'));

L=128;
Ps=[0:L]/L;

H2=pds_shomega(Ps,2);

H3=pds_shomega(Ps,3);

H4=pds_shomega(Ps,4);

H5=pds_shomega(Ps,5);

H10=pds_shomega(Ps,10);

H24=pds_shomega(Ps,24);


plot(Ps,H2/2,Ps,H3/3,Ps,H4/4,Ps,H5/5,Ps,H10/10,Ps,H24/24);grid

ylabel('Taxa media da entropia conjunta');
xlabel('Probabilidad \rho com a fonte escondida');

legend('H(2)/2','H(3)/3','H(4)/4','H(5)/5','H(10)/10','H(24)/24',2);
print('grafico.png','-dpng')
