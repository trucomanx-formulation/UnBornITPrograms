addpath(genpath('/home/fernando/lib/octmat/pdsit'));

X=load('data10.txt');
N=10
ITER=length(X(:,1))
Ps=X(:,1)';
Perrfonte=zeros(1,ITER);

Perrfonte=pds_bersbceo(Ps,N);

semilogy(Ps,X(:,2),'-*',Ps,Perrfonte,'-.'), grid on
axis([0 0.5 1.0e-6 1]);
xlabel('Ps');
ylabel('Perro');
legend('Perro','Teorico Perro',2);


