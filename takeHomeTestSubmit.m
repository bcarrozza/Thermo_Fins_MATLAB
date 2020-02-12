% close, 
clc, clear, format shortg, format compact
tbl=@Table;
%Table Headers: T rho cp k alpha mu nu Pr
Tbl=readtable('TableA15.txt');
%to call 'k' @T=10degC from table, use:
% Tbl.k(Tbl.T==10)
kk=@(x) x+273; cc=@(x) x-273;
Tw= 70;  TwK=kk(Tw); %Twall:
Ta= 20;  TaK=kk(Ta); %Tambient:
H = 0.0254*2; W = 6*2.54/100; t = .003;
L = linspace(t,3,201);
dT = TwK-TaK; TfK = (TaK+TwK)/2;
g = 9.81; B = 1/TfK;

Tf=cc(TfK);
fRaL=@(L)(g*B*dT*L.^3)/(tbl(Tf,'nu')*tbl(Tf,'alpha'));
fQo=@(L) 0.59*tbl(Tf,'k')*W*dT*(fRaL(L)).^(1/4);
fQcmax=@(L)(1+0.36*(fRaL(L)).^(1/4).*H./L).*fQo(L);
ratioQs=@(L)(1+0.36*(fRaL(L)).^(1/4).*H./L);

% figure(3)
% subplot(2,1,1)
% plot(H./L,fQcmax(L),'bo',H./L,fQcmax(L),'k')
% title('Q_{c(max)} vs H/L'), grid
% xlabel('H/L'), ylabel('Q_{dot(max)} (W)')
% axis([0 max(H./L) 0 max(fQcmax(L))])

% subplot(2,1,2)
% plot(L,ratioQs(L),'bo',L,ratioQs(L),'k')
% title('Q_{c(max)}/Q_o vs L'), grid
% xlabel('L (m)'), ylabel('Q_{c(max)}/Q_o (W)')
% axis([0 max(L) 0 max(ratioQs(L))])

%%%%%%%%%
L=.0254*12;
%%%%%%%%%

sop=3.15*L.*(fRaL(L)).^(-1/4);

Qcmax=fQcmax(L);

nSpaces = round((W-t)/(sop+t));
t=(W - nSpaces*sop)/(nSpaces+1);
nFins = nSpaces+1;

x = linspace(-pi/100,2*pi,501);
x1=x*W/(2*pi);
duty=t/(sop+t)*100;
y = H/2.*(square((nFins-1+t/(sop+t))*x,duty)+1);
y1=(y>0);y2=(y1==0);
Z=L*y1-.01*y2;
[Z]=meshgrid(Z);
[x1,y]=meshgrid(x1,y);

Z(5:end-5,:)=Z(5:end-5,:).*(Z(5:end-5,:)~=L);
C=ones(size(Z));

figure(2)
s=surf(x1,y,Z,C);
title('Simple Rendering of Heat Sink (m)')
xlabel('W'),ylabel('H')
zlabel('L','Rotation',0)
axis ij
axis equal
axis([-pi/1e4 W -.0001 H+t 0 L])
RaL=fRaL(L);
% t=round(t,4);
sop=round(sop,4);

table(Tw,Ta,Tf,L,H,W,t,RaL,sop,nFins)

function yi = Table(knownTemp,desiredProperty)
Tbl=readtable('TableA15.txt');
x=Tbl.T;
Y=desiredProperty;
xq=knownTemp;
switch Y
    case 'rho'
        Y=Tbl.rho;
    case 'cp'
        Y=Tbl.cp;
    case 'k'
        Y=Tbl.k;
    case 'alpha'
        Y=Tbl.alpha;
    case 'mu'
        Y=Tbl.mu;
    case 'nu'
        Y=Tbl.nu;
    case 'Pr'
        Y=Tbl.Pr;
    otherwise
        fprintf('\n     The header option does not exist\n\n')
end
yi = interp1(x,Y,xq);
end

