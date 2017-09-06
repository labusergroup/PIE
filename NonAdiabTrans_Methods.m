% working ver, 2017-1-11

% Non-adibatic transtion calculation comparing 3 methods, manually put
% population on thw lowest adiabatic state at t=0 and see it transits.
% (1) surface hoping, cf. Tully, J. Chem. Phys. 93, 1061 (1990); http://dx.doi.org/10.1063/1.459170 
% and http://doi.wiley.com/10.1002/wcms.1158
% (2) cf tannor (9.100)
% (3) solve TDSE, using Euler method
% (4) sovle TDSE, using ode, see another code PIE_2level.m

% ------------------------------------------------------------
%% atomic unit , 
clear

hbar = 1; 
auTime = 2.42e-17; %[s]
Hartree = 4.35974417e-18; % = auEnergy, [J]
auDipole = 8.47835326e-30; % C m
Epsilon = 1/4/pi; % Ft/m
auField = 5.1422065e+11; % V/m
NStates = 17;   %number of state

load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\BIM_MRCI_2015\DataPESandTDM.mat')

% Energy
for j=1:NStates-1
    w(j,j) = PESRFC(j,5); %first N ionic states in au.
end
% TDM (y direction) 16x16
TDM = zeros(NStates-1);
TDM = TDMyRFC;
mu = TDM - diag(diag(TDM));

% pulse parameters
EIR = 1*.0237; % 0.0237au ~ 20 TW/cm^2
tIR = 10e-15/1.6651/auTime; % Gaussian FWHM_int
hbar_ev = 6.58211928e-16; %eV s


% simulation
dt = 1e-2; %[au]
clear a b U H_ad Omega
tL = 0; %ini time
xt = tL:dt:round(tIR*4);
a = zeros(2,length(xt)); % adiabatic/dressed states' coeff, tully
b = zeros(2,length(xt)); % temp coeff
c = zeros(2,length(xt)); % adiabatic/dressed states' coeff, messiah
d = zeros(2,length(xt)); % temp coeff
Phase = 0;

ft = zeros(2,length(xt)); % a la tannor, adiab
f0 = zeros(2,length(xt)); % a la tannor, diab
fH = zeros(2,length(xt)); %Heisenberg
Ft = zeros(2,length(xt)); % with adiabatic approx
F0 = zeros(2,length(xt)); % with adiabatic approx

g0 = zeros(2,length(xt)); % using bare H basis
gt = zeros(2,length(xt)); % dressed

% putting poplation on lower dressed state at t=0
a(1,-tL/dt+1) = 1;
c(1,-tL/dt+1) = 1;
ft(1,-tL/dt+1) = 1;
fH(1,-tL/dt+1) = 1;
Ft(1,-tL/dt+1) = 1;
% a(2,-tL/dt-round(0.6486/(auTime*1e15))) = 1; % putting poplation on lower dressed state at t=-0.6486

% 2-level
w2 = zeros(2); w2(1,1) = w(1,1); w2(2,2) = w(5,5);
mu2 = zeros(2); mu2(1,2) = mu(1,4); mu2(2,1) = mu(4,1);

%set carrier frequency with possible detuning
% omegaIR = 1.59/hbar_ev*auTime; %[rad/au.]
% omegaIR = (+0.0 + 1.1923)/hbar_ev*auTime;
omegaIR = w2(2,2)-w2(1,1)+ .3/hbar_ev*auTime; 

% flags
flgTully = 1;
flgMessiah = 0;
flgTannor = 1;
flgSchr = 1;

for j = 1:length(xt)
    yE(j) = EIR*exp(-1*(xt(j).^2)/2/tIR^2).*cos(omegaIR*xt(j)); %field
    H_di(:,:,j) = w2-yE(j)*mu2; %diabatic Hamiltonian
    [U(:,:,j), H_ad(:,:,j)]=eig(H_di(:,:,j)); % [V D] = eig(A) => AV = VD
    if j > -tL/dt+1
%         going from j-1 to j
        Omega(j-1) = H_ad(2,2,j-1)-H_ad(1,1,j-1);
        if j>=3
            Phase = Phase + Omega(j-2)*dt;
        end
%--------------------------------------------------------------------------------------
if flgTully == 1
% b(t-1) temporarily stores the amplitude of a(t-1), from which a(t) is calculated  
% eq(8) in Tully, 
% CORRECT simultation
        b(1,j-1) = a(1,j-1) - 1i*Omega(j-1)*a(1,j-1)*dt...
            - a(1,j-1)*sum((U(:,:,j-1)*[1;0]).*(U(:,:,j)*[1;0]-U(:,:,j-1)*[1;0]))...
            - a(2,j-1)*sum((U(:,:,j-1)*[1;0]).*(U(:,:,j)*[0;1]-U(:,:,j-1)*[0;1]));
        b(2,j-1) = a(2,j-1) - 1i*0*a(2,j-1)*dt...
            - a(2,j-1)*sum((U(:,:,j-1)*[0;1]).*(U(:,:,j)*[0;1]-U(:,:,j-1)*[0;1]))...
            - a(1,j-1)*sum((U(:,:,j-1)*[0;1]).*(U(:,:,j)*[1;0]-U(:,:,j-1)*[1;0]));
% CORRECT transformation        
%         a(:,j) = b(:,j-1); 
        a(:,j) = U(:,:,j-1)*(U(:,:,j)\b(:,j-1)); 
        
% test Hermitian, const coupling
        %             b(1,j-1) = a(1,j-1) - 1i*0.03*a(1,j-1)*dt - 1i*a(2,j-1)*dt;
        %             b(2,j-1) = a(2,j-1) - 1i*0*a(2,j-1)*dt - 1i*a(1,j-1)*dt;
% with renormalization
        % a(1,j) = a(1,j)/sqrt(sum(abs(a(:,j)).^2));
        % a(2,j) = a(2,j)/sqrt(sum(abs(a(:,j)).^2));
        %             LtoU(j) = abs(sum((U(:,:,j-1)*[1;0]).*(U(:,:,j)*[0;1]-U(:,:,j-1)*[0;1])));
        %             UtoL(j) = abs(sum((U(:,:,j-1)*[0;1]).*(U(:,:,j)*[1;0]-U(:,:,j-1)*[1;0])));
        %             LtoU2(j) = abs(sum((U(:,:,j)*[0;1]-U(:,:,j-1)*[0;1]).^2));
        %             UtoL2(j) = abs(sum((U(:,:,j)*[1;0]-U(:,:,j-1)*[1;0]).^2));
        %             MeanLap = 0.5*(sum((U(:,:,j-1)*[1;0]).*(U(:,:,j)*[0;1]-U(:,:,j-1)*[0;1]))+sum((U(:,:,j-1)*[0;1]).*(U(:,:,j)*[1;0]-U(:,:,j-1)*[1;0])));
        %             a(1,j) = a(1,j-1)*exp(-1i*Omega(j)*dt) - a(2,j-1)*MeanLap; % eq(8) in Tully
        %             a(2,j) = a(2,j-1)*exp(-1i*0*dt) - a(1,j-1)*MeanLap;
end
%--------------------------------------------------------------------------------------
if flgMessiah == 1
% Messiah p755
% CORRECT simultation
        d(1,j-1) = c(1,j-1)...
            - c(1,j-1)*sum((U(:,:,j-1)*[0;1]).*(U(:,:,j)*[1;0]-U(:,:,j-1)*[1;0]))*exp(1i*Phase)*dt...
            + c(2,j-1)*sum((U(:,:,j-1)*[1;0]).*(U(:,:,j)*[0;1]-U(:,:,j-1)*[0;1]))*exp(-1i*Phase)*dt;
        d(2,j-1) = c(2,j-1)...
            + c(1,j-1)*sum((U(:,:,j-1)*[0;1]).*(U(:,:,j)*[1;0]-U(:,:,j-1)*[1;0]))*exp(1i*Phase)*dt...
            - c(2,j-1)*sum((U(:,:,j-1)*[1;0]).*(U(:,:,j)*[0;1]-U(:,:,j-1)*[0;1]))*exp(-1i*Phase)*dt;
            
% CORRECT transformation
        c(:,j) = d(:,j-1);
%         c(:,j) = U(:,:,j-1)*(U(:,:,j)\d(:,j-1));
end
%--------------------------------------------------------------------------------------
if flgTannor == 1
% a la Tannor p200
        ft(:,j) = ft(:,j-1) - 1i*H_ad(:,:,j-1)*ft(:,j-1)*dt - U(:,:,j-1)\(U(:,:,j)-U(:,:,j-1))*ft(:,j-1); % equ 9.100 tannor
        f0(:,j) = U(:,:,j)*ft(:,j);
        
        Ft(:,j) = Ft(:,j-1) - 1i*H_ad(:,:,j-1)*Ft(:,j-1)*dt; % with adiabatic approx
        F0(:,j) = U(:,:,j)*Ft(:,j);
        
%         fH(:,j) = fH(:,j-1) - 1i*H_ad(:,:,j-1)*fH(:,j-1)*dt; % adiabatic approx + transform of basis = Heisenberg
%         fH(:,j) = U(:,:,j-1)*(U(:,:,j)\fH(:,j)); 
end
%--------------------------------------------------------------------------------------
if flgSchr == 1
% using bare H basis, solve schrodinger eq
        if j == -tL/dt+2
            gt(:,j-1) = a(:,j-1);  % adibatic basis
            g0(:,j-1) = U(:,:,j-1)*gt(:,j-1); % diabatic basis
            f0(:,j-1) = U(:,:,j-1)*ft(:,j-1); % diabatic basis
            F0(:,j-1) = U(:,:,j-1)*Ft(:,j-1); % diabatic basis
        end
        g0(:,j) = g0(:,j-1) - 1i*H_di(:,:,j-1)*g0(:,j-1)*dt; % solve TDSE in diabatic basis
        gt(:,j) = U(:,:,j)\g0(:,j); % transform to adiabatic
end
%-------------------------------------------
    end
end

if flgTully == 1
    dfigure; hold on % surf hopping
    plot(auTime*1e15*xt,abs(a(1,:)).^2);
    plot(auTime*1e15*xt,abs(a(2,:)).^2);
    plot(auTime*1e15*xt,abs(a(1,:)).^2+abs(a(2,:)).^2);
    % ylim([0 1])
    title('surface hopping');
end
if flgMessiah == 1   
    dfigure; hold on % messiah
    plot(auTime*1e15*xt,abs(c(1,:)).^2);
    plot(auTime*1e15*xt,abs(c(2,:)).^2);
    plot(auTime*1e15*xt,abs(c(1,:)).^2+abs(c(2,:)).^2);
    % ylim([0 1])
    title('messiah');
end
if flgTannor == 1
    dfigure; hold on %tannor
    plot(auTime*1e15*xt,abs(ft(1,:)).^2);
    plot(auTime*1e15*xt,abs(ft(2,:)).^2);
    plot(auTime*1e15*xt,abs(ft(1,:)).^2+abs(ft(2,:)).^2); ylim([0 1])
    title('a la tannor')
    dfigure; hold on %tannor
    plot(auTime*1e15*xt,abs(f0(1,:)).^2);
    plot(auTime*1e15*xt,abs(f0(2,:)).^2);
    plot(auTime*1e15*xt,abs(f0(1,:)).^2+abs(f0(2,:)).^2); ylim([0 1])
    title('a la tannor, diab')
    dfigure; hold on %tannor
    plot(auTime*1e15*xt,abs(Ft(1,:)).^2);
    plot(auTime*1e15*xt,abs(Ft(2,:)).^2);
    plot(auTime*1e15*xt,abs(Ft(1,:)).^2+abs(Ft(2,:)).^2); ylim([0 1])
    title('a la tannor, adiabatic approx')
    dfigure; hold on %tannor
    plot(auTime*1e15*xt,abs(F0(1,:)).^2);
    plot(auTime*1e15*xt,abs(F0(2,:)).^2);
    plot(auTime*1e15*xt,abs(F0(1,:)).^2+abs(F0(2,:)).^2); ylim([0 1])
    title('a la tannor, adaibatic approx, diab')
    % dfigure; hold on %tannor
    % plot(auTime*1e15*xt,abs(fH(1,:)).^2);
    % plot(auTime*1e15*xt,abs(fH(2,:)).^2);
    % plot(auTime*1e15*xt,abs(fH(1,:)).^2+abs(fH(2,:)).^2); ylim([0 1])
    % title('a la tannor, Heisenberg')
end
if flgSchr == 1
    dfigure; hold on %Schrodinger, adiab
    plot(auTime*1e15*xt,abs(gt(1,:)).^2);
    plot(auTime*1e15*xt,abs(gt(2,:)).^2);
    plot(auTime*1e15*xt,abs(gt(1,:)).^2+abs(gt(2,:)).^2); ylim([0 1])
    title('Schroedinger Adiabatic');
    dfigure; hold on %Schrodinger, diab
    plot(auTime*1e15*xt,abs(g0(1,:)).^2);
    plot(auTime*1e15*xt,abs(g0(2,:)).^2);
    plot(auTime*1e15*xt,abs(g0(1,:)).^2+abs(g0(2,:)).^2); ylim([0 1])
    title('Schroedinger in Bare basis');
end
dfigure; hold on
plot(auTime*1e15*xt, yE);plot(auTime*1e15*xt, zeros(length(yE),1))
title('E Field');
dfigure; hold on
plot(auTime*1e15*xt, squeeze(H_ad(1,1,:))*Hartree/1.6e-19)
plot(auTime*1e15*xt, squeeze(H_ad(2,2,:))*Hartree/1.6e-19)
title('adib states energy')
