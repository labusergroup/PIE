% - See README.md for general desciptions.

% - Modified on 2017-9-5, 2 states from BIM, SI unit
% - This code simulates non-adiabatic transition among ionic states

function [ttt, w, mu, t, a,  Pop,Pop_c, omegaIR, U, H_ad, yE] = PIE_2level(varargin)

%%
%Define parameters
tic

global hbar hbar_ev auDipole w mu EX EIR omegaX omegaIR tIR tX Power Basis NStates a0

Power = 12; %power factor used to select unit, in this case, "Power=12" gives unit of ps

hbar = 1.054571726e-34*10^Power; %J ps
hbar_ev = 6.58211928e-16*10^Power; %eV ps
Hartree = 4.35974417e-18/1.6021765e-19; % = auEnergy, eV
auDipole = 8.47835326e-30; % C m
Epsilon = 8.8541878176e-12; % F/m
auField = 5.1422065e11; % V/m
NStates = 3;   %number of state
auTime = 2.42e-17*10^Power; %[ps]

gap = 2; %eV
%pulse parameters 
EX = 0*0.00217*auField; %V/m, peak int = EX^2*Epsilon*3e8/2 [W/m^2]
if nargin == 0
    EIR = 1*.0237*auField; % 0.0237au ~ 20 TW/cm^2
    omegaIR = (-0.7 + gap)/hbar_ev; % [rad/ps]
% omegaIR = (+0.0 + 1.1922)/hbar_ev;
tIR = 10e-3/1.6651; % Gaussian FWHM_int
else
    EIR = varargin{1};
    omegaIR = varargin{2};
    tIR = varargin{3};
end



% omegaX = 10.37/hbar_ev;
% omegaX = (10.37+0.7062)/hbar_ev;
% sqrt(9.69/2/1.6e-19/240e10^2*4*9.11e-31*(3e8/780e-9*2*pi)^2)
% IR_MinI =10e16;     %W/m^2, 0.0169 auField
% IR_MaxI = 100e16;    %W/m^2, 0.0534 auField
% TickY_factor = 3;
% Ymax=3*TickY_factor + 1;
% omegaIR = 0.05843*Hartree/hbar_ev; %radian/ps, 3e8/(omegaIR/2/pi*1e12) [m]


% 'tIR' and 'tX' are characteristic times for IR and VUV pulses,depending on the pulse shape. E.g.,
% exp(-1*(t.^2)/2/tX^2) --> FWHM_field=2*sqrt(2*log(2))*tX=2.35482*tX, FWHM_int=2sqrt(ln2)tX=1.6651*tX [ps];
% (cos(pi*t/tX))^2.*(t>-tX/2 & t<tX/2) -->FWHM_field=0.5*tX,FWHM_int=2*acos(0.5^0.25)/pi*tX=0.364*tX [ps]
% tIR = 30e-3/2*pi/acos((.5)^0.25); %ps, cutoff time for cos^2 waveform
% tX = 20e-3/2*pi/acos((.5)^0.25); %20e-3/2*pi/acos((.5)^0.25); %ps in cos^2(pi*t/tX) -->FWHM_int=20 fs

tX = tIR/7^.5;
% tX = 0.1e-3/1.6651;

Basis = 'a';
% Basis = 'd';

%%
% Population transfer calculation
% adiabatic basis with scan between RFC and MIN
if Basis == 'a' % adiabatic basis
    display 'adiabatic'
    clear PopE
    load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\BIM_MRCI_2015\DataPESandTDM.mat')
    
    w = zeros(NStates);
    w(1,1) = -9.69/hbar_ev;
    w(2,2) = 0;
%     w(3,3) = PESRFC(5,5)*Hartree/hbar_ev; %[radian/ps]
    w(3,3) = gap/hbar_ev; %[radian/ps]
    
%     omegaIR = (w(3,3)-w(2,2));
    
    TDM = zeros(17);
    TDM(2:end,2:end) = TDMyRFC; % + (TDMyMIN - TDMyRFC)/Np*p;
    
    mu = zeros(NStates);
    mu(2,3) = TDM(2,5)*auDipole;
    mu(3,2) = TDM(5,2)*auDipole;
    
    omegaX = w(2,2) - w(1,1); %irrelevant for Ex=0
    
    H_di0 =  w*hbar-EIR*mu; %[J]
    [U0, H_ad0]=eig(H_di0);
    c0 = zeros(NStates,1); a0 = zeros(NStates,1);
    c0(2) = 1;
    a0 = U0*c0;
%     a0 = [0;1;0];
    [t,a] = SolveTDSE_IwPE(); %solve TDSE
    clear ProbAmp Pop;
    aTilde = zeros(size(a)); %to get the amplitude in the schroedinger picture (from interaction pic)
    for j=1:size(a,2);
        aTilde(:,j) = a(:,j).*exp(-1i*w(j,j)*t);
    end
    %             ProbAmp = a';
    Pop = abs(aTilde).^2;  %population
    Pop(:,end+1) = t;
    
    yE = EIR*exp(-1*(t.^2)/2/tIR^2).*cos(omegaIR*t);
    for j =1:length(yE)
        H_di(:,:,j) = w*hbar-yE(j)*mu; %[J]
        [U(:,:,j), H_ad(:,:,j)]=eig(H_di(:,:,j)); % [V D] = eig(A) => AV = VD\
        c(:,j) = U(:,:,j)\aTilde(j,:).'; % transforming to adiabitic in Schroedinger pic, note to use .' instead of '
    end
    Pop_c = (abs(c).^2).';
    
% %     plots
if nargin == 0
%     dfigure; hold on
%     plot(t,EIR*exp(-1*(t.^2)/2/tIR^2).*cos(omegaIR*t));
%     plot(t,yE);
%     title('E field')
%     dfigure; hold on
%     plot(t,squeeze(H_ad(1,1,:))/1.6e-19);
%     plot(t,squeeze(H_ad(2,2,:))/1.6e-19);
%     plot(t,squeeze(H_ad(3,3,:))/1.6e-19);
%     title('adiabatic energy')
%     dfigure; hold on
%     plot(t,Pop(:,2)); plot(t,Pop(:,3));  plot(t,Pop(:,2)+Pop(:,3));
%     title('diabatic states population')
    dfigure; hold on
     plot(t,Normalize(yE-min(yE)));
% plot(t,yE>0); ylim([-0.05 1])
    plot(t,Pop_c(:,2)); plot(t,Pop_c(:,3));plot(t,Pop_c(:,2)+Pop_c(:,3));
    title('adiabatic states population')
    plot(t,Pop(:,2)); plot(t,Pop(:,3));  plot(t,Pop(:,2)+Pop(:,3));
end
ttt=1-Pop_c(end,2);
end
%%
% diabatic basis
if Basis == 'd' % Energies and SOC in diabatic basis, directly from Tamas' calculation
    display 'diabatic, do nothing'
end

toc
end

function [t,a] = SolveTDSE_IwPE()
global  tIR tX NStates w mu a0
options = odeset('RelTol',1e-6,'AbsTol',1e-7); % ode solving accuracy
% a0 = zeros(NStates,1)';
% a0(1) = 1;
% a0(2) = 1/sqrt(2);
% a0(3) = 1/sqrt(2);
% [t,a] = ode23(@TDSE_IwPE,[0 0.025],a0,options); % set solver, time interval, initial condition.
[t,a] = ode23(@TDSE_IwPE,[0 tIR*4],a0,options); % set solver, time interval, initial condition.
end

function da = TDSE_IwPE(t, a)
global hbar w mu EX EIR omegaX omegaIR tIR tX Basis NStates
% t in ps
% da = zeros(NStates,1);

MU = zeros(NStates);
for j=1:NStates
    %     approx, cos^2
%     MU(1,j) = 1/2/hbar*mu(1,j)*(EX*cos(pi*t/tX).^2.*(t>-tX/2 & t<tX/2)*exp(+1i*omegaX*t)+EIR*cos(pi*t/tIR).^2.*(t>-tIR/2 & t<tIR/2)*2*cos(omegaIR*t))*exp(-1i*(w(j,j)-w(1,1))*t);
%     MU(j,1) = 1/2/hbar*mu(j,1)*(EX*cos(pi*t/tX).^2.*(t>-tX/2 & t<tX/2)*exp(-1i*omegaX*t)+EIR*cos(pi*t/tIR).^2.*(t>-tIR/2 & t<tIR/2)*2*cos(omegaIR*t))*exp(-1i*(w(1,1)-w(j,j))*t);
%     approx, Gaussian
    MU(1,j) = 1/2/hbar*mu(1,j)*(EX*exp(-1*(t.^2)/2/tX^2)*exp(+1i*omegaX*t)+EIR*exp(-1*(t.^2)/2/tIR^2)*2*cos(omegaIR*t))*exp(-1i*(w(j,j)-w(1,1))*t);
    MU(j,1) = 1/2/hbar*mu(j,1)*(EX*exp(-1*(t.^2)/2/tX^2)*exp(-1i*omegaX*t)+EIR*exp(-1*(t.^2)/2/tIR^2)*2*cos(omegaIR*t))*exp(-1i*(w(1,1)-w(j,j))*t);
    %     approx, Gaussian, picket UV
%     MU(1,j) = 1/2/hbar*mu(1,j)*(EX*exp(-1*(t.^2)/2/tIR^2)*exp(+1i*omegaX*t)*(cos(omegaIR*t-pi/2)>0.98)+EIR*exp(-1*(t.^2)/2/tIR^2)*2*cos(omegaIR*t))*exp(-1i*(w(j,j)-w(1,1))*t);
%     MU(j,1) = 1/2/hbar*mu(j,1)*(EX*exp(-1*(t.^2)/2/tIR^2)*exp(-1i*omegaX*t)*(cos(omegaIR*t-pi/2)>0.98)+EIR*exp(-1*(t.^2)/2/tIR^2)*2*cos(omegaIR*t))*exp(-1i*(w(1,1)-w(j,j))*t);
%      MU(1,j) = 1/2/hbar*mu(1,j)*(EX*exp(-1*((t-.000651).^2)/2/tX^2)*exp(+1i*omegaX*t)+EIR*exp(-1*(t.^2)/2/tIR^2)*2*cos(omegaIR*t))*exp(-1i*(w(j,j)-w(1,1))*t);
%     MU(j,1) = 1/2/hbar*mu(j,1)*(EX*exp(-1*((t-.000651).^2)/2/tX^2)*exp(-1i*omegaX*t)+EIR*exp(-1*(t.^2)/2/tIR^2)*2*cos(omegaIR*t))*exp(-1i*(w(1,1)-w(j,j))*t);
   
    %     no approx, cos^2
%     MU(1,j) = 1/2/hbar*mu(1,j)*(EX*cos(pi*t/tX).^2.*(t>-tX/2 & t<tX/2)*2*cos(omegaX*t)+EIR*cos(pi*t/tIR).^2.*(t>-tIR/2 & t<tIR/2)*2*cos(omegaIR*t))*exp(-1i*(w(j,j)-w(1,1))*t);
%     MU(j,1) = 1/2/hbar*mu(j,1)*(EX*cos(pi*t/tX).^2.*(t>-tX/2 & t<tX/2)*2*cos(omegaX*t)+EIR*cos(pi*t/tIR).^2.*(t>-tIR/2 & t<tIR/2)*2*cos(omegaIR*t))*exp(-1i*(w(1,1)-w(j,j))*t);
end
if Basis == 'd';
    for k=2:NStates
        for j =2:NStates
            if j~=k
%                 MU(k,j) = 1/hbar*mu(k,j)*EIR*cos(pi*t/tIR).^2.*(t>-tIR/2 & t<tIR/2)*cos(omegaIR*t)*exp(-1i*(w(j,j)-w(k,k))*t) - SO(k,j)*exp(-1i*(w(j,j)-w(k,k))*t);
                %cos^2
%                 MU(k,j) = 1/hbar*mu(k,j)*EIR*cos(pi*t/tIR).^2.*(t>-tIR/2 & t<tIR/2)*cos(omegaIR*t)*exp(-1i*(w(j,j)-w(k,k))*t) - w(k,j)*exp(-1i*(w(j,j)-w(k,k))*t);
                % Gaussian
                MU(k,j) = 1/hbar*mu(k,j)*EIR*exp(-1*(t.^2)/2/tIR^2)*cos(omegaIR*t)*exp(-1i*(w(j,j)-w(k,k))*t) - w(k,j)*exp(-1i*(w(j,j)-w(k,k))*t);
            end
        end
    end
elseif Basis == 'a'     % SO adiabatic basis, no SOC H, approx, matrix form
    for k=2:NStates
        for j =2:NStates
            if j~=k
%                 MU(k,j) = 1/hbar*mu(k,j)*EIR*cos(pi*t/tIR).^2.*(t>-tIR/2 & t<tIR/2)*cos(omegaIR*t)*exp(-1i*(w(j,j)-w(k,k))*t);
                MU(k,j) = 1/hbar*mu(k,j)*EIR*exp(-1*(t.^2)/2/tIR^2)*cos(omegaIR*t)*exp(-1i*(w(j,j)-w(k,k))*t);               
            end
        end
    end
end
da = 1i*MU*a;

end

