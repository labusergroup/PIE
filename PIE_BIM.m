% - See README.md for general desciptions.

% - This code simulates the AC Stark effect ionic states
% - Modified on 2017-9-3, 11 states for BIM.
% - ONLY in diabatic basis.
% - SI units
% - Two pulses are used: Low freq and high int "IR" induces Stark shift. High freq and low int "X" induces ionization.

function [w, mu, t, a, ProbAmp, Pop] = PIE_BIM()

%% Define parameters
% tic
global hbar hbar_ev auDipole w mu EX EIR omegaX omegaIR tIR tX Power Basis NStates

% constants
Power = 12; %power factor used to select unit, in this case, "Power=12" gives unit of ps

hbar = 1.054571726e-34*10^Power; %J ps
hbar_ev = 6.58211928e-16*10^Power; %eV ps
Hartree = 4.35974417e-18/1.6021765e-19; % = auEnergy, eV
auDipole = 8.47835326e-30; % C m
Epsilon = 8.8541878176e-12; % F/m
auField = 5.1422065e11; % V/m
NStates = 11;   %number of state

% pulse parameters 
EX = 0.00217*auField; %V/m, peak int = EX^2*Epsilon*3e8/2 [W/m^2]
EIR = 0*0.0237*auField; % 0.0237au ~ 20 TW/cm^2
omegaIR = 0.05843*Hartree/hbar_ev; %radian/ps, 3e8/(omegaIR/2/pi*1e12) [m]

% 'tIR' and 'tX' are characteristic times for IR and VUV pulses,depending on the pulse shape. E.g.,
% exp(-1*(t.^2)/2/tX^2) --> FWHM_field=2*sqrt(2*log(2))*tX=2.35482*tX, FWHM_int=2sqrt(ln2)tX=1.6651*tX [ps];
% (cos(pi*t/tX))^2.*(t>-tX/2 & t<tX/2) -->FWHM_field=0.5*tX,FWHM_int=2*acos(0.5^0.25)/pi*tX=0.364*tX [ps]
tIR = 5e-3/2*pi/acos((.5)^0.25); %ps, cutoff time for cos^2 waveform
tX = tIR/7^.5;

% Choose between adiabatic and diabatic bases. 
% Basis = 'a';
Basis = 'd';

%% Main simulation

%% adiabatic basis
if Basis == 'a' % adiabatic basis
    display 'adiabatic basis'
%     energy levels, Omega w, use 72th row in PES.dat (S0 min)
    PES = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\PES.dat');
    w_di = zeros(11);
    for j = 1:11
        w_di(j,j) = PES(72,j+1)*Hartree/hbar_ev; %radian/ps
    end
%     spin-orbital coupling
    SO = zeros(11);
    SOreal = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\S-O.dat');
    SOimagine = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\S-O_cmplx.dat');
    for k = 1:11
        for j = 1:11
            SO(k,j) = Hartree/hbar_ev*(SOreal(781+k,j)+1i*SOimagine(781+k,j)); %radian/ps
        end
    end
    w_SO = w_di+SO;
    [v, w_ad]=eig(w_SO);
    u = inv(v); % use either u*a or v\a to transform the final state vector 'a' to adiabatic basis.
    w = w_ad; % energies in the adiabatic basis. [rad/ps]
    
%     TDM from dipole2.dat, 11*72=>781
    dipole2 = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\dipole2.dat');
    mu_di = zeros(11);
    for k = 1:11
        for j = 1:11
            mu_di(k,j) = dipole2(781+k,j); %au
        end
    end
    mu = v\mu_di*v*auDipole;
    clear PopE
    omegaX = w_ad(2,2) - w_ad(1,1); % set to be on resonance
    [t,a] = SolveTDSE_PIE(); %solve TDSE
    clear ProbAmp Pop;
    ProbAmp = a';
    Pop = abs(ProbAmp').^2;  %population
    Pop(:,end+1) = t;
end

%% diabatic basis
if Basis == 'd' % Energies and SOC in diabatic basis, directly from Tamas' calculation
    display 'diabatic basis'
        
    %energy levels, Omega w, use 72th row in PES.dat
    PES = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\PES.dat');
    w_di = zeros(11);
    for j = 1:11
        w_di(j,j) = PES(72,j+1)*Hartree/hbar_ev; %radian/ps
    end

    %spin-orbital coupling
    SO = zeros(11);
    SOreal = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\S-O.dat');
    SOimagine = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\S-O_cmplx.dat');
    for k = 1:11
        for j = 1:11
            SO(k,j) = Hartree/hbar_ev*(SOreal(781+k,j)+1i*SOimagine(781+k,j)); %radian/ps
        end
    end 
    w_di = w_di+SO;   
    [v, w_ad]=eig(w_di);
    u = inv(v); % use either 'u*a' or 'v\a' to transform the final state vector 'a' to adiabatic basis.
    w = w_di; % energies in the adiabatic basis. [rad/ps]
    
    %TDM from dipole2.dat, 11*71=>781
    dipole2 = load('D:\group2\code\Halomethane_Stark_Shift_Codes\BIM data\dipole2.dat');
    mu_di = zeros(11);
    for k = 1:11
        for j = 1:11
            mu_di(k,j) = dipole2(781+k,j); %au
        end
    end
    mu = mu_di*auDipole;
    omegaX = w_ad(2,2) - w_ad(1,1); % set to be on resonance
    [t,a] = SolveTDSE_PIE(); %solve TDSE
    
    for j = 1:NStates
        a(:,j) = a(:,j).*exp(-1i*w_di(j,j)*t);
    end
    clear ProbAmp Pop;
%     IMPORTANT: when making the transformation, use .' instead of '
    ProbAmp = v\a.'; % to get population in the adiabatic basis
    Pop = abs(ProbAmp.').^2;  %population
    Pop(:,end+1) = t;
end

% toc
end

function [t,a] = SolveTDSE_PIE()
global  tIR tX NStates w mu
options = odeset('RelTol',1e-6,'AbsTol',1e-7); % ode solving accuracy
a0 = zeros(NStates,1)';
a0(1) = 1;
% a0(2) = 1/sqrt(2);
% a0(3) = 1/sqrt(2);
[t,a] = ode23(@TDSE_PIE,[-1.5*max(tX,tIR) 1.5*max(tX,tIR)],a0,options); % set solver, time interval, initial condition.
end

function da = TDSE_PIE(t, a)
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
