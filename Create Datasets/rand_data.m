clear all;
close all;
clc;

%% Variables declaration
Fc = 865*10^6;              %carrier frequency
c = 299792458;
lambda = c/Fc;
TagRate = 40*10^3; %tag rate/ backscatter link frequency (BLF) ???
Tpri = 1/(2*TagRate); 

Fs = 10*10^6; %sampling rate
Ts = 1/Fs;

Nant = 3;   %number of antennas
d = lambda/2;
X_Tx = [0.465  0];
X_Rx(1,:) = [d/2  0];
X_Rx(2,:) = [-d/2  0];
X_Rx(3,:) = [-3*d/2  0];
Rx_midp = sum(X_Rx)/Nant;
d0 = 1;     %reference distance
d_cr = zeros(1,Nant);   d_tr = zeros(1,Nant);

for i=1:Nant
    d_cr(i) = norm(X_Tx - X_Rx(i,:));
end

v_ct = 2;
v_cr = 2;
v_tr = 2;
k_ct = 20;
k_cr = 20;
k_tr = 15;
sigma_CT2 = 1;
sigma_CR2 = 1;
sigma_TR2 = 1;
G0 = 0;
G1 = 1;
As = 0.6047 + 1j*0.5042;
s = 0.6;    %tag's scattering efficiency (0.1 realistic)

L = Fs/TagRate; %oversampling factor

m2_preamble = [1 0 1 0 ... %miller 2 preamble
               1 0 1 0 ...
               1 0 1 0 ...
               1 0 1 0 ...
               1 0 1 0 ...
               1 0 0 1 ...
               0 1 0 1 ...
               0 1 1 0 ...
               1 0 0 1 ...
               0 1 1 0].'; 

Nb = 128; %payload
N = Nb + length(m2_preamble)/4; %total bits

Nt1 = round(10*Tpri/Ts); %Nominal cw duration
Ntotal = round(Nt1 + L*N*4); %total duration

t = (0:Ts:Ntotal*Ts - Ts)'; %time vector
Ttotal = Ntotal*Ts;

Pc_dBm = 30;                % power of carrier wave in dBm
Pc = 10.^(Pc_dBm/10)/1000;  % 30 dBm = 1 Watt
count = 0;
iter = 10^6;
data = zeros(iter, 6);
for o=1:iter    %parfor
    tic
    X_tag = [(-4+8*rand) 4*rand];    
    d_ct = norm(X_Tx - X_tag);
    d_tr = zeros(1,3);
    for i=1:Nant
        d_tr(i) = norm(X_tag - X_Rx(i,:)); 
    end
    L_ct = (lambda/(4*pi*d0))^2*(d0./d_ct).^v_ct;
    L_cr = zeros(1,Nant);   L_tr = zeros(1,Nant);
    for i=1:Nant
        L_cr(i) = (lambda/(4*pi*d0))^2*(d0./d_cr(i)).^v_cr;
        L_tr(i) = (lambda/(4*pi*d0))^2*(d0./d_tr(i)).^v_tr;
    end

    %Rice channel fading assumption 
    h_ct = (sqrt(sigma_CT2*k_ct/(k_ct + 1)) + sqrt(sigma_CT2/(2*(k_ct + 1)))*(randn(1,1) + 1j*randn(1,1)));
    %induced phase of the propagation channel
    h_ct = abs(h_ct)*exp(1i*2*pi*d_ct/lambda);
    h_cr = zeros(1,Nant);   h_tr = zeros(1,Nant);   h = zeros(1,Nant);
    phase_noise = deg2rad(25)*randn(Nant,1);
    for i=1:Nant
        h_cr(i) = (sqrt(sigma_CR2*k_cr/(k_cr + 1)) + sqrt(sigma_CR2/(2*(k_cr + 1)))*(randn(1,1) + 1j*randn(1,1)));
        h_cr(i) = abs(h_cr(i))*exp(1i*2*pi*d_cr(i)/lambda);
        h_tr(i) = (sqrt(sigma_TR2*k_tr/(k_tr + 1)) + sqrt(sigma_TR2/(2*(k_tr + 1)))*(randn(1,1) + 1j*randn(1,1)));
        h_tr(i) = abs(h_tr(i))*exp(1i*2*pi*d_tr(i)/lambda);
        h(i) = h_ct*h_tr(i)*exp(1i*phase_noise(i));
    end

    %% Tag's signal construction 
    x = round((sign(randn(Nb, 1))+1)/2); %bits
    x0h = [1 0 1 0].';
    x0l = [0 1 0 1].';
    x1h = [1 0 0 1].';
    x1l = [0 1 1 0].';
    s1 = [+1 -1 -1 +1].';
    s3 = [+1 -1 +1 -1].';
    c = zeros(4*Nb, 1); %miller2 bits

    if x(1) == 1
        c(1:4) = x1h;
        prevsymbol = 1;
        prevstate = 1;
    else
        c(1:4) = x0h;
        prevsymbol = 0;
        prevstate = 1;
    end

    for i=2:Nb
        if(x(i) == 0)
            if(prevstate == 1)
                c((i-1)*4 + 1:i*4) = x0l;
                prevstate = 0;
            else
                c((i-1)*4 + 1:i*4) = x0h;
                prevstate = 1;
            end
            prevsymbol = 0;
        else
            if(prevsymbol == 1)
                if(prevstate == 1)
                    c((i-1)*4 + 1:i*4) = x1l;
                    prevstate = 0;
                else
                    c((i-1)*4 + 1:i*4) = x1h;
                    prevstate = 1;
                end
            else
                if(prevstate == 1)
                    c((i-1)*4 + 1:i*4) = x1h;
                    prevstate = 1;
                else
                    c((i-1)*4 + 1:i*4) = x1l;
                    prevstate = 0;
                end
            end
            prevsymbol = 1;
        end
    end

    c = [m2_preamble; c]; %appended preamble

    u_c = zeros(1, 4*L*N); %upsampled

    for i=1:4*N
        u_c((1 + (i-1)*L):(i*L)) = c(i)*ones(1,L);
    end

    %% Noise Signal
    kb = 1.380649e-23;            % Boltzmann's constant
    Tth = 300;                    % Reader temperature @ 300K=27C\
    N0 = kb*Tth;
    W = 2e6;
    sigma_n2 = N0*W/2;            % 2sigma_n^2 = N0*W
    sigma2 = 2*sigma_n2*Ts;
    e_bit = Pc*s*L_ct.*L_tr*Tpri/2;
    SNR = e_bit/(2*sigma2);
    SNR_dB = 10*log10(SNR);
    n = zeros(Nant,Ntotal);
    for i=1:Nant
        n(i,:) = sqrt(sigma_n2)*(randn(1,Ntotal) + 1j*randn(1,Ntotal));
    end
    e_bit = Pc*s*L_ct.*L_tr(1)*Tpri/2;
    SNR = e_bit/(2*sigma2);
    SNR_dB = 10*log10(SNR);

    %% Signal Model
    m_dc = zeros(1,Nant);
    m_tag = zeros(1,Nant);
    for i=1:Nant
        m_dc(i) = (sqrt(L_cr(i))*h_cr(i) + (As-G0)*sqrt(s*L_ct*L_tr(i))*h(i));
        m_tag(i) = (G0-G1)*sqrt(s*L_ct*L_tr(i))*h(i);
    end

    ys = zeros(Nant,Ntotal);
    for i=1:Nant
            ys(i,:) = sqrt(2*Pc).*(m_dc(i)*ones(1, Ntotal) + [zeros(1, Nt1) m_tag(i)*u_c]) + n(i,:);
        %^^ recieved samples vector
    end

    %% Estimations
    theta0_est = zeros(1,Nant); dc_est = zeros(1,Nant);
    Y = zeros(1,length(ys)-Nt1);
    for i=1:Nant
        theta0_est(i) = mod(mean(angle(ys(1:Nt1))),2*pi);
        %dc_est = mean(Y_step(1:Nt1).*exp(-1j*angle(Y_step(1:Nt1))));
        dc_est(i) = abs(mean(ys(i,1:Nt1)));
        Y = ys(i,Nt1+1:end).*exp(-1j*(theta0_est(i))) - dc_est(i); %for |dc_est|
    end

    %% Divide signal during preamble bits
    %Preamble bits have equal number of zeros and ones
    y_preamble_0 = zeros(Nant, L*length(m2_preamble)/2);
    y_0_count = 1;
    y_preamble_1 = zeros(Nant, L*length(m2_preamble)/2);
    y_1_count = 1;
    for i=1:length(m2_preamble)
        if m2_preamble(i) == 0
            y_preamble_0(:, (1 + (y_0_count-1)*L):(y_0_count*L) ) = ys(:, (1+Nt1+(i-1)*L):(Nt1+(i*L)) );
            y_0_count = y_0_count + 1;
        else
            y_preamble_1(:, (1 + (y_1_count-1)*L):(y_1_count*L) ) = ys(:, (1+Nt1+(i-1)*L):(Nt1+(i*L)) );
            y_1_count = y_1_count + 1;
        end
    end

    y_diff = y_preamble_0 - y_preamble_1;

    %% Phase at the Rx antenna
    %Calculate phase on the i-th Rx antenna
    phi_meas = zeros(Nant,1);
    for i=1:Nant
        phi_meas(i) = mod(atan2d(mean(imag(y_diff(i,:))),mean(real(y_diff(i,:)))),360);
        % if phi_meas(i)<0    %convert phase measurements to [0,360] from [-180,180]
        %     phi_meas(i) = phi_meas(i)+360;
        % end
    end
    
    %% Phase Difference
    D_phi_12 = mod(phi_meas(1) - phi_meas(2) ,360);
    D_phi_23 = mod(phi_meas(2) - phi_meas(3) ,360);
    
    %% DoA
    doa = zeros(Nant,1);
    for i=1:Nant
        %create pairs (1-2, 2-3, 1-3)
        if i<Nant
            j=i; k=i+1;
        else
            j=1; k=3;
        end
        if d_tr(j)>d_tr(k)
            if phi_meas(j)>phi_meas(k)
                doa(i) = asind( lambda*(phi_meas(j)-phi_meas(k))/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
            else
                doa(i) = asind( lambda*(phi_meas(j)-phi_meas(k)+360)/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
            end
        elseif d_tr(j)<d_tr(k)
            if phi_meas(j)>phi_meas(k)
                doa(i) = asind( -lambda*(phi_meas(j)-phi_meas(k)-360)/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
            else
                doa(i) = asind( -lambda*(phi_meas(j)-phi_meas(k))/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
            end
        %if equal take the non complex solution
        else
            if phi_meas(j)>phi_meas(k)
                doa(i) = asind( lambda*(phi_meas(j)-phi_meas(k))/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
                if ~isreal(doa(i))
                    doa(i) = asind( -lambda*(phi_meas(j)-phi_meas(k)-360)/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
                end
            else
                doa(i) = asind( lambda*(phi_meas(j)-phi_meas(k)+360)/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
                if ~isreal(doa(i))
                    doa(i) = asind( -lambda*(phi_meas(j)-phi_meas(k))/(360*norm(X_Rx(j,:)-X_Rx(k,:))) );
                end
            end  
        end
    end

    data(o,:) = [X_tag(1) X_tag(2) doa(1) doa(2) D_phi_12 D_phi_23];
    toc
end

save rand_data_iid_noise data
