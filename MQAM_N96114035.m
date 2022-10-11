clc ; clear all ; close all;
N = 1000; %iterative number
d = 1; % distance
M = 64; %M-ary QAM % 4 16 64

k = log2(M);
snr_in_dB = 0 : 1 : 20;
snr = 10 .^ (snr_in_dB ./ 10); % signal-to-noise ratio
%BER_theo = (3/2)*erfc(sqrt(0.1*(10.^(snr_in_dB/10))));
%x = sqrt(3*k*snr/(M-1));
%BER_theo = (4/k)*(1- (1/sqrt(M)))*(1/2)*erfc(x/sqrt(2));
%BER = zeros(1 , length(snr_in_dB));
L = sqrt(M);
SER_theo = 4 * ( (L - 1) / L) * qfunc(sqrt(( 3 * k * snr) / (M - 1)));
BER_theo = SER_theo / k;


four_mapping = [-d d ; d d ;
                -d -d ; d -d]; % 4 QAM mapping

four_gray = [0 0 ; 0 1 ;
             1 0 ; 1 1];

four_E = [sqrt(2) ; sqrt(2);
          sqrt(2) ; sqrt(2)];

sixteen_mapping = [-3*d 3*d ; -d 3*d ; d 3*d ; 3*d 3*d ; 
                   -3*d d ; -d d ; d d ; 3*d d ;
                   -3*d -d ; -d -d ; d -d ; 3*d -d ;
                   -3*d -3*d ; -d -3*d ; d -3*d ; 3*d -3*d]; % 16 QAM mapping

sixteen_gray = [0 0 0 0 ; 0 0 0 1 ; 0 0 1 1 ; 0 0 1 0;
                0 1 0 0 ; 0 1 0 1 ; 0 1 1 1 ; 0 1 1 0;
                1 1 0 0 ; 1 1 0 1 ; 1 1 1 1 ; 1 1 1 0;
                1 0 0 0 ; 1 0 0 1 ; 1 0 1 1 ; 1 0 1 0];

sixteen_E = [sqrt(18) ; sqrt(10) ; sqrt(10) ; sqrt(18) ;
             sqrt(10) ; sqrt(2) ; sqrt(2) ; sqrt(10) ;
             sqrt(10) ; sqrt(2) ; sqrt(2) ; sqrt(10) ;
             sqrt(18) ; sqrt(10) ; sqrt(10) ; sqrt(18) ;];

sixtyfour_mapping = [-7*d 7*d ; -5*d 7*d ; -3*d 7*d ; -d 7*d ; d 7*d ; 3*d 7*d ; 5*d 7*d ; 7*d 7*d ;
                     -7*d 5*d ; -5*d 5*d ; -3*d 5*d ; -d 5*d ; d 5*d ; 3*d 5*d ; 5*d 5*d ; 7*d 5*d ;
                     -7*d 3*d ; -5*d 3*d ; -3*d 3*d ; -d 3*d ; d 3*d ; 3*d 3*d ; 5*d 3*d ; 7*d 3*d ;
                     -7*d d ; -5*d d ; -3*d d ; -d d ; d d ; 3*d d ; 5*d d ; 7*d d ;
                     -7*d -d ; -5*d -d ; -3*d -d ; -d -d ; d -d ; 3*d -d ; 5*d -d ; 7*d -d ;
                     -7*d -3*d ; -5*d -3*d ; -3*d -3*d ; -d -3*d ; d -3*d ; 3*d -3*d ; 5*d -3*d ; 7*d -3*d ;
                     -7*d -5*d ; -5*d -5*d ; -3*d -5*d ; -d -5*d ; d -5*d ; 3*d -5*d ; 5*d -5*d ; 7*d -5*d ;
                     -7*d -7*d ; -5*d -7*d ; -3*d -7*d ; -d -7*d ; d -7*d ; 3*d -7*d ; 5*d -7*d ; 7*d -7*d ]; % 64 QAM mapping

sixtyfour_gray = [1 0 0 1 0 0 ; 1 0 0 1 1 0 ; 1 0 1 1 1 0 ; 1 0 1 1 0 0 ; 0 0 1 1 0 0 ; 0 0 1 1 1 0 ; 0 0 0 1 1 0 ; 0 0 0 1 0 0;
                  1 0 0 1 0 1 ; 1 0 0 1 1 1 ; 1 0 1 1 1 1 ; 1 0 1 1 0 1 ; 0 0 1 1 0 1 ; 0 0 1 1 1 1 ; 0 0 0 1 1 1 ; 0 0 0 1 0 1;
                  1 0 0 0 0 1 ; 1 0 0 0 1 1 ; 1 0 1 0 1 1 ; 1 0 1 0 0 1 ; 0 0 1 0 0 1 ; 0 0 1 0 1 1 ; 0 0 0 0 1 1 ; 0 0 0 0 0 1;
                  1 0 0 0 0 0 ; 1 0 0 0 1 0 ; 1 0 1 0 1 0 ; 1 0 1 0 0 0 ; 0 0 1 0 0 0 ; 0 0 1 0 1 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 0;
                  1 1 0 0 0 0 ; 1 1 0 0 1 0 ; 1 1 1 0 1 0 ; 1 1 1 0 0 0 ; 0 1 1 0 0 0 ; 0 1 1 0 1 0 ; 0 1 0 0 1 0 ; 0 1 0 0 0 0;
                  1 1 0 0 0 1 ; 1 1 0 0 1 1 ; 1 1 1 0 1 1 ; 1 1 1 0 0 1 ; 0 1 1 0 0 1 ; 0 1 1 0 1 1 ; 0 1 0 0 1 1 ; 0 1 0 0 0 1;
                  1 1 0 1 0 1 ; 1 1 0 1 1 1 ; 1 1 1 1 1 1 ; 1 1 1 1 0 1 ; 0 1 1 1 0 1 ; 0 1 1 1 1 1 ; 0 1 0 1 1 1 ; 0 1 0 1 0 1;
                  1 1 0 1 0 0 ; 1 1 0 1 1 0 ; 1 1 1 1 1 0 ; 1 1 1 1 0 0 ; 0 1 1 1 0 0 ; 0 1 1 1 1 0 ; 0 1 0 1 1 0 ; 0 1 0 1 0 0];

sixtyfour_E = [sqrt(98) ; sqrt(74) ; sqrt(58) ; sqrt(50) ; sqrt(50) ; sqrt(58) ; sqrt(74) ; sqrt(98);
               sqrt(74) ; sqrt(50) ; sqrt(34) ; sqrt(26) ; sqrt(26) ; sqrt(34) ; sqrt(50) ; sqrt(74);
               sqrt(58) ; sqrt(34) ; sqrt(18) ; sqrt(10) ; sqrt(10) ; sqrt(18) ; sqrt(34) ; sqrt(58);
               sqrt(50) ; sqrt(26) ; sqrt(10) ; sqrt(2) ; sqrt(2) ; sqrt(10) ; sqrt(26) ; sqrt(50);
               sqrt(50) ; sqrt(26) ; sqrt(10) ; sqrt(2) ; sqrt(2) ; sqrt(10) ; sqrt(26) ; sqrt(50);
               sqrt(58) ; sqrt(34) ; sqrt(18) ; sqrt(10) ; sqrt(10) ; sqrt(18) ; sqrt(34) ; sqrt(58);
               sqrt(74) ; sqrt(50) ; sqrt(34) ; sqrt(26) ; sqrt(26) ; sqrt(34) ; sqrt(50) ; sqrt(74);
               sqrt(98) ; sqrt(74) ; sqrt(58) ; sqrt(50) ; sqrt(50) ; sqrt(58) ; sqrt(74) ; sqrt(98)];

if M == 4 
    for t = 1 : length(snr_in_dB) 
        symbolerror = 0;
        numofbiterror = 0;
        numofsymbolerror = 0;    

        for i = 1 : N % do N iterations
            a = randi(M); %random number from 1 to M
            sgma = sqrt(four_E(a) ./ (2 * k * snr));
            n(1) = normrnd(0 , sgma(t)); 
            n(2) = normrnd(0 , sgma(t));
            decis_index = 0; %index of the maximum distance
            r = four_mapping(a , :) + n; %received signal after adding gaussion noise
            for l = 1 : M %find maximum distance    
                c(l) = (r(1) - four_mapping(l,1))^2 + (r(2) - four_mapping(l,2))^2;
            end 
            
            [min_c decis_index] = min(c);
 
            if (four_gray(decis_index , 1) ~= four_gray(a , 1))
                numofbiterror = numofbiterror + 1;
            elseif (four_gray(decis_index , 2) ~= four_gray(a , 2))
                numofbiterror = numofbiterror + 1;
            end
            if (decis_index ~= a)
                numofsymbolerror = numofsymbolerror + 1;
            end

        end
        BER(t) = numofbiterror / (k * N);
        SER(t) = numofsymbolerror / N; 
    end
end

if M == 16
    for t = 1 : length(snr_in_dB)
        symbolerror = 0;
        numofbiterror = 0;
        numofsymbolerror = 0;    

        for i = 1 : N
            a = randi(M);
            sgma = sqrt(sixteen_E(a) ./ (2 * k * snr));
            n(1) = normrnd(0,sgma(t)); 
            n(2) = normrnd(0,sgma(t));
            decis_index = 0;
            r = sixteen_mapping(a,:) + n; %add noise to received signal
            for l = 1 : M %find maximum distance
                c(l) = (r(1) - sixteen_mapping(l,1))^2 + (r(2) - sixteen_mapping(l,2))^2;
            end 
            
            [min_c decis_index] = min(c);
 
            if (sixteen_gray(decis_index,1) ~= sixteen_gray(a,1))
                numofbiterror = numofbiterror + 1;
            elseif (sixteen_gray(decis_index,2) ~= sixteen_gray(a,2))
                numofbiterror = numofbiterror + 1;
            elseif (sixteen_gray(decis_index,3) ~= sixteen_gray(a,3))
                numofbiterror = numofbiterror + 1;
            elseif (sixteen_gray(decis_index,4) ~= sixteen_gray(a,4))
                numofbiterror = numofbiterror + 1;
            end
            if (decis_index ~= a)
                numofsymbolerror = numofsymbolerror + 1;
            end
            BER(t) = numofbiterror / (k * N);
            SER(t) = numofsymbolerror / N;
        end

            
    end
end

if M == 64 
    for t = 1 : length(snr_in_dB)
        symbolerror = 0;
        numofbiterror = 0;
        numofsymbolerror = 0;
            
        for i = 1 : N
            a = randi(M);
            sgma = sqrt(sixtyfour_E(a) ./ (2 * k * snr));
            n(1) = normrnd(0 , sgma(t)); 
            n(2) = normrnd(0 , sgma(t));
            decis_index = 0;
            r = sixtyfour_mapping(a , :) + n; %add noise to received signal
            for l = 1 : M %find maximum distance
                c(l) = (r(1) - sixtyfour_mapping(l,1))^2 + (r(2) - sixtyfour_mapping(l,2))^2;
            end 
            
            [min_c decis_index] = min(c);
 
            if (sixtyfour_gray(decis_index,1) ~= sixtyfour_gray(a,1))
                numofbiterror = numofbiterror + 1;
                symbolerror = 1;
            elseif (sixtyfour_gray(decis_index,2) ~= sixtyfour_gray(a,2))
                numofbiterror = numofbiterror + 1;
                symbolerror = 1;
            elseif (sixtyfour_gray(decis_index,3) ~= sixtyfour_gray(a,3))
                numofbiterror = numofbiterror + 1;
                symbolerror = 1;
            elseif (sixtyfour_gray(decis_index,4) ~= sixtyfour_gray(a,4))
                numofbiterror = numofbiterror + 1;
                symbolerror = 1;
            elseif (sixtyfour_gray(decis_index,5) ~= sixtyfour_gray(a,5))
                numofbiterror = numofbiterror + 1;
                symbolerror = 1;
            elseif (sixtyfour_gray(decis_index,6) ~= sixtyfour_gray(a,6))
                numofbiterror = numofbiterror + 1;
                symbolerror = 1;
            end
            if (symbolerror == 1)
                numofsymbolerror = numofsymbolerror + 1;
            end

        end
        
        BER(t) = numofbiterror / (k * N);
        SER(t) = numofsymbolerror / N;
    end
end

semilogy(snr_in_dB , BER , '*' , snr_in_dB , SER , 'o' , snr_in_dB , BER_theo)
axis([0 20 1e-4 1])
xlabel('SNR[dB]')                                    
ylabel('Bit Error Rate'); 

