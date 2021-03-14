close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_ADC = 3.2e6;
DEC = 1;

Fs = F_ADC/DEC;
Ts = 1/Fs;

Resolution = 0.5; % in Hz
N_F = Fs/Resolution;
F_axis = -Fs/2:Fs/N_F:Fs/2-Fs/N_F;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FrontEndSampleRate=F_ADC;
EU_FM_band_start=87.5e6; 
EU_FM_band_stop=108e6;
FM_BW=EU_FM_band_stop-EU_FM_band_start;
FM_BW_steps= FM_BW/FrontEndSampleRate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i= EU_FM_band_start+Fs/2
index=0;

%%
fi = fopen('spirosfifo', 'rb');

t_sampling = 0.1;      % seconds
%t_sampling = 0.05 ;      % seconds

N_samples = round(Fs*t_sampling);

t = 0:Ts:t_sampling-Ts;

counter = 0;
packets = 0;

HIST_SIZE = 10;    
 
fprintf('Lock the SDR to F=%3d Hz\n',i) 
input('Press ''Enter'' to continue...','s');

while i  < EU_FM_band_stop +Fs/2 
     
    index=index+1;
    ONE_band_temp_max_mval=0;
    while(1)
    

            x = fread(fi, 2*N_samples, 'float32');  % get samples (*2 for I-Q)
            x = x(1:2:end) + j*x(2:2:end);          % deinterleaving
            counter = counter + 1;
            
                
                packets = packets + 1;
              
                   % fft
                    x_fft = fftshift(fft(x, N_F));
                    % cfo estimate
                    [mval mpos] = max(abs(x_fft));
                    %mval
                    DF_est = F_axis(mpos);
               
              if  mval > ONE_band_temp_max_mval
                        ONE_band_temp_max_mval=mval;
                        max_freqs(1, index)= i+DF_est;  
                       
                       dBm_mval=10*log10((abs(mval).^2)*Ts);
                       max_freqs(2, index)=dBm_mval;
              
              end
          
                
                fprintf('Packets=%d Indexi=%d\n',packets, i) 

                %F_sense_hist(ii,packets) = F_sensor_est;
                
                 % fft
                x_fft = fftshift(fft(x, N_F));
                
                
                if 1
                    %plot
                  figure(1);
                    subplot(2, 1, 1);
                    plot(abs(x).^2);
                   % hold on;
                    %plot( imag(x_fft), 'g--');
                    %hold off;
                    drawnow;
                    subplot(2, 1, 2);
                  semilogy(F_axis, abs(x_fft).^2);
                  grid on;
                  axis tight;
                  drawnow;
                end
              
        
        if(mod(packets,HIST_SIZE) ==0)
            break; 
        end
  
   
 
    end
      
        max_freqs
        i=i+Fs;
        fprintf('Lock the SDR to F=%3d Hz\n',i) 
        input('Press ''Enter'' to continue...','s');
end

[M,I]= max(max_freqs(2,:));
maxfreq=max_freqs(1,I)
save('Powerfull_FM_stations_table.mat','maxfreq')

% evgale  95799857


