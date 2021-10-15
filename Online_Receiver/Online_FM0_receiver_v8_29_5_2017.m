%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis  
clc; 
close all; 
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pause(10)
%%RTL SDR parameters
GAIN=-15;
F_ADC = 1e6;
DEC = 1;
Fs = F_ADC/DEC;
Ts = 1/Fs;

%% Sympol parameters
Tsymbol = 0.990e-3;   %408e-6; % 4.95e-3
% vazo to length to mikroterou palmou

%0.990e-3; gia 500 bps
%500e-6; gia 1 kbps
%202e-6 gia 2 kbps (try 198 )

Tbit=Tsymbol*2;
over = round(Tsymbol/Ts);
newover = 10;
%% PAcket parameters  
preamble_length=10;                % NoFM0_prample=[1 0 1 0 1 0 1 1 1 1];
id_length=2;                       % NoFM0_ID=[0 1];
util_length=2;                     % NoFM0_util=[0 1];
codeword_length=12;                % NoFM0_DATA=[0 0 1 1 1 1 0 0 0 1 0 1];
dummybit=1;

total_packet_length=id_length+preamble_length+util_length+codeword_length+dummybit;
total_packet_duration=total_packet_length*Tbit;
preamble_duration=preamble_length*Tbit;

% Create the preamble in FM0 format.
preamble_bits=[1 1 0 1 0 0 1 0 1 1 0 1 0 0 1 1 0 0 1 1];
preamble = preamble_bits;       % Variable containing the preamble
%try 2*preamble_bits-1  
fixedata=[0 0 1 1 1 1 0 0 0 1 0 1];
fixedpacketdata=[0 1 0 1 0 0 1 1 1 1 0 0 0 1 0 1];


%% Sigmal Prosesing  Variables
Resolution = 1;   % in Hz
N_F = Fs/Resolution;
F_axis = -Fs/2:Fs/N_F:Fs/2-Fs/N_F;

framelength=3;
t_sampling = framelength*total_packet_duration;     % Sampling time frame (seconds).
N_samples = round(Fs*t_sampling);



fi = fopen('spirosfifo', 'rb');

t = 0:Ts:t_sampling-Ts;


HIST_SIZE =600;
dataset=[];

% Debug Print variables
DEBUG_en1=1;
DEBUG_en2=1;
%% Decoder variables
correct_packets=0;
nopacket_ind=0;
nodroped_packets=0;
cut_packets=0;
error_packets=0;
negative_starts=0;
negative_starts2=0;
packet_yess=0;
droped_packets=0;
pos=1;
FLIPPED=0;
packets = 1;
counter=0;
BER_int=0;
%% Decoder variables
bits_FM0_2sd_wayB=[];
decision_bits_B=[];
BER_int=[];
infomatr=[];

%% Orthogonal pulces for detection

D1_ups=zeros(1,newover*2);
D1_ups(1:newover)=1; 
D1_ups(newover+1:newover*2)=-1;

D2_ups=zeros(1,newover*2);
D2_ups(1:newover)=-1; 
D2_ups(newover+1:newover*2)=1;



while (1)

            x = fread(fi, 2*N_samples, 'float32');  % get samples (*2 for I-Q)
            x = x(1:2:end) + j*x(2:2:end);          % deinterleaving
             
            counter = counter + 1;
         %dataset= [dataset ; x];  
      if ~mod(counter, 2)   %delay every two pakcets || ===> capture_____delay(duration=packet)____capture_____......  
    
      
           packets = packets + 1;
           fprintf('Packet=%d|\n',packets)
           
           
      %% Absolute operation removes the unknown CFO
      abstream=abs(x).^2;
      %% Matched filtering
      matcheds=ones(round(Tsymbol/Ts),1);
      dataconv=conv(abstream,matcheds);      
      
           %% Signal Smoothing with  1-D median filtering
            %dataconv1= medfilt1(dataconv); 
            dataconv1= (dataconv); 
            %% Downsample same prosedure
            total_env_ds = dataconv1(1:over/newover:end);
            
            %% Time sync
            total_envelope = total_env_ds(newover+1:end-newover+1); 
            
            total_envelope=total_envelope-mean(total_env_ds);
            
            
            % Calculate the energy in the packet. If the energy is less than a
            % threshold, discard packet.
             energyther= sum(total_envelope.^2)
            if(energyther <= 5) %200 for 1ms 
                  droped_packets = droped_packets + 1;
                  disp 'Drop packet';
               continue;

            else
                nodroped_packets = nodroped_packets + 1;
                packet_yess=1;
            end
             
            %% Flipe the packet  if its nessesarry 
            %%  Position estimetion of packet with  packet's energy synchronization
            for k=1:1: length(total_envelope)- (total_packet_length*2*newover)+1
                    energy_synq(k)=sum(abs(total_envelope(k : k+total_packet_length*2*newover-1- newover)).^2);          
            end
            % find the starting point of packet
            [energy_sinq_max  energy_sinq_ind]=max(energy_synq);
            
          
             
             % exo frame_window= 3*packets_length  
             % tha anoixneuso an to paketo einai pano  i kato apo tin mesi (frame_length/2)
             % *******1) an einai pano apo tin mesi tou frame_window (part=2) tote
             % elenxo an o meso oros  enos window (diarkias total_packet_length/2)  prin tin axi tou paketou einai
             % megaliteros apo ton meso oro ollou to frame_window -> an einai tote kano flip to paketo 
             % *******2) an einai kato apo tin mesi tou frame_window (part=1) tote
             % elenxo an o meso oros  enos window (diarkias total_packet_length/2)  meta to telos tou paketou einai
             % megaliteros apo ton meso oro ollou to frame_window -> an einai tote kano flip to paketo 
             
            if (energy_sinq_ind*framelength/length(total_envelope) >= 1.5)
                part=2;
                
                if  mean(total_envelope) <= mean(total_envelope(energy_sinq_ind-floor(total_packet_length*2*newover/4): energy_sinq_ind-1)) 
                   FLIPPED=1;
                   total_envelope=-total_envelope;
                end
                
            elseif  (energy_sinq_ind*framelength/length(total_envelope) <1.5)
                part=1;
                Packetend=energy_sinq_ind+ floor(total_packet_length*2*newover+1);
               
                if  mean(total_envelope) <  mean(total_envelope(Packetend : Packetend +total_packet_length*2*newover/4))
                   FLIPPED=1;
                   total_envelope=-total_envelope;
                end
            end
            
%             if(mean(total_envelope)>0) %200 for 1ms 
%                   disp 'Drop packet';
%                   droped_packets = droped_packets + 1;
%                continue;
%             end 
            %total_envelope=abs(total_envelope);
            %% Print Plots

           if DEBUG_en1==1;
                time_axis= 0:Ts:Ts*length(abstream)-Ts;        %same as xaxis_m= (1: length(abstream))*Ts Captured signal time axis.
                     % fft
                     x_fft = fftshift(fft(x, N_F));
                 figure(1);
                  subplot(2, 1, 1);
                    plot(time_axis,abstream);
                    title('Time Domain')  
                    xlabel('Time (Sec)');
                  
                   subplot(2, 1, 2);
                    semilogy(F_axis/1000000, abs(x_fft).^2);
                    title('Frequency Domain')  
                    xlabel('Frequency (MHz)');
                   drawnow;
                          
           end 
          
         
             if DEBUG_en2==1;
               
                 figure(2);
                    subplot(2, 1, 1);
                    plot(dataconv);
                    title('Matching Filtering')
                    subplot(2, 1, 2);
                    plot(total_envelope);
                    title('FLIPPED DOWNSAMPLED')
                    drawnow;                
             end 
             
             
            
            
            %% dc zero offser
             %% Assume symbol synchronization, which can be implemented using correlation with a sequence of known bits in the preamble       
              % comparison of the detected preamble bits with the a priori known bit sequence
              %convert the header to a time series for the specific sampling frequency and bit duration. 
            
            %% create the preamble neover format
            preample_neover=upsample(preamble, newover);
           
            
            %% Sync via ENERGY
             for k=1:1: length(total_envelope)- (total_packet_length*2*newover)+1
                    energy_synq(k)=sum(abs(total_envelope(k : k+total_packet_length*2*newover-1)).^2);          
             end
            [energy_sinq_max  energy_sinq_ind]=max(energy_synq);
            sumxor=0;
            
           if energy_sinq_ind-2*newover<=0
               negative_starts2=negative_starts2+1;
               negative_starts=negative_starts+1;
               disp 'Negative start_2';
                continue;
           end 
            
            %% Sync via preamble correlation
            corrsync_out = xcorr(preample_neover, total_envelope);
              
            [m ind] = max(corrsync_out);
             %notice that correlation produces a 1x(2L-1) vector, so index must be shifted.
             %the following operation points to the "start" of the packet.
            start = length(total_envelope)-ind + 1;
            

             if(start <= 0)
                negative_starts = negative_starts + 1;
                disp 'Negative start_1';
                continue;
            end
            
              %% Check if the detected packet is cut in the middle.
            if start+(total_packet_length*2)*newover > length(total_envelope)
                cut_packets = cut_packets + 1;
                disp 'Packet cut in the middle!';
                continue;
            end    
            %% move  the signal to 0
            preamble_thress2=mean(total_envelope(start: start+length(preample_neover)));
            signal_to_zero2=total_envelope-preamble_thress2;
            
            %if (start+total_packet_length*2*newover <= length(signal_to_zero2))
            %o assos den eimai sigouros oti xreisimeuei !!!!!!!!!!!!!!!
                shifted_sync_signal_B=signal_to_zero2(start+length(preample_neover)-newover-1: start+total_packet_length*2*newover);
            %end
            
              
                  for xi=1:newover*2: length(shifted_sync_signal_B)-newover*2
                   
                        sample2=shifted_sync_signal_B(xi: xi+newover*2-1);
                        
                        sumD1_ups=sum(D1_ups.*sample2');
                        sumD2_ups=sum(D2_ups.*sample2');
                         if (sumD1_ups > sumD2_ups)
                            bits_FM0_2sd_wayB=[bits_FM0_2sd_wayB, 1];    
                         else
                            bits_FM0_2sd_wayB=[bits_FM0_2sd_wayB, 0];  
                        end  
                    end
                   
              jim=1;
            for indx=2:1:length(bits_FM0_2sd_wayB)
               if bits_FM0_2sd_wayB(indx) == bits_FM0_2sd_wayB(indx-1)
                     decision_bits_B(jim)=0;
               else
                     decision_bits_B(jim)=1;
                end
                   jim=jim+1;
                
            end
            
            id_est_B = decision_bits_B(1: id_length);
            sensor_id_est_B = decision_bits_B(id_length + 1: id_length + util_length);
            data_bits_es_B = decision_bits_B(id_length+util_length + 1:end);
               
            if  isequal(decision_bits_B, fixedpacketdata)
                        disp 'Packet Correct !!!!!!!!!!!!!!!!!!!!!!!!!';
                       correct_packets=correct_packets+1;
                       
            elseif ~isequal(decision_bits_B, fixedpacketdata) 
                 disp 'Packet WRONGGGG';
                error_packets=error_packets+1;
                BER_sum(error_packets) = sum(xor(decision_bits_B,fixedpacketdata));
           
            end
            
            
            
              
          % pause(1);       
         bits_FM0_2sd_wayB=[];
         decision_bits_B=[];
         

      end
        
      if(mod(packets,HIST_SIZE) ==0)
          
          
         
           infomatr(1)=correct_packets;
           infomatr(2)=error_packets;
           infomatr(3)=negative_starts;
           infomatr(4)=cut_packets;
           
           
         
            
            
            fprintf('Corecct Packets=%d|Packet Error=%d\n',correct_packets, error_packets) 
            fprintf('Negative Starts=%d|Cut Packets=%d\n', negative_starts, cut_packets)
            fprintf('Negative Starts2=%d\n', negative_starts2)
            fprintf('Drop Packets=%d\n', droped_packets)
           if (error_packets ~= 0)
               
           infomatr(5)=error_packets /(correct_packets+error_packets);
           infomatr(6)= sum(BER_sum);
           infomatr(7)= sum(BER_sum)/((correct_packets+error_packets)*length(fixedpacketdata));
            %PER is the number of incorrectly received data packets divided by the total number of received packets.
            fprintf('Packet Error Rate=%d\n', infomatr(5)) 
            fprintf('Bit error rate (BER)=%d\n', infomatr(7))
           else
              infomatr(5)=0;
              infomatr(6)=0;
              infomatr(7)= 0;
            fprintf('Packet Error Rate=%d\n',0) 
            fprintf('Bit error rate (BER)=%d\n',0)
              
          end
            
%         rootname = 'FM_amb_1ms_15G_95_8_MHz_300_pack_meas_025m_'; % Root filename
%         data = clock;
%         filename = [rootname, num2str(data(4)) num2str(data(5))...
%         num2str(data(2)) num2str(data(3)) num2str(data(1))];
%         save(filename,'dataset')
        
           %save('FM_1ms_15G_1000_pack_meters','dataset')
           save('results_FM_amb_1ms_15G_95_8_MHz_100_pack_meas_demo','infomatr')
             return;
      end
        
end




