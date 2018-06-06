function dat = ft_synchronize_trigger_for_trial(dat, fs, info_Flag)
%% synchronize the EEG and ECoG data with N_Trigger trigger signals,
%      taking all the data between the first and last trigger. 
%      Trigger information: about 0.08 - 0.1 S; heigith: 600 - 700 unites
%   by Vincent
%%
N_Trigger = 6;
TRIG_width_TH = fs*0.07;
TRIG_height_TH = 50; % the difference between up and down edge
[n_channel, n_sig_length] = size(dat);
ch_trigger = n_channel;
trigger_pos_List = zeros(1, N_Trigger);
tmp = zeros(1, N_Trigger);
trigger_pre  = dat(ch_trigger,1:n_sig_length-1);
trigger_post = dat(ch_trigger,2:end);
trigger_sig = trigger_pre - trigger_post;
trigger_sig(1:10)=0; trigger_sig(end-9:end)=0;
for iTrigger = 1 : N_Trigger
    tmp_trig_pos_1 = find(trigger_sig==min(trigger_sig),1);
    tmp_trig_pos_2 = tmp_trig_pos_1 + TRIG_width_TH*0.8 + find(trigger_sig(tmp_trig_pos_1+TRIG_width_TH*0.8:tmp_trig_pos_1+TRIG_width_TH*1.2) ...
        == max(trigger_sig(tmp_trig_pos_1+TRIG_width_TH*0.8:tmp_trig_pos_1+TRIG_width_TH*1.2)), 1);
    if ((tmp_trig_pos_1 - tmp_trig_pos_2) > TRIG_width_TH) ...
        &&  (abs(trigger_sig(tmp_trig_pos_1) + trigger_sig(tmp_trig_pos_2)) < TRIG_height_TH)
        tmp(iTrigger)  = tmp_trig_pos_2 - tmp_trig_pos_1;
    end
    trigger_pos_List(iTrigger) = tmp_trig_pos_1;
    trigger_sig(trigger_pos_List(iTrigger)) = 0;
end
% trigger_pos_List;
trigger_pos = [min(trigger_pos_List), max(trigger_pos_List)];
if info_Flag ==1
    disp(['data are taken from ',num2str(trigger_pos(1)),' to ',num2str(trigger_pos(2)),' and totally ',num2str(trigger_pos(2)-trigger_pos(1)) ,' points...'])
end
dat = dat(:,trigger_pos(1)+1:trigger_pos(2));
