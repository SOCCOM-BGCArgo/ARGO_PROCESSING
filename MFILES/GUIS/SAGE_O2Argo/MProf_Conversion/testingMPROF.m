f_info.WMO = '5904474'; 
f_info.fn = '5904474_Mprof.nc';
f_info.dac_path = '/ifremer/argo/dac/5904474/'
f_info.local_path = 'C:\Users\tmaurer\Documents\ARGO\EXTERNAL_DACs\DAC_files\5904474\';
f_info.dac = 'aoml';

% mprof_out = mprof2mat(f_info);
tf_odv = mprofmat2ODV(f_info, f_info.local_path)