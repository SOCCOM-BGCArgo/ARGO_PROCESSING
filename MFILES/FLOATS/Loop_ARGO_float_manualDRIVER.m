
% Loop_ARGO_float_manualDRIVER.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT TO RUN Loop_ARGO_float.m IN MANUAL MODE ("ALL").  THIS SCRIPT ALSO
% SERVES AS A LOG HISTORY OF SETTINGS USED IN CALLS TO THE FUNCTION.  IE 
% A HISTORY OF FLOAT SUBSETS USED IN THE CALL, ETC.
%
% TANYA MAURER
% MBARI
% 06/05/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 05/25/2017, TM: processed all files for floats 0037, 9634 (for testing on
% Sterna).
% Loop_ARGO_float('all','0037|9634',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

%05/25/2017, TM: processed all files for floats 9254, 9275 (for testing on
% Sterna).
% Loop_ARGO_float('all','9254|9275',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 06/12/17, TM: reprocessed NAVIS floats 0568 and 0570
% (complete msg file for cycle 40 finally came throguh for both floats,
% along with cycle 41.) R/D files for #40 now exist on GDAC for these floats.
% Also reprocessed float 7622; was an error in generation of CFG file.
% Loop_ARGO_float('all','0568|0570',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','7622',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 06/13/17, TM: reprocessed APEX floats 12361 and 12388 from PIPERS cruise
% in order to incorporate O2 Gain.  (Initial QC, 6 cycles for each float)
% Loop_ARGO_float('all','12361|12388',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 06/19/17, TM: reprocessed APEX floats 12573 cycle 9 (had an incomplete msg file)
% Also reprocess APEX floats 12371 and 12549 (incomplete msg files now
% updated)
% Loop_ARGO_float('all','12573',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12371|12549',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 06/19/17, TM: reprocessed APEX floats 9101 9761 12363 12378
%Loop_ARGO_float('all','9101|9761|12363|12378',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 06/21/17, TM: reprocessed APEX floats 9652 for cycle 52 (052.isus was
% late
% Loop_ARGO_float('all','9652',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 06/29/17, TM: reprocessed APEX float 9659 (WMO 5904843) that had crazy O2
% value.  Testing mods to change O2 vals > 99990 to 99990
% 07/03/17, TM: reprocessed 9659 once again to fix DOXY Adjusted data parameter
% crazy bad values.
% Loop_ARGO_float('all','9659',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 07/17/17, TM: reprocessed APEX float 12543 (WMO 5904980). Was getting
% matchtest fail for BR transfer because we had an incomplete msg file for 018 (and
% mat file) on our end (only ~15KB).  Full file finally came through.  Need
% to reprocess, then BRtransfer routine should pick it up next iteration.
% Loop_ARGO_float('all','12543',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 07/18/17, TM: Reprocessed data for NAVIS float 0569 to incorporate new pH
% pressure coefficients derived by J.P.
% Also, tested updated processing code that now includes 'bad sensor list'.
% Loop_ARGO_float('all','9018',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0569',{'8514HAWAII' '8501CALCURRENT'})
% Also, reprocessed ALL floats at the end of the day!
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 07/19/17, TM: reprocessed float 0569 to pick up notes for header in ODV
% file.  Reprocessed floats with updated QC adjs.
% Loop_ARGO_float('all','0569',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12371|0691|0692|12381|12541|12390|12549|12559|0567|12543|12545',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 07/20/17, TM: reprocessed float 12551, 9274, 9659 to pick up notes bad O2 sensors. 
% Loop_ARGO_float('all','9274|9659',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 07/24/17, TM: reprocessed 4 floats for updated QC. 
% Loop_ARGO_float('all','9750|12575|12573',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','8501',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9646',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

% 07/27/17, TM: 4 float QC updates
% Loop_ARGO_float('all','12382|9642|9762|9634',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

%07/29/2017, TM: Process 12380 float 1st profile
%Loop_ARGO_float('all','12380',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

%07/31/2017, TM: Pick up salinity flagging for float 9630
% Loop_ARGO_float('all','9630',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12545',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

%8/01/2017, TM: Reprocess to pick up cycle 58 that came in full
% Also, reprocess floats with bad salinity to propagate changes in code
% that auto-set other paramters to 8
%Loop_ARGO_float('all','9668',{'8514HAWAII' '8501CALCURRENT'})
%Loop_ARGO_float('all','9630|0506|0568|9018|',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------

%8/02/2017, TM: Reprocess to pick up QC changes
% Loop_ARGO_float('all','9632|9660|9265|9744',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9630|9744',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9650|9637|9602',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
% 
% %8/03/2017, TM: Reprocess to pick up QC changes
% Loop_ARGO_float('all','9757|9749|9662|9655|9657|9652',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
%
%8/11/2017, JP Deleted cal file & reprocessed so WMO could found
% Loop_ARGO_float('all','12380',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
%
%8/14/2017, TM reprocess 9766, typo in NO3 cal file so no NO3 data was
%being processed.  Reprocess 12366 to pick up complete cycle 17 (for
%BR-transfer).
% Loop_ARGO_float('all','9766',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12366',{'8514HAWAII' '8501CALCURRENT'})

%--------------------------------------------------------------------------
%
%8/14/2017, JP reprocessed 8486 & 7672 using udated LIPHR. 8486 may have a
%stry light issue with nitrate. Hawaii - no surf nitrat so corrected to 0
%at surface so off at depth
% Loop_ARGO_float('all','8486|7672',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
%
%8/16/2017, TM reprocess 9274 hawaii after fixing bad sensor list bug
%(related to case insensitivity in file names).
% Also, some QC reprocessing for select floats.
% Loop_ARGO_float('all','9274',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9668|8514|9125|9260',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9313|9254|9101|9095',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
%
%8/17/2017, TM reprocess 9099 to pick up ad-hoc salinity adjustment
%affecting MLR in SAGE
% Loop_ARGO_float('all','9099',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
%
%8/22/2017, TM reprocess 12372, and updated O2 gain for 0691
% Loop_ARGO_float('all','12372',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0691',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
%
%8/23/2017, TM reprocess numerous P16 floats after QC updates
% Loop_ARGO_float('all','9031|9092|9091|7613|7567|7557|6091',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9092|9091|12552|12558',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9018',{'8514HAWAII' '8501CALCURRENT'})
%--------------------------------------------------------------------------
%
%8/24/2017, TM reprocess floats after QC updates
% Loop_ARGO_float('all','9631|9094',{'8514HAWAII' '8501CALCURRENT'})

%8/24/2017, JP cal files got corrupted. Rebuilt files & reprocessed floats
% updated QC adjustments whil I was at it O2 gain ws OK
% Loop_ARGO_float('all','0506|0507',{'8514HAWAII' '8501CALCURRENT'})

%8/24/2017, 8486Hawaii - pH pressure coefficents have changed since lab calibration. new
% coeficients were determined by optimization against HOT - ALOHA bottle pH
% data and GLODAPV2 pH data > 1000m. A new config file was uploaded to
% Sterna. Nitrate data was also problematic WL offset was set to 208.5 in
% the natrate cal file and then nitrate was corrected to the surface (0 µM)
% the udated isus cal file and QC adjustment file were uploaded to sterna
%Loop_ARGO_float('all','8486',{'8514HAWAII' '8501CALCURRENT'})

%8/24/2017, 0412Hawaii - txt files messed up with odd data. Corrected pH to
%LIPHR at 900m and HOT -ALOHA bottle data, no nitrate only 33 good profiles
%this is a float without a WMO number
% Loop_ARGO_float('all','0412',{'8514HAWAII' '8501CALCURRENT'})

%8/29/2017, 0412Hawaii - QC update for some NAVIS floats
%Loop_ARGO_float('all','0566|0571|0565',{'8514HAWAII' '8501CALCURRENT'})

%8/30/2017, Reprocessed a bunch of Hawaii floats. Gain adjustments for all
%floats with nitrate sensors have been set = 1. The WaveLength offset
%varaible was used instead to tune the floats. In general nitrate was
%corrected at the surface (=0) and then the wavelength offset was was tuned
%so that float nitrate matched the LINR or HOT bottle data at depth (ideally both). pH was
%corrected to LIPHR and HOT bottle data at 900m. (LIPHR and HOT bottle data
%diverge below this depth a bit per Yui's investigation. A new set of pH
%pressure coefficents were determined for 8486
% floats_str = '7593|9274|8514Hawaii|7672|7622|6966|6401|6403|6891|5145|0412|8486|8497';
% Loop_ARGO_float('all',floats_str,{'8501CALCURRENT'})

%8/30/2017, TM, Reprocessing a number of NAVIS floats for updates to QC.
%Many of these had not been revisited in many months.  O2 gains were
%changed on many to reflect bottle data, with weight given to surface
%samples. 
%Loop_ARGO_float('all','0566|0571|0565|0570|0510|0507|0511|0037|0508|0509',{'8514HAWAII' '8501CALCURRENT'})

%8/31/2017, TM, Modified bad sensor list for 9750 pH, rerun.  QC update
%0564 (dead float but updated data available).
% Loop_ARGO_float('all','9750',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0564',{'8514HAWAII' '8501CALCURRENT'})

%09/01/2017, JP,Updated all Station P floats, set gain back to 1 for 5143
%and tuned wavelength offset = 202 to matsh surface and deep - A big change
% Loop_ARGO_float('all','5143|6400|6881|6972|7601|7641',{'8514HAWAII' '8501CALCURRENT'})

%09/12/2017, JP,Updated 7550 to recover late *.isus data for cycle 237 and to
%see if this helps Uday's nitrate data problem for this float per email 
% forwarded from Ken on 9/12/17
%and tuned wavelength offset = 202 to match surface and deep - A big change
% Loop_ARGO_float('all','7550',{'8514HAWAII' '8501CALCURRENT'})

%09/12/2017, TM,Updated two more floast to recover late .isus file.
%and tuned wavelength offset = 202 to matsh surface and deep - A big change
% Loop_ARGO_float('all','12396|12537',{'8514HAWAII' '8501CALCURRENT'})

%09/17/2017, TM,Reprocess 12729.  Dock test files still in directory
% Loop_ARGO_float('all','12729',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12744',{'8514HAWAII' '8501CALCURRENT'})

%09/17/2017, TM, Added 12729 and 0508 to bad sensor list (pH).  Most
%getting autoflagged, but a couple squeaking through.
% Loop_ARGO_float('all','12729|0508',{'8514HAWAII' '8501CALCURRENT'})

%09/20/17, TM, reprocess 12729, giving error (old msg files still existed
%in mat dir).  Reprocess 12723, isus file was absent when processing
%occurred.
% Loop_ARGO_float('all','12729|12723',{'8514HAWAII' '8501CALCURRENT'})

%09/25/17, TM, reprocess 19662.  matchtest error in BRtransfer (cycle 60).
%Reprocess to generate full .mat file
% Loop_ARGO_float('all','9662',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12369',{'8514HAWAII' '8501CALCURRENT'})

% %09/26/17, TM, Reprocess 12380, 12372, initial QC
% Loop_ARGO_float('all','12380|12372',{'8514HAWAII' '8501CALCURRENT'})
% 
% %09/27/17, TM, Reprocess 12540, and others with nowmo stuck.
% % Reprocess 12540, wrong DF#, per Hans.
% Loop_ARGO_float('all','12540',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12723',{'8514HAWAII' '8501CALCURRENT'})

% Reprocess Arctic floats - [O2] was greater than max allowable range. Set
% range limit to 550 instead of 450
% Loop_ARGO_float('all','7564ARCTIC|7596ARCTIC',{'8514HAWAII' '8501CALCURRENT'})

% %10/02/2017, TM, Reprocess 0568, mismatch in BR transfer
% reprocess 7641, missing isus files
% Loop_ARGO_float('all','0568',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','7641',{'8514HAWAII' '8501CALCURRENT'})

% %10/09/2017, TM, Reprocess 12369, update to O2 flagging for crazy large
% values
% Loop_ARGO_float('all','12369',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9630',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12733',{'8514HAWAII' '8501CALCURRENT'})

% %10/10/2017, TM, Reprocess floats with missing .isus files (these three
% have isus plumbed but no longer reporting files)
% Loop_ARGO_float('all','0569',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','6091',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','7622',{'8514HAWAII' '8501CALCURRENT'})
% Revisited QC on handful of newer floats (starting to go through SOCCOM
% floatQC again, newest floats first.)
% Loop_ARGO_float('all','12372|12380|12361|12388|12542|12386|12366',{'8514HAWAII' '8501CALCURRENT'})

% %10/16/2017, TM, Reprocess 12733, now has WMO#; reprocess 9757 (problem
% with ArgoTransfer)
% Loop_ARGO_float('all','12733|9757',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0570',{'8514HAWAII' '8501CALCURRENT'})

% %10/23/2017, TM, Reprocess 0570, BR transfer issue;
% Reprocess numerous floats for updated QC
% Loop_ARGO_float('all','0570',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12744|8204|12396|12366|12558|12541|12549|12559|12543|12545|8501|12552|12382',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9642|9766|9762|9634',{'8514HAWAII' '8501CALCURRENT'})

% %10/24/2017, TM, 
% Reprocess numerous floats for updated QC
% Loop_ARGO_float('all','9602|9637|9650|9600|9631|9744|9265|9660|9632',{'8514HAWAII' '8501CALCURRENT'})

% %10/25/2017, JP, 
% Reprocess 3 floats because there were FLBB cal mistakes
% Loop_ARGO_float('all','7652|9274|9646',{'8514HAWAII' '8501CALCURRENT'})
% Added 9632 to bad sensor list for pH.  Reprocess additional floats for
% updated floatQC.
% Loop_ARGO_float('all','9632|9757|9645|9749|9662',{'8514HAWAII' '8501CALCURRENT'})

% % %10/30/2017, TM, 
% % Reprocess 9094, problem with interpolated position values
% Loop_ARGO_float('all','9094',{'8514HAWAII' '8501CALCURRENT'})
% % Process new floats 12717EqPacE and 12728EqPacE
% Loop_ARGO_float('all','12717|12728',{'8514HAWAII' '8501CALCURRENT'})

%10/31/17 TM  Reprocess for updated qc
% Loop_ARGO_float('all','9094|12742|12729|0570',{'8514HAWAII' '8501CALCURRENT'})

%11/6/17 TM  Reprocess 12717 after updates to pH diagnostic range checks
%(they were loosened)
% Loop_ARGO_float('all','12717',{'8514HAWAII' '8501CALCURRENT'})

%11/7/17 TM  Reprocess 0570, BR mismatch error
% Loop_ARGO_float('all','0570',{'8514HAWAII' '8501CALCURRENT'})

%11/13/17 TM  Reprocess 12728, .dura file time lag --> reproccess to pick
%up pH data in full for cycle 3.
% Loop_ARGO_float('all','12728',{'8514HAWAII' '8501CALCURRENT'})

% %11/20/17 TM  12573 generating matchtest error.  Try to reprocess.
% Loop_ARGO_float('all','12573',{'8514HAWAII' '8501CALCURRENT'})

%11/20/17 JP  UPDATE SOME BBP QUALITY FLAGS.
% Loop_ARGO_float('all','9602|7614|7552|12573|0569|6968|6967 ',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','6967',{'8514HAWAII' '8501CALCURRENT'})

%11/26/17 TM UPDATE 0569, BR TRANSFER ERROR
% Loop_ARGO_float('all','0569',{'8514HAWAII' '8501CALCURRENT'})

%12/4/17 TM UPDATE 9655, BAD PH SENSOR IDENTIFIED
% Loop_ARGO_float('all','9655',{'8514HAWAII' '8501CALCURRENT'})

%12/11-12/17 TM QC UPDATES
% Loop_ARGO_float('all','0507|0567|7613|7652',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9657|9652|9646|9668',{'8514HAWAII' '8501CALCURRENT'})
%12369 and 12370 FLBB were switched (per Dana Swift).  After corrected,
%need reprocess (12370 not yet deployed).
% Loop_ARGO_float('all','12369',{'8514HAWAII' '8501CALCURRENT'})

%12/13/17 TM QC UPDATES
% Loop_ARGO_float('all','9096|9260|9101|9254|9313',{'8514HAWAII' '8501CALCURRENT'})

%12/14/17 TM QC UPDATES
% Loop_ARGO_float('all','12366|12379|12559',{'8514HAWAII' '8501CALCURRENT'})

%12/16/17
% Loop_ARGO_float('all','0690',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0566',{'8514HAWAII' '8501CALCURRENT'})

%12/20/17 Initial QC for a handful of floats.
% Loop_ARGO_float('all','12745|12733|12540|12723',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12369',{'8514HAWAII' '8501CALCURRENT'})

%12/21/17 Initial QC for a handful of floats.
% Loop_ARGO_float('all','9275|9752|9750|12390|12537',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9666',{'8514HAWAII' '8501CALCURRENT'})

%01/03/18 Test complete reprocess with CHL modifications.  DO NOT COPY TO
%NETWORK
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})

%01/04/18 Float-specific updates to BSL and QClists
% Loop_ARGO_float('all','0510|0565|0571|0691|0692|9091|9092|9095|9125|12361|12369|12573|12575|0276|9750|12542|12372|12380',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0565',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12372|12380|12542',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12537',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12386',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12388|12744',{'8514HAWAII' '8501CALCURRENT'})
% % complete reprocess again to pickup finalized headers.
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})
% 
% %01/08/18 Float-specific updates to QC (out from under ice)
% Loop_ARGO_float('all','7614|12541',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12551',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','7614',{'8514HAWAII' '8501CALCURRENT'})
% 
% %01/16/18 BRfile error.
% Loop_ARGO_float('all','12542',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12551',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0688',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0690',{'8514HAWAII' '8501CALCURRENT'})
% 
% %01/22/18 error with 12736 processing (missing .isus and .dura on first
% %cycle)
% Loop_ARGO_float('all','12736',{'8514HAWAII' '8501CALCURRENT'})
% 
% %01/23/18 reprocess certain floats after activation of their alternate
% %directories
% Loop_ARGO_float('all','12381|12361|12366|12379|12559',{'8514HAWAII' '8501CALCURRENT'})
% 
% %01/25/18 reprocess 12757 only partial file before.
% Loop_ARGO_float('all','12757',{'8514HAWAII' '8501CALCURRENT'})
% %01/25/18 reprocess 12757 only partial file before.
% Loop_ARGO_float('all','9101',{'8514HAWAII' '8501CALCURRENT'})
% 
% % 01/29/18
% Loop_ARGO_float('all','12381',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0068',{'8514HAWAII' '8501CALCURRENT'})
% 
% % 01/31/18 Reprocess to QC O2 outlier
% Loop_ARGO_float('all','9642',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12709',{'8514HAWAII' '8501CALCURRENT'})
% 
% %02/01/2018 Reprocessed all that had empty alternates after Dana
% %initialized
% Loop_ARGO_float('all','9666|12363|12390|12537|12551|12575|12717|12543|8374|9750|0567|0569|0692',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0690',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12549|12378',{'8514HAWAII' '8501CALCURRENT'})


% %02/01/2018 Reprocessed 7564ARCTIC - add 3 bad bbp profiles to the bad sensor
% %list: 44,45,106 -jp
% Loop_ARGO_float('all','7564',{'8514HAWAII' '8501CALCURRENT'})

%02/02/2018 Reprocessed 12540SOOCN Adjusted NO3 & pH corrections 980-1420 depth range corrected to LIR's
%list: 44,45,106 -jp
% Loop_ARGO_float('all','12540',{'8514HAWAII' '8501CALCURRENT'})

%02/02/2018 Reprocessed 12396SOOCN Adjusted NO3 & pH corrections 980-1420 depth range corrected to LIR's
% Nitrate really drifting may need to mark bad soon
%list: 44,45,106 -jp
% Loop_ARGO_float('all','12396',{'8514HAWAII' '8501CALCURRENT'})

%02/05/18 reprocess 0037 for BR transfer update test
% Also: reprocess all files to propagate <param>_DATA_MODE variables into
% historical Argo mat files.
% Loop_ARGO_float('all','0037',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})

%02/06/18
% Loop_ARGO_float('all','9101',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})

%02/07/18 FLOAT QC....
% Loop_ARGO_float('all','9646',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9313',{'8514HAWAII' '8501CALCURRENT'})
% 
% %02/12/2018
% Loop_ARGO_float('all','8501',{'8514HAWAII' '8501CALCURRENT'})
% %float QC
% Loop_ARGO_float('all','12371|9666|12381|12549|12734|12755',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12558|12545|12552',{'8514HAWAII' '8501CALCURRENT'})


% 02/20/2018 reprocess 5 floats which had incorrect FLBB calibrations -jp
%Loop_ARGO_float('all','8501CALCURRENT',{'8514HAWAII'})
%Loop_ARGO_float('all','7652|9091|9650|9652',{'8514HAWAII' '8501CALCURRENT'})
%Loop_ARGO_float('all','9091',{'8514HAWAII' '8501CALCURRENT'})

% 02/20/2018 More BB cal fixes - right cal but wrong serial # -jp
% Loop_ARGO_float('all','12730|12782',{'8514HAWAII' '8501CALCURRENT'})

% 02/26/2018 
% Loop_ARGO_float('all','12363|12542',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','8204|12543|8501',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9757',{'8514HAWAII' '8501CALCURRENT'})

% 02/28/2018 QC & bad sensor updates -jp
% Loop_ARGO_float('all','12549|12541',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12717',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12757',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9275|9762|9744|12369|9018|9091',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12575',{'8514HAWAII' '8501CALCURRENT'})

% 03/01/2018 rerun to remove correction -jp
% qc updates -tm
% Loop_ARGO_float('all','9275|9766',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12390',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9637',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12382|9642|9766|9762|9634',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12730|12757',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9668',{'8514HAWAII' '8501CALCURRENT'})

% 03/05/2018 Initial QC & other updated QC
% Global reprocessing this evening prior to snapshot creation for March
% 2018
% Loop_ARGO_float('all','12781|12736',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9660|9265',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9602|9650',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12396',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})
% 
% Loop_ARGO_float('all','9600|12366|9752',{'8514HAWAII' '8501CALCURRENT'})


% 03/12/2018 Initial QC 
% Loop_ARGO_float('all','12769',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12370|12702',{'8514HAWAII' '8501CALCURRENT'})

% 03/19/2018 Update QC 
% Loop_ARGO_float('all','12545',{'8514HAWAII' '8501CALCURRENT'})

% 03/28/18 Reprocess to fix ArgoB transfer error
% Loop_ARGO_float('all','9092',{'8514HAWAII' '8501CALCURRENT'})

%04/03/18
% Loop_ARGO_float('all','12363',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9660',{'8514HAWAII' '8501CALCURRENT'})
% 
% Loop_ARGO_float('all','0571',{'8514HAWAII' '8501CALCURRENT'})

%04/11/18  Initial QC for recently deployed floats.
% Loop_ARGO_float('all','12779|12784|12782',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12709|12741|12748|12370',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12782',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0690|0693',{'8514HAWAII' '8501CALCURRENT'})

%04/14/18  Initial QC for recently deployed floats.
% Loop_ARGO_float('all','0569|12552',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12396',{'8514HAWAII' '8501CALCURRENT'})

%04/23/18  Test addition of "GMT" tag to file header section.
% Loop_ARGO_float('all','9650',{'8514HAWAII' '8501CALCURRENT'})
% %Reprocess ALL to pick up (1) updated BBP coeffs, (2) GMT header tag
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT' '7597HAWAII' '7593HAWAII'})

%04/24/18  Reprocess 
% Loop_ARGO_float('all','8501CAL','')
% Loop_ARGO_float('all','12386',{'8514HAWAII' '8501CALCURRENT'})

%04/26/18  Reprocess bad pH sensors identified.  :-(
% Loop_ARGO_float('all','9602|9666',{'8514HAWAII' '8501CALCURRENT'})

%04/30/18  Reprocess bad pH sensors 
% Loop_ARGO_float('all','9666|12371|12379',{'8514HAWAII' '8501CALCURRENT'})
% reprocess 9096 with updated wavelength range for fit (max wl set to 238).
% Loop_ARGO_float('all','9096',{'8514HAWAII' '8501CALCURRENT'})

%5/1/2018 Float QC
% Loop_ARGO_float('all','9662|9749|9757|9631|0037',{'8514HAWAII' '8501CALCURRENT'})

%5/7/2018 Reprocess 6967SOATLANTIC and 7558ETNP; gdac rejection for CHLA
%(datamode was set to 'A' when there was no chla data for certain cycles)
% Loop_ARGO_float('all','6967',{'8514HAWAII' '8501CALCURRENT'})

%5/10/2018
% Loop_ARGO_float('all','9096',{'8514HAWAII' '8501CALCURRENT'})

%5/14/2018 QC updates
% Loop_ARGO_float('all','12742|7613|12729',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9094|0570|0507|7652|9657|9652',{'8514HAWAII' '8501CALCURRENT'})

% 5/15/2018 QC updates
% Loop_ARGO_float('all','12366|12379',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12366|12379|12559|0566|12745|12733|12723|9752',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12701',{'8514HAWAII' '8501CALCURRENT'})

% 5/16/2018 QC updates
% Loop_ARGO_float('all','9125|12575',{'8514HAWAII' '8501CALCURRENT'})

% 5/17/2018 QC updates
% Loop_ARGO_float('all','12361',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0691|0510|0565|0692|12573',{'8514HAWAII' '8501CALCURRENT'})

% 5/21/2018 QC updates
% Loop_ARGO_float('all','9650|12769|12545|12370|12779|12784',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12782',{'8514HAWAII' '8501CALCURRENT'})

% 5/23/2018 Added 12396 to bad sensor list for NO3; 6960 bad profile
% registration
% Loop_ARGO_float('all','12396',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','6960',{'8514HAWAII' '8501CALCURRENT'})

% %5/24/18 profile 51 added to BSL for 0569
% Loop_ARGO_float('all','0569',{'8514HAWAII' '8501CALCURRENT'})

%6/4/2018 New QC for newer problem floats (missing profiles...)
% Loop_ARGO_float('all','9766',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12363|0688',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0688|12702|12363',{'8514HAWAII' '8501CALCURRENT'})

% 6/5/2018 Reprocess 12779 after adding to BSL for pH 
%          Global reprocess to implement spike test!
% Loop_ARGO_float('all','12779',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12728|12717',{'8514HAWAII' '8501CALCURRENT'})

% 6/6/2018 Reprocess 0688, 0690, 0691, 0692 after bug fix in SageO2 for
% calculating gains based on WOA.  Only affects these floats.  All other
% NAVIS cal'd to bottle.
% 12361 update to QC while testing sage GLT version
%  Loop_ARGO_float('all','0688|0690|0691|0692',{'8514HAWAII' '8501CALCURRENT'})
%  Loop_ARGO_float('all','12361',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12748',{'8514HAWAII' '8501CALCURRENT'})

 % June 8,2018 Rerunning floats to pic up bad qc flag for backscatter
 %Loop_ARGO_float('all','8501|0691|12371|12541',{'8514HAWAII' '8501CALCURRENT'})
%  Loop_ARGO_float('all','8501|0691',{'8514HAWAII' '8501CALCURRENT'})

% 6/18/18 Updates to float QC for 0571 and 9749
% Loop_ARGO_float('all','0571|9749',{'8514HAWAII' '8501CALCURRENT'})

% 6/28/2018
% Loop_ARGO_float('all','12398|12758|12754|12787',{'8514HAWAII' '8501CALCURRENT'})

%7/6/18 TM, 9101 PH sensor bad profiles toward end of record.
% Loop_ARGO_float('all','9101',{'8514HAWAII' '8501CALCURRENT'})


%7/24-25/18 TM, float qc
% Loop_ARGO_float('all','12741',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12545|12575|12573|12543',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12549|12381',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12542|12371',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12723|12729|12386',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12575',{'8514HAWAII' '8501CALCURRENT'})

%8/6/2018 TM, float QC
% Loop_ARGO_float('all','9655|9031|12382',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9018|9313|9602',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9632',{'8514HAWAII' '8501CALCURRENT'})

%8/7/2018 TM, float QC
% Loop_ARGO_float('all','9632|9744',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9645',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0567|9254|9750',{'8514HAWAII' '8501CALCURRENT'})

%8/8/2018 TM, float QC
% Loop_ARGO_float('all','9095|9744',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9095',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9091',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9101',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12369|12372|12537|9091|9101',{'8514HAWAII' '8501CALCURRENT'})

%8/9/2018 TM, float QC
% Loop_ARGO_float('all','12380',{'8514HAWAII' '8501CALCURRENT'}) %added to BSL for large nitrate drift at cycle 33
% Loop_ARGO_float('all','12388|12744',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','9666',{'8514HAWAII' '8501CALCURRENT'}) %manual flagging of some bad salinity (cycles 27, 51, 63, 64).  Not 
% Loop_ARGO_float('all','12755',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','12755|9666|9646|12540|7614|12744|12388',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','12734|12558|12541|8204|8501',{'8514HAWAII' '8501CALCURRENT'}) 

%8/13/18 TM
% Loop_ARGO_float('all','12543',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0569|12757',{'8514HAWAII' '8501CALCURRENT'}) %BR file transfer hiccup.
% Loop_ARGO_float('all','9101|12543|9652',{'8514HAWAII' '8501CALCURRENT'})

%8/15/18 TM
% Loop_ARGO_float('all','9092',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9275|12390|9637|9642|9762|9634|12730|12757|12781|12736',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9668|9660',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9265',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0690',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0949',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12709',{'8514HAWAII' '8501CALCURRENT'})

%8/23/18 TM
% Loop_ARGO_float('all','12783|12788',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12744|12769',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9749',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0948|0949',{'8514HAWAII' '8501CALCURRENT'})

%8/27/18 TM
% Loop_ARGO_float('all','0569',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12784',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12741',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','12744|12784',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9631',{'8514HAWAII' '8501CALCURRENT'})

%8/28/18 TM
% Loop_ARGO_float('all','0948',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9637',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9757',{'8514HAWAII' '8501CALCURRENT'})

%09/04/2018
% Loop_ARGO_float('all','0948|0949',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9662|9757',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0037',{'8514HAWAII' '8501CALCURRENT'})

%09/5/2018
% Loop_ARGO_float('all','12779',{'8514HAWAII' '8501CALCURRENT'}) %update QC.  CSUMB float, "polar otter"

%09/6/2018
% Loop_ARGO_float('all','0569|0510|0566|0571|0507',{'8514HAWAII' '8501CALCURRENT'}) %NAVIS float QC
% Loop_ARGO_float('all','9101|9092',{'8514HAWAII' '8501CALCURRENT'}) %9101 added to BSL.  pH backup battery issue -flagging

%09/10/2018 TM
% Loop_ARGO_float('all','9662',{'8514HAWAII' '8501CALCURRENT'}) %BSL 93 onward for pH
% Loop_ARGO_float('all','0690',{'8514HAWAII' '8501CALCURRENT'}) %updated f(P) correction (JP)

%09/10/2018 TM
% Loop_ARGO_float('all','9662',{'8514HAWAII' '8501CALCURRENT'}) %manual QC, cycles 86 and 58
% Loop_ARGO_float('all','9092',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9275',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9757',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','9637',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0569',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0948|0949',{'8514HAWAII' '8501CALCURRENT'})
% Loop_ARGO_float('all','0565',{'8514HAWAII' '8501CALCURRENT'})
% %Lots more navis QC, will get updated with global reprocess
% %Global reprocess prior to sept2018 snapshot
% Loop_ARGO_float('all','',{'8514HAWAII' '8501CALCURRENT'})
% 
% 
% Loop_ARGO_float('all','12702',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','9092',{'8514HAWAII' '8501CALCURRENT'}) 
% 
% Loop_ARGO_float('all','0507',{'8514HAWAII' '850floatvi1CALCURRENT'}) 

%10/01/2018
% Loop_ARGO_float('all','7641',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','8374',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','0276',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','0949',{'8514HAWAII' '8501CALCURRENT'}) 

%10/08/2018
% Loop_ARGO_float('all','12573',{'8514HAWAII' '8501CALCURRENT'}) 

%10/09/2018  Changed flagging from 'bad' to 'questionable' for pH and NO3
%data for floats with corresponding bad O2 values.
% Loop_ARGO_float('all','0565|9274|9659|12551|12573',{'8514HAWAII' '8501CALCURRENT'}) 

%10/11/2018  Modify the SF for FDOM channel on 0948 and 0949.  Swap them
%back to original manufacturer scale factors, and then divide the SF for
%0949 by 2 (per Ken, see email 10.10.18)
% Loop_ARGO_float('all','0948|0949',{'8514HAWAII' '8501CALCURRENT'}) 

%10/30/2018 12736 and 12781 added to bad sensor list for pH
% Loop_ARGO_float('all','12736|12781',{'8514HAWAII' '8501CALCURRENT'}) 

%11/01/2018.  GLOBAL REPROCESS TO INCORPORATE UPDATES TO CANYON AND LIR
%ALGORITHMS!!! (1) UPGRADE TO CANYON-B FROM CANYON.  (2) UPGRADE TO FIX BUG
%IN TRAINING DATASET FOR LIRS (SEE EMAILS FROM BRENDAN CARTER ON 10/31/18).
% ALKALINITY ESTIMATES AND ASSOCIATED CARBON PARAMETERS IN ODV FILES WILL
% BE UPDATED!  *** LIR update was done first
% Loop_ARGO_float('all','',{' '})

%11/05/18 reprocess 0948
% Loop_ARGO_float('all','0948',{' '})
% Loop_ARGO_float('all','12380',{' '})

% 11/9/2018
%reprocess floats with Aandera optode SVU coeffs.  47 floats.
% Loop_ARGO_float('all','9125|9668|12371|9657|12552|9761|9757|9766|12542|12723|12755|12778|12781|12742|12734|12757|12747|12779|12784|12769|12758|12709|12702|12754|12768|12736|12730|12741|12787|12733|12729|12748|12745|12782|12744|12701|12700|12717|12728|12792|12783|12780|8474|12696|12883|12788|12881',{'8514HAWAII' '8501CALCURRENT'}) 
% Loop_ARGO_float('all','12542|12701',{'8514HAWAII' '8501CALCURRENT'})

%11/15/18
% Loop_ARGO_float('all','9642',{' '}) %reprocess 9642, bad optode.

%12/17/18  GLOBAL REPROCESS TO INCORPORATE MANY UPDATES!  O2 GAIN DRIFT,
%ISCHANGE, UPDATED QC-MATRIX MATH, GUI INTERFACE ENHANCEMENTS.  TESTING
%LOCALLY ON STERNA FIRST!  EXCLUDING A FEW FLOATS THAT NEED LAST MINUTE
%UPDATES.  WILL CHECK PROCESSING SUCCESS TOMORROW AND FINISH LEFTOVER FLOATS.
% Loop_ARGO_float('all','',{'8374HOT' '7622HAWAII' '7657SOOCN' '7620SOOCN' '7619SOOCN' '7552SOOCN' '6968SOOCN'})

%12/31/18.  Some QC updates to floats that popped from under ice.  9657,
%12783, 12361, 12388, 12381

% 1/2/19 TM testing cals for some of the recent navis floats.

% 1/10/19  Update code for processing 12768 (grab position info from
% % 000.msg since no position info ever for this float.
% Loop_ARGO_float('all','12768',{' '}) 

%1/14/18 bad pH added to BSL
% Loop_ARGO_float('all','12398|12758|12386',{' '})
% Loop_ARGO_float('all','0690',{' '}) 

% 1/16/19: global reprocess to move from using phcalc_jp to phcalc
% (reflection of what is in Argo documentation.  difference in these two
% pieces of code amounts to difference in calculated pH on the order of 0.1
% millipH)
% Loop_ARGO_float('all','',{' '})

%1/22/19, reprocess after alternate dirs were activated.
% Loop_ARGO_float('all','12768|12787|0886',{' '}) 



% %2/4/19; rsync hiccup reprocessings
% Loop_ARGO_float('all','9762|12381|12549',{' '}) 
% 
% %2/7/19; rsync issues
% Loop_ARGO_float('all','12549|12881',{' '})
% %2/11/19 qc update on bsl for 12549, update qc for 9634, marked prof 1 for
% %pH questionable for 0566
% Loop_ARGO_float('all','12549',{' '})
% Loop_ARGO_float('all','9634',{' '})
% Loop_ARGO_float('all','0566',{' '})
% %12702 cycles 36, 37 bad salinity
% Loop_ARGO_float('all','12702',{' '})
% %initial qc for recently deployed floats:
% Loop_ARGO_float('all','12708|12889|12727|12711|0886|0689',{' '})

%2/12/19 9265 removed from bad list.  Only cycle 85 remains bad.

% 2/13/19 rerun 0412 after deleting files on sirocco to see if qf >34 go
% away.  pH failures on 9650 and 9099
% Loop_ARGO_float('all','0412Hawaii',{' '})
% Loop_ARGO_float('all','9650',{' '})
% Loop_ARGO_float('all','9099',{' '})
%initial QC for new floats
% Loop_ARGO_float('all','12714|12886|11090',{' '})

%2/25/19 Reprocess; Btransfer error.
% Loop_ARGO_float('all','12889|12782',{' '})
%Float QC
% Loop_ARGO_float('all','7652|0037|9313|9125|9668|9646|9655|9662|9757|12575|9637|9744',{' '})

%2/26/19-2/28/19 Float QC
% Loop_ARGO_float('all','0949|0571|9660|9630|9762|9752|12545|12559|9766|12537|12541|12558|12552|9750|8501|12371|12542|12388|12363|12380|12396|8204|12729',{' '})
% Loop_ARGO_float('all','12782',{' '})
% Loop_ARGO_float('all','12723|12742|12540|12733|12744|12734|12755|12730|12781|12736|12757|12784|12769|12709|12741|12748|12370|12787|12696|12881',{' '})

%3/4/19 Float QC
% Loop_ARGO_float('all','0569',{' '})
% Loop_ARGO_float('all','12883|12700|12747|12778|0569|0691|0692|0688|12745|12369|9600|9652|9602|0510|12386|12754|12543|9631|9275|9749|9657|12361|12758|12768|12701|12702|12381',{' '})
% Loop_ARGO_float('all','12381',{' '})
% Loop_ARGO_float('all','12749|12739|12688',{' '}) %initial QC
% Loop_ARGO_float('all','12702|12543|9099',{' '}) %floatQC

%3/10/19 Global reprocess prior to March, 2019 snapshot
% Loop_ARGO_float('all','',{' '})

%4/3/19 reprocess to pick up files in limbo from Dana's server issue (might
%not all be resolved?)
% Loop_ARGO_float('all','12788|12741|9757|12696|12727|12709|9634|12757|12700',{' '})

%4/10/19 More file transfer issues (Dana fixing on his end).
% Loop_ARGO_float('all','12727|12709|9634|12757|12700',{' '})

%4/15/19  BR transfer matchtest errors.  Likely due to updated msg file
%with more data.  Reprocess first, then try BR transfer.
% Loop_ARGO_float('all','12708|12712',{' '})

%4/30/19 Initial float qc for ANDREXII floats
% Loop_ARGO_float('all','12707|12752|12786|12884|12879|12880',{' '})

%06/03/19 Initial QC for I06S floats
% Loop_ARGO_float('all','12888|12878|12892|12885|12882|0888|0889',{' '})

%6/4-6/19 QC updates
% Loop_ARGO_float('all','12549|9650|12708|12889|12727|12711|0689|0886|0566|9094|0690|0568|12714|12886|11090|12398',{' '})
% Loop_ARGO_float('all','6091|9645|7613|9092|9018|7614|12782|7652|0037|9646|9655',{' '})
% Loop_ARGO_float('all','9094|6968|7619|7620|7657|9092|9125|9744|9668|9662|9757|12575|9637|9660',{' '})
% Loop_ARGO_float('all','9762|9752|12545|12559|9766|12541|9750|8501',{' '})

%6/11-14/19 float qc
% Loop_ARGO_float('all','12371|12542|12388|12363|12396|12380|8204|12729|12757',{' '})
% Loop_ARGO_float('all','12723|12742|12540|12733|12744|12734|12755|12730|12781|12736|12784',{' '})
% Loop_ARGO_float('all','12757|12769|12709|12748|12741|12787|12696|12881',{' '})
% Loop_ARGO_float('all','12778|0569|0691|0692|0688|12754|12369|9600|9652|9602',{' '})
% Loop_ARGO_float('all','',{' '}) % Global reprocess to implement the
% update to Calc_O2_4ARGO (coefficient updates to match exactly what is in the Argo processing manual).
% Loop_ARGO_float('all','0510|12386|12543|9631|9275|9657|12361|12758|12768|12701',{' '})
% Loop_ARGO_float('all','0569|0692|12702|12381|12749|12739|12688|12537|12366|12779|9313',{' '})
% Loop_ARGO_float('all','12712|12631|12652|12370|12883|12752|12879|12707|12884',{' '})
% Loop_ARGO_float('all','12786|12880|9265|12552|12558|12747|12700|12888|12878|12885|12882',{' '})

%JP 06/19/19 float QC for NON SOCCOM floats
% Loop_ARGO_float('all','8474EQPACW|8374HOT|12728EQPACE|12792EQPACW|12788EQPACW|12783EQPACW|0948STNP2',{' '})

% 6/24/19 TM - Added questionable profiles to BSL for pH sensors exhibiting
% variable "pump offset"
% Loop_ARGO_float('all','12714|12768|12787|12381|12379|12739|12880|12742',{' '})
% a few fixes from JP:
% Loop_ARGO_float('all','12701|12559|9313',{' '})
% Some salinity and nitrate manual flagging, TM:
% Loop_ARGO_float('all','7613|9655|9762|12543|12545|12714|12757',{' '})

%7/8/19 Adjust NO3 wavelength offset from 210 to 208.5 for 17307HOT and
%17267HOT
%Loop_ARGO_float('all','17267|17307',{' '})

%08/28/19 JP Reprocessing NAVIS float with updated QC adjustments
%TM additional floatQC
% float_list_str = ['0037SOOCN|0507SOOCN|0510SOOCN|0565SOOCN|0566SOOCN|', ...
%     '0571SOOCN|0688SOOCN|0689SOOCN|0690SOOCN|0691SOOCN|0692SOOCN|',...
%     '0886SOOCN|0888SOOCN|0889SOOCN|0948STNP|0948STNP2|0949STNP'];
% Loop_ARGO_float('all',float_list_str,{' '})
% Loop_ARGO_float('all','12386|12543|9275|12361|12758',{' '})
% Loop_ARGO_float('all','112714|9631|9657|12537|12779|9749',{' '})

% 9/3/2019 TM, 
%Updated manual flags for 12778 (cycle 6 at depth)
% Loop_ARGO_float('all','12778',{' '})
%Float QC for soccom and non-soccom floats:
%Loop_ARGO_float('all','12472|12631|12752|12884|12786|12780|12783|12788|12792|12728|8474',{' '})

%09/06/19 JP Reprocessing APEX floats which had last qc date of 06/11/19 with updated QC adjustments
%  float_list_str = ['12396|12729|8204|12388|12723|12742|12540|12733|12736|12784|12781|', ...
%      '12696|12881|12741|12709|12769|12778'];
% Loop_ARGO_float('all',float_list_str,{' '})

%09/26/19 Final updates before global reprocess for snapshot
%  float_list_str = ['17267|12747|9744|12708|17307|12551|9752|12885|12573'];
%  Loop_ARGO_float('all',float_list_str,{' '})
 
 %09/26/19 Global reprocess
 %Loop_ARGO_float('all','','')
 % GLOBAL REPROCESS
 
 %09/26/19 Final updates before global reprocess for snapshot

%  % 10/04/19 JP Reprocess to fix one aspect of GPS date bug
% Loop_ARGO_float('all','9634',{' '})

 % 10/08/19 JP 17534EqPacE deployed
%Loop_ARGO_float('all','17534',{' '})

 % 10/31/19 JP Running 9634 after deleting all instances on SIROCCO and
 % coding a fix for GPS weeknumber rollover bug
% Loop_ARGO_float('all','9634',{' '})

% 11/04/19 New EqPac floats added 
% Loop_ARGO_float('all','17350|17862|18608',{' '})

% 11/05/19 JP Rerun 12888 after bias cor calc fix to improve pCO2 data
% coverage
% Loop_ARGO_float('all','12888',{' '})

% 11/05/19 JP Rerun 12880 to pick up bad qf fro BBP cycle 21
% coverage
% Loop_ARGO_float('all','12880',{' '})

%11/18/19 first QC for some equPac floats and a few othe QC updates
%  float_list_str = '12889|17267|17307|12888|17862|17350|17534|18608|12892';
%  Loop_ARGO_float('all',float_list_str,{' '})

 %11/20/19 rerunning to see if it fixes QC flag bugs in cycle 1 for NO3
 %after deleting text files on sirocco
 %Loop_ARGO_float('all','6967',{' '})
 
  %12/05/19 rerunning to includen cycle 1 in QC adj. (adjustments were starting on cycle 2 & cuase br file creation to crash) 
 %after deleting text files on sirocco
 %Loop_ARGO_float('all','17307',{' '})
 
 %12/05/19 rerunning floats that Raphaëlle noticed had bad profiles
 %after deleting text files on sirocco
 %Loop_ARGO_float('all','12755|6975|6967',{' '})
 
  %12/05/19 1st SOCCOM deployment for the season!
 %after deleting text files on sirocco
 %Loop_ARGO_float('all','18721',{' '})
 
 %12/11/19 Rerun both float now that lab k0's are available
 %after deleting text files on sirocco
 %Loop_ARGO_float('all','18320|18721',{' '})

 %12/13/2019 JP, 0569 was breaking masg parsing code due to partial hex line
 %of cp data. re-did hex parsing code to better deal with partial hex data
 %lines
  %Loop_ARGO_float('all','0569',{' '})

 %12/17/19 12892 popped up rerun to check pH after QC - is funky
 %after deleting text files on sirocco
 %Loop_ARGO_float('all','12892',{' '})
 
  %12/19/19 12749 popped up rerun to check pH after QC - is funky
 %after deleting text files on sirocco
 %Loop_ARGO_float('all','12749',{' '})
  %Loop_ARGO_float('all','0689',{' '})
 
 % *** GLOBAL REPROCESS ***
 % 12/20/19 Global reprocess after final QC of all floats by JP
%  Loop_ARGO_float('all','','')
 
 % 01/02/20 1st Marai deployment for the season!
%  Loop_ARGO_float('all','18299',{' '})
 
 % 01/03/20 12739 under ice but redue. txt file cycles 1-13,35, but more
 % availoable
 %Loop_ARGO_float('all','12739|12688',{' '})
 
 % 01/02/20 1st Marai deployment for the season!
%Loop_ARGO_float('all','12733',{' '})

% 01/09/2020 JP updating some changes made to bad sensor list for pH
% sensors based on high calulated pCO2 values
%Loop_ARGO_float('all','12549|12537|12388|12739',{' '})

% 01/09/2020 JP updating some changes made to bad sensor list for pH
% sensors based on high calulated pCO2 values
%Loop_ARGO_float('all','12786|9265|12386|12708|12781|12787|12727',{' '})
%Loop_ARGO_float('all','12786|12892',{' '})

% 01/13/2020 1st Investigator float
%Loop_ARGO_float('all','18852',{' '})

% 01/14/2020 Updating under ice floats
%Loop_ARGO_float('all','12559|12381|11090|0568|12366',{' '})
%Loop_ARGO_float('all','12688|12701|12714|12758|12787|12878',{' '})

% 01/21/2020 Forcing processing for 18864. NO3 sensor not activated so
% there will be no *.isus msg files for cycle 1 & 2
%Loop_ARGO_float('all','12559|12381|11090|0568|12366',{' '})
% Loop_ARGO_float('all','18864',{' '})

% 01/22/2020 Reprocessing for 18417. k2 was missiing a negative sign
%Loop_ARGO_float('all','18417',{' '})

% 01/30/2020 force processing of 1st porfiles 18643, 18829
%Loop_ARGO_float('all','18829',{' '})

% 02/13/2020 force processing of mostly complete msg files
%Loop_ARGO_float('all','12886',{' '})

% 02/27/2020 threw a processing error
% Loop_ARGO_float('all','0889',{' '})

% 03/03/2020 ph went bad starting 12778 cycle 50
%Loop_ARGO_float('all','12778',{' '})


% 03/09/2020 Initializing float deployed off mointerey today
%Loop_ARGO_float('all','18340',{' '})

% 03/10/2020 Reprocess after 1st QC of CUSTARD cruise floats
%Loop_ARGO_float('all','18242|18320|18721|18545|18771|18098',{' '})

% 03/12/2020 Reprocess after 1st QC of I7S cruise floats
%Loop_ARGO_float('all','18299|18994|18821|18739|17898|18013|18864',{' '})

% 03/15/2020 Forcing 1st processing od 18081
%Loop_ARGO_float('all','18081',{' '})

% 03/15/2020 Forcing 1st processing od 18081
%Loop_ARGO_float('all','18097',{' '})

%03/21/2020 isus files for 18340 showing power failure & no data. Added
% 18340 to the if statement  in Process_APEX_ floatnfor floats which will
% never have usable isus files
%Loop_ARGO_float('all','18340',{' '})

%03/22/2020 run TPOS float which have just been deployed
%Loop_ARGO_float('all','18114|18709|18820',{' '})

%03/29/2020 initial QC for 3 of the Thwaites cruise floats
%Loop_ARGO_float('all','18643|18829|18169',{' '})
% Loop_ARGO_float('all','0507SOOCN',{' '})

%03/29/2020 TM, reprocess after flagging of bad S,O profile at cycle 128
%Loop_ARGO_float('all','0691SOOCN',{' '})

%04/07/2020 JP, reprocess after change pH on BSL and removing files from
%sirocco
%Loop_ARGO_float('all','12559SOOCN',{' '})
% Loop_ARGO_float('all','9766SOOCN',{' '})% updated WL offset =209

%04/08/2020 TM, reprocess to repopulate missing NO3 on cycles 100-101
% Loop_ARGO_float('all','12543SOOCN',{' '})
% Loop_ARGO_float('all','0566SOOCN',{' '})
% Loop_ARGO_float('all','0510SOOCN|0566SOOCN|12543SOOCN|0068ROSSSEA|9615SOOCN|9762SOOCN|9313SOOCN|7613SOOCN',{' '})

%04/09/2020 JP, deleted *cal.mat & reprocessed to pick up change in NO3 WL
%Offset  to 209
%Loop_ARGO_float('all','12733SOOCN',{' '})
% Loop_ARGO_float('all','12788EqPacW',{' '})

% 04/13/2020 TM initiated global reprocessing prior to Spring snapshot.
%Loop_ARGO_float('all','',{' '})

% 04/16/2020 JP updated QC corrections for 18340CalCurrenet, stong pH pump
% offset
% Loop_ARGO_float('all','18340',{' '})

% 04/17/2020 TM bug in sage - if qc matrix starts on cycle other than 1,
% missing values added to QC file for early cycles.  Manually forcing qc
% node of 1 for select floats in this category.  Need to implement a
% software fix to prevent this.
%  Loop_ARGO_float('all','12369SOOCN',{' '})
%  Loop_ARGO_float('all','12370SOOCN',{' '})
%  Loop_ARGO_float('all','12747SOOCN',{' '})
%  Loop_ARGO_float('all','17307HOT',{' '})
%  Loop_ARGO_float('all','12783EQPACW|17534EQPACE',{' '})
%  Loop_ARGO_float('all','18417SOOCN',{' '})
%   Loop_ARGO_float('all','7672HAWAII|8514HAWAII|7618CALCURRENT|6960ETNP|6391BERMUDA|8486HAWAII',{' '})
% % Also, 12781 missing cycles 74 and 77 but msg files are there.  reprocess
% %  Loop_ARGO_float('all','12781SOOCN',{' '})
% 

%4/19/20 reprocess 12733, BR bug because older msg file came through in
%full
% Loop_ARGO_float('all','12733SOOCN',{' '})

%4/24/20 modified qc for float 17862eqpace.  Bad profile not getting
%flagged, and issue with qc adjustment on first profile (LIR), flagged
%cycle 1 S as bad
% Loop_ARGO_float('all','17862EqPacE',{' '})

%4/28/20 reprocess 5426drakepass after receiving updated oxygen cal coeffs
%from dana swift.  Also cycle 185 sal added to bad sensor list.
% Loop_ARGO_float('all','5426',{' '})

%5/4/2020 TM reprocess floats 9600 and 9650 after updating sage qc matrix
%(there was a NaN in matrix entry for NO3, causing empty NO3 profile...)
% Also, additional surface spikes from Seth
%  Loop_ARGO_float('all','9600|9650',{' '})
%  Loop_ARGO_float('all','7613|9313',{' '})
%  Loop_ARGO_float('all','9662|9766',{' '})

%5/11/2020 TM, reprocess 9634 after fix to code dealing with gps bug
% reprocess 12542 after fix to code to deal with flbb dropout
% reprocess floats 12768 and 12783 after fix to code to deal with missing
% lat/lon on the first profile.
%  Loop_ARGO_float('all','12542',{' '})
%  Loop_ARGO_float('all','9634',{' '})
%  Loop_ARGO_float('all','12783|12768',{' '})

% 6/11/2020 TM, testing Process_APEX_float.  

% 6/23/20 TM, there was a hiccup in the 8am processing.  Rerun select
% floats from that run.
%Loop_ARGO_float('all','9645|7567|7597|9260|12363|12390|12398|12472|12542|12549|12551|12631|12652|12708|12712|12745|12755|12758|12786|12879|12888|',{' '})
%Loop_ARGO_float('all','12792',{' '})
%Loop_ARGO_float('all','18643',{' '})

% 6/29/20 TM, QC updates
%Loop_ARGO_float('all','12701',{' '})
% no qc update needed on these, but changed qclist date to reflect they had
% been reviewed: 7619, 7620, 12712, 9313, 12380

% 6/30/20 TM, QC updates
% Loop_ARGO_float('all','9650',{' '})
% Loop_ARGO_float('all','6091',{' '})


%7/6/20
% Loop_ARGO_float('all','9602|9265|7614|8204|9125',{' '})
%7/7/20
% Loop_ARGO_float('all','9637|9634|9632|9630|12372|12396',{' '})
% Loop_ARGO_float('all','9762',{' '}) %salinity flagging
% Loop_ARGO_float('all','12372',{' '}) %salinity flagging and update NO3 on BSL

%7/9/20 Float QC
% Loop_ARGO_float('all','9645|9646|9655|9657|9662|9668',{' '}) 
% Loop_ARGO_float('all','9642|12379|12537|12552|12575|12733|12736|12778|18417|18601|18796',{' '}) 

%7/13/20 Float QC
% Loop_ARGO_float('all','12575|9744|9750|9752',{' '}) 
% Loop_ARGO_float('all','9757|9762|9766|11090',{' '}) 
% Loop_ARGO_float('all','12366|12381|12386|12388',{' '}) 
%7/14/20 Test bug fix JP implemented for NAVIS parser 
% Loop_ARGO_float('all','0571',{' '}) 
%7/15/20 Bad profile, detected by NO3 CANYB test audit! 
% Loop_ARGO_float('all','7550',{' '}) 
% Loop_ARGO_float('all','6381',{' '}) 
 %7/16/20 Float QC.  Lots of pump offset issure on pH.
%  Loop_ARGO_float('all','18771|18299|18098|17898|18739|18821|18852|18013|18864',{' '}) 

%7/16/20
% Loop_ARGO_float('all','9660',{' '}) %flagged some bad profiles, ctd
% failing?  this float stuck surface. 
% 
%7/23/20 TM SAGE QC
% Loop_ARGO_float('all','18862|18994|18829|18820|18796|18721|18709|18643|18608|18545|18320|18169|18161',{' '}) %sage qc
%7/24/20 TM SAGE QC
% Loop_ARGO_float('all','12883|12882|12881|12880|18862|12787|12878',{' '}) %sage qc
% Loop_ARGO_float('all','12741|12742|12749',{' '}) %sage qc

%7/26/20 SAGE QC
% %7/27/20
% Loop_ARGO_float('all','12631|12688|12696|12700|12702|12709|12711|12714|12754|12757|12758|12769|12778',{' '}) %sage qc
% Loop_ARGO_float('all','12757|12575|12779|12781|12784|12886|12888|12889|12892|18081|18069|18097|18110|12369|12370',{' '}) %sage qc
% Loop_ARGO_float('all','12537|12540|12543|12545|12549|12552|12558|12559|12573',{' '})

% %7/28/20
% % Loop_ARGO_float('all','11090',{' '}) %bad pH cycle 47
% Loop_ARGO_float('all','12542',{' '}) %bad salt, 122-125, failed?
% Loop_ARGO_float('all','12542|12747|18417|12768|12551|12744|18321|18677|18431|18764',{' '}) 
% %7/30/20 flags
% Loop_ARGO_float('all','18242|9749',{' '})
% %7/31/20
% Loop_ARGO_float('all','9631',{' '}) %sage QC; pump offset, pH qf=3
% Loop_ARGO_float('all','12744',{' '}) 
% Loop_ARGO_float('all','9642|12382',{' '}) 
% 
% 
% %8/2/20
% Loop_ARGO_float('all','0568|9091|9101|9096|6091|12575|9637',{' '}) 
% Loop_ARGO_float('all','12372',{' '}) %nitrate cycle 1-26 qc=3.  bsose anomaly and poor diagnostics; fails cycle 27-

%8/4/20
% sage QC
% Loop_ARGO_float('all','17862|17534|17350|17307|17267',{' '})

%8/10/20; 8/11/20
% Sage QC
% Loop_ARGO_float('all','9260|9761|9092|9091|9031|9631|0889|12361|0689|12398',{' '})
% Loop_ARGO_float('all','7614',{' '})
% 
% %8/12/20 sage qc, non-soccom
% Loop_ARGO_float('all','18601|18820|18709|18114|18608|17862|17350|17534|12792|12788|12783|12780|12728|12717|12712|12652|12631|12472|8474|8374',{' '})

 % *** GLOBAL REPROCESS ***
 % 8/26/20 Global reprocess after fleetwide QC 
 % upgrades to processing include: (1) QC header description of NPQ-data
 % flagging scheme, and (2) modification to pH ref temp from 2degC to
 % Temp@1500m (for application of TCOR)
 %Loop_ARGO_float('all','','')
 
 % 09/25/20 % JP up dating float with correct no3 cal file
%  Loop_ARGO_float('all','12758',{' '}) 

 % *** GLOBAL REPROCESS ***
 % 10/8/20 Clean global reprocess -- local reprocess & file transfer error to chem?
%  Loop_ARGO_float('all','','')


%10/26/20 TM EQPAC,CalCurrent recent deployments - initial QC
% Loop_ARGO_float('all','19139',{' '}) 
% Loop_ARGO_float('all','19412|19090|19017|19644',{' '}) 


%10/27/20 TM full reprocess to incorporate change to INFO.gps held in each
%  Loop_ARGO_float('all','','')

%11/17/20 BBP QC
% Loop_ARGO_float('all','7552|12885|0691|12744|12892|0037|12363|9645|12573|9602|0276|0068|12379|9096|0949',{' '})
% Loop_ARGO_float('all','0688|0888|0949|0689|0507|0568|0566|0569|12783|0690|9645',{' '})

%12/13/20 TP full reprocess to incorporate change to UTF-8 ASCII encoding.
%   for all generated text files .Matlab 2020 fopen writes now default to
%   UTF-8 vs system default prior. UTF-8 allows for better cross platform sharing
% Loop_ARGO_float('all','','')

% %12/14/20 TM float QC
% Loop_ARGO_float('all','18821|18852|18739|17898|18098',{' '})

%12/17/20 JP float QC
% flt_list_str = ['19875CalCurrent|19644EqPacE|19017EqPacE|19412EqPacW|19090EqPacW|', ...
%     '19139CalCurrent|18601EqPacE|18764SOOCN|18431SOOCN|18677SOOCN|18069SOOCN|18110SOOCN|',...
%     '18796SOOCN|18169SOOCN|18829SOOCN|18643SOOCN|18417SOOCN|18864SOOCN'];
% Loop_ARGO_float('all',flt_list_str,{' '})

%12/18/20 TM float QC
%Loop_ARGO_float('all','17350|17862|18608|17534|17267|17307|0889|0888|18601',{' '})

%12/18/20 JP float QC - All pH QC = 8, but I don't think they should be
%trying sirocco flush & reprocess small dura file size because of many
%missing pH samples + small bug in Ik testing set QC for all remaning good
% pH to bad
%Loop_ARGO_float('all','12888',{' '})

%12/19/20 JP float QC batch #2 - The problem kids!
% flt_list_str = ['12882SOOCN|12888SOOCN|12885SOOCN|12707SOOCN|12688SOOCN|',...
%     '0886SOOCN|12889SOOCN|12778SOOCN|12747SOOCN|12700SOOCN|12881SOOCN|', ...
%     '12696SOOCN|12792EQPACW|12780EQPACW'];
% Loop_ARGO_float('all',flt_list_str,{' '})

% %12/21/20 TM happy winter solstice.  :-)  Global reprocess prior to soccom
% %snapshot creation.  Update includes switch to UTF-8 encoding.
% Loop_ARGO_float('all','','')

%12/29/20 TM Odd processing hiccup on the 28th.  Select files on the
%network do not reflect the UTF-8 encoding update.  But, these were not
%present on sterna.  Someone processing data on another machine with old
%code and inadvertently copying to sterna??
% Loop_ARGO_float('all','','')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SYSTEM TRANSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WE ARE NOW USING WMO FOR ALL FILE NAMING.  WE ELIMINATED THE
% BASIN-IDENTIFIERS FROM THE INTERNAL MBARI FLOAT IDS.  ALL FLOATS WERE
% REPROCESSED TO BE INTEGRATED INTO THE NEW SYSTEM ON 3/13/21. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4/29/21 Reprocess 1113 and 1115; testing implementation of the suna temp
% modification (using temperature from HR ctd data, 1 measurement deeper
% than no3 sample for processing no3 data (due to the lag in the two
% sensors based on position on the float)). 
%Loop_ARGO_float('all','un1115|un1113',{' '})

% 4/29/21 change wavelength offset from 210 to 211 due to float no3 reading
% negative at sfc
% Loop_ARGO_float('all','ua19875',{' '})
% NON-soccom TPOS/CALCURRENT 19xxx floats
% Loop_ARGO_float('all','ua19139|ua19090|ua19412|ua19644|ua19017|ua19875|ua19063|ua19050',{' '})
% pH failing (38-39) on 18161, added to BSL; 19875 pH QC to 980 due to pump
% offset, so added to BSL 1- as questionable.
% Loop_ARGO_float('all','ua18161|ua19875',{' '})

%4/30/21 Global reprocess prior to snapshot archive!
% Loop_ARGO_float('all','','')

% 5/3/21 A few more floats with updated QC!! somehow missed initial qc on
% 17328.
Loop_ARGO_float('all','ua17328|ua17307|ua12892|ua12758|ua12786|ua12787|ua12878',{' '})

%5/17/21 Wrong Fdom & chla dark counts from Ryan's API -- he was notified, and cals fixed on our end.  Reprocess:
%Loop_ARGO_float('all','wn1203',{' '})



