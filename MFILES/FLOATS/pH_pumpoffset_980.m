% pH_pumpoffset_980.m
%
% This is a list of floats that have been QC'd to 980m  We are storing this
% list in processing in order to automatically:
% a) mark PH_IN_SITU_TOTAL_ADJUSTED_QC below 980m as QUESTIONABLE
% b) inflate PH_IN_SITU_TOTAL_ADJUSTED_ERROR
% This is a temporary solution to recover data with pump offset issue.  The
% next step will be to implement a more comprehensive adjustment scheme
% that 'fixes' the whole profile, at which point we will remove this
% temporary solution.  
% File created 8/25/2023
% Last updated 8/25/2023
%
% 11/29/23 2903869  added for NG
% 12/12/23 5906583  added by TM; variable but some cycles >0.01 offset

% pH_pumpoffset_980_floats = [5905991;
% 5906003];

pH_pumpoffset_980_floats = [5905991;
5906003;
5906000;
5905636;
5906033;
5906031;
5906551;
5906207;
5906215;
5906314;
5906495;
5906295;
5906317;
5906546;
5906313;
5906493;
5906488;
5906499;
5906498;
5906496;
5906294;
5906494;
5906562;
5906497;
5906550;
5906526;
5906468;
5906523;
5906554;
5906553;
5906542;
5906531;
5906556;
5906519;
5906521;
5906529;
5906552;
5906517;
5906525;
5906470;
5906559;
5906548;
5906560;
5906561;
5906577;
5906574;
4903587;
5906570;
5906581;
2903869;
2903472; % JP ua21267 15-20 milli pH
5907055; % JP ua21286 50 milli pH
5907051; % JP ua21535 10 milli pH
5906583; % TM ua21302 > 10 milli pH, variable
3902556; % NG un1516 NAVIS OFFSET!!!!!!!  Checked again 6/13 - still offset!
4903748; % JP ua21939 03/07/24 
2903867; % NG ua20265 03/21/24
6990582; % NG ua22751 03/21/24
6990583;
1902650; % LG ua20075 04/26/24
1902494; % NG un1557 NAVIS!? 04/30/24
5906534; %TM ua20528, gulf of alaska. 5/14/24.  Was previously qc'd to 980 but on BSL as '3'.  I think ok to incorporate here.
2903862; %TM ua21285
5906312; %LG, TM ua19118 failed optode and pH pump offset; was previously qc'd to 980 but never added to this list... 5/16/24
5905099; % TM 5/29/24
6990587; % LG 6/7/24 
6990588; % LG 6/7/24 Inital qc for float, really bad pump offset +.03
7901108; %LG 6/7/24 Initial qc for float +.02
7901107; %LG 6/10/24 Initial qc for float, bad pump offset in first profiles but may be getting better
5906514; % JP 06/10/24 offset at high pH profile curvature
5906471; % NG 06/12/24 
2903474; % NG 6/13/24
5906568; % JP 06/14/24 ua20039. stong pump offset but goes away
5906557; % NG 6/17/2024
5906228; % JP 06/18/24 ua18097. stong pump offsets cycles 10-30 but absent on both sides but goes away
7901103; % LG 07/17/2024 un1571, noticeable pump offset for NAVIS float, found during initial QC
1902489; % LG 07/17/2024 UN wn1529. another pump offset from a NAVIS float that appears to be getting better 
5906469; % NG 07/31/2024
2903470; % NG 07/31/2024
5906340; % LG 08/1/2024 ua19364, adding pump offset float from DMQC that should've been added before (my bad)
4903757; % NG 08/09/2024 Logan caught in RTQC, pump offset appeared in most recent cycles
5906489; % LG 08/19/2024 Caught with RTQC, I originally noted this as having a pump offset but never corrected to 950m. 
7901106; % LG 08/27/2024 Floats caught with RTQC, 7901106 and 5906502 just need to be added here, 5906209 re-qc'd to 950m
5906209;
5906502; 
3902259; % SB 08/27/2024 Initial QC 
3902331; % SB 08/27/2024 Initial QC 
7902101; % SB 08/28/2024 Initial QC
4903755; % TM 08/29/2024 added; caught with RTQC; was second dmqc, pump offset was previously noted on first few profiles.
7902103; % SB 08/30/2024 Initial QC
7902104; % SB 08/30/2024 Initial QC
1902643; % LG 08/30/2024 Float caught with RTQC, develops a pump offset in more recent profiles but the float is still young
7902105; % SB 09/03/2024 Initial QC
7902106; % SB 09/05/2024 Initial QC
3902258; % SB 09/09/2024 Initial QC
3902260; % SB 09/09/2024 Initial QC
3902333; % SB 09/09/2024 Initial QC
7902108; % SB 09/13/2024 Initial QC
7902109; % SB 09/30/2024 Initial QC
7902110; % SB 10/01/2024 Initial QC
3902559; % LG 10/17/2024 RTQC catch
2903866]; %LG 10/17/2024 RTQC catch on NAVIS float that is still young