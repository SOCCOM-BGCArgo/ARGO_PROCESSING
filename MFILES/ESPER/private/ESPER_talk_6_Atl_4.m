function [Y,Xf,Af] = ESPER_talk_6_Atl_4(X,~,~)
%ESPER_TALK_6_ATL_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:46.
% 
% [Y] = ESPER_talk_6_Atl_4(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 7xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-0.2178];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0479365703301392];
x1_step1.ymin = -1;

% Layer 1
b1 = [-8.8156875656465932;-2.3665145824083851;9.9213484895279898;-0.68139016755844473;6.3618057366886651;0.52052794901336907;-0.36877892564673764;3.4804447297790877;0.83086256649334478;0.90197985749733223;1.1878240938837954;-5.7380922031529016;1.2915199412788252;0.078942536648681033;0.31541199238079809;2.5466164969719105;-4.5454533698953412;1.4520663757966179;-6.9229382005897993;1.2398094218146241;2.9399385799452311;-0.031650499273749917;-1.0579171112577666;3.3873315142064713;-3.6104217976671227;1.4576877237045469;0.78122634234252886;2.5859760096180118;3.5262547586385624;0.074972333439992864];
IW1_1 = [20.791888207640287 -2.1795352540108475 -3.0177000291379894 -4.2840561503583912 -3.1245849816134585 26.051028615549516 17.994278720616119;-0.31298280983965743 1.3680393315487651 -1.7078531279716846 -0.12186770470837928 0.19008021060558358 2.0174596467194812 -0.91921895387657038;2.0352061680598443 -1.3629067500111387 -0.69308614828442039 0.007812398514670079 -17.83651464280095 1.9108875576695661 -0.57977247639449714;2.1261076095809468 -2.0410805086265364 -0.24127678718239345 0.22970840009831664 -0.68256914888637876 1.2931957472058699 0.1794463485793544;1.4309110354545094 -1.2086683542946817 -2.1160082874845663 -0.044090013917402873 -4.4859895824287292 1.4434438854474398 1.7072416078123003;-0.88167124688993925 0.77614420083059865 0.20781526792095242 -0.33637116278756313 -0.11996162003195163 -1.0225734200432277 -0.43943319330263231;0.18309977008095812 -0.260551899274032 -0.66408093170783122 0.27634387353727669 2.1773802267781495 0.051554505501858666 -0.12056672624228208;-0.74661204769819534 0.22430132749530732 -0.50341309548312874 0.40547749818399037 0.77317036969387998 2.6597292824414875 -0.69798265594115005;0.15111589447189522 0.204804698991818 0.7019517956022665 0.65046721723392087 0.15749569661235968 1.0924452368664725 -0.04907804085715007;-0.17926623450820506 4.8515812828822069 -0.42694349332890619 -0.2130942268943353 -0.83711197875793042 -0.84843132879173611 -0.18106980766581998;1.5181600665288539 -0.61744054655418246 -0.81977288822923411 0.28517249664381211 -1.4932829286407676 1.2092309778553199 0.46621232318150774;0.47117821621702444 -1.7057172265552665 0.42607196397040437 -0.21846191833632181 5.2440480413020554 -1.346480940811164 -0.52786718762103857;-0.21037447570388623 1.4874974022668983 0.43457414169424602 0.22719571834382951 1.1732678142995325 0.96024702557895181 0.64944999415793669;0.36530645972882991 -1.3243425461554192 0.092570285496842694 0.0081699168928939993 -2.1209006967861512 0.17326959810957573 -0.41778651365406927;-0.34824798390492701 2.6378424052721861 -1.5483030097099513 0.16397453162442591 0.91584478025625105 2.875884281766508 0.68913138561398124;0.79308564446559027 -1.0752867750484825 0.907716331082558 -0.20492224229326281 -4.1641329965615714 1.2223243925939378 1.2043840550270251;1.1419911065754353 -3.2612916708963318 1.9247451693760458 -0.92598100019328511 0.029573681652569191 -1.265097284726094 -0.71949859594170684;-0.19061627723692381 0.46587461563813748 -0.26961852450846713 0.26015309431063854 -2.4210499014296674 -0.45882105371533399 -0.57135308247763794;0.51161149533296246 -1.5783629836024697 0.65560071333911563 -0.28311835287343046 6.8904069562824022 -1.7379561188598667 -0.64506035790257799;-0.53567339221430876 -0.3593363786222753 0.5896608328840407 -0.02885382055878254 0.80390253903802555 2.3497207849829334 -0.50172760555361329;0.1986728064657779 -0.78589050455048626 0.15193063409342639 -0.66877899295298882 -3.071785512561779 2.0729953211410508 -0.59911120187702083;-0.89536484913891001 -2.8899586728415567 -2.1280761667291217 2.471578032755914 -1.2733032745000896 -0.086117480404345317 2.6945222017157269;0.89855036909846775 -0.89064832699370855 -0.20514293385633786 0.33210224526462534 0.79411101655973393 0.95428766576613167 0.41934152801422653;-1.7385827150045328 2.289513628835127 0.62711811863970623 -0.45128357688326387 -1.784828733499968 -0.82396501442389836 -0.32534774746291695;-0.053173615936253617 0.025025319051786551 -1.4288679081908666 0.64137581051552384 1.2446664245474552 -0.40027590721564787 -1.3666255481207392;-0.91638682835946972 0.4085247906138369 1.4963436629684923 -0.59883249647884162 -2.8325847833945077 -0.49225751988414618 -0.52841971751508066;-1.3689681134914486 1.3799556131435082 0.43664457903119197 -0.44818881629271368 0.61800182569156947 -1.1967286857843065 -0.49855166528391603;0.16398891428992163 -2.2690258159485777 -3.3860336951337398 -4.8100857258748286 -10.350541794494545 3.6663788134797222 -2.8540124151465815;-1.4216380800549737 3.0945891857355505 -0.55092926465833159 -0.28215340180283788 4.269330778228686 0.20280122991756652 -0.40358570398681171;-1.2732095190972839 -0.57192750619371946 0.79470443110373823 1.075162982252504 -9.7025401516690479 -5.3736377834646625 1.1727754507118482];

% Layer 2
b2 = [-2.3596367975413273;-9.6053175621052986;-1.5938974051838566;-0.30566819712812876;-0.17336192853346616;11.249853315420099;-2.4794971216399184;-1.8203662906527478;2.6033087804731059;2.968053595700292];
LW2_1 = [-4.8455373125026018 5.163325597396665 3.4594031093332522 -4.7854530326413522 7.2070056877222264 -6.4167752428573372 -1.8854811007322025 -13.11952443878886 0.44239312681694948 -12.398399889284327 0.56801356436008332 5.0849969854281563 8.5578330275122152 -1.8655371234437415 -1.9330081860743638 -4.5863288341801267 4.8426628141250623 7.8996857396624049 -6.0819529355991797 12.076229630222953 -3.3394305374159181 -0.44268077605140188 -23.190568173205978 -8.294270763695577 1.1716386213050978 -5.0074883988046457 4.9719741870653831 -0.081513981698313934 -4.5368651516600496 -1.2151147811292813;4.9181432232938382 -23.23693552416006 7.1169239804456392 -0.5407952619050671 -14.183051213552529 -9.0968803949023069 1.1557768421344383 1.6441205473242821 3.8926282594648587 8.6972957336232906 -1.5780268606789569 0.25905433971433678 2.3187575459333125 2.8544306983967536 0.41941393782914971 -1.6896198562886477 4.9358332078245386 -2.3473367745208757 4.0135598443619145 1.4156930736369007 0.10471609168424748 3.056941544764463 -3.9702896285595157 -3.3878603972369814 -3.5958460410428605 -7.8886908482605174 -0.85475667270194577 5.791046827222174 17.879404473083767 5.5489789178276787;-0.001993579251720791 0.49546846358032548 0.045862474356130153 -1.5755299899717141 -0.14675476628193307 0.46635364679736685 0.55819211557926984 -0.85306674916520575 0.44744338777571385 -0.35269733473065934 1.7915614779884832 0.18826382977209621 -1.2623858259035847 -1.2412253892492211 0.19052378042327836 0.68360905985684095 -0.42624172978051589 1.333469399984311 0.13669967538162964 0.57558194523656836 -0.55417773725956554 -0.000698950806951304 2.8509993729675971 -2.2091982401120469 -2.3197814417036668 2.3251551209620649 3.8129894258957728 -0.058890783574209793 0.9006236539264183 -0.050633050846756768;-0.83745391125691548 -2.0597310353176566 -3.754376131505992 -3.3113153296776883 -1.8005062358389958 1.245819434701825 2.0297280011063767 1.5189180140683889 2.4342706472210369 1.2206070822945805 0.32858380158885631 -1.1411072255291923 -3.7380850937933854 3.6326729535725262 -1.5993985576892276 1.6953189618956599 -0.94906606464287746 2.0311849643936704 0.31885698856074424 0.97422049155217405 -0.92972896108909908 0.93372125861649324 1.5461328694356336 0.40306976526186999 -3.7081973136243009 -1.556371279403993 -0.28039487993625806 0.094791171352555714 -0.73944562710450168 3.2556086345187785;0.33668433382306717 -1.9945793743352831 -1.5708864954815225 0.78430339003914507 4.3476335416360063 3.873880918868581 -2.8307640723418919 -2.8225237638148997 -2.4003497818526123 1.5045056719993395 1.2622534645150252 -0.70027317922476195 -4.2517577053432447 -0.71557423037701073 2.1940518352297582 1.7853682957571169 -0.054156352379070155 -0.69234903696230909 -2.4102200265047062 2.8498241604326591 3.0806990896490558 -0.18084874455779115 -2.4206056224094521 -5.2497514443934739 2.2309274940604977 -0.57840229542711519 0.39782104529209 2.938920076048837 -3.1705239547719284 0.79953372662103062;-4.7116718520521106 18.100454452884275 1.6463746444586447 0.98289623953420535 2.4097572069044224 3.4068903386728602 -1.0368772685978964 -4.9634719124933238 -4.4116332315399553 -4.3437117331970496 3.4996926605521903 -1.9610091052340559 0.46490170632165478 -1.7218876443559266 -4.6861988515727031 -1.7984513312329493 2.7797999144689349 -3.499036050173022 -1.1987976566618501 -2.6354951713169061 -0.93652009630113509 2.0211626022496647 0.12824231788269219 4.6073750718943476 -1.17870047492232 5.6810883742941085 -0.71411813122959467 -8.2539736947392441 -3.3083675331840614 3.3089267661706994;3.3436075137541157 -8.9722375450983645 -0.47630404278114186 0.23816917269420138 -0.13005898130181995 2.7782430789218573 -1.1480456029518886 -0.079218328777078573 0.31778079104214224 0.13641188426406642 0.08421380113988014 6.1856435782766539 0.77786015958086729 -2.2479874849164267 2.0406824425826358 -2.8866352056851707 0.89282024978246943 1.043087386097584 0.3692619774585929 0.22734230130018526 1.8691559112765177 -0.95635288879877201 1.1578842214210012 -2.2312098541577363 2.6053205710306959 0.64132581681787815 3.8458014795647024 2.9805399953712586 -0.99066559193278525 -2.2814227839440249;7.006705659185271 0.15929129948673199 7.6161747294236699 -4.4425808445964181 -3.5241931210142332 15.436815354915018 16.899852881366972 0.99662836035859836 -3.9058837950782799 -5.0464088358239954 -6.7178168433049716 -2.8190810189654765 7.4797548250702199 -19.817618522574975 0.67537943647148635 -13.717359839250282 6.4576021250819364 -2.0763731076488985 -12.166345820715978 0.087782067223910576 5.4189442705038582 3.1874266527124342 13.562663976406578 5.7328785751931131 26.090022822905542 -5.7317970238529101 -13.73954948577563 6.8886046957123295 -7.2170655559958456 -35.286423582441969;0.0021586851913024349 -0.45580437633022025 -0.064041681338355472 1.488658980652948 0.25883292444346384 2.5029307713804849 -0.85414224092157776 0.77405961453077132 -0.50542698657629437 0.2790404714213276 -1.721499843710987 -0.56945279425459794 1.267970622370334 1.3555315897132421 -0.2178388066459383 -0.62632078180916062 0.44467427014963806 -1.028514709582169 0.071592581787155718 -0.58831068749345505 0.6281922096092768 -0.00021124484808112979 -0.064254373084398131 2.4556810778892237 2.0309106096526723 -2.336040773145045 -4.2653662695674068 0.061523973760483873 -2.0645199078121665 0.022257619155937176;-3.4719057675376832 6.8722110458162051 0.0090700440762239345 1.2997280124226109 1.2315526846845142 -0.80799460661129985 -16.508695859961158 -0.030505506568411386 5.9360557556384501 -0.42473512297119714 -0.73574505544088664 -0.084327030092801292 -4.8913833453081415 0.850397157229412 -0.24901300550146402 0.54938091951563151 -1.7930887688494286 -2.2932446111409672 0.30394015717714939 -1.2693943637661369 -0.5166994761926117 1.1367814979250221 0.1046729459279086 -1.2169211324206468 -5.768930047906383 -0.26395040983303164 3.0132461435541593 -1.7925030162970581 1.7255203084405137 -4.2511646807375518];

% Layer 3
b3 = 0.39786769713485726;
LW3_2 = [-0.0021663577083553594 -0.0028530785346562392 -2.0599983950724181 -0.0054026432135379958 0.0036518556197891352 -0.0010812129881373413 -0.0030256327760650767 0.0010845114160319416 -2.0825882961765685 -0.0070312458268387982];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00111018595614765;
y1_step1.xoffset = 1025.5;

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
  X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
  Q = size(X{1},2); % samples/series
else
  Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS

    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = tansig_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Layer 3
    a3 = repmat(b3,1,Q) + LW3_2*a2;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a3,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(3,0);

% Format Output Arguments
if ~isCellX
  Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end
