function [Y,Xf,Af] = ESPER_talk_7_Other_4(X,~,~)
%ESPER_TALK_7_OTHER_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:47.
% 
% [Y] = ESPER_talk_7_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-133.803853032615];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.7024987077283102;0.72068021638320123;-1.2483133166064808;7.2115299985719687;-1.26919183754053;2.4661385060894565;3.7276196529122574;0.75119877355367837;0.26109433383516334;0.54135375589207424;0.52175932860705743;0.17319183651296549;0.16107101289693576;0.88712516899057947;-0.8429934755213927;-1.3654204955336651;0.57294120539428095;4.9837045453307658;0.57251305223116811;-0.014959453577475209;0.44243671535467949;2.738360402059917;1.3342154892002529;3.1272998011220383;1.6308976986419601;1.7408911013683936;-1.97379483794065;0.97129367829484858;-1.3425215232094141;-4.233388806209339];
IW1_1 = [0.82144708443501158 0.9106792547898066 -0.54122684713923275 0.04874788312071323 2.3896166487767401 -1.2375996438145009 -3.1738260941909515;-0.13017759991811531 0.055430947335479402 -1.4196268622464594 0.50683233614768664 0.44246218680207694 -0.79191508502374863 -0.62753565054023808;2.0524229086246231 0.81794678741059834 1.6713352032055599 -0.59079248583065902 -2.276764926509633 -0.53616535326833492 1.2873687391964534;-0.60974599781618199 -0.43499551868274072 3.5646819560242493 -0.5495147501865667 -3.086988028027446 6.8241842420929677 0.42190438136115621;-0.14165654778799583 0.18251251701205809 0.20219623984966367 -0.13788503453040113 2.3770885878303383 0.67486594206973871 1.2748283988299784;-0.83353199436821312 0.95782251295853826 1.0549478079386019 1.1389551384388359 1.5363136283132439 0.26451670418430345 0.7258638801461561;-4.8392926750144447 -0.21908113424591424 2.8804122798218152 0.7074118662587634 -1.7698599006067643 0.38475278566552873 -0.033766683523729787;0.77606692725076332 -0.1825170715887848 -0.44265078534431562 -0.56628960711824317 0.79897878577321113 -1.5571751397486944 0.18258431272267239;-0.07271051655964797 0.49543098882951941 1.8077328554158238 0.48092389350099846 2.851762529813668 -2.8741948358587557 -0.97498641101654882;0.6744763692711323 -1.2527725150373754 -0.28213084565951391 0.058239346033441684 -1.1894768554048576 0.43646192323462246 -0.96503084784632687;0.01880172660091927 0.5019658092991216 0.62209452994590475 0.30216292829110586 1.2661370302524673 -1.8429528258899366 -0.36784575814164494;0.42386497734689071 0.56970645744248161 1.2884528491649789 0.78257555192364836 -1.901626242666719 0.40850367904788099 0.17751559551366858;-0.041051060599106994 -0.63083574650816032 0.74522606085960019 0.629082464116967 1.774526624560804 0.15221406973979756 -0.66432724892279948;0.47963785696348393 1.1823922077519518 -0.025950318324115654 -0.24145285936339053 -0.036908232951122083 1.5244890029657392 -1.1139710171929365;-0.13821503857983278 0.22589468307801489 0.81272920558940931 -0.11984674712025829 1.410360823211092 0.34037850605047915 -0.36968678829120949;0.50310574387700746 0.24984255719481449 -1.5970241332295549 -0.37155154305332944 0.93191578816134435 -1.0820294300628972 1.114757567748514;-0.21420664096681141 0.27657453042107144 0.014576906543302367 -0.1979468813359985 -1.431885499776427 0.28718564394496188 0.20622435330407748;-0.82078784498221735 0.043461223334739237 -1.2324965886302248 -1.1298627301786057 -3.183431661622818 4.2219175812961547 -1.8759915411214241;-0.2728064767298003 -0.27819408197802392 2.1811648885996018 -0.28651014127273011 -1.1031114633081991 0.49974559045477235 3.2210372638725042;-0.080442119505108114 -0.084094370915000702 0.15192047948743911 -0.28160038195188003 -1.5342767492507587 -1.4613169437348521 0.053871287554837874;0.46054300812539867 0.37210111710817267 0.47940322509595462 1.5638998745556831 1.1026223902533112 -2.9566443412585501 -0.79118282667528472;1.4074008836852319 1.4571483642353904 -0.77782671959843563 1.3144014041597836 -2.1352250689583152 -0.35862731472869142 1.4256350784959433;0.37888531269274672 2.0579609380112349 0.90228917756390503 0.36986254005878133 -0.41479883816991542 0.22043005342236452 -0.1142608890939238;2.1731968906752832 2.2652500165609561 1.8921338075653589 0.28769597258730184 0.067567352318302798 -0.90384780421948041 -0.32642908680285837;1.1325541349128627 1.4860238138654085 1.2712512395936317 1.6096725348047201 -2.402163305707762 0.79048986324463 2.7221344566656427;-0.02614512583322131 -1.4662881486602677 -2.4615734252381816 -1.5645463479155002 -0.24245074530174801 -1.6870362944825934 -0.61238110881608032;-1.8180832708513608 0.68809946916339171 0.16302197567252213 -1.1372028766523785 -1.8268313352553363 1.7657301593442771 0.034235113830572608;0.14338416952427477 -0.20161090442557553 -0.16887310617444398 -0.2133400049235008 -2.5071657663933227 -0.43986283806789955 -1.1020951442461016;-0.11241395057026822 -0.22560513414605832 -0.42136003851556725 -0.67542140509379012 -0.36459424095375492 1.0399944316676824 -0.32769608263679251;0.52871865721483224 0.74433320966174132 1.3620334016340683 0.54655158075102317 7.6983896532645479 -3.2746906844683283 -0.33641117700591344];

% Layer 2
b2 = [1.2390766293314359;-0.63685984909809656;1.2393523618408306;3.2842247619531175;-1.0147310683764317;-5.8718857464898813;0.83508085218258499;3.1336065537158846;1.6264205952766726;2.9150388365897344];
LW2_1 = [-0.05571361139309336 -0.017162084153468204 -0.18583625952527466 0.10537033214380755 0.35996668629478951 -0.091021543669906155 -0.028702391168314315 0.20600703051151917 0.16812528060066062 -0.046261889772814861 -0.88777989956904213 0.2762254254389559 -0.16515947997633418 -0.055603326608631584 0.048157922665335577 -0.068014855649281614 -0.062390223241766143 -0.052919026557614199 -0.1776020495251317 -1.0115825995237429 -0.11848333678517178 -0.11137561813996816 -0.039579864987840933 0.06217638237267735 0.045116992870921517 -0.23369666854690385 0.17206112087769121 0.57689041032849941 -1.6274134149511976 -0.32072247783145025;-0.72386829887180359 0.39777825241686915 -0.31863922440648057 -0.96670398040592131 0.99107423624283175 0.092341397052784588 -0.12658046357408209 -0.68596340943236389 0.89192526545427919 0.49583615643781065 0.86222120829857829 0.61061926070961436 1.0197972753848916 -0.048444216392947312 -2.0119748344624639 -0.47376662238389106 1.4940305783744194 0.28952005507047573 -0.15547532989578361 -1.5529900365007443 0.13042848055296874 0.84597871976274674 0.18230727572309988 0.030105289107238319 -0.42008694449585954 0.032150859142865723 -0.23267315726324861 -0.66664019915792139 2.3747869929499079 0.84945744824005531;-0.040412359700246078 -0.072791766705450484 -0.18569366925265315 0.083073147781253728 -0.083191771031231668 0.021285124126696316 -0.037552972829537568 0.25911064534602396 0.050392247129901389 -0.0019248482856147082 -0.2977214815667506 0.27390562740805791 -0.22454477092786287 -0.046633614617071259 0.15629575346626109 -0.037775231895033998 -0.35321967111360741 -0.021963837384042628 -0.14014468133815997 -0.57719821753300848 -0.055232620217271239 -0.096431840124389115 0.0017883277827175058 0.060041260915041444 0.03330547529429595 -0.2400085589505872 0.076674590312995061 -0.12321942565319917 -0.62802078786292903 -0.39447627952599457;-0.65680893088546521 -1.6895656696717904 0.21371111663877793 -0.23463859657162803 -7.5179415267713035 -3.4673703057152516 -0.51117051865694352 4.3951680230816619 -1.4113328945755494 0.68921745541500512 0.77189165437771856 -1.3676881984456359 -2.0169374083773546 -1.4679583470188522 4.1473146033443351 0.52549902206710153 2.1189488588966241 -1.0567220127104597 -0.78807069867090862 -0.10077961625359043 1.3574603208328513 -0.42919031875891334 0.36203378219014798 -2.5493168568363926 0.065271783277486145 -0.14496862350778877 2.1599516492708717 -6.5975700659190704 -1.3686109139744183 -1.5527923985722323;-0.31679005064067933 -3.9362662680912814 0.085298901150684434 0.55133096006272142 -5.8888201898146253 -0.59448696503544218 -0.34157026459319384 6.3374648090902097 -1.8197546474093156 0.21346967170092546 4.4887898305060059 1.1082070974317648 3.9317100877329922 -0.045354138438227845 -2.8051973391923943 1.7009630309974293 5.6457065506460475 0.36582312677433704 -0.32267185972127327 -0.084880885667429073 0.96968298019113086 0.90447801728807065 0.97267903063805872 -2.3925893929117343 -1.2105283689256745 0.18443858471971664 1.3963580067472987 -7.3767255594353323 7.1197345822177782 -1.3665706933230344;-0.73573867395892734 -3.4188430346483734 -0.080377148952034377 0.75744821189449807 -5.8273622919426673 1.0763924908960976 -0.16727946393732834 4.9319864648661991 -0.28041667175821777 2.7508857681967469 1.5833430704521443 0.039171151615388665 -1.8547430108627554 0.90347725359356545 1.0495179790070051 1.0536537260322199 0.67913734165132167 2.0920873701681497 -0.76170900142251918 -3.2299719199321193 2.6872652179293977 -0.11176656647577379 1.2690123858281424 5.2334368483764324 -0.68089939164875446 -2.411616646542158 0.42127329914994427 -6.7053566029932581 4.5161406900569903 0.85499366879656769;-0.63927997594547281 0.67814477154032826 0.65059961912509956 -0.46014055569506696 -0.25897815989719347 0.22289270680516601 -0.16155066275667951 -1.6197526810578506 -0.72971690169203074 -1.2674692476838394 1.1233823293963952 -1.5736209890290553 0.26860356193759694 -0.52026940739859884 0.19417404117875006 -0.69476185214363717 -4.0223220138930147 0.08273041481455197 0.46070961923934101 0.17159927454145724 -0.19461982042394324 0.68490906957814934 0.8933018353802975 -0.5795039075424786 -0.037354557922939423 -0.8691133313777385 -0.88473125483714332 0.57583773868118637 -0.058856609680308619 -0.41564023659271415;0.48839156787225924 -4.1037012512081965 -0.14749271235258618 -1.4104690608468442 -5.4039442565574829 0.17429923236541855 -0.42723912522922902 -2.284611988103685 -4.2276958733360841 1.5028558065649578 9.3286745448087558 -5.3149058360691175 -3.6393243619667519 1.5382473327382882 -3.6287927774500965 -5.6148001317779395 -1.2611989005143589 3.124141869637953 1.261861164127446 -0.74340551638251917 1.527158389576734 0.78968050638193232 -1.3392037506362122 -0.28302996108933587 -1.6904932996376998 -3.4872022828095868 -0.60684357594550353 -9.6584114067719309 5.8958247755472781 2.898066713286457;0.5048871187913746 0.082617260735422443 0.27545594303621329 0.46627909415575641 -0.18515371409204848 -0.47793275469380042 0.075552344512381756 0.64379438698646452 0.11202128448536905 -0.25918047345552447 -1.1189918573367656 -0.95342818305750132 -0.6606078341223246 0.14630320139226227 1.8041085390699543 0.49790796627842521 -0.46419157093701241 -0.55311266191428143 0.16693305343417419 1.1315063321362282 0.033310226963958467 -0.49035627597583614 -0.46560633632150389 -0.012908518205370244 0.37288825194368813 -0.072633701668856729 0.32801803043297229 0.56057577790164792 -0.74851826961278745 -0.23825286827078324;-0.96544677846975802 -1.0081100626464885 -0.35253278526350934 -1.4505513012125095 3.9478552720136215 0.41087806367602325 -0.69237155597108391 -3.4446385074961152 -0.03022735055718876 -0.28816198708307483 4.71729648046297 1.6970398605059878 -1.9237427024633798 1.0521965498387889 -1.6969145959521377 -2.7259895317492111 -4.6319540389019638 -2.3049020143232926 0.39026940866643267 1.3188720490339871 -0.094375684871624427 0.14134047357232424 1.0726149158761369 -1.8035704035121412 -0.22454274100285901 0.6147527423481236 -1.2616934560088815 1.8649599592505535 2.4223645691970019 -3.048293382250955];

% Layer 3
b3 = 0.40517337788316576;
LW3_2 = [-7.3514679994456635 0.32622719679246637 7.2433168082698538 -0.094931000138801569 0.064119663034760979 -0.04283527315270453 0.12525350033538085 0.042216364891449192 0.54029175680114971 -0.048067114230349865];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00492610837438424;
y1_step1.xoffset = 2075.6;

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
