function [Y,Xf,Af] = ESPER_DIC_13_Other_1(X,~,~)
%ESPER_DIC_13_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:14.
% 
% [Y] = ESPER_DIC_13_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-0.9;-133.803853032615];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0419287211740042;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [7.4156061349624398105;-0.6829174374700790473;1.6875769722300684172;2.6706674510376040566;-0.21772576956084033206;1.2972148560440124321;-1.448008348663527789;0.025106555887079790718;0.011283135092246297099;6.187368097391281907;1.9512294155964819975;0.21590459128253236543;-2.1097460653582396084;0.094426827802485485686;-1.8848954909606632224;1.2894293823148639255;-4.156923538268360474;-1.4941400296365106026;1.2865013617622180764;0.319923819195712833;2.1522186794212667671;9.6132556827117650045;-0.25430761615462649239;-1.2947959284886503273;3.5447992593984567478;8.0759968911667989744;-2.2980154179289171701;1.5485549472269388716;0.7026449778218282205;-1.3739358808726083172;-0.10362288759624815881;-6.4334374447528510643;9.3292842804638613785;-1.3341794060366181096;0.074896224090341245727;2.923389751514460233;0.43188171211769504909;3.3485677569401448572;-1.7480301457599709547;3.2736787775053621807];
IW1_1 = [1.6564826164097552486 -0.26891647363282494831 -1.2090635377176581766 -6.0374853236607872731 -1.170160327970040326 20.832340765117347559 -24.615239615503501369;-0.13747798167920519119 -0.073944059474996576276 -1.6806131106181330104 -0.22845743127199172329 0.98837283422221633433 -1.087464034176516936 0.35720083461310625506;0.22502567142579338322 -0.096780226417464712552 0.031824785591747575819 -2.7054909079223405932 -9.5281420200960678102 -1.7683108242558842615 -0.93099947215945100254;-25.68253148522942908 -28.395901301377513448 -8.0243972438269182135 -13.823119172952173628 27.623318184532180908 -1.9651347022746590287 5.0418152950293535497;-0.014355642544736100033 0.29103936223811716077 -2.7752445023296949955 0.05843712304949880143 -1.2918356212705115471 -1.1937770012633022176 -0.10959582756693904415;-2.3620233752376429415 -0.26738018648499245478 -3.9292626293791936831 -0.096744559988566597442 -0.64826352454108848722 2.1143112047603667136 -1.4844026554685176844;1.3599369478656828658 -2.9260638639190195143 6.1615064536459644984 2.0797547249788572721 6.6894818785628720192 -3.745630382100041178 3.2585834729463143411;0.12858989673889284489 -0.13505125561108058396 1.2739158134595340854 -0.69880948279604870876 -0.51471598740597523935 1.3540238898815506463 -1.0461666953645312184;0.19459058834051093378 -0.58948866631384300163 -0.14038118264774329669 0.085756260075696352785 -0.12081207322610505728 0.91623131344933372855 0.84030618434310755394;3.9379318787677179969 -3.4365999261224144945 3.6902126754149717769 4.2218174343771845969 1.3667058083159338011 -13.128055977568134693 6.1770833746956004262;-1.7411405901305863075 -0.88570192145945847439 0.73661347470453086839 -0.25449053423051426481 -5.0649352898208253038 -3.1287531416179890087 3.3156573448836286211;-0.13753557417624215353 0.92292462042972089886 -0.18443047549805896557 0.0084580261111027697785 0.25497930203487978451 1.2086545203045560015 -0.71063579733711634834;0.032631964849473733825 0.0031129576068776882575 -4.0536628422626588275 -0.42604604091461784066 -0.43430939886664199356 -0.3941991482859264484 0.90816478655086796401;0.13790140205193740508 -0.13325496633212749464 1.3450519261232005697 -0.69044448309948180942 -0.54784438158759951598 1.3446349579361758231 -1.1007003442704206275;0.001637920817739173425 0.01011018364132335555 -3.6973055440727904397 -0.36880563600675042579 -0.50170643572193418525 -0.50641538261310459212 0.9144218635153522623;-2.498496504762480086 -0.25720238655193250477 -4.046151674689477673 -0.17308615228771756511 -0.59962299916244399256 2.3548533979607344513 -1.7206082768678372386;-0.74304103419353229043 1.1658121413680251344 0.45922033699301934018 -1.5473677672439500874 0.3677799095284105535 -9.9866519132071598364 10.932059003091490368;-0.10944323955901145984 0.034811322674890579632 -3.0672034575766042863 -0.26930000412821336608 -0.64504359947292111421 -1.1831824637365702113 1.1758540397724241977;0.075989484649253286053 -0.046431453975435230119 -2.6571177250119624169 1.7946380036560380944 -0.81746636539749539807 -0.47044466306380133824 0.62410212654448726788;-0.89619902728468736619 -1.1541429803242682084 2.9729081565371591189 1.2592238060611962158 4.9261470226016061957 -2.212739735221684878 -1.033417013793943573;2.4402483150754532559 -0.37039333575725275072 5.8823119697636512271 1.6533399900403100702 -9.2238692134806878897 0.046005309733025474828 -0.68475766296981488424;-0.031870155370993802357 -0.27272577598434771806 0.74173763911035062169 6.9580754484034557095 2.0260529416788557455 3.5713270935333465417 1.0020819169799832782;0.057898958023350002855 0.19128024105056534721 -2.2999751071054248541 -0.71234981571278166435 1.2568343438103737331 0.039977877608874667903 0.24949030378544220787;-0.087851707938847140023 0.050822644935441772807 2.6935533284617210192 -1.8434378782316152012 0.9646228687187712314 0.4102823926430948509 -0.64569468609277880233;3.7895923412101821803 -3.9303602251571145842 2.670689348661831275 2.1569214427304261861 14.350058952732965167 -0.34453973921575281603 -2.4508758478371821532;1.6318009845434866545 -2.3574144539405534893 -2.7925188491299808291 5.1427156888073817242 -5.4687559122537630074 9.1726050326065475815 -9.7084876504625903237;4.3117339941511954038 6.6893085171160251434 -9.9504167419892954882 0.69437570919880642339 -21.779410971539963526 -14.505304472539256366 11.759943569495975524;2.7844611624348969414 1.5363808195209904817 -0.5994138580174719344 0.32761880315788360907 -1.7119974588975499508 1.2323693386614671397 1.119945063073784608;0.14611285068801907472 0.13892251189251694776 1.895221631899300041 0.24382936085312822971 -1.0382997296624572048 1.1406665375624001069 -0.3835696385620580573;-0.57282151474386944212 2.327135259182344651 -3.0636591824567909192 -0.70289816946996031 0.87770458057942291052 2.5646649231919558254 0.17252421047456667491;-0.04841781437120203091 0.12290737365125867231 -0.97100325436632872744 -0.30112883154602426705 0.6631039526517619187 -0.16312843650277589203 0.10487046305315642081;-2.8256657325851590734 -4.1098332493991582837 -2.5862629908371790322 -1.5521048124246377142 15.094013665696317972 3.0098720372785678911 -7.1590013071885723051;-0.022204285825032407664 -0.28378415279895469192 0.74586018358632222292 6.7367988055745131604 2.0508459544481176806 3.4650714530325097407 0.94643321038424721436;-0.54722263368445189791 -0.094398318542841294643 1.5590623089102195742 0.5179139346670905697 -1.8152389437346958534 0.38703608110819576194 -0.295725488615780574;0.1149265580668048653 -0.32155548587873417077 0.79602195409858778063 0.5679833126798456977 -0.48009118778725212717 0.7110344115279790822 -0.792488068913749788;1.0071820719474462447 0.17656152612416828318 -0.52815406508052575862 -0.44624778092833289334 0.68730012319967515033 0.54084304163187368886 -0.83382834015985318743;0.3335218205032450256 -1.6413967879785551496 -0.83305893130224073495 0.4258882473273259972 -2.4009795835032257116 0.73014357752647573996 2.312398707780922269;0.97629980801756366393 0.17635472675870991766 -0.035748124780548215518 3.1611479036696876754 4.1873737625716538346 1.6255235230045443906 -3.1334429723792034395;-1.756976687713961649 -1.2932857319677473207 0.99995519398022314483 0.89976795231217332471 0.90531500640866247043 2.7143480372481545437 -2.7811949535159246061;2.6983860286042116172 4.5600785418211851052 2.9780362156625175274 0.75430215981669712644 -4.0937387025698122045 -1.0569215468870709618 6.1518113093884698728];

% Layer 2
b2 = 2.6938032193830596661;
LW2_1 = [-0.036701166031028388537 4.9964952137544953459 -0.05890035309007064096 -0.0077873781573290406025 -0.21811675090383758624 1.1248869785790629461 -0.019466722802074645088 6.8481741285146648934 0.26458382204910940727 0.006174354015593856021 -0.039083916718934706702 0.29748861581455088299 -4.0049119059991333813 -6.6409806253782859287 6.0069321887102749002 -1.0056166636512373636 0.018958720574181818497 -1.7787541275281057995 4.6566922824975138795 0.02456462971332224951 -0.043937657812700071924 -2.0086532866119539698 0.74990163540986221058 4.3896905812188498075 -0.01438270983863055956 0.016811709761430877752 -0.0096216864385291651918 0.058755677560207872356 3.7178792960237294629 -0.040778276966837154194 -4.222063691925273865 -0.022891536694839433042 2.149196213484469542 -0.98073516422999762465 -0.80768241611477575859 -2.6604952725113926171 -0.080903244627902645703 0.067891217651871671457 0.055541049452047026869 -0.022821030831962390434];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00314743155299998;
y1_step1.xoffset = 1774.26765577007;

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
    a2 = repmat(b2,1,Q) + LW2_1*a1;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a2,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

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
