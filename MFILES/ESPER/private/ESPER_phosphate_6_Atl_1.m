function [Y,Xf,Af] = ESPER_phosphate_6_Atl_1(X,~,~)
%ESPER_PHOSPHATE_6_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:28.
% 
% [Y] = ESPER_phosphate_6_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.28];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.0470499670650231];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.2130927943271061498;-3.3088427249891756432;3.9868832972150252125;1.8224507740536415312;3.148049917273767484;0.63533570622784163451;1.5853602617640996097;0.15649794727674251615;-1.4635139157431564794;-2.034554445569831671;0.51047100937912315288;-0.53149570733336271022;-0.28502894702778641678;0.52131522559544885009;4.6473038939338549014;-3.138313965710167075;1.3955595219741827595;0.93164265199334128642;2.0799059885716548735;-2.110889234145852722;0.77909486428275065428;-2.5633001990154089533;4.9601080004445918803;-3.772595426640771521;-0.66235594861238933362;-1.2120939524274869648;-1.2159557023497951089;0.7888238007820951303;-0.34663506575041391766;-1.0102930478235325751;1.2228425161772953356;1.5777682839932627434;0.39944640455220931097;5.1554801113822277259;-1.1939510390839016551;-5.7228904185822626616;0.021156354647600437296;0.25749026446687944825;2.5559649730577413962;1.0583160221113852018];
IW1_1 = [-0.33656469842366348866 0.83172386478228443885 0.86101265773068325071 0.60340339830876643656 -1.673849301558703484 -0.64740234381468053648 1.6584160045421507679;1.1272039042018102695 -0.7755609651452466391 -1.2441012311776840082 -1.2813296392484974096 -1.1937181025576946158 -1.0905438533791658173 -1.2390793257179641973;0.88912679858549026335 -1.5813374007229756835 0.67982536421286843265 1.8091748789281199006 -1.6157475933747589281 1.0071349725898024463 -1.7189148833540459638;0.029029879431224748371 -0.014940247226235311984 -0.85924744700060029512 0.59254109404583821252 -0.9946373623633092409 -0.59600157996465208399 0.47250142277316381501;-1.1633949638128986415 1.0282445048195254511 1.4202789866473854996 0.6450252591305331018 1.1985352520999876802 1.5153714505950801605 1.2496032812821136204;-0.61602227342814219746 -2.4171262520253917749 -0.59445744734050931513 2.1042803758338690656 -1.0348874407980444001 -2.8355906257272724957 0.15359216853493767863;-0.94359625450911555067 0.23560213717839484371 -0.76784014926412691437 -1.2628373916864945592 -0.9211012814036959595 1.2929914243384761363 -0.35055980852928231295;-0.7877685085855324898 1.0236338503860300531 -0.13264710998016901455 -1.4366098679775318736 -4.710680958616291214 -1.9437264177706294888 1.6784037441143737635;0.66459693213535087075 -1.4875818910250562244 0.44474827894721097321 0.36736338044276500847 0.32505131049559748968 0.27285520078469033267 2.2062729494277695252;-0.047965539972626300425 0.24675830406289195329 -0.84490412049857210874 -0.7166178813410173376 3.6553380336694694108 0.70319478728922713007 -1.0338754063743329503;-0.35395736803911931334 0.043328605753345822427 -0.77454733620028282104 -0.76250784045210007189 -0.15971460502450154917 0.82939669753815858932 -0.93065537948715792993;-0.012327730805897975375 0.22628371947116857066 0.66274062642794584477 0.50259450231851487967 0.65046262317177927503 -0.5854089772943561476 1.1285760430861577763;-0.83178734877845550422 0.82488790632869846142 -3.7178031092637282562 -0.12444406005821907757 -0.4091797505619804709 -2.0648607760201747041 -0.33912908538287245275;-1.4760279293763198805 1.6069506755261113717 1.4597978915149480983 -0.58668008676412186642 -0.84413981434706530926 0.094955985344755272637 -0.63381765291046499033;-1.2347520549651127553 0.19766166794770012594 0.16585416374310529841 -0.27924168354017792071 4.2195317261939173648 0.49660279015805397496 5.6904646604649089525;1.0360126661677380167 -0.6931958606343385787 -1.0826379871353792872 1.0372523979949439532 1.0823952931024951596 0.5303090525840602476 0.0077318924768670109146;-0.67632527262436437088 0.36455572093696336777 0.87915291064358513928 -1.167168951924198117 1.9740331554223431443 0.51830451250080700198 -0.0046356879599862305802;-0.40287961686211248846 0.087687147402812554153 -1.322457299817431986 0.13170618730097588545 1.1671789242159531508 -1.5118791398316837604 0.39086236171955213159;-0.30940540036414410086 -0.17212693323034039539 0.73285086289159417916 0.73784770075694205982 -3.3056995909805957368 -0.64529736429095996275 1.0885949554754252855;1.6323466303166962454 -0.47182842376106848281 -1.8720546047344768237 1.1824901485599301765 1.3681722423437066016 0.32606161437787856139 -1.5332593142393562413;0.034262523776842357393 -0.55630333429154066316 -0.1529239506616610067 -0.33638310268835003303 -0.58159465461430037703 0.37738822730373305925 1.4111286467026795677;-0.77942734221530651784 0.080157816873580325923 1.2424466848234483329 0.2446339786142383288 3.7967047737371086846 -0.94686038639888325807 1.3259318339410006704;-0.75704726306775049771 -0.57741048713034481832 -0.29384624323367369847 0.60710452397437264516 -7.5607034403176465176 0.40422802829167764349 0.19803527466549084024;0.0052036208810277160067 -1.364538590294571696 1.282679631605585957 -1.1813222309022448631 1.3655521373987584699 -0.37025241656860030659 -0.89155685753938240801;-0.68543230142725308784 0.68903800517239943701 -3.6072457402327242626 -0.14347226414479444245 -0.10765020591866737709 -2.0852686917371330289 -0.39521578656909134297;0.83866982452620775224 -0.7656675874536127413 0.22380678855556512818 -0.2539775318414077554 6.581791598714581859 0.62141452201628866181 0.52818308578841310919;0.16714606242349391851 -0.25790956091752181845 0.32944988270527353347 -0.10193691331306821424 0.077295155779419227882 0.53293687050937987859 -0.39973898598500201551;0.74967114051581595913 -2.2782170164059452055 -2.676054393759028649 -1.4310293667655156202 -0.63755117508671410587 2.0355388576676474166 3.9101982619529125351;-1.7953075711174781848 -1.1724824000979501282 -2.0827909955960866029 -0.95489375969744727524 -3.7224132778454714909 -1.4253667583345397407 -2.0383307576105877956;-0.015254268441129243891 0.71179134902802199125 -0.98103587299540240263 1.6900577609294875625 4.7031883355527739354 -0.82648552121784601976 -0.53694796603928973688;1.8078142574511528728 -0.1342782443989105623 -3.7767429261843368238 -0.47521727945097763657 1.3827896405245290179 -4.4774249549231175394 -3.9399243208490228874;2.7102668164173064902 0.37929925037651257558 -4.2338885581357565968 -1.6449519844339668406 0.81947353295192693068 -5.696910605986678533 -4.79846912993028063;0.36316257687272995858 -0.20668970281217555196 1.129055228208689865 0.51412030694993382518 -5.3465111133530491827 -0.28419722500105787377 -1.2007844591381591215;0.19181141632568854405 -0.51808138542858639397 -4.519150715656563122 0.98729558598926603352 1.1215916298043508981 2.2190059155787125178 -3.3715739049862762933;0.5662540721792459486 0.66422544540494232734 0.20788642949356911371 1.4997516949948768161 6.8355032340789447787 -0.61105771917414009042 0.82098955700436793137;0.34079103796837389995 1.1391981536980071432 6.6609373202005768277 -0.67593800901223577426 -4.3901829287191302953 -3.0462279870347859223 2.135095093087471696;0.52258802456181829843 1.1034494243543704339 -0.1804097146733617818 -0.12374552200643115785 0.64966611243624938776 -0.2373209825108033022 0.1448580764047924796;-0.16888209158225750417 -0.72341946175523508256 -0.757064613653302243 -0.49987708389564511346 -1.1035341104721116956 0.27889866164261212234 -1.5510090456929965796;2.9401054749333410676 -0.049390685103687215207 1.4302956711708463278 3.6786733998403939871 7.1194738503549794828 2.7077320611466517697 -2.4000549941977831203;1.5512901205179099318 -0.57953799937169692136 1.58391518711915813 2.5399873487091548796 0.80422915486712720501 0.57852815700411708288 -0.94094070710991506434];

% Layer 2
b2 = 0.16293792423602015962;
LW2_1 = [1.1766770891072666494 -0.8826355855559742114 -0.4748476005136448852 1.1715413938380436765 -0.79349917491374566403 0.055783027328931741873 0.53224670475337410025 -0.5489198506625275753 0.24710422930001041397 2.7070467431220057719 -3.5462718515869093849 -5.3924319054086522485 -1.6397066219982223689 0.21855098543362758146 0.31849551868827513657 -1.674590342002167942 -1.5314337560595423504 0.72918649489647968664 2.2867331653360070121 -0.092592440030501987902 0.4533105095954682584 0.51817024560448743387 0.40089871088713119329 0.15799507448747926452 1.7577655966722403047 0.28960907138036162678 2.9002069707177748903 -0.044755819178954110715 -0.083104526299005831991 -0.2153853835779719339 3.6479647376052581897 -3.0239812123361304863 0.67214421677773306119 -0.3959169392613651195 0.50392990921326918929 -0.15090158730202068993 0.71230147583191905092 -0.87911702586877149557 -0.44204083202372779215 -0.16061715608120202825];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.604229607250755;
y1_step1.xoffset = -0.03;

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
