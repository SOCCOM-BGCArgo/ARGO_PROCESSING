function [Y,Xf,Af] = ESPER_phosphate_9_Other_1(X,~,~)
%ESPER_PHOSPHATE_9_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:29.
% 
% [Y] = ESPER_phosphate_9_Other_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 8xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-0.9;-138.904784718693;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0415843642790311;0.0042772329208737;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.2075306100839298118;2.100060060168057241;-9.7157260612445668357;9.337678704056912693;8.6497594492530947008;8.2786108269749743016;-0.69271348719893421908;0.81301924898519217777;0.73957400256014571571;0.28381966524636226268;-0.86266895605261351143;1.3267266593666762198;-0.16207220579562453544;2.9663445253900881937;-2.0097381657670907451;-1.2286231523797566112;-1.6698130929550853985;2.0503545864777055208;-1.9779016520249075306;0.73368989649076932213;-1.6456160741990863983;0.19307226972504973328;-0.99478720589563562893;-1.0615059002810511846;0.49589511067647945586;-0.52751761788435602707;-0.53863028008980951;-5.4868890592236265746;-0.073630352337826332487;-0.73043567604973036467;0.09422462791465220644;-1.6240290205397673873;3.2980939465771097474;-3.3167645203861448877;-0.26534255078999841038;-1.1672536087798288662;-2.631140665421018987;-1.6662368837216705586;-3.5069191626367337555;8.8429765630585706049];
IW1_1 = [-0.045421649942938673028 -0.1223471547707112167 2.8466928617492035514 0.82984729548975466606 -0.84364809301426912924 0.75719228008119299211 0.99094724287045554689 -1.3175755628291172528;-0.94115304327517723681 -1.3336613890382533754 -0.28356289260045214329 0.43596583014144019153 1.1691115866292589587 0.31899453712877373635 1.8935394173563044884 -0.98928024177904261283;4.8542330591120324002 10.965642026896176731 -0.36954449118028265264 -0.52073807878993572107 5.2217485568083112923 0.24704612009330326194 0.47833642388901209674 0.6591548860338395377;1.0231455122248012124 -9.5901324431111376612 -1.3765028809274530985 -0.27929375329839950526 0.3898802885755558445 -0.6444495234853120369 0.79198395087974871043 1.6001891723060970296;-0.68717575546199971548 -4.3142957143104574413 -4.6761568514428208942 0.64073457975476211246 1.6715345901103537685 -0.29200402470910469255 -0.72743844825650527675 3.031465098591702656;-1.0231720440298643915 -3.2513220859669780261 -4.4832683260090764321 1.0845605071879960146 2.1379560318963402032 0.24970378364216347378 -1.180895465495649832 2.9554232035255800781;1.3339144418274100357 0.60004762152090984628 -0.42423148267234289666 0.19369262035861170568 -0.23939353198678195445 -0.20143795947976503302 0.44552523398188931258 -0.8411292206894755985;-0.30964988903816309529 -0.30486480465365695114 -0.72116989246130569491 0.75203485718159857054 1.2385080792935845295 -1.2618999252823686064 2.0302148490211462395 0.03873971071589235593;-0.30002011107955339764 0.30593257247741900695 -0.010492830434756777319 0.11816677259118432697 -0.68257706337032275634 -0.41773360507836182487 0.26485478855159405587 -0.13901807975556815244;-0.15927554819251471541 -0.15998390071949650237 -1.6889572523862008691 0.3165806003977868488 -1.9877973737019789535 0.22860160975596807709 -0.17717223032626752688 -1.3621638036903513935;0.40651641366741975148 -0.5799124298266622457 -0.19374738075849268215 -0.037472596086861958975 1.1497603351775047553 0.47460390676269059718 -0.19963123624600007111 0.065593122887909652285;4.9171926805841978947 0.25984802165840520383 4.4477271999094893928 -0.65635722884594338566 3.5225530190262555053 0.44962622691909376327 -2.125276989630478397 -3.2502522201595214213;-0.014487774121350412496 -0.43991661726980951785 -1.276017395514017938 0.49901621278150171124 0.34610794070677225465 -0.63788883827555475925 1.3379280741037753444 -0.5585069245988638631;9.6777187111310798429 3.5004807651803369062 0.034690459195496847755 -0.11762834669056052017 1.7776645105508723788 1.4954780299102898855 -1.5249586266251535793 -0.76359259585712824236;-0.23622170670994718789 -0.0012605596291266189281 -3.2841778312110885274 0.63466653110706228169 -1.6908931423238673375 -1.2351487092414981106 0.63557667109800564109 -1.9099599954874726393;-0.087798235241981117949 0.06375486416073948559 0.91807800710184905846 -0.75526505623366169129 -0.57647806691618208585 0.17586557419741699615 1.122493496732043905 -0.98138778779556579135;0.042691481270219870658 0.10224696905743753395 -0.28305957900915001968 -1.269075418040946035 0.63066214755851013329 -0.088478398682131878528 -1.7616946585621910959 0.66128508273317665456;0.22069642786404644608 -0.017655140101064649427 3.3535888826106536875 -0.6323218494598499273 1.8152156417657483445 1.2547413729452385756 -0.7060205196216047252 1.8872303690721650149;0.07139341944786907479 0.13416127752445125565 -0.28214940535131249888 -0.94691389051440633118 -0.16904625180502200199 -0.68173278846861573577 -0.6843133808611042701 -0.27691237543345370575;-0.30379434823001155719 1.8579083315632389795 0.36227989592696813181 -0.13138220617880380647 -0.79176161587869897662 0.18730704963690386622 0.3226716006009911486 0.22364408023162229977;0.10703761975537018358 0.084394834500796447885 -0.22468544925890063246 -1.3798831484989071416 0.86719301317283092878 0.052473754739121018786 -2.0950769390474985876 1.1307112230451790147;0.56727503349709906821 0.9365734939081763466 -0.057028589300281144758 0.036344118025063223976 -3.4308110410568870563 -1.1488717553737366472 -0.81320103923978115379 1.4476973239503143631;0.35904980989861617902 0.49665857512186928346 -0.61577504414292028034 0.21934432408286386185 -1.2836307363988899155 -0.30351975619010368446 1.3031017647630658285 -0.61219859046853553064;-0.51416499900706802428 -0.2074993552575182254 -0.46347537580289871473 -0.70633013953644108174 1.0284797097736424476 -1.6319718483375367235 0.28523010736929116371 -1.6112123809768077454;-0.20314361232211761155 0.077478037512504455142 1.4346149026353456346 -0.53251549645845452918 -1.5490088808130038078 0.4026400483383653639 -0.40682626399384236304 1.3409538713432709311;0.20451809539802903726 0.98257054165278767677 -0.46906276183053974282 0.15653212854471457161 0.86905240024501217988 -0.20992265984926461453 0.4325303210381074015 -0.66994265709528244557;0.25817399824728604063 -0.0079996012902773926817 -1.3669102734732652138 0.54124224878642657188 0.9782788802521076299 -0.33513900407898938871 0.36960695166512952392 -1.3059434682267900829;-6.5308151576670310234 -7.0088254534459197842 5.3537720065209937204 -0.90209922481529503102 1.2098667322261245083 -4.1979438710159548975 2.5375092357646003016 0.81545441278239905181;-0.060716492383346490347 0.27236568158396357031 1.1348534129241056867 -0.19537367162001925691 -0.65534176620198125729 -0.090706380359659583013 -0.56479688726699794898 0.79762935089603281824;-5.16063161965900008 1.5322429930977570223 1.5635508310551389233 0.95004002844428314667 0.19220522805852441328 -0.64629122269834549108 1.2355439978453670591 -2.8278738071365405737;0.1237186934563438423 -0.40534386394985405522 -1.2139522681468573051 0.23572156525329562848 0.50319417054693482694 -0.0023770466694112120723 0.69766132780938849844 -0.65162078513746835018;0.032533808367902088132 -0.56020291733290494829 1.3191756242925924969 -0.40369517016327577386 -0.039856334347345229174 -4.0382961932119094683 3.8502334156511661689 1.8870363166224313201;1.0303537922160201923 -1.3416361515898986845 -2.5871066050496045463 0.084952291765437060178 1.8768752856564567821 0.9894808485017170252 -0.48600986697710718554 -0.2975967252527818574;0.037122460430567626177 0.12084641422448419434 -3.0422198222655980615 -0.82419324126682691478 0.24063003286412915083 -0.82118686815617059072 -0.91316771700962562353 1.2715225094471001022;-0.72226554486000804634 0.36727513059331967638 -0.53453879649768332172 1.0038697437032728654 -1.7966225985103219109 0.30505679225650372821 0.89454341752189181491 -0.31858578778686513777;-0.055565185315987537795 -0.15337079880811230193 1.4355748688465175267 -0.85057011138198401756 -2.8263685746574367386 0.76849873589575479205 -0.09110062300987771533 -0.56060230935338972724;-1.0906748118241289536 -1.150590771553070546 -0.2421152130643920708 0.26523709400178308027 1.5758297513925387801 0.0090129133876040018469 -0.40462178096800999194 -0.81034715219866104707;-0.29254997234228719005 -0.52739089517639425875 2.5066852120018521255 -1.4757526394229547595 2.4303462432123250636 2.2503212212478946874 -1.4699353170170341709 -1.1877805852290392874;-1.025660142607673686 0.57441622327305175144 3.2603054676702161174 -0.79514508457811816644 3.5313515357718232579 -0.37049741442175526673 0.36285137600538175384 -0.10108160040063386698;-0.11819642288210219228 -0.13731989346073827929 0.34931517060798211727 1.4998716426619134179 0.99614150747454532908 4.2066678245444419915 -0.23807747530052564278 1.5463569811809492194];

% Layer 2
b2 = -0.48156231566931195776;
LW2_1 = [-3.9781497549899684785 -0.35855917224833022683 -0.030123417277862021141 -0.055021520194069847842 0.18655102755790875513 -0.2457388378269884377 -0.3822403515283877562 -0.22598940097289943463 -5.0770763801940796967 -0.34417445795736839598 -3.2073384310967956168 -0.014526934916775544396 0.74135680488147615019 -0.026328918037493272142 1.9098097698391736099 0.4436898813711440126 2.8004479713713830158 1.9163455319003219746 -1.3067686470779582564 -0.19054761887197410597 -2.3259266244576557447 0.1514991850147806296 0.77648320755130517679 -0.074084963451738400364 2.2965002996068490937 -0.64660714560767251768 2.3243802467751137897 0.010432486856716339726 -3.0801317820917173407 -0.026826380077238084348 -3.0582125685280225902 0.11135220251509117539 -0.20914730734101596998 -3.4104709482406168064 0.20627180796629848714 -0.6170540944001652317 -0.77415249448394096721 0.13025372878394059217 -0.14651284742791043092 1.8303400381343868375];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.568038194888224;
y1_step1.xoffset = -0.04979;

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
