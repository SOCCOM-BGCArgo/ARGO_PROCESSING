function [Y,Xf,Af] = ESPER_oxygen_8_Other_3(X,~,~)
%ESPER_OXYGEN_8_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:46.
% 
% [Y] = ESPER_oxygen_8_Other_3(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987];
x1_step1.ymin = -1;

% Layer 1
b1 = [-6.5640761751368925658;-1.7069965223579386393;-6.6228466439734603455;9.1354310881969436053;0.88535933105661446252;-2.4161668845922901205;0.69567843912397431261;0.0074367803413638820798;-2.261696935125404373;1.6182893342987512142;1.2299761716895003527;1.3271994605513999232;-0.51227058144505710491;0.49915874011550664946;0.31224626766431085034;-0.26384481308419194256;-0.57500799767461063805;-1.3478714780602223389;-6.5401118053399542873;3.1590855109599447204;2.2479883035884742348;-0.9574531899634904919;2.0247833480421144614;-1.1608172895655399781;-0.51227721087174349357];
IW1_1 = [0.83730030381899345215 0.71968765036163662163 1.7510356948080250117 -4.6394967917814957659 4.6785856631608160328 -2.1492036499515028147;0.37088838747479679903 0.11364761134379486318 -0.67116983146486841605 1.162983737516085192 6.0221588049363337092 -1.9232698570513537817;0.87615679513040012338 -6.8643778762075733368 0.39284835197779466309 -0.55940601180453308761 0.20171982882752378941 1.3417472614086243254;-0.71561814968030357242 -0.55255416735906537351 -0.69357585072670790982 9.6075285806200607652 1.4084499318584735583 -1.2739193146414513613;-1.388641861923238352 -0.52636952016347016148 0.25973678413842676838 -0.76214322499287989743 -1.1306719116196388697 -0.10796498182420978529;-1.0493481729452391882 -1.7929334007569182052 0.79563490526487212762 0.57510566528547424348 -3.0431942064104666734 2.4187172682610147412;-1.4950097576731040672 -0.438204219738162426 -0.62444433067531945891 -0.4770589794586739818 -0.9624259163253708893 0.14137251444335369177;-4.4295116727809835666 -1.5364787718156780461 2.0067705993215336768 -0.51065708336313209159 3.6102964801975510056 1.355467030413469498;-0.16672745452687959267 -0.62357666849390935848 0.83732653693018810337 -1.3989495268853056054 2.2495923633328787083 -3.7506021158062292642;-0.74219185236849938292 0.8520025183207364794 0.063869052812488558146 0.81386737709498135906 -4.8278986456265204907 4.6448808440671625064;-0.10469058530339062085 -0.40441263808961691595 2.1963743078163151523 -1.0701558363387957407 -0.93149193886411718601 2.2431116048832389076;0.061064397905487767526 0.11216850642309236785 0.22324336788488410854 -0.10880909575678317236 -1.8573472738785818947 1.361847842708412859;-0.49305847711694461566 1.2790349072503033323 0.82633278247539132533 0.0745812074966921007 0.58320320321291507781 -0.73184819318274030309;-1.419470799588518517 0.19250113084417905251 -2.1601939306319963841 0.23353444021683672505 -0.3476463343998865696 0.54699035864463785295;-0.37997384825213825765 1.0446904099273028876 0.89221145621695330252 1.5152847525784365335 1.4568126870343665047 0.27266036879918958968;-0.99159154804735649691 -0.68624916406174896188 -0.70193380395325732923 1.1834275750847100817 0.95031613705281448023 -0.14976977221403561802;0.40953718112708958898 -0.44067719284720924255 -1.6441499386567417496 -0.5413010000198408278 -1.4171082882664935987 0.58146873587769964953;0.29751223273970700545 -1.7080548228692635782 0.12077437812971195119 -0.62530814359440090033 -0.34184308375229066623 1.0973807093379768229;0.32472526353813496724 -6.6435279646697464884 0.58736566848042082167 -0.36083678165821481132 0.80769027181164509965 1.3437263504891936439;0.39836208225508296099 0.20169994002222488283 -1.9373206961514992663 3.4854484101304645804 0.7457089281342040854 -1.592388496677342502;0.087207840156809435239 -0.27123236870631139972 -4.1233950175230686952 0.68346906285442399209 3.7894288202367758522 -2.2038255613942681066;-0.0049075449829458547141 0.71597682431558529181 1.6418255683610589735 -0.2502429115272667981 -2.70377712113777946 0.72491428392763468391;0.13122263584200491437 0.17447342492947681691 0.48496993558543166625 -0.23896515140425952928 1.6846077001591381883 2.1566500368705194823;-1.3810740200733304661 -0.56573064484206481151 0.73632081585856801187 -0.54933200015787819748 -6.3810071246761213004 3.6649278698170957291;0.160506347554496831 -0.67770071793024144746 -1.1993617239972007571 0.76431237823612008242 1.219791846040680916 2.7312400427965353167];

% Layer 2
b2 = [0.46863139026260253184;-0.3768213272070142783;11.708899280531079867;1.0275501831151767629;-5.1600230345327631198;1.1771502901007495634;3.7303644844533683234;-1.3294689591065864676;-4.603583111667454375;-0.29246586774044597279;3.1891717952137015324;-4.691359987139836818;-0.31616496456830583961;-5.0718637520006177866;-7.4652516085025926174];
LW2_1 = [0.80986409432405515041 0.19080739467740023829 -1.615829520764072269 -0.065268769806914275766 0.29811993508718548362 -0.33696616082442470974 -0.11767027386760202379 0.048621593645321772859 0.13440770995953921219 -0.76406575370615814613 -0.03268123991336834705 5.334974077870573872 0.71937093423786424573 -0.49207705790837519277 0.15341025957799164781 0.01888202907462609853 0.8873579068958541427 -0.0019142331434560603924 1.6976331885622542028 0.31758939695204152498 0.15494670550923292374 -0.30377744623881369757 -3.5967557830119285533 0.64559453303324221274 0.22875524345180175967;-0.37242670635638208676 -1.6641505598836239255 1.3321026957034249882 0.46591925625225882879 -1.0795871424748169076 1.3378033780454998691 0.85834772599418818206 -0.34453701933363173326 -0.28556365816100198751 0.5313576378596340799 0.62353460499112667303 -5.7509702019993325095 0.04979599175679380596 -0.56646237436878799087 0.94100931640558627134 0.51010186191590600746 0.024733983377095709877 -1.014858600859008364 -1.3254336692766550687 -0.23775985339633373683 -0.10718253425370179166 -0.86588856747940212699 4.0666239978786666853 -0.1540480003601143455 -0.47830213273144794872;0.69516227441383471319 3.8270044496656536737 -0.28140891325490435815 0.87709573191898893008 -4.1787709606561938003 0.71716597448030339201 4.8990108825526403891 0.098815854458214502642 -1.9157512217279399902 -0.079295515080828687937 -0.39552361421181730927 -5.1779004443534235591 -0.61571566100074304195 -2.284519112948331987 -0.087151017538101172022 0.14583235037752545482 -1.6668775572299947552 -0.23250525539336994174 -0.45619464809645954295 -0.33282858488945993303 0.24725412137979341232 2.463831378758682078 -1.6146033805817745499 -0.78352046347295156714 -0.38870422954401362414;2.8327052201718525559 -1.519632440887651148 -4.6710108818464810554 1.3260869532968968532 -4.5789809590904173575 -2.2031436551334726559 4.6653577296666108865 -0.37981914885840695195 -0.26566649412733994717 -0.58627006789378122154 0.79581696581337524332 6.4792318204595549247 0.13717194235056540119 -1.8717111996092017367 0.42892754542381572591 0.62619747922628377079 1.3337827006515046158 3.1505895106249277582 2.0554806746007052709 -0.40351409214866629238 2.0486607007428592198 1.4965805602847246369 -4.6630691505750210268 0.27624234770165262409 1.4857239706660547718;-2.1061570253156016719 -0.40002146880731720424 4.825012324298606714 0.90360058632698048164 -3.735331198984399137 -0.63848602300452250713 1.1325133148910593306 0.80377194322243239455 -0.92273007821216823654 -0.97249145773915601776 0.18554185366322140105 -0.4131641717795031532 0.90589659087749418998 0.16414686985587503831 -1.4380761383015652211 0.89543749623363888546 -1.5513253097239123601 1.284322604273339774 -6.3508285669797279738 -0.39626732019714944188 -0.69928108636975094825 1.7577057023182574014 2.7162508763879089635 -1.0974738134792807465 0.49775201897265064499;-1.890066454430523013 -0.29505372259336370977 -0.88468431725453466896 0.82771627918591250683 -3.6879814882585701952 2.2019707695841663764 4.1963717035975065528 0.024508173824734250795 0.5961193464704630518 -0.71768914461299826435 0.11040382877809140627 -1.3225682522074453384 0.30866722276962710314 -1.3980189895554566704 2.2832240324245067775 1.2532243353937986274 2.5332752306502346329 -0.19874217310934932934 1.1001334831825415606 -0.84180360508666751684 -0.45904815802663184732 -1.8368983552032449591 4.5611091775082526212 0.4997230465484928108 -0.24726256048832040424;-0.30491165902252370667 -1.7014341006587339056 -1.9126713277515290379 0.68369846088863261979 -3.6331972958059828294 -0.35137651422319271521 3.830389288802880543 -0.023915261418220611078 -0.32512328764031767081 0.071686644379910999247 0.025071094407927572295 -5.009919718984670034 -0.42860311884850249564 -1.0430208643235805432 -0.81725577532083915067 -0.30344985224685755432 -0.8741415520289210761 -0.543687987495818148 2.9700590852066315328 0.22793387438717513382 -1.5345926802454938898 0.55720254959801818195 -0.64808580041552865669 -0.90622138279404496064 0.50288174352775172338;-1.1279120096323014266 -0.13455991457755031582 1.6385767091100096327 0.092466429032613234318 1.0407643834789905846 0.15132846423859086515 -1.2091375033618785295 -0.059544900287288399232 -0.15277137782830566004 0.82551169441972338525 0.11415618097996327318 -5.3234310505892663556 -1.1801510157938088774 1.0475476842392090759 -0.20949295063674114314 -0.22270387597474430863 -1.1434612830068637823 0.17964624621260069648 -1.7496966318139708374 -0.35964504854419054869 -0.43678990876372658647 0.38180045382725402714 3.6242728999123485778 -0.90157535549449652201 -0.2431672919346536399;-2.6159540255696289712 2.0061345328967057355 -2.7957413195093501201 1.4405734477585301345 4.2939011387520436003 -1.2058717231718683038 -2.6445621873635829502 0.70175912512871196913 0.66574429283427005721 0.065428644121760665131 1.4442597524521125774 3.2788896538520044821 0.51739603313320259304 0.011712235804841727646 2.3116049580438473576 -1.3816519335195398543 4.5439379752867266049 -1.5745879487081257952 2.2699496466995578992 -0.40450386225784518146 -1.7221014617555123749 -3.5606421297334134302 -1.7696096480201997103 -0.029523775667979347548 0.91227306723966361535;0.93813896191077594811 2.850724626544338669 -1.1529360092226998979 0.55851567012810332535 4.4890928140825030468 0.048760058255833606256 -3.5671664958780753274 0.73182542392756511074 -0.0032852779026950901936 -0.38465957471964151715 -0.34412319713376249286 -1.1386040273957243052 -2.4533900097765850568 1.8281533572475598959 -0.36378193104009809433 -1.3609421804762993347 -1.7439869285544260524 0.094925445668592567383 1.1850278980130641848 1.1069391526554477334 -2.4960544948564020729 -1.3125252321545686218 3.5460778792549740857 -0.97571228427510292569 -0.40021997694009248114;1.8489647754353499831 0.50205611280573392285 -2.0102924302163622805 0.53877167783264612311 -1.0216801373995632041 1.3358974906848053443 1.3403382355729673847 -0.024566551385750780984 0.59934613570498418422 -0.30253868601298111729 -1.3095096514429325207 -3.1560225355232551792 -0.47131839641117323492 -0.25695296014532031936 0.4661290852420413855 0.26702331403351542427 -1.41560829681862832 0.53888630584144858737 3.2668167887297334673 -1.3725659953228019905 1.1802075531253186647 -0.016771389327689803506 3.0925172782603200794 -0.979632980527568753 -1.6780902472121019819;-4.0579204309168144249 -0.48435972226293821929 -5.0284465734156489347 -0.65033598427791905472 3.7328994780329951553 -0.19797863904234697108 -4.3841659679379363368 0.45742071137164985739 0.095748951910976251578 0.99873626963896255937 -1.1073292530010494517 3.4790583175735632082 -1.2096267617061455812 1.3302455477572083797 1.0549940887482884833 -0.42251969873968198144 -0.3154001625280831167 0.65437331322351444118 2.6600774173491661401 0.97044299479730500302 -1.0216784528105298246 -2.6296951012715812368 -0.23141393226356010704 0.58212692811820299088 3.1477566846481481733;0.063552830850222508108 -0.58994556077921866954 4.6458509750471463917 -0.016712367912539122933 1.4069277567965305931 0.16393359715989797021 -1.0505809421455751806 -0.18886125096380115829 -0.9339122208654082602 -0.43744954092935522949 0.12643962285183510397 -3.3857065870357305748 0.70821316167405279351 -0.34589004630558739528 -0.62219136685003606502 0.11620261599076839309 -0.33476204448949536951 -0.35925436054129317176 -5.7758058678644790263 -1.1550525630926828935 -0.86561097824754684726 0.19750265826942797753 0.81770012792790147671 -0.97471908490851943974 -0.09203480886681313089;-3.6790620809195933028 1.6884281998878578612 -3.8401242950142080446 1.782812389083885396 2.3257992256496535788 -1.1540345047841671366 0.50355252898344626811 0.80059500542617378915 0.034438992993809813825 1.6283966109601777372 1.8823070145941684927 1.268631723566308489 1.3021277929926178274 -1.6823015453543435083 1.8915419934669452928 -1.6313265198516246901 6.2151360793524483839 -4.1222130182201652815 3.547711003975280164 -2.007311510962777934 -2.005754913980255516 -5.3004143626632638231 -2.0497350033071355391 -2.0675032865062017429 1.9918968981068920776;2.1162220105787139524 0.63369383008981994276 -1.2519229870261840976 0.67984830967704223958 -2.1740184037623238034 -4.5936071855158173705 1.8728737078871715571 -0.18741381710487151069 1.3126092620815623668 -1.3205631150479089087 -0.097131377288053460939 -1.1348695161199127135 1.4515829715787655907 -0.71567957545551386112 1.0703615551717360344 0.71050659301122620626 1.6280620069627140811 0.62398789743995253243 2.3659065677191355448 3.2787778945562746635 1.5969779895544888859 0.26201979030031935514 0.61608029719771761901 1.8983696831944738559 -2.1153564109230136303];

% Layer 3
b3 = 0.052826337135872744599;
LW3_2 = [-3.3512773508583646453 -0.47865309672908468164 -0.38117118510918274188 0.21371665474501816417 -0.61310893308984026895 -0.33797589682534240341 0.46595103072128701838 -2.617676767897991752 0.47487748582949529386 0.30900551894188282853 -0.50584672684526776898 0.20666401361588485064 0.5582258291681487794 -0.27315691333593761225 0.52972662473353970825];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00454442172233583;
y1_step1.xoffset = 0.3;

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
