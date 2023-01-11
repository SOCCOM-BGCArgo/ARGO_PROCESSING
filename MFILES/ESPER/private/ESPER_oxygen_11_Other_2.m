function [Y,Xf,Af] = ESPER_oxygen_11_Other_2(X,~,~)
%ESPER_OXYGEN_11_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:47.
% 
% [Y] = ESPER_oxygen_11_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-0.9;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0415843642790311;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [-4.5471121035350998696;-1.336268723292824312;2.9486255923845905436;-1.9312343150948991788;1.8120351470623146106;-3.8829590814653922592;-1.6409325285111013848;-1.3287189220320514416;0.65193351192545201389;0.13538696935052782488;-2.1237787641015635209;-0.30222903135263667895;0.55061895144434014338;-0.11241796212342726469;0.082193912388390660828;-4.1758277109772254931;-0.99947398673248410983;4.4018772668640719914;2.7403296825102025913;-1.9775849375113332407];
IW1_1 = [1.4446880880071066855 2.2555303521422902513 3.6259033758060237318 -0.52478347167744587853 -1.1333300488883693014 1.1578727284428043109 -0.96198602445110259485;-0.33366160653252219248 0.22759454112980570528 -0.28378794529947581893 -0.023585337295942050001 -3.4362401326822369363 0.027255198522637053005 1.5055932784635230881;-0.35131233110032183964 0.90982636946191053351 -0.75994662272063107089 1.9447933268788979433 -1.3584892871623304167 -0.11831281019714805225 -0.20785129950980088531;1.1024854379951145589 0.52679242185649022279 -2.5337385914652506003 -0.033740161216926981003 -1.1387079459252285041 -0.39169176437580965899 -0.052121438382883071572;-0.86670314925303471032 -0.27802716936787136293 2.5248010139658045858 0.056068510454155798262 1.918621870448443234 0.41454303314828544469 0.080341354508410484092;1.7141736190441629883 -1.3871932486712217347 1.5928766143075632211 -2.2888901503675112004 -0.93149383487474579635 0.090584482322423692802 0.22426700801703977617;1.1106063747903098804 0.49566449921667266043 1.0664562220999334041 -1.4019901758049782003 0.46847371145787836122 0.42630027020836430696 0.91070481790473090111;-0.097004934789832214515 0.59656009890040662125 -3.698567933725391832 1.0003148225453803999 -2.1880318685994990524 0.023436219700113107567 0.040214976568679491886;0.41686854255695648952 0.29043081602595316459 -1.6329653975028350921 0.81150596212936265683 -0.73369430726865636494 -1.4782407252750511262 0.40492604297809164926;-0.35294675459103869164 -0.13484231328177342935 0.82279594353589313727 0.82821986926571167231 -0.22775395671563611333 1.5716248885690020654 0.068252388913777592427;-0.38861975098468121148 0.10501339552621508089 2.0364337933142842019 -2.5331580867617096153 -0.33795330905164533419 0.39478348690552717004 0.5522056480462399497;0.26402640833092705375 -0.14757907219030894086 0.021248307485021582924 -0.35495575820792135691 -2.3647055261353484568 -0.97071987149993244071 -0.014441450781271467432;0.39764233537086957071 0.011984035191081487959 -2.7061007022157550495 0.65220274599869521026 -0.38165226469171359414 -2.6963939135506205247 1.2866681285383396727;-1.1101886235996052399 0.48202497282814710244 0.93349492842605508169 -0.88189376314176570748 0.1975528403072201733 -0.5220680357102668312 0.42576653363327565671;0.45225458946122376602 -0.2146437973601905691 0.25574605702536679086 -0.10541633273025680095 -3.6119815687732153009 0.50853392562461741111 -1.6773381111302505531;-0.17824166101953778374 -0.48199712397257671315 -0.48812644341812000404 -2.0042105618649972421 7.2419636451137847644 0.96532222006095436129 -2.6624089573642106821;0.33330771746686266299 0.20992322133497012815 0.71226365263701085784 -1.0345208949324373737 3.6746102540695555483 0.8764442817938785657 0.89207646055308276178;-0.24563402107185919321 0.74401366068437113643 -1.0927013547933197568 2.6492526540806924373 2.0297708089852113389 0.66883898051586787581 0.74813319800320521313;0.17018002141900212543 1.003888432736493419 -0.6456072917113595544 1.366851133412472219 0.22672836234282184309 1.0645993680354342015 -0.13359235366297500525;0.01668118584054806286 0.22604856549500090823 -2.4146879725101166869 -0.37439485259334726042 -2.2688568246054341415 -0.041970706863130494146 0.21166484715774369629];

% Layer 2
b2 = [-2.2338612866961788939;-3.1962881418396382749;-2.1390170162261865627;-2.5691067280264787698;-2.3286528693815826685;-9.9406632643883305889;-2.4356531586895653341;-1.6829426598777359469;1.7276377896430930026;1.7923552387130148311;-1.5160474249419195747;1.8210090571795496395;-0.39221073576398868266;-0.53002409913324888269;0.58921841413629438833;-0.92034568569868180354;-3.8336072862526195948;-1.884008851429420206;5.5645182132664530172;1.6340542806681512822];
LW2_1 = [-0.1379775361714705717 -1.1557777750641071712 -2.0100468527289736542 0.68519048840037355408 -0.16476232087283038119 0.19535160917964927707 1.614077896344582097 1.440810898348939828 -0.43030986606375093428 1.0044081515999152554 1.0711677181569894124 4.1571167349213133946 -0.33072285924459038364 2.4520056493213830962 -1.6103852672530305501 -0.29758680141914278794 -0.20027261335056634195 1.7603549361941079887 -0.10666422432132863229 -0.96147811127155680122;-0.61164656140228368386 0.087161019597556066429 -0.3024073804460587489 -0.068021682286131424089 2.4577457487629543742 -0.40687550158906810038 0.43587840407809996801 -0.037541611347564507462 -0.57956928903899129768 -0.75549307238240215323 0.26263976090839818989 -0.70534210747695702537 -0.23045419980979978125 0.045552149953519915915 -0.47650982217535464702 0.21317849250859557952 -0.25226629832837432765 -1.1018629977773961581 0.47447654215290979707 -0.14600824242258675367;-0.036508920501035567463 -0.27210454574549103324 1.5820202097929521035 -1.7218315879094368004 -2.6454404207329500487 0.39226335936514550307 0.092678834478773175487 -0.034963578799324435131 0.14155462678956964861 0.4498284350567502976 -0.67357175514727773979 -0.065367332011679865622 0.40833886873214625401 0.604955275270150028 0.056546472053997612339 0.016409481472178903683 0.93621369586025404708 0.96137534032598959133 0.27350081715325591514 -0.83969439884292351728;0.22088783794815716743 -0.21289892142172697764 -0.25026497714274853479 -0.79834130115967294472 -2.6734702199574096326 0.83427927413750080365 0.46745223272733449349 -1.451970535719907085 -0.31658692791605741457 -0.92433500422265646446 1.8426262347617492487 -2.4118596142853441222 0.12775580437300726455 0.14316257919110844665 1.0795704870892426364 -0.088596399305059275742 -1.8067957572031034807 0.46157014469745227814 -0.10160754541734426382 1.5814430900872378238;-0.20405001903124500573 -3.0610326065653272742 2.6837534919698620151 0.014860466644864614091 0.61094799263185950444 2.3376651313314558145 0.20865517086919554979 -1.6549711953808248222 1.895397300534662266 0.51246024632380415831 1.1428258527192032279 -3.280373326174607751 -0.09431033728863867649 0.28602499582293622327 -0.21886617408773836879 0.0010238420268877241031 -1.4030327072215666995 0.7829719607970418993 -1.7199536055027315751 0.88610497910441843317;-1.8341112182240377759 -2.2462605456469777465 -3.5828244672174802687 2.9462801000984755184 2.0646202338395083231 -1.1015122456870722889 0.82579780070702202277 -5.1561301805122612407 0.11174042430045603991 -3.1934093505403859936 -1.405728498281676675 -4.8714740802291478516 0.63425615230859666571 0.45582882549821990281 0.69147820346931188062 0.56769292250163294611 -5.9750083009870476047 1.2565966269427739643 -0.13315065110067855758 6.2779999896392437719;-0.074722723613359520733 -1.3873787976641456243 -2.113487667646868573 1.2857324746647993408 0.84894060407051508133 -0.08480554391999231556 1.6215798601072191065 -0.29780828714492973974 -0.53025324905411197474 0.92150029307692882252 0.96999778064637776076 4.0013118132935625582 -0.21070193026559849714 2.8663848977388499684 -1.7808175102079666363 -0.45140179743218677988 -0.30502097252491505719 1.65321433730844225 -0.16722953007543286952 1.4708280044021135868;-2.3794178441313071914 -0.14782152937681350413 6.0507010132967291582 0.37256109692778804909 -2.4208806278381840293 1.8137969348399241021 2.2131350470694823862 -0.24166280225003758386 -0.68866691595481488442 -0.24964461628681750005 1.4376737900181468799 -3.2510308453415781571 0.75414359068438074907 -0.58582623872035888279 -1.2359060239664989567 0.5377402215698438015 -2.809683508110361938 -2.293701735231085781 0.11673486208788243224 -0.84643826376050101956;1.1113620404764184002 -1.9431457673010490428 -3.2834865906217611808 2.1621652520157050859 2.2805732451825648788 -3.049975903715874459 -0.26182028601120282962 0.71652367687894857884 -0.23343267354137514635 1.8733132366832503379 1.1092884689508175633 3.2767922935518409844 0.070564283803905933357 0.27341470903754228061 -2.0092686255581706511 0.41837963007054196352 -0.46357305747552585151 -2.8838909672998163103 2.7691387308552704027 0.72610172989393140242;-1.0355035574123079289 1.3006242761062118962 -4.5029173252505998448 0.25792020183045694104 -0.077137133813771310842 -2.5227305152973520386 1.5617990721473951066 -0.82486566016314954997 -0.89111216703423035135 -1.1109969116983795523 0.75449016867232576722 -0.84695269425677166542 0.83967102034123830556 0.72674114639620723821 0.21904611035286003506 0.36915385723625571135 -1.2386973624954153994 -2.524670014935028739 1.2516380221211611534 -0.69841441581482488488;-0.97232236292500118235 -0.7731156042164292197 0.57842403204031556108 -1.1415041903054525818 -1.3679859116264161401 0.054988286549584451668 0.69062748139165686112 0.76492187291064428489 0.2265148396963054811 -0.77700255970112808157 -0.88812630197081254835 -1.1992076299445570964 -0.27912625549916614265 0.40847027452653583701 0.059981563039415118632 0.40238689919288700425 0.4884709035607722738 -0.16885016218978635139 0.044042563623829057706 -1.5255948771790814877;0.49648707562420452133 0.28878214159130510819 -1.6885953573103618108 -1.1194994156264799212 -1.867632764475174012 -0.27322057956552964209 -0.55183727647723301857 0.56308651710517831912 0.074623694786053990091 0.17251522882650291435 0.61068820645655308432 1.3126328367442989808 -1.0532857050832546975 0.085883065088140669863 -0.029508133439333131265 -0.12206738328941876826 0.75662346820418413529 1.5777350618948058347 0.020263908901884866015 -0.26480087917299849343;-0.27951441710183599909 -0.6056002596111031977 -0.13249789338153336593 0.29733478914576888519 0.54099287922892513869 -0.99648956975489499932 -0.014520984111153271284 0.34847004121595687165 0.74587492815740352725 -0.42051321450504969279 -0.66978562523854789923 -0.35276765347555616437 -0.22868848481121620275 0.51534623164025972208 0.14697607498870179921 0.38616119045516028585 0.92829492970613325387 -1.925173010661782147 0.89291546587583880967 -0.53806991079610932971;-0.88065808231123565708 1.7307269742844937355 -0.71471365262608954616 1.7609066779693025584 3.1877771628903448509 -1.8222241912781671491 -0.97502656644581742196 -0.741506686257518699 2.221210500882905059 0.88122001437015029524 0.73565008229074324575 -1.5636429958062973977 -0.39194714071564079605 -0.23922830061362246323 2.1099936550393971757 -0.13615455316849098266 1.0580710470214078622 -1.222893408230173673 0.18492682009077748839 0.93329698113150860728;-0.56013343530974801787 0.51727841208935754214 -1.4985373957710705906 0.95416440354839215932 1.6635079327941160177 -0.9556035314694100169 0.74823483165286697183 -0.12345772140804470163 -0.68786423171262767173 0.28265400622508884965 0.47412046564473331323 0.17326513183907737981 0.2209905870108282433 -0.51938484589986622098 -0.21503928006161054931 0.041080712755219836452 -1.0474549688361023403 0.07784197428945451025 -0.02302942859791032304 0.86069394632039619708;0.6694279221707659655 0.72547805892160965868 -1.197314251548548647 -2.1374871941795010599 -3.903158819536073576 0.32675226573301385802 -0.54368464702809737865 -0.23250178291033193889 -0.3975874585086823787 -1.6922027208456495995 0.94409429798376076626 -3.1106096519819659285 -0.21559710225379510873 -0.56983559906676350959 0.0090494110502771368809 -0.20873859706520989099 -0.96084851698188822322 1.6135371680403558514 0.34629843069232740804 -2.0475909465367942985;-4.4586641463938176955 5.1443962396678895033 -4.6169154104218179668 -3.2986717222564521101 -1.7232343263867322403 -8.7819196602670359653 5.0168914026658981697 1.0221201421863104208 -2.4813582659220712046 -1.2271401655204905445 -0.46430352877982039583 -1.5984906142181491973 1.6265000375109626152 0.59469884948906115429 1.5173610945028916941 -0.71609550145024614132 -4.1730899042612632499 -1.9870812249172722375 -1.2728963213807393817 -0.53105569723341872201;-1.3207375713226445768 1.1037424694241033674 -2.486988734565117376 -0.090087781336208372363 -0.7896521063166873633 -2.0256000697314169301 1.7178269484264439004 0.25471516714702246986 0.13040748987258440894 1.1183019083611047684 -0.045712405053163950031 0.90618244165555617009 -0.34207336047584618877 0.65273194596841310755 1.0535161564770729736 -0.003027529551289220161 -0.66106191778930512015 -0.35254515622492449411 1.5380400831430731667 -0.50195860127538050044;0.92670564431647328707 -2.3990240234574304878 -5.3812371660455475819 10.419390001271491286 8.1801137131773149491 -9.0373655302809421386 -0.99943267332012641901 1.4892837041549473742 0.77085107535332142437 2.1387114508185733897 1.8802515924683960158 7.8373516407302670217 -1.4560008927092713549 -0.33220007906906418027 -3.7630745322760752458 -1.3266449703660536663 0.046753546185306069094 -5.5175847956743799827 3.1450055365196103452 0.82252439027356560874;-0.84673577974836278237 3.3401965731467551102 0.88367626835802182139 -2.2076216789362961279 -2.9282066767049519385 1.8528634043331027037 2.166822059893538821 0.45724127366606454048 -2.6687722732078578147 -1.176052349837653388 1.8502742987639122507 -2.4461297839688009859 1.9950004298341112996 0.034664635998901384417 1.0214774096266427961 0.0057630631385899161456 -3.8007949530788258663 0.86992542704475706028 -1.2994736741404606661 0.098471486861423457304];

% Layer 3
b3 = -3.6607346072335977816;
LW3_2 = [-10.080081189487792415 1.2548112299897913324 -1.0085030894703963522 -5.0325689876410315549 -0.50641956404765342548 0.38331273443123231592 9.9429552082120338241 -0.16046078779697897643 0.39081579992905957077 -0.6558692767458822459 1.201532295050711685 -0.76168156532139030812 -1.551958272839417674 0.36481496488605402728 -1.7496540420894546308 -0.38904712740449154396 -0.11901051940648170457 0.74955175398535822673 0.14798688078765454001 0.37011907183300174484];

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