function [Y,Xf,Af] = ESPER_nitrate_8_Atl_4(X,~,~)
%ESPER_NITRATE_8_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:35.
% 
% [Y] = ESPER_nitrate_8_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417];
x1_step1.ymin = -1;

% Layer 1
b1 = [7.732135502987394382;4.4909257988245423476;6.3287957715390401603;-1.9286583724153230079;-1.5364212283076801313;2.6552059745060514473;-5.8444455124214673702;-3.3559460257002968753;-4.4248600100554522641;1.1125950799717261308;8.3832559873881127288;-2.329909397346181521;-8.2777639888112979349;-1.516429211792922116;2.3448033074779011109;4.8879811818067917528;-1.2734827091671021559;-4.5469895245099252179;-1.264724328513843199;-5.4228227702276337041;-0.92530182710177721894;-1.9675952855454739154;5.266748901548602646;2.2421559048329284813;3.0352140003545478031;14.491231934583462149;9.384295277769789223;4.3724229916624066661;5.8003591178217606128;-0.61034318460663883332];
IW1_1 = [-2.6461405623293092759 3.2373473257247473711 -1.2001488208561998938 -0.50196040451081291867 -2.7085560634037735284 0.59002369283343025774;-0.69192526157832912048 1.7609332268839614777 -1.3086263196732106362 0.88969184387476418774 -0.0470913092895341584 0.68841691538230587177;-2.2399385094052615131 1.3226060971893618401 1.5474795306041222442 -0.71386035747004139562 -6.3054601159623562268 -0.62418881231015543865;1.2983612920513354361 -0.85160845992417388839 -1.3676102169189268132 0.70028682247696460994 0.44624342212192474877 -0.10691629934134712354;-0.26571380799444049847 0.16142729653071544171 4.2477386532988585799 -0.729280317164121894 1.9991729654409882144 -0.26258369043246476693;0.25162743502063167345 -1.5102955396901394103 1.0635728410115188325 -1.337646513577366969 -5.5697464940128131872 1.1245707840709442404;-0.33150597246954766062 1.3601285056353438208 0.27443151020610206947 -2.4246773299607067287 9.2895005599669708829 1.2135545118614501359;0.14228539512154725344 0.19061237721487134245 0.19213911835466226941 0.62936842746603793852 4.5084342225171276652 -2.2406187706561904172;-0.18709079236742764962 0.48376511641850417078 -0.67277601065705305849 1.3204665163134821082 5.9430983049974814847 -0.30304809940692101211;0.41988859953730417462 0.25976107583160323866 -0.6823506376094311543 -1.3962250565229559207 -3.218248565715549514 -0.32693933754654314061;0.055112228350638059471 1.2209595370303423234 -2.1232840540551167763 -1.0637380329856473971 -8.5055824788208980181 2.2377260577495619032;0.29042117596973449078 0.073361718767460851187 -3.8226115746800091166 -0.59185023620838805858 6.0716572994054276435 -0.53173797523027743761;1.7781261268302155365 1.2261121347416401584 1.196175854408331185 -3.3496945914710827275 5.691984256621463345 -1.680586282621569838;-1.4367053765623392092 -2.9028739133773333592 0.91744569814472054503 -0.16097348859990703862 -0.056845767127531379725 0.24109595991134710746;0.4022287146376391509 -0.61740089580544921244 0.99320428620867717839 -1.3324762942818559885 -4.2933487017997657986 -1.1239720061684106422;0.44637125200634869548 -0.83508986722322442464 -2.0664238385563309386 0.19942039387290338981 -8.5560773101910356786 -0.24776393612108352227;-1.5971888106629457127 1.555025513233818657 -2.4569789257311009933 -3.4687104997911513138 -1.3134249356403546738 0.86163125989308886332;-0.46949617247666125541 -0.17405714779678344017 4.6254314826049052201 0.73045102776093906449 3.0338305765566708772 0.35513168268917644221;0.83362186191287235282 0.13192339666255092911 0.96912453474572424206 -1.1572464381028526592 -1.1475980164513885917 -0.81808082346535382356;-0.46387190379788306727 0.015753895083176309988 -1.8769354043106452057 1.2254260154602867328 7.1911010961739192382 -2.1818041553181384273;-0.77397701400226093238 -0.74922367899806308067 -0.32208360201347641105 1.2203732502062392129 1.0583205743059547199 0.49117655780659685849;-0.25685538074801672614 -0.18192801568990221717 0.3475697098919240946 -1.0095363678721334111 3.6102195305811721759 0.38237234189514146543;0.37083666242940915847 -1.0446563518161717798 -1.398711926192481636 1.0815177966165467893 -8.1103633474157881267 0.11084038634310725002;2.0367304290535552802 1.5315207281650489168 -2.1699998361678276559 -0.29318504004759871062 -1.3596755027985227837 -0.14675510105393776983;2.585237235531463007 1.842333010799179327 1.0162375224548825248 -1.7359373640751658563 3.1660972290811852936 6.4689905708891819458;-0.015286472474123485865 0.20461032719351007758 -1.0449257061835641736 18.141758236716853503 5.5894989324373840489 -2.2212831059282152069;-0.34303848521212698675 0.20760104915912069856 -8.570940449033148667 2.1267945422871039618 -4.3530154177496278578 -4.0603465248232746276;0.75963210374673717507 1.5120443064225010943 -1.9209320044074968958 1.8391449212687625003 2.2533700659736108918 1.6333236460907041288;-0.18918511939793208887 1.2414027969053602618 -2.5962487226049502098 2.1341556501268112989 0.46224367864093540703 1.9566181254496857189;0.0051407073503125181452 -0.14530549884988994136 -0.57959899094086764215 -0.10585134368144227812 -0.40269123820779895162 -0.86633851348617330057];

% Layer 2
b2 = [-6.0796911713504151464;2.7387392722492682751;4.4655645816648199542;-3.9183982972464064254;-0.40799482452173269076;7.0066584005563612436;5.8135742633003513191;2.5674469862078685978;-2.3474546721918096992;-4.5121717921911228544];
LW2_1 = [-3.4475691263605057202 5.5792308768006435571 -0.54926551100120679649 8.389659479798952546 0.6793118793814089118 9.7281301504988668682 -3.1798209720050385485 4.1355379161407057609 -5.5093847313538111266 -2.4906044173261467911 2.1563863998569336999 1.9257677810792719342 6.1844736138862588248 6.072597739124264038 -8.1636984615128813658 -1.7011347685974698152 3.3968278940232656815 0.90302845011363774841 8.9766106050932918947 7.4847851062659387367 -5.6613333510473760057 1.5450115254404888088 -1.1845688661611548653 -0.45796317836780831012 -1.437166130035161693 -1.162244437862921087 4.129179582037499685 4.3044920636512324918 9.3930740358003976809 -6.434040908711412321;1.6004836268995314974 0.86443051475398358097 -1.3105866089819988574 -1.498133898318823265 0.29247966333888480994 0.16009547063889636997 -1.0367531268026424218 -0.64176016347747977253 -1.7789642471807112223 0.16472902756051052697 0.31361953206581644205 2.0515271558691723364 -0.21519670567253804117 0.55286125727992396772 -1.2396581106875748102 -0.46008220394505333628 0.73027741004571866679 2.0725122905241932791 1.4994447889253617046 -0.92712197771563698012 -1.5592609481375960812 -1.6137315991375817958 1.3842500659242078509 -0.27713272210324990663 -0.38304304291644530878 -4.8977933572614125168 0.46808035849161694353 1.716078396152203478 -1.0030317497550687023 -0.22592995600915749543;2.6799399986730438883 1.6324499552177140771 -0.46055508265227873554 0.3845451161671734619 3.0601600214524333765 -0.92063299419492117259 1.6612746747242006506 -1.876346593947630037 -0.73589216343984042989 -2.0278846561914196656 -1.9893229485158090419 -1.9150134824217710516 -1.0383183428749573896 -0.77335313318136478244 -0.8436805876844931662 2.6447539675505975509 2.6574348827365508185 -1.0384428044901286015 1.5518699692168917981 2.7136652439142929794 0.97553349025157609642 2.6345724097757448057 -0.70286787268608408308 0.70509053741572635943 1.005740023159670038 -1.4120432312240829287 -0.16096163569630667078 -0.014186975966314161607 -0.88911734737262437811 5.6204941461179718232;-2.540689263779033702 3.6229595018754245928 -1.9493365497955670307 3.1034546901769934912 -3.3279159864492804388 3.4122438624515916139 -1.9901692369523773518 3.6657998549570032232 8.3529611501852603084 6.3890777637670757017 5.5534932235656961907 13.24301165123877766 0.96520694354531844095 1.5296838065067395185 7.5063654802267603117 -1.9700515287655580021 1.4409699423851387579 10.86374217873397896 2.9875659522361037723 1.5503049134528494868 4.6471725452120091759 3.0124725287649174454 -2.1413621830392735035 -4.0779715162570013831 -0.66370054340138129323 1.3574370358139289383 0.19235696891086734306 5.7700744065009352823 -4.7330206199026960689 -8.2526905731352169227;0.23927904672324137869 -1.1916598943295437785 -0.091533083814908322573 -0.62136911339765710238 0.058751036586836162534 0.35098603729797422934 -1.2953815984932961491 2.7987869111045822201 5.4526012453476111475 2.822014095467147321 -0.90940461772335523527 1.2839233982898776709 -0.66444067749368096454 -0.31148631414880650903 5.6824138521575333627 5.5693927402972605023 0.27745541664230216616 1.7356291203282705471 -3.4310759955564678059 -0.28642646338578003906 3.4605989612289724811 -0.02262898448026488471 -8.0906979676908221677 0.39252696840724765615 0.84403580698499724999 0.93005646889481530248 0.10609808960859759552 -1.3707152695966069444 1.0391970878855143479 -3.3570516726203964986;-0.20618363343261092435 -3.5426449617582602691 -6.2812317532121477726 -9.5894336669208186663 -0.40319963153432569758 1.7377997804477853805 -1.134285479567030297 -1.7418947227831134139 1.1130355614527145391 2.2971731462481481678 0.98094189189864644529 1.5148731988839756468 1.1725542023478952025 0.88671365386049316726 0.20640360440163280087 1.9721481121113209056 -1.3585525173598322901 1.492435297676563799 -1.1488422508262157162 4.2173939346787117088 0.11272169075291255114 -1.5565477733522499193 -4.5538772499257129667 -0.55033813132597930728 -0.019562900569640301329 -1.8483066943827488604 -0.80804258431421538411 0.90697771680457184296 2.2306757069796931958 -3.7112641420503242884;0.84448755564081290359 -1.9826012324456518954 1.5435456569193442 2.8257643767907238974 1.5718821884628950691 1.8042416743254152678 0.16696762115702956653 -3.1758128583152527291 -5.0345079683299314155 -4.1577692571188409332 2.2364875745284229502 -0.72968531003648295918 0.28701177272397226936 0.71737639416870868647 -7.1912618085494033338 0.67655868565466170672 -1.2351209354095391557 -1.3582424160818460379 3.0722854272757702532 5.3293859569845611901 -0.61108971269301592955 -0.29591826622573164096 -2.1178777128089465442 -0.4254429116285613488 0.13296268789099724428 -0.2209389342708000803 -0.096161189767295929265 1.7212879341264777366 0.83849242253385614099 -2.7206950110994250203;1.7567857091150786975 0.66903014999680887609 -1.3319048734568494918 -1.6646939149363408905 0.23886180281704885719 0.17689944019285602783 -1.0319904769492698193 -0.69516431274237333593 -2.0009146022187103142 0.12151520947437161047 0.29405161623547809713 1.9412438277132380016 -0.24188970757582084481 0.52167221724115309112 -1.4493435805507566272 -0.90427829438283757746 0.75618094612921427089 2.0796092408178741273 1.4402529937285344364 -0.8922131221667011225 -1.6196311289006022349 -1.547833948520236147 1.5717575367394507335 -0.2018119284188421847 -0.35999954705400483146 -4.9399247101326793086 0.4573003745388125596 1.6378219777211346919 -0.89133750553181001486 0.13111790903807823971;-4.9968733022581028891 -0.50619470993251658442 -0.56542973898323645177 1.0192340813384221931 -1.3011760580505649099 0.4499503672940913801 -0.36009654449019579125 1.5513078613232944569 1.5593101025909301693 1.663673304434958089 2.1239953652071745793 0.48530685708240844933 -0.024718813110817899198 0.31103623384158024523 0.49741231584017769896 -1.2678553817788120739 -1.266007499410729098 -1.1406326338273458099 1.0000542383873414209 -0.017941147048778888218 0.69259172268773716574 2.4999947759144687254 1.2378218462501067165 -1.8660165054940125451 -0.15755382417839747244 1.0527750228675460953 0.12273703336516897378 -0.65472073042684442168 3.2255930263693279869 -3.6593457562445177444;-3.8046722589703820816 4.2325277145384818667 -0.47242675433945718888 1.3048755951958734833 -1.1819515345154387731 0.60845820582206877614 -0.55640153359684973999 0.8023401884002021589 2.0719450660547593657 1.0474683112816198793 1.5312932603512618979 2.7218210620217764273 0.057042204471795983356 0.56367352944248694246 0.83699163067650705905 -1.123276073849729606 -0.94599486489892081753 2.0966081524871165165 1.7250125186204805505 0.044766607781583742098 -0.14420636431032121672 2.1690074328061195175 0.71233451183865681422 -1.6117303344921682928 -0.19792742111815198047 1.0292592568363034822 0.40590035638993515743 1.1825295012933676464 -0.12864460294388282979 -2.3999857770102726029];

% Layer 3
b3 = -1.5747637757402586089;
LW3_2 = [-0.034701834694803673775 4.3589546788096473762 -0.59044534886496091719 0.33077013938508265101 0.38478755982734214047 -0.56729995407278455133 0.90389995493324548637 -4.2980041284640053689 1.8742496418086767029 -3.3950020603802402164];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0470499670650231;
y1_step1.xoffset = -0.28;

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
