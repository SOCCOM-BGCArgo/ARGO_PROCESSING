function [Y,Xf,Af] = ESPER_DIC_9_Other_3(X,~,~)
%ESPER_DIC_9_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:13.
% 
% [Y] = ESPER_DIC_9_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-0.9;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0419287211740042;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [5.4099338162050596424;-0.47768120905633071338;11.369023029132936387;-0.11916665094780917733;0.035375991944639140208;0.42325254076221208166;0.30421542061255690825;-0.40011735647271212679;-2.0765976187607408932;3.7099235090565190198;-0.17414076436129000136;0.94375153805045575162;0.95586290460020462145;-6.029015600144391307;-0.046210274091794831253;-0.079664896278142105501;-0.68626913912638409787;-2.3073356529673398896;1.4591810539189420393;-0.077605119541036077235;-0.87033490785631917319;1.5765674027032292592;-1.1582738775585725843;-0.3468452517426479087;0.16559184069830198971];
IW1_1 = [-1.7795488298175041741 -2.1593968883598684272 2.364817693110131458 -1.5841788439301218716 -2.3272599386750365547 0.28083333453976089267 0.11177922097821042791 2.0682241396542280221;-0.22964440746853109765 -0.32970256598842562212 0.30347295420147668565 0.3831907332524537857 1.5152004451337244273 0.64495307802032508793 -0.088912872715215071628 -1.154329178722748761;-0.2588001892024360262 -0.12126614534682736424 -0.016119600456046618425 4.6738681443686189709 1.5900024815233440023 1.3595548878481824584 1.5345219362749069614 4.6940764767483003084;0.21799480988388111879 0.24582307291140281791 -0.30781997578444697528 -0.38206439282745857655 2.1031655487837230822 0.39305504851331846128 -0.17256733448123767749 -0.37486112299868695752;-0.11232771807113109408 -0.045027468944748710866 1.0803338222455705075 -0.066681303020858520569 -1.1962752345171505652 -0.20938807878058834921 0.080075320928082185934 -0.97106774598450162905;-0.18009156926748121075 -0.10828839284359183248 0.052919768761231084087 0.38960589500097730919 0.056503440902613839836 1.1074464850337075728 -0.45649705906199039562 -0.61608155837275113687;0.57196250491174027708 -0.80141925421693560594 0.038655532981887764932 0.1007720696585664294 -0.53381402823680268455 -0.62171977274620360276 0.86181156805214065297 -0.70851402889471504221;0.13503936594374463254 0.50547973755026498477 -1.0963528470468886233 -0.41872134640881542023 3.1254455611041858099 0.51690723012820072313 -0.47368485839770263501 -0.010721159571402235494;1.2139250210788379292 -0.13527792009870601153 2.0389717417376997943 -0.31256574726115238239 1.2815720449079177623 -0.16884794192104177202 -0.39032345372160037744 -0.98011277626188342271;-0.13730946165093352307 -0.4467138838080874641 -5.1546130792254079012 7.2324340454266957678 2.1799171407383690635 1.2746752766252889 1.4646614053981903414 -2.665103249772958538;-0.0048539894305627950127 0.085561022528890320893 -0.4040402050462139627 -0.054106466707822056805 0.097839769710568555428 0.16988352821516100599 -0.43406516246196652498 0.45279406692546825264;-2.0195586755378598376 0.66951688970406297852 -1.4223347152011844852 0.76722117302772319114 -1.4240629620527418542 1.0902971119737845385 -0.26107204844059384019 0.98441164438863315578;0.32160847061855191154 0.70243501029852695172 -1.6839089298557385899 -0.31293508213771303872 -0.58507874615557708964 1.2402953913405931896 -2.235359245810959905 -0.86149434600472474255;0.31614712884262252679 0.80514913772238783096 4.4482825212267309567 0.038466658590067821788 1.149011493580351484 0.59867824975762795603 -1.127117954383437981 -3.0276918734073237083;1.2154042867562357966 -0.49775492654030872863 -0.92737461285615407913 -0.012041746423945480157 -1.3470765693026081067 -0.69386992658483637175 0.93157355804625940898 -0.50629232751739272622;0.0016304944260265761234 0.044291162343532940182 0.42080246518334230021 -0.037519302020102750228 -1.0094719928432751566 -0.55515753240115828238 0.48011290452836713483 -0.76807500497707814624;-0.029105088823316829089 0.90551886801517533776 -2.5596443915045723472 -0.36071761930573048582 1.819255289198734582 0.11152601684594815046 -0.6545013513335767863 0.33368589940703868857;-0.43644382287416155775 -1.1591126560178468541 5.323492228394076875 -2.6403626064958798381 -1.7456092175303239422 4.0144415212220003752 -6.6861233991375108587 3.3181362939115062183;-0.15397974762296684492 0.092477633184982629966 2.6394790927914986156 0.65483338492145815035 3.4007928957207433918 -0.99001271854557337804 2.3802255847022459889 -1.3841842334844809415;-0.020991474684154305674 0.1435992079049786263 1.8951605953038035413 0.046453878164327083722 1.4876948212566980523 -0.20062236450256687448 0.65924408764295039376 0.59492647901287964984;-0.26038737475830781021 -0.093437977798301347088 0.82166712322475987662 -0.26089569193674572167 0.80326074645190803647 0.5733942439295007798 -0.8886500367291690905 -0.52062107396881129873;1.2102501766762299162 -1.9250378454992285171 -1.6715216392532394885 -1.3772984788172313309 2.1374291248845636026 -1.3991321343477676198 -0.44521083059301919116 0.13367923346356963554;-0.1326941448624658626 -0.066658315940161291802 0.77234231086066473804 0.59821157501435395787 1.9091584755037216414 0.6374286708656751177 0.58315045988574709668 -1.8618907390153658366;-0.12731040734764104139 -0.18467435031008089252 -0.061776534576598304827 0.28907429083277730308 -1.2061180204484722722 -0.11274786601019265442 0.57559823594714565331 0.31783043011800027333;0.06453622798840732977 -0.13791066557708331342 -0.75493485506288948717 0.13156097615466996742 -0.40460259267637355141 2.2625482107573016144 -2.4991622804586257622 -1.0583675578080404822];

% Layer 2
b2 = [-0.4990240709469591085;-0.74444275305882101623;1.4264503796820757753;1.8942280575237375739;4.5832807772745720598;0.43367145103723891975;4.7108252260674703393;-5.2860921097453639206;-4.3328345294445300695;2.9862408959401087749;1.3050589285023066477;2.6567653947492884114;-1.331352262346797044;-1.2874338880963218656;0.25212990642225902027];
LW2_1 = [1.3607036239944167999 2.5206595808518739155 -2.0196260601291378123 -0.45434339696490555482 1.9048048247706066416 -2.1002115741264053561 -2.3463470284667966581 1.1513355462466159373 -1.763000142160022321 0.79663378723340472209 2.1251132187971566445 0.56731566813921607828 1.2668228787934896484 -2.0020635017486583251 1.4516872067728214812 0.85288726670407333241 0.049661430729046086852 -0.2698510888642728367 0.60229657171554951134 0.84163940334184339509 1.5247863867180666553 0.5230888684766495933 -0.70477186990265572319 3.4096424186788114596 -1.5183124315926652681;4.1714369700639819527 3.57111201301708725 7.9133904332513393953 3.1726373975416737316 -0.047939239282663401176 -4.8993800114782652955 -1.2254600420431074337 -0.68576708899319616108 -0.13520794697277652152 -0.50939170949350109918 -2.1705929947932869339 -0.38181649307188569864 -0.8276089014621372808 5.2566621428537017024 -0.22000479554912275026 1.7135582653534686681 -5.4809874679684771337 3.4000922375559317068 0.42144161551513498187 -6.5128812981115267178 2.4006039772777030272 -2.2259060238976084989 -0.64374446469882151689 0.69375282182674091747 2.0956589730908805613;-0.15334942140896395957 1.0209590921001359121 -0.39691716608740129635 2.4483750218579110047 -0.80648860278055545248 -0.15516898166496900791 0.54483015838031501055 -2.4731009018301941538 0.91751850797585110531 0.398370467015151275 3.1344501552093788632 -0.30415687030656457512 0.84209734275526515468 -0.79214333621176891143 -0.099853233813006744835 1.3818171047805802854 0.53704259166737844389 -0.029084438082523061953 0.075392498131799795535 0.4563934796404393035 0.56199171009959980605 0.065078279487547074922 -0.45274424525760226246 1.0233027155679241815 -0.38226355199447142885;-0.039200923176010239735 -1.4316294333352925783 -0.98431329333552253225 -0.010771636108322499639 -4.1301894786180586294 2.5127319736596072097 -0.84628584470824819164 1.1235117873681326195 0.38306664346130658849 -0.33612742991098581991 -2.0167792773249808569 -0.053668536229237437385 -0.77999675282607772164 -1.6944107715940006642 0.90272609733741748794 2.6370126663760706798 -1.0324169425307221371 -0.13776085569776699846 -0.54805592202080788411 -0.21226741330581178069 3.4083664508582107722 0.31239724671327961936 0.83322006371207124609 1.3624879081590755003 -0.66560672774810458741;-3.6062698194499542836 -3.0363973436484714696 -2.5455326824755708515 1.4792574055535891997 -0.12243663963282871987 5.6816614277194004146 0.81693585246533717559 -1.2130957003008742134 4.4393473086434411812 1.4199854261198263661 -3.8143706904834044913 0.44602222342758529594 -0.38669643375058365686 -4.4260573547901742586 -7.6781340349044544524 -2.4909213493667428097 -4.6046384342156843061 1.3988219569524729113 -0.0046695307109290542028 1.1646345711905345954 -3.4812432812228193413 3.8231939365801799191 -1.0403835645348777916 0.032055474136516409756 1.7511618360469032307;-0.043674337922021452818 -0.035136548137723308538 0.11989412340300582838 0.1438738477414547301 0.2264226713273652658 0.12320664863419815427 -0.048569054409136259987 -0.020235452148122660965 0.074746685462074538564 -0.04077567740461507384 -0.93826231319482145121 0.019604927142829513265 -0.1333707952845633038 0.1099613869338222083 -0.1376758002755890109 -0.61369316536600782008 0.024901752940129077973 0.027462063602275568058 -0.03163857620547458338 -0.097195513130672706037 -0.44669031397514474824 -0.083876078261884035037 -0.14658802023910694201 -0.059530438941216512472 0.091359138599285308335;-2.4829438551720617312 -3.7476594878308691072 -0.15310705862631698859 1.3262700432863216982 2.6155168577075178504 4.217316415871408708 -0.17173168973021307737 -2.0010241611023804609 -1.9415435167509846881 1.0856772552671958287 -4.8814123736764578609 -1.5502778292467034316 0.077711488170834286393 4.3002760435005020767 -0.64493061906667004113 -0.67787244005715896478 2.3389141176549053114 0.13393270349046920709 -1.1106844415001739801 -0.3723388710667733803 -0.88264839347666657066 -0.61791661076293324939 -0.38805523162967203721 -0.31912226899455914397 1.7882495513399989662;1.2376200249348838778 -2.6781219220179157681 -1.2141214094774852583 3.0067491415629632456 -6.6870255513024048355 5.2378185583723624319 -0.12690824900922240137 -1.0358438049285920712 1.6236528595613444192 -1.4674564244512100508 -5.1197004790924163231 -0.60571463545041848153 2.2178149189105811345 -1.4345264818317120525 -0.29838801405643922138 5.4270006857142165657 1.0238812227856357584 -3.1414731781522933218 -1.050995730160616537 -4.3307775924602420758 5.0231380287636477888 0.44215150327135233077 0.6341816401233728584 5.8053385038109670546 -3.4209523386526741184;0.57833219448596140388 -1.8502117353534088817 0.75625000787770368937 1.2198496065280135525 0.075086930571726162587 -2.6374136279094115842 -0.71131103639530979255 1.4520910748457074302 -1.0139660655952011936 0.44245014897316792268 -2.49823904422325338 0.098068193073576989671 0.23283073508572449395 -5.2839956565352146356 0.61971274697958689259 -0.29294005550038926255 0.13465313866879433102 0.18182859089912406514 0.16191934607153082859 0.28988376655452824382 1.7158698581311142561 0.020972613780272861272 1.4776782565360757626 5.3299289662902422648 0.46107074131777842085;3.1574374573443271252 -6.226607719891857684 -0.31981547519020769865 0.86484815716536456165 -0.48117805724775347676 -4.2216906146934771726 0.70379155746254140613 -3.1033918891508762705 -0.44456147893649583525 0.93471113207421807623 6.4958343327607490636 0.83083111262619402382 -0.26740154643129199208 0.7127465083952811975 -0.0025986836444696266059 -2.1298570455749907815 -0.67093332697905982975 0.025146011515626898181 -2.2251828813566487675 0.8188881795828294452 -2.6698142659503742991 -0.41748341961435631031 5.7543835081302141532 -5.5809660916061725189 2.0939948278457625186;-3.5650620177172491765 1.4441391331021975475 -0.01811009730518539293 -1.9902751859974876503 -1.7778020390526358963 1.9325912813738830298 1.9499762240032563199 -3.3344274609807946774 0.59565702487774263574 1.0248197077498581997 4.2570145345425123296 1.1803344738473240394 2.7352291813985276647 -1.4449344803548460625 3.7780478776010557596 -0.9135841402050808302 -0.52837144890252107921 -0.59433754597978505352 0.54025301661752023197 -1.5172681744683060856 2.3504262311987687006 -0.26858882241358328358 -0.892143119575980581 -0.95062298770934405834 -1.1241041323296780963;1.8725688777530449247 -1.451219119033134497 2.3089167638664735982 0.48726678987826171907 0.30207883322392836689 2.8571221127818478003 3.2418178529898740159 -1.2713803061239419279 -2.3571090438854787052 -1.5417337514102635865 5.7435676960891015597 -1.9597360185869530103 -3.1786547452268969494 4.6029380736965643806 -1.3760013168559390717 3.9116787627706282038 0.9372078617693089253 -1.7003885497816983996 -0.90784290216508922722 -0.30648243203060365758 1.2739641422521839065 -0.97263096263072879033 -0.88714661074243506 -2.7138598097887300931 0.93523788937320040393;-4.9878392015580708119 0.047786028364433678906 0.12090039141488029273 -1.9057023425728434329 6.8145306831573266848 4.3907676691486932086 -1.2389489649464915377 -1.7235255401570930278 4.1722798798271796272 -0.33405095088466180142 1.4802100947440526291 0.82082019608310652359 -3.5527016912466602427 1.3015035445537332937 6.0371700296732324631 4.7212629919403754641 -1.8497148543837409473 -2.3464818863861331977 -5.2752268745099062386 10.92473217998514734 -1.5481544477453645925 -0.99485903933370056773 -4.7500872930838626829 -3.8126451447990730159 2.545096315413888366;3.4061510326055217668 -8.8226682717433959624 -1.1554372561038521727 5.2212625842359194905 -0.09931681120825525122 -0.90062497850657519116 1.9499250306780875341 -0.7996649206578595459 -0.90076072872023282834 2.7776046114974524137 4.6923061788609823353 -1.1025592656779172884 -5.0406862983615869922 -1.6793848202637666756 -2.2128686063524769523 -1.4921951840182454685 0.69497608245416808703 1.0952125261907748044 -0.24950985985508353227 -3.738659630247983312 -4.2239048882069640811 -2.3759732756942955056 6.4526340120412593038 -1.6086955258852262407 1.498503754761588258;1.1596894777136956822 -3.8156850152425936251 2.5350604849914462058 -6.8275461398035837135 1.4403258475857894094 -2.4659892738646620458 4.6703712783520598961 6.4659297727579003734 -0.57878986500266571547 0.55266260208399531351 6.7181500240667642032 6.0673276070838975116 0.62513738784124806003 -2.5215960038501497031 -4.1235821682842335179 0.47961814894480242488 -3.5517986168885942178 0.46397451832077246658 2.6091229746605106854 -5.116701151726309682 2.8452659641561934833 5.8756798666577250145 3.0547016722665061295 2.6724139539687716827 -0.26945986222931561649];

% Layer 3
b3 = -1.6214752698269874198;
LW3_2 = [0.056282390037665205873 -0.0090247281767907427258 0.26747544941166939703 0.16144189610681694558 -0.025577171027753274613 2.8322049487292848546 -0.042252975973983149249 -0.22313510551843149243 0.13462026305620963806 -0.055344033894003254181 0.012038609451964778566 0.024116458125881368169 0.02487535547579281317 -0.011217439014245197959 -0.055829340404110958107];

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
