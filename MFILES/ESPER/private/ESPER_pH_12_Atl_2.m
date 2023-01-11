function [Y,Xf,Af] = ESPER_pH_12_Atl_2(X,~,~)
%ESPER_PH_12_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:24.
% 
% [Y] = ESPER_pH_12_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999025224415;-0.999999994516886;-40;0;11.29;-0.12];
x1_step1.gain = [1.00000048738803;1.0003169936794;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.3851462013201234491;-3.8357694970138136448;0.39372689724556847457;2.4538544541156359635;0.40166228708501472422;-0.67411025946747837878;2.3661790280338546033;0.98776880753777540622;0.74205756293546420999;-0.55885691443233487252;2.213588205628543637;1.3394006312402773062;0.8499909057078384178;1.1625443594446978057;-1.9000301998882991317;3.6759966966295856139;3.0655498051644913815;-0.53142840139772451469;-3.2964515102341875874;-4.5914013180263415492];
IW1_1 = [-0.36091822898725739455 0.57159139572617945202 0.69176318173123441113 -1.9346560509719057386 -6.6982052689974169724 -0.84530660421003134264;1.6719931849248166333 1.9384258360236736873 2.632588165465647112 0.12416916497362588678 -1.0652497044064801468 -2.2570703848133053526;0.26659150908346485043 -2.6970392864802787969 -0.47381225156404149779 -1.2661566461232867109 0.23152459248625570365 0.012870993752917038455;-0.88005628531169155515 1.1182180633366372025 -1.5907573061465156705 -0.3038110616863498592 -2.1518908094688122823 0.74288054703442396587;1.6402722674364267341 0.81864235926610917637 -2.5444096738146320114 -1.1881095995103396756 1.6584760913553251438 2.0921217231720663499;0.44299683667955419564 -0.02574336037685593867 -1.1392340076481664024 -1.1786469622777975275 -0.41327284863186014574 -0.50523059900607181216;-0.1746825815807822313 -0.42006690816486835782 -1.4181618411959933113 2.4828290460965796171 -1.3649379186469003944 0.21462689793579958764;-1.1788274756570238377 0.99823630536706586813 1.9849476247039865573 -2.1941790684686983681 -2.1459060517612855712 -1.0753151549961750266;-0.66704704019237437418 1.2359362293140665479 -1.5399584942260575726 0.5640936438825570054 -1.3020462677645230354 -0.0062664466620600814059;0.55891123445987267448 -0.39215118124233605545 -1.824937751875267189 -0.58354406162199257135 -0.30089731381902817819 0.75004573128841323104;1.096928358686110716 0.36086761137443323921 1.9427930793411503085 3.1631099956731718592 1.4804170505413931647 0.67703682142285004364;0.025892237489316766164 0.70831426818798781664 0.74159729192843704126 0.13301903466612105542 0.87137024529006124052 1.1931763867476541385;0.31723462402926094761 1.4852358344397866841 -0.9412882632333195998 0.048952019407809133078 2.1756624667835446907 1.4846709501515835683;2.0428419541482516486 -0.63433414017870815016 -0.28618599854801290983 -0.61615530363199455532 0.6545628504128698788 -1.5210187509614290224;-0.71283687727285471869 -1.0431362983719894011 1.3569955687392516452 0.50847029565607315682 1.0123815793072441505 0.9089246407471917033;2.8887470290865806533 -1.9100047598199330068 -2.5565300612638122146 -1.7781104690048954442 0.63357084104528371693 1.2604248717403880864;1.8477380658116822065 -1.487312781024866748 -1.1808305040521582807 3.1332364860930090877 0.74738064566738748873 2.2752564501168390798;-0.96187306076132250698 -0.9255377532118754802 0.4537940932887571388 0.57569203128497348221 -3.1691892619501422779 -1.6554905096451626445;-0.63548354309987808985 -0.14876870867170241541 2.7537978030917051875 -1.351336713269935208 -0.31086767752661581543 0.28183458889668910885;-0.8661512326390791161 -0.037998149713032898511 -1.2588402391917186129 -2.2432442095844393748 -1.0464240400887629967 -3.6327789976897157942];

% Layer 2
b2 = [1.7833233985883834016;-2.0219792901406030694;-0.37663946217053673182;1.8050147775188631893;-1.839909909182387171;-0.96589589294890954019;1.1592492277523462896;0.067521234308141026004;0.018195572374044927111;0.31118710493501150482;-0.66113424272638277035;0.90879825229829469535;-0.74224007165113659035;1.9733226090847009182;-1.2095635412030703026;1.9710961476243746837;-0.76351637758142076873;-1.7342628917695246393;-1.422356864044169722;1.5891916089626916264];
LW2_1 = [0.10926718516807829018 1.0013752824448611811 0.40606852335801163489 -0.12783767655507799477 -0.66355716114596874355 -0.96912602238613476668 1.6044337665044552388 1.3026201167173403483 0.3126656202612645985 1.7130519676948565699 0.46224999924744009094 0.14114520939560248891 0.83865361938381033369 0.61830860208533755529 -0.21960748130067658801 -0.36068170856618825271 -0.62432125708347874227 0.23225286497648575668 -0.73549848110390658729 -0.50411964139105180926;0.2640834019913112507 -0.97316636662408506631 -0.12841363737453581617 -0.59128695842012923833 0.07172437473182304013 -0.52984688609242958091 0.63483579828986580473 0.21517399030739456367 -0.3194415124670066386 0.2544645182231789482 -0.40827944742168170711 -1.0263806962654498101 0.29720619727890229544 -0.05149927259671967239 -0.93584038676698100279 0.65923048207761703043 0.60749693182241037981 -0.99100629482231983314 0.81864753314664362449 -0.36177786402238731878;1.5990828149591875817 -0.61794127033983059505 -0.23764362656515930494 0.92104347883806925168 -0.74435103490500076351 0.18182585376729065718 -0.20413918451472243265 -3.3608500975723716309 0.50127151298199967844 -0.80966690053460432797 -0.78236500517223994677 -1.3162282031144598715 0.18098404956454391956 1.0638145637526392751 -1.1653495883647966291 0.76673748271263131304 -0.82779027076829780896 0.5130714251761416822 1.0384481681060679925 0.75205920570456430418;-1.7137959127106114909 0.11770855884587597129 0.70898194943769798915 0.27751768647619406982 0.78295392796394569501 0.51379890515055270939 -0.44599277571822587563 0.90646563333758467529 -0.65016176638297429324 1.8787978008337820057 -0.03916998451848582169 0.1703859969353083148 1.2197118004531664059 0.51000463127816830777 -0.54441187084475994862 0.77994546683360999051 -0.050997093338073361257 -0.60437667244259474231 1.1334693516491975984 -0.011244348656203474734;0.93747507297900267709 1.3069560899867456172 -0.16121449068007487337 0.22204676542095394542 0.33620313888058295237 -0.11027360542755156037 -0.363657695868907338 1.1922966411786481977 -0.47754774454658899252 -1.9083804885295698739 -0.53577862556436528507 -1.332861814086644392 -0.57762295730717594378 -1.578729584613081105 0.62333855202413701591 -0.63063193041939624717 -0.20048607713089916516 0.15285139603475014769 -1.0760254681082817019 -0.25560683633538178183;0.63479436131475408533 0.42936775259132287097 0.68717171501717411175 -0.72100570148383902325 -0.88133626405431664264 -0.46199452554282782746 0.85807139488595574051 0.16975886596746539436 -0.0647509438941580584 0.68176873779833135902 0.047732671517547002604 0.23905051617401390507 0.57795265226362380062 -0.25133191715753377959 -0.24657242018130257066 0.37444108614676030466 0.50975479514084898192 -0.20917923597193158192 -0.4947652947488955788 1.1250987768578688542;-0.44856117036848802693 0.64867398330549985719 0.41870328123150230004 -0.33157841900418144609 -0.65165437053579355897 2.9717064881711992008 -1.4583118993229851057 -1.8311206270424666975 0.75915675642755442443 -1.3379808885693624809 0.43461936339339396485 -0.27858794804719233529 0.22283129982243110589 -0.24016087319191453675 0.49992193444455412799 -0.45332673132828948948 0.041915178112017420498 0.28019887993156739059 1.1818084436378513136 1.2237042387977814517;0.092537908718827624477 0.025537267457021246148 0.9602430568899499308 -0.81735319645001813793 -1.261307910000539545 -0.60645110617991782043 -0.38346467525380034669 0.18779278579070193822 0.1997291447355757088 0.75241058221913181825 -0.18803681763708435248 -0.53611508674851904299 1.9292695620013251556 -0.29155224485019298752 -0.62909018738004340854 -0.232826048701404692 1.1721359049158441845 0.37748851223404988886 0.16822992142421125394 -0.04365448076024002616;-0.47980617271856174089 -0.15458665098103827051 -0.66707249162037696166 0.40360854909322441131 0.79170862100071481837 0.81156551532442167129 -0.62486653819376380969 -0.33215096934742838197 0.39057830181403790926 -0.94944111332396352587 0.2979999482209633177 0.85559868272927008803 -0.72321102178935103844 0.79569123059851454283 0.87014985588944526729 -0.30736113615867671767 -0.3075136395690331037 0.51869831702960556008 -0.70198762571778594577 0.14686700049438472515;-3.3148646659098703893 -0.098128327754730765298 -0.24703116554440665564 1.4765113686099644585 1.3429297281590792945 -0.32413934918703873889 -1.700684177905276151 0.3604701318193395676 0.82286347372378843001 -1.5695071216263547509 0.81349794840328559786 -2.4793355609724381239 -1.0263894254759549085 0.79425970025560344023 -0.72627671698501283792 -0.20912411817052931062 0.018020764577422774405 -1.4883498092923987954 -0.83567470486063000212 0.47513142013499531924;-0.26642782292834044089 0.0181529492776048261 0.4278352410241864634 -0.16554767311343945857 -0.64546734680119677297 0.14901474289694441322 0.17397398805914460196 0.1629645715617660251 0.065009132821391840684 0.29543792643536387343 -0.036073052556370888144 0.34284107531515989331 0.92645262755661861309 -0.085637549586159766868 0.052302675046621546762 0.55692173174106684996 0.73025636803598026336 0.067378065911703088897 -0.25801707365067355626 -0.1844873067350838447;0.79172038164002056515 0.54423446339305536767 -0.49925517105939554163 1.5601213856987867779 0.13557312750251215849 -1.0217681032511569672 1.1154938774008660296 -2.2376268279142887963 -0.020837283220736770722 -0.4114215323157065507 -1.0699014043990178191 -1.3505086973163140662 -0.32209224857770757389 -0.30051663111164256259 2.4220963355819398188 -1.647678439084035773 -0.87793879638832361678 -0.22906838050903419179 0.45293720786259150657 -0.29164461204036851161;-1.1054398940383571848 -0.56040837676296562098 0.021652284643085505794 -0.87491721580867465313 0.17916158158804162959 1.1747974137058370125 0.74087677495713122333 -0.047727555503649814195 0.72673410241298685985 -0.1728379610898608254 0.4486997560945937269 0.90860940705535708783 -0.054403638679776093745 0.22815935867263520409 0.15584991083406260493 -0.62469713402554749582 -0.95288084801826533621 -0.26186117134180270094 -0.31091147108220351658 1.1239932160646173287;-0.034493841467025367076 -0.16000298486255762342 0.20618532376017828511 -0.94060823792111203367 -0.85482367256117974819 -0.40836923716601758461 -0.2198682901213062646 -1.185320113923150398 0.55638694365065488601 -0.4368171397030533587 -0.1326105240135704022 -0.17337861782756128215 0.28216840414784438718 0.76898466796523634592 -1.0600048686139447174 -1.0654298697460766743 0.5600319262379964691 0.71923534772238117974 0.34703325410684110475 1.3133325575077019476;0.83127841411470693167 -0.96994515601809661032 0.92632018017381101416 0.19165396168344403516 0.3019231229167302355 1.0001500767016002769 1.0050819249021012425 2.0283838667471756345 -0.213327834845683878 -0.35985703991801082369 -0.10239952092833333464 -0.51126010979531177902 0.29206964458566070375 0.10012027595665255741 0.6616373833510070801 -0.80809296097238803025 -1.1957269820023046503 -0.010561096008411686997 1.0445221619933340307 0.12982697418084743202;0.098383081286082627037 0.057275649100736049368 0.17815363014223056792 -0.42451764897522537012 -0.7489679351953031361 -0.32425225271602864385 0.32807061790838804027 -0.13695938539074367646 0.16477153947402919942 0.049079245241399137301 -0.12266436624126048072 -0.95965195605989117578 0.2903872886330038372 0.24466554149570501897 -0.60371972563191311778 -0.60304688231770309237 -0.5028540912197518864 0.42721421995777442149 -0.76157482884971716697 1.1527398448595689295;-1.4037218219733049729 0.10604894862573802095 -1.2932615373060836195 -0.85561848700293952064 0.61216193657130479888 0.95742180468852666753 -0.10180850392673222726 -0.47037126684897467843 0.89600777130884456945 0.56379521114151576189 0.40498327122340116446 0.23121053215441134121 -0.93051204576945290015 0.13410728813953939897 -0.0025918771830353692226 -0.30267422990154507767 2.5040562332073417195 1.050481497176093626 0.30428071681304119656 0.56616628627354170789;-1.0368663153856072512 -0.78444677445887378564 0.81539336901447145856 0.16508546197657800181 0.20330770712381038279 1.0986805823281642347 0.87545198394951440068 0.22514061413572880599 0.084923900505529628657 -0.51907323007861627673 -0.26885028410493139184 0.60166704659115266551 0.020080995247146477684 -0.31660218331488398924 -0.24627715938284133745 0.50248232102673218158 -1.3721699441273200737 0.040852775507435587499 0.53658106204892119795 -0.55880590896372306364;-0.12039325932353342474 0.30474354364311273535 0.28697151445625040944 0.0078004962160669027299 -0.92302835866117516339 -0.08651244081444564582 0.32934611970201110864 -0.45884389285112558987 0.84634256747969505152 -1.6263849456314618536 -0.30260711082097913049 -0.76571385671477232293 0.52209406198492669393 -0.16278095344489124097 0.20866310961101280919 0.89565030590883887296 1.7349615206256898681 -1.2907846213489060361 0.84364163072068820703 0.92155425335984453294;0.30769732813565148222 -0.68425295786773654427 -0.28685558488126089438 -0.43708113145310828473 -0.08477911836896616371 0.50889871771487826813 0.58219414548556647215 0.11615541793573039286 -0.095920251316898172411 0.51426439216361341522 -0.74095109646517409718 -0.38962640672027432087 0.23802575168438144804 1.2000074579549839804 0.045349451058910937584 0.27722860132391191623 -0.49520585594441829391 -0.21335727138569865469 0.70424846886385128464 0.36176996277659517753];

% Layer 3
b3 = 2.1647814139381176268;
LW3_2 = [1.0591161240871886928 1.8783238675306535015 -0.26758751536298774898 1.1872075484436568349 -2.1791335091329631091 2.1770811841750785653 1.0952357676470492187 1.5062485936333358527 3.1647061714599096227 -1.0235060608280499572 -3.2302227223185155935 1.7714620114280443808 -1.5013826722054965934 -0.99536400568563243585 -0.86161930617222615147 1.2626924259114880034 0.92306997668915768962 1.5411399546254269666 -0.83813019802164367977 -2.5363482397841572435];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2.03830178109291;
y1_step1.xoffset = 7.47181919175367;

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
