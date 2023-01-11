function [Y,Xf,Af] = ESPER_nitrate_5_Other_2(X,~,~)
%ESPER_NITRATE_5_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:33.
% 
% [Y] = ESPER_nitrate_5_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-0.04979;-138.904784718693];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.568038194888224;0.0042772329208737];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.6935251428300794885;-1.1250589470672307524;-2.0142217818375249827;-0.99765736050813857982;0.87770271089154416888;2.1314699580203666507;-2.3581536792165964656;-0.39879102600274712875;-0.234701809402590722;-0.36392722852697700908;0.60327411727026003518;-0.21057811006736001924;2.2772424484376081466;0.97904693978643320129;1.7706659276615248633;-2.1718320163104412046;1.2857369119446735439;-1.5312609164797317263;-0.08465643439780715096;-2.0255582476165656658];
IW1_1 = [0.19603711606351234553 0.17902360066595590649 -0.70783451556427146034 -0.2139286052775181135 2.1466450492038942066 -2.5147215797079929267 -1.2705622580184428028 -0.026961882372448859463;0.20474239462842622195 0.47091213922773600409 -0.82371770825714640818 -0.12798143859359212282 0.7604087114961923799 -1.3918002492574059747 -0.21423615996125244809 0.12654767175143202818;0.62411708787507047536 0.58210785262907727144 1.5163864993013693727 -0.13037951971231981352 -0.96515260665309676735 -1.1072742464491287517 -1.0913184309779746517 0.0636908051160247346;-0.28116577312675422595 -1.1224096724264047698 -1.0791659321028881369 0.039759415070173553808 -0.14898420025113118248 -0.24092603282492205952 0.11726909624535990229 0.16856869364426402269;-0.30172512104480830875 0.2649537646131332469 1.1802544709632849518 -0.17649141030988546586 2.2684421726567189914 0.25992368545297722671 -1.158205672305748557 0.14561057596438586881;0.076035217142469857521 0.21099935930242655346 -0.0086232994562085475698 0.089836020450960271666 0.49911307101357488269 -0.67371974690294822441 1.6346016216091405759 0.4675795237143779981;0.057389687663565086551 -0.18042262182387694414 -0.14865330662991249611 -0.29199153636479202589 1.8036298105035499795 0.82587741066883180174 1.3966826833363850113 0.66122520827176145541;0.37623345094407711908 0.00063320694356633182038 0.027551180441727985554 0.36483401687253369516 -0.44587876218532734862 2.4150853402003007453 0.5101640656380782568 0.40204012499148317206;0.014402794834532384183 0.171967186249704862 -0.39556967373412660338 -0.031363160721497577998 1.3729066438194363897 -0.87543118353159776746 0.29754622250213313972 0.58802310149873115996;0.2292814159925403994 3.7212655733843074497e-05 0.25166133891379732779 0.13789571328824057339 -0.59252632448652275965 2.2796730957630102488 0.53858329669919569049 0.80383307004650472738;0.60304411484837838753 -0.61091865075297047127 -0.61620154588449149724 0.056564068712053439114 -0.76447268535359735875 1.1681099990187036042 -1.6963396061872761944 1.0463894038396492903;-0.16790415157531013435 -0.11025279238525330372 0.91078468278457969109 -0.089130245313763900583 -1.1770691112737492467 -0.22086594444111576352 -1.2558967882798817151 0.93828243721621340523;0.37376298338816532629 0.79534595868193513901 -1.318011119631697925 1.0068135637486153033 0.58533195430764717315 -0.27823673402987547654 -0.023780578434968682733 -1.2893727828172136807;-0.060950915320524795549 -0.22646381292497769211 -0.77375628855732936362 0.16213666086039824066 -0.24617537577853823572 1.7347732320569728692 1.3415732066583523352 -0.90927080457399900304;0.62988709133356379777 0.69220906301824669526 -1.1280137589950396393 0.27784566920996878192 0.24477250735284170058 -0.17603519468446163221 -0.57020074901450579219 -0.20372057237456606127;-0.039490360586357374506 -0.25975434135471442421 0.83602901344757585278 -0.42252926779161242843 -1.2977237300609465454 -1.1505100557170659936 -0.41649333384170228456 -0.94631612256283070828;0.71611910071000095535 0.26804617768463012695 -0.5630652863657951146 0.0091601913132232787307 -0.28203803861340426273 -0.47251389844684654751 1.5983222392133187295 -2.2861076658430969921;0.16196825020911523318 -0.081571422370904825283 -1.8428264046684761368 -0.79776378535400560565 -0.31417842104780818735 0.35358722810683501647 -1.8653400317249237972 -0.72435176879463492039;0.57223658472212013404 0.0076256090269128948625 1.8578940031231809815 -0.058417941843605769581 -2.3281954910954190119 -2.0171013083899422824 -0.86638292623478774246 0.41391743122794871779;-0.15533851483701380247 -0.77328176591042585564 -0.2801536018287191232 -0.044481889485160706899 -2.2842056459859891859 -1.8048738567525979271 -0.22712327071151797164 0.76607199871203479624];

% Layer 2
b2 = [2.4091959816741659317;-5.6422176108777417625;-0.49323177960549008381;-5.6231822656038481512;-0.73262047511094541274;2.9732846805355315212;-4.2227584220159322115;-1.0710537374725344062;0.20106319060990257563;1.6009124257417732817;1.2973936393387819344;-0.78522645135543978068;-6.4873868975678643523;2.1168235163905353424;1.1206210299512775475;0.17712957989873118203;-1.4044906653343129932;-0.44224225450069526433;-1.151146134205634608;-0.01758887708280624726];
LW2_1 = [0.14809371584701525038 0.70115006089345599527 0.36974923649118912605 0.8866295396256463901 0.092616605325337372778 0.43392591373532496402 0.88490995605559386039 1.4627750932509477888 0.095056397797346051104 -0.92401612083373707573 -0.1171983072479413357 0.94358853906737416573 0.62104003516375683169 0.20708248967270953544 -2.1302439784970319536 -0.68633634090604123124 -0.28685685301466601338 0.46920716429608211318 0.52911053053403855895 -0.18896416930625684061;4.31530006731143434 -4.6816114113962479237 -0.74275141424801527368 -0.15502289032806823577 -2.3122209991223336978 2.103320140544622685 -2.4276285256502228016 1.2815106891769254549 -0.14821270295162361497 -2.3356580155989439262 -1.3442270864952823306 1.1806631159798066655 -4.956222575118275131 1.8656723770167171761 4.193699152247474693 0.80299412406290293642 0.72871900846389259243 1.6563637042358845264 -1.2700654310870145292 -0.87266912650427586051;-0.41010027855808928177 -0.95827619638341265951 0.96693098713444358694 0.17385525798711767842 -1.1152349421414358144 -0.7217402838356125061 0.81244285862887810623 -1.2129924619217160586 -0.88279321323820159062 0.96049203343793987386 -0.28158933161855553839 -0.36444106822234978438 -1.1683202969830128737 -1.2653846286886982764 0.10804568447582020563 -1.2526496663761552686 0.72250385585128884536 0.33810365319852853094 -1.1931039603751358857 -0.13218252844967870452;4.5953666307831673166 -5.2354365095273900366 -0.076217065968642902507 -0.16832592208665508871 -2.3751973995357440472 2.0712465261100176939 -2.953750800115322761 1.574404596488663266 0.16584991137280319373 -2.6774973028596993352 -1.4449351799754126535 1.1565147218982745336 -5.5929602754603262582 1.8671199889736751487 5.1540147407233591892 1.0028617756417943507 0.74680703611429333222 1.7455115135573413543 -1.3212395847600344467 -0.88495476662428829151;0.41262067367657562222 -1.1874415155237252328 -0.1587316307131579618 0.0073884804191644195409 1.2784521725897375344 0.64347091606523421081 1.2128413457606939119 0.15456765923559040443 1.2593235787774992307 0.23781219799249625479 0.86858449904444501311 -1.424579021177422522 0.45891308470760788785 0.4453956538878539595 -0.57161225335485787991 -0.38197029761645695523 -0.10940726824171140286 -0.68603134938592313929 -0.081226368019910538676 -0.9580765506466959458;-0.51793190682768353028 0.24626465205667993907 1.2378288190992980944 0.3055242362683038837 -0.11212635265234298376 0.45872816075466976704 2.3831692493746770012 0.57999691023475041796 0.51088864414934920255 -1.0100076333754304247 -0.14707726148451674697 0.028380332208249932358 0.75810294476945738751 0.66360398710038381864 -0.27879305291876288475 -0.0032306321295729988063 -0.40049786396773162656 0.37132354720119331581 0.41915482573000084621 0.36932282491053436413;0.55664154311829128208 1.1003299580776684241 -1.2456312703994205648 -0.091740916257446486481 -1.0706801108641150755 -0.85902218479269787377 -1.3835810504856409509 -3.2687631078394492334 -1.9136921374582920308 2.2422655214625240028 -1.1171002786803410078 2.1836919792495668169 -0.21062453905029415724 -0.33095956646586915939 -0.41371880607509586847 -2.3607941862007510991 1.1361762045406400556 0.69875615336793617249 0.069101563067468688151 -0.21444545680176038682;0.5575547168941089593 0.47766689789187516402 -0.30248379097509325852 -0.3289271105709977161 0.93008975127391502191 0.23930591350440622023 -1.8548868212197366923 0.085907284809889544741 0.83389476628701886263 -0.1027367006888494777 0.090428623262463037058 0.042797063071051660854 -0.057986333420711959341 0.98815786946846961847 0.74010759842064011149 0.40684791422124694238 -0.2996956989374581215 -0.031342268884802446105 0.78125930902971363334 0.60096461372124232447;0.48207324755449548759 -0.11137375830451133607 0.10363997408726519089 0.72014260423893827756 -0.9239431141077421028 -0.24918670874751563948 -2.1280257022105693032 0.062564858242168980285 0.17635787573477179624 -0.23764074035356835823 0.65110947481850167673 0.46529119083916359312 0.03110848371628837164 -1.2413009002076969001 -0.89983123894979621937 0.74015752223932707032 0.88211015671978276753 -0.86485018997289264675 -1.0450628992999959621 -0.85293264523530287757;-1.3367400286949637778 1.6448478118581604246 1.8496598229830853999 0.37786891346175638384 -0.2705129032862538252 1.3214349953869302112 1.6605337285881827114 0.38691893901826590607 -0.51106528049472577102 -1.1260107436655630853 0.0090380427702474869234 0.2064860122446362245 0.55049594785331645852 0.28928336577809554697 -2.5979368914322500395 -1.4306624913251038134 -0.6130857668178300024 0.32869424099138322415 0.71447422970933238684 0.29222798991217718845;-0.58552787120145055599 0.36840094993022171055 -0.29544120679192065193 -0.41411746346846917133 0.64446555840143926464 -0.75662987184783758732 -0.98093642592436269823 -0.4077147843152936213 -2.8374225475698997556 1.5048737257109152665 -0.70940576578300174049 -3.0993211605711765699 1.2140224100747771807 -1.9046330347481288392 -4.0418676894754339557 -1.7936781197890381279 -0.16156616617050306717 -0.14045207212241472838 -0.181657468805763056 0.51745676668091178385;1.311652326200504648 -0.51636199908909086176 -0.54578273278180877703 0.26206056098720653669 0.57546172815129781331 -0.72332791484076763844 -0.79431094167823101948 2.7009961471373835629 1.5257629341502898512 -2.6664337255510290881 0.58183506981231836175 -0.36487161564904430078 -1.6780410973936172869 0.76667573377674969048 3.4991594808530028793 0.5163617232299869686 -0.48098178120827989046 -0.82339783477819095125 -0.29385708264128862499 0.35919788877806035909;0.055753140320011565501 -0.074817053973812136514 -0.65565092423025495005 -0.35642644416039165334 0.56380835976744625437 -0.71370691524521212568 -5.9345731773795717956 -2.1902576447022932271 -0.39303394208853242864 1.7838540954054613596 0.20328443028001813575 -1.3047042875388998873 -0.99846108095434349394 -0.3009402585844737632 0.58681493904061332589 0.0069792168592526159995 -0.064311320408974151364 -0.08114703563311714829 0.058134828546778816172 0.18005890647538849714;0.30331817269347532973 -0.19171821645607015072 -1.4010223699195329949 -1.0484438004495739705 -2.9141056083307597113 1.519718566692807471 -1.6000384885397380685 5.2664546039047230863 -3.5893377117634153173 -0.66097568617762192389 -0.34275438768788457633 2.2386290195341986831 1.67582336241606078 -2.3087784081522872981 -4.8124389960825428858 -2.6080445491903372357 1.4391722403367408489 -0.072742293726474269167 -1.4486713637489891937 -1.0624503137595000535;-1.1998190968092383724 0.51526343844563071706 0.25548209585054965087 -0.53346110091205001247 -0.34533343740249289322 0.71823774787511540652 2.0751874334198352479 -1.8985841107654344295 -1.3505061100308122235 1.9707938165709244505 -0.6078299290421841361 0.35586536573986210508 0.99499512606095608103 -0.33168133074270889082 -1.9369247253553683041 -0.51722502472684184216 0.32940615919909332909 0.84681064198053546122 0.47021030932470858632 0.22947405889880614915;0.12522225883281640213 0.23379581771887450459 -1.412613967770447676 -0.5333191528100650336 -1.0334669570773424851 -0.67770026415658701779 -0.19384356560574533823 -2.1234157005858684819 -1.5315892571437701442 1.7987655588179367605 -0.70577926918971511938 0.54683656158377724754 0.79812846949760485327 -1.5913779300720072296 -2.4258027279386840114 -0.15920329298352764447 0.73115795564275165574 0.3478887960985423411 -0.50180122954154571246 0.039211404373204396256;-1.6443978425858569725 0.38847993245110751248 -0.1027593268828428108 -0.4723009137603257912 -1.253229923414109992 0.4126740058804385658 -3.4448238150832195537 -2.7429556862224866443 -1.4709186791003425743 3.1712264517545616371 -0.11671303883905412291 -0.41353271648275063965 0.86077043223580729236 -1.2594178102940039032 -1.3758350939843755345 0.52232022816242329455 0.42741082001186286865 0.089099523294026247711 -0.43094011850234986705 0.63975043793549324356;0.45320909140218212929 -0.56337876103088113666 -0.23004558173089839634 -0.76030439493172663479 0.52312061196306269828 -1.6402038049057408031 1.1366277222079292919 0.51015408241766835307 0.12124857532805850824 -1.1174029708646027181 0.41537499356977497511 0.5932061243867047251 -1.6331059324074952421 -0.57971519333719412348 2.7321437230363025961 0.31990221878047447879 0.66245512078683077561 -0.65881258228321348813 -0.4719744468731570386 0.078671381652828553666;0.94025003463183165486 -0.43738882761193631499 -1.311562224873386473 -0.15284125817590901719 0.21155012122927077667 -0.2946875627992673885 -0.61192026660618070721 -0.006723072033929534605 -0.22044013777023793721 0.43237302322960258216 0.17522784163192742235 0.33514977734353984529 -0.56482673588189458158 -0.41316775677307066017 0.13773962352671079157 0.24973105860855540628 0.29441330484795591449 -0.39554575990109774875 -0.44871424204381127332 -0.67693779324137048903;-1.7321210619770304184 -0.65163892049667948569 0.20235412727475560191 0.69625642670156551173 -0.54986184397074833186 -0.34109575831240968036 1.8596883795067000733 -0.84940372495746296888 0.15664281043047070807 0.54465926091900918049 0.13231836722119225125 0.14110698578258815128 0.62605465559897044425 -0.762588473909511122 -1.653412378022703022 -0.25926565196440587435 0.35007994942785153336 -0.20349959674731252224 -0.36947539805125600987 -1.3697094124455662101];

% Layer 3
b3 = 2.1830150723638515053;
LW3_2 = [-0.49292766885696376322 4.2622046046937525077 2.3359529804883560899 -4.0944895505283378512 -0.32412669024751394886 -3.914853500654057239 -0.34586197683422692517 2.2671366830021915106 -0.61211310941702679056 0.39437710662700936615 -0.29807723985260548005 -0.9707923945350767081 -1.0525585798810872529 -0.12562536735175178215 -1.4587017202628485002 0.5260542136834224225 -1.6116017774965549769 -0.53455807277390077203 -3.0269566629004462932 0.89020817646725269867];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0415843642790311;
y1_step1.xoffset = -0.9;

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