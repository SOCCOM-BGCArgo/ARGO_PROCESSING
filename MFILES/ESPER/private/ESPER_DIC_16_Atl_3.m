function [Y,Xf,Af] = ESPER_DIC_16_Atl_3(X,~,~)
%ESPER_DIC_16_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:15.
% 
% [Y] = ESPER_DIC_16_Atl_3(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552];
x1_step1.ymin = -1;

% Layer 1
b1 = [-4.614629891693062369;4.7340351421833926082;-3.7937304376191556265;-2.6347376134710955853;-1.9067906632944617851;-3.1056108355754266093;-1.7746386562740510406;-3.4950192959066033183;-3.2627436238545577396;2.8917148349686647713;1.3784727536728962338;0.69650861483514359129;-1.4867857495673861834;-1.0231031500718639204;1.8183270289848749712;-0.48439766057937810873;-5.1324023006393870716;0.46442933763548399639;6.0858412269305732778;-2.3281970389161501522;-34.929075410279651237;4.9842537184748367807;-4.9637790793556870739;-4.1967064060835337358;-0.80082274342530357369];
IW1_1 = [0.18420669765279773511 0.38946046797965327979 -1.4463267280080138733 -0.74116576341938944505 6.5318769128589373452;0.35904510174047071791 0.2919371492104168575 1.3292614621476539494 1.9562539212619420859 -6.5463594948018135611;-0.40133526627869725578 0.34375674256935906392 -0.48120580890672987051 -1.5986030610532453444 4.0955934865555008173;0.92022854452229563282 -2.0612511721563424594 -0.78807027216131098069 0.69193824339347553831 -0.81596532118200093553;-0.27042618432525517003 0.64936903117771671567 -0.15068867683956213566 -2.350831279989712197 -0.16898953604698363495;-0.24991922204537497798 0.034375765857054624286 -0.17238601228994918091 -1.3406951698749394186 3.4962656115369403764;0.081611495093214150853 -0.28133903971897900043 -1.4311262256999925935 -1.6283565297025448526 -1.8337294871049092393;0.19324798106093749039 -0.66680738093729341465 2.0827186001731678999 1.5591418344832022047 4.1013863851271805672;0.27116902140834886481 -0.34158649401539081802 1.9938869157422491885 2.8837608931974534343 8.1729177853020029687;0.30574863634273319324 -0.63579622653957179601 0.57875721127352564022 1.5514863455396095127 -2.5264276418132176616;-0.14665797902872101943 -0.29553178737192065118 -0.54390853251857207429 -0.97364683442429988958 -2.2344121630325504135;-0.041286449334363817254 0.58056800948586129163 0.091534292076131840155 1.4247260932746526052 0.69340734079350430097;-0.18403124029951223317 0.02910933434534677261 -0.071020464615889358773 -0.63055285518760439967 1.1989033097591859089;-0.24316588196671215516 0.38686824840844824536 -0.91542073069120533102 -2.7500618229402982173 -1.0980806674516425669;0.30339991631247748005 -0.75632795351132942052 0.0072453049817143280395 1.3328770894876269537 -0.40644284816479475664;0.28151388056428444662 -0.71315566434650057293 -1.107738986803876724 -0.24550579907118919021 -2.8552851855851275609;0.22140956512896245445 -0.25398269411452079058 1.4132336344886682333 0.075526753218883718732 7.542287202276799718;0.041705703351706632787 -0.35025461394271739657 0.45441910148532188574 2.7283274829958505947 2.4380832467367801186;-0.15962348972168477479 -0.40345359815465503539 -0.58529626637015441659 -0.37855291052255429785 -9.6579571814211284675;0.13121816533903768454 0.21427782993597221939 -1.3696295270102463171 -1.3897801132615865694 2.8765127019463534452;-0.4843535005598595955 -0.032685355005495887737 -1.2502075575608244673 -34.525867720032543673 -0.86637657797390010828;-0.11040301803497279409 -0.37567424698941093153 -0.48021038419678147013 -0.34641847107276907414 -7.7848038916738220294;1.8994755851665712232 -3.1325647945130810079 -3.0759294258594391813 -0.525188897668363186 -0.90389015361223090661;-1.1067221503531214921 -1.7782698211551606349 -0.55265034888521191725 -4.0240526817724751041 -0.19201728204459739269;-0.98975683598777774996 1.0277192389696907515 -3.7701380080415733964 -0.39051654146629477538 4.3272770577266257419];

% Layer 2
b2 = [-6.6565543390846890759;2.7575822792624924418;3.8888803887186953467;-20.450974397735091515;-13.47560895761091615;-18.820611075494781517;-0.15555967523895575244;21.25368131306814945;1.73814852652700913;-6.2432022683705241661;9.9416828631936464689;-21.56606571838691977;-0.23400587901079600295;46.271413832727176896;-13.11800508012008315];
LW2_1 = [6.0237998386187996402 -21.026718216053232879 -10.513814714663253014 -2.8312160304616562634 6.4279755723293607872 35.709716962344344893 -31.338322675258730499 8.3068622109980054802 3.8567552053558209479 -28.496697669606181336 -18.642331574826211948 -15.170157697905017002 13.954186939447616567 -12.798551880450270701 4.2278588865628412563 -28.475715353902554483 -21.308139228397632081 -10.138921407431068289 8.4573494865191420899 1.6677934058363601544 -8.3739509412304613534 0.24600382142270280084 17.398028718644209079 30.317193440148336947 -7.9137106412986559789;0.40543688707427943951 -3.7450979157878361292 -20.454989745042198024 -5.139830878586938212 -6.9695274107989799717 -3.0810008081720199336 7.417860559167576362 -4.2248754406686952478 -9.8537203782839437594 -7.537372887040358016 -1.6097757201102012026 -12.546430864715672371 8.2232729429936544108 -4.7058207800072660021 -4.9727581022931648391 1.0759587877957272894 -2.3309748422362703835 -6.5499882281839330034 -11.100724640639212382 -4.6115924168182713672 -12.354056062484534806 -3.9496667542978882715 -0.82674336132076242656 4.6685558822364843223 2.8865733896629852673;-0.34923898097004979491 -12.03489435448227951 2.5321342149172485492 -26.079903992419712466 20.062291582602568951 -5.6475440124670361897 5.9810719543299537548 -0.030084125044686920586 -5.928915242786377604 -4.3191877193827856729 3.1607176751829988426 -3.1687432510658002371 -7.2358303340185630148 -7.7632689135170842576 8.2230472910999452552 -12.836478556967987785 -1.4486606817643525869 -9.2309879613335201753 10.866052329315868619 -4.12112397006247555 -8.3120124334042664316 17.097671527156467874 7.4473280956899197136 -16.248221649175423664 -18.604117447172775712;2.6176465005865363622 -12.67042733654133535 -14.188418709874079937 22.836315432524944669 9.7878647460603342978 -9.0537408279253615007 -14.020024305343000037 24.822632460847405156 -3.5333592675978446884 -13.803445822240499652 -2.6140073778871708576 -1.1012071748619551848 -21.976032809867611206 22.54976412007685127 13.561913001610431451 20.655406895158009206 12.490504072076829445 -6.4995443786984132828 -2.6203728979549056888 1.8396810110373806868 -18.72662021133801602 3.7196542347667276296 -2.1922007151029223593 -12.540671231665555752 10.7927138457330134;3.9631709449853107508 21.372092453422190772 6.077294423956632663 -8.2455436710408562817 17.527094583499998492 -5.8123624926017480874 -0.12162303169948474046 16.797558807335665421 2.8457520422965147944 -18.638566460795200186 12.653035525294269448 -2.6602737691623952188 -8.662092340181210659 -6.9670939294725924285 4.7772326932978304015 -3.4131009687944473541 12.600061594622049554 -7.3212089828703152961 1.2387751909943498507 -9.6564357198944907879 5.1180676020717790209 6.6856906066814882905 11.99916768946527057 2.5952172288715735782 0.50812222302052412193;-4.4470484205048661508 1.0807601243500248867 -2.6358816790330847901 11.468124800726547718 -6.2007130388230571327 -9.1579768257995528558 -24.39136164975316845 -2.3858208433554652927 -9.2942953119317976984 13.228845039283381624 15.46432030362773169 -0.85693775188463927606 -8.3395019372585483808 16.815526896747226004 -3.3766131233827181291 15.614060191011663647 -4.9175513521554181651 5.3703636103449374772 3.4082221979303275461 7.6006191003911771631 -9.5936140917143681861 -36.517235438226812505 -2.5748255196758886854 7.1796310759993904327 -4.3447957087065693926;0.35000653021436917145 1.6823823384978360274 5.559751382759403171 6.7977881753248077246 3.5134894643977170325 -4.128914650501807948 3.6706444032729121929 -0.96895201315238654871 -2.8508006131248708748 23.427895742982890681 4.4384836386399735275 8.1202252505368655022 20.085893139471426849 8.3892244574564536208 2.2285345424472864373 -6.3492298701745557565 0.91509279203125726276 0.90620910400316501132 -4.4177450316025659305 2.7325673030969341859 2.0996326949906216619 0.76990407335239585684 -1.6232301619694289307 -1.8696908165135546209 -0.30809195050087800904;-0.13659091902758502268 -0.31862683554796433771 12.197528442273231519 -2.1781814961809238262 -10.8393908231653473 -10.871399906947434744 7.4788173013673526768 1.3932072804699575652 4.8428557670886416275 -0.57637978513028964578 -0.32827356226807502537 -6.8931766083983543325 -4.497895357942247152 -5.7926831924739428104 -10.976663117634572586 -12.363785750833253374 -2.3530449241064368415 -6.1784854869775074349 6.6818182835907320438 0.70063690314244264634 21.866062771038052404 -5.8075822592457342353 1.6991227590720314744 -0.1015598498761822871 -0.31345641660677936446;-9.0714705339649217564 3.0038603166999209293 0.28648458494792622409 5.8299857561649988114 2.0775266407809613867 2.0892933931431127981 9.5538701509226662267 -0.86938437192147666188 -1.2568209471289812917 2.8839740863410532334 15.937186943716051601 3.8975049613798002923 -7.3329072978937865912 5.2321706722756777808 2.7967234489714591916 -10.210897841326685764 -0.75230300575089936732 7.9612796813566699683 10.35120349239874038 8.05867231315521515 8.2638233810014600778 -21.107931841868943224 -0.85828760944885751361 -0.28115774369451934644 -1.032311026228748263;-7.0897669662370734756 7.0792196516728287392 -11.528534833768132373 -0.95096319998636824611 0.96369781104217155931 8.1649305482098384346 -10.117820217929068605 1.884841485347821477 5.6817655830437265863 -29.675688612674473887 5.2823092939512532595 -6.3453628398524388743 -8.3891063294112040438 -0.48552124821174996727 13.813599212608062672 4.8096376464372552206 -4.8519840351925402899 9.6436402489958439332 12.366249816389986549 9.8672409582082227786 0.5191755603449031975 -14.304166843783697516 1.4867730625640767084 0.26985136789630786991 0.13811460624273674824;-7.8958792719955530615 7.4429184464996369641 20.573222305799102827 8.5755135849023869099 9.4485520357713213002 -16.482925286634188211 12.610639885197594268 0.495089045305201092 4.749654575134941048 13.391775445301135505 9.7000774736963677469 -9.4265164169605597522 -7.6269580710997502138 3.0869778457665146831 6.0690658231698533243 18.877473775039948833 -0.55937124553930028803 15.240139083667202158 8.1835054037204475463 7.6558799714437393291 15.807744895744376024 -12.118564805691676156 -21.251495171801540351 -1.806401828509403007 0.25142516383969615301;2.198989489889915383 -3.4100520465126473724 -34.001589202909585197 -9.6824637429385465026 23.917298834609216129 -12.376941680847679095 5.405561380168863117 10.854729676092171431 27.256131909161990023 10.24656943491903327 3.0842393240938301702 13.050121488552569815 -11.360703417475958688 2.4654545319768037537 0.76255846032339891227 6.4749351332014812854 -0.49452424877356715438 6.7171168184932597001 13.325750806966395245 -9.8581313663076013398 1.3115617615853421007 -1.2717949597469497824 19.17317701452256884 27.440110879666974597 -23.752153868708258955;-0.45239607207960580437 -0.83443888468987537088 11.968330496971040233 1.5730541515826343968 -6.6947885151818011806 2.5891146780842104747 -16.608831432815442497 -1.6127835119088547167 -1.6257517723561709655 6.3972430843318797145 -7.3075404699842323453 -4.4398895837788776575 -6.2094134503183520124 6.1979917012142680122 -0.3713923331104763137 9.3749306771287468365 -5.5138683998131270414 5.5104614911759943752 -9.1421442657150926436 -3.6956239810517312883 9.1061265650010358996 9.0104632891280473928 -1.5231011082255059286 -1.7725845846982901666 -2.0835504644511395611;-3.659526127769769932 10.408222283727790014 3.5060928502965906794 -1.6370823447885458624 47.837837361528229962 -37.3437774273890426 11.805985102005422505 9.4854390283611866863 -12.953634321723054512 -20.945341375891288749 20.961573045623072886 0.93898287812322356327 28.159838208402870663 -19.973536315753733561 9.2911671041836623175 17.655570135818106081 -1.9243823448995924785 15.405793711040887572 -13.035531468054958992 11.191985948663601746 -12.097986937818333786 1.3412212945555508359 -3.176157793264303919 9.3635369724354795551 1.0412382002441078477;0.16357510904109298644 1.0204436033903727665 -10.614370246630672412 18.277704873544777797 -1.5021090290774041787 -7.3982149759229827168 19.818608606677102557 6.0494979061244844587 -2.5942168194192309905 -11.360457627563883776 5.8232209010391446924 0.040622764062252424144 -2.1107938864358994913 11.730719607346925315 11.599163796773012791 -22.73100170127780828 13.54755747067707361 10.41888890179403937 12.400946430789463903 12.355958443758551013 -6.0099791843170677907 -6.7938700603210264717 -3.9026137188904224118 -1.6982786044469095899 1.9784430417662652868];

% Layer 3
b3 = -0.060016218739372878543;
LW3_2 = [0.00031930704244243938517 -0.019202200435670703504 0.0031253389900237839479 0.017580720394178059318 -0.71165649462950253401 -0.013024094466380326668 -0.067450181514558865681 -0.08470188913568472977 -0.073975278993313187126 0.047116863592763083557 -0.078272768262338393574 -0.0036042879739664728661 0.042798469401135047496 -0.0097577244418074766258 0.028900852197925402942];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00115315130997576;
y1_step1.xoffset = 678.286531932822;

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
