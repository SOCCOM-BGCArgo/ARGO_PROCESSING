function [Y,Xf,Af] = ESPER_pH_5_Other_3(X,~,~)
%ESPER_PH_5_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:22.
% 
% [Y] = ESPER_pH_5_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-0.9;-133.803853032615];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0588888685539285;0.0422119037568594;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.3251719755542481494;-3.2801678014062574107;1.4854396622594625565;1.7762908491904270125;1.3039707875337633869;-0.25299580491626444401;-0.62657554784531432901;0.44593020874754007998;-1.00444663211489571;-0.37849018322430688643;-0.25714610445643748537;0.52799663217580006247;-0.052195351114229786293;0.20256821143692546605;1.0351017519139227296;0.1050194579468192585;-0.3754639110339026975;-1.3220392942768155464;0.37312947566854426373;-0.40166451675862024917;-1.2662570150613110531;-0.54991440869107355294;2.1078687198370369416;-2.3692035275187204313;2.2649040467798515763];
IW1_1 = [0.088793804421749653044 0.70609051211243401003 -1.2952311965300340457 -0.1320710633997931871 -0.029051863051660981735 -1.1689528122772432361 -1.1829690641818213503 -0.83316327873607476384;-0.52436627093036669578 -1.3114160520196056492 0.53557765502271081903 -0.043499361438851899542 -1.7288783928939621681 -2.4886166400959912615 -1.4057586264812138577 -0.1026229939621524645;0.2609859001433553205 -0.016352145386928064463 0.9621415867795184429 -0.40996051320895621206 -1.019692533257874123 0.625151118322677668 -0.47703529909894348826 0.014748310015256118136;-0.58507688605175767194 0.46148688594099585014 -2.3643534217342745407 0.3860420713569199247 -0.14958329748061544029 -0.13753947526860291206 0.73667787118127614487 0.63905979689741332361;-0.82322290290572130989 -0.22329280710303736623 -0.56276539114938128883 0.62354629711852149754 -0.4479188924860715848 1.179993554825314872 -0.12423352468151099837 -1.4577978104211657406;1.0772330119959532269 -0.41811427228649244991 -0.48532236501251452276 -0.12177962752102873545 0.19923360944622095503 0.87362558312484606127 -0.051683289906207595388 -0.012165998734600887848;0.40897369634891128642 0.44041285431785243132 -0.48809920449380667096 0.16824360290260023665 0.37089739342418398538 -0.044824028336950275109 0.55622912758517140119 0.36779024122967979293;-0.84252504008886486275 0.19709931279534959048 0.7160051610902985475 0.077687171006111443994 -0.2699389077670036885 -0.57833981127881695095 0.12191466940334896862 0.19533047839294864767;0.1757784209618498128 -0.6900179405976900604 -0.32228805782376795941 -0.27160319195838394934 0.3739558487848225754 -0.83095593739838580216 -0.63801579619901283014 -0.42003629464870989585;0.15126165964371185435 -0.23591137608600912245 1.1312486132620651613 0.2952934441686415612 0.97191275387138820196 0.7578557125770521985 0.25013270446311214323 -0.50375755993002357958;-0.40001584538685769354 -0.066262941681557732676 0.93147595448979803034 -0.74570573722484112533 -1.3677359723877744457 2.0082063939742611502 1.0008014334407819845 -0.60515657783610388165;-0.19896950982281594311 -0.45713365290464913482 -0.26906871385912600747 -0.28095829220329510134 -0.34614092573217936089 0.45247127569845058126 0.009494811893251860202 -0.51591031526738517776;0.342643823283058091 -0.14579109513553595301 -0.49931571929890933292 0.1448249438213231044 0.82530235194196721604 -0.42263674570494441651 0.22532423836730069988 0.33216804368713426232;0.060659177840929383307 -0.27437477317908021135 -1.6342712458164301292 -0.067698258552855591907 0.78918519183428192942 0.66483659076814316524 -0.23001940196686818108 0.4729060446068870105;-0.80380565689909566629 -0.072701738972641469227 -1.276008552768903126 -0.15463493097487607253 -1.8064033291642169754 1.281023914411651754 -0.88714572095200650015 0.090949921731471844444;0.13494573968069201175 0.34309095801241223933 1.1923932380546808929 -0.3772767313646405718 -0.2495570550983160496 -0.13252775697485313922 -0.5998794522945231078 0.45927510152007128275;-0.24041233649490076352 -0.35728506354544453671 -0.30635793137849576251 0.27588589467506141606 -0.03096347482220227279 -0.50154626771718857636 0.20442547375435615731 -0.20798063086991677717;-0.4341786345752334797 1.1791904213108037602 0.15924976681742003537 -0.62026836527270678978 -0.7605695314423238651 -0.4752916281133784393 0.1843004239724566895 0.31596370728003342343;1.0062954905567753627 -0.2893623554978859147 0.15617743432139727378 -0.070674225800074944037 -1.7153323013400563202 -0.14558178408525401792 -0.25971180071098343145 -0.70041364989086674075;0.21240493974719623638 0.15814388536994450352 0.75040889314140946098 0.26421171181375113779 -0.18361882668124573925 -0.0067030530644124446116 -0.24389025766100019821 0.65762871303431769121;-1.3633565078017906291 -0.13208548813158815793 -1.3452282550394254468 0.36645138475818622714 -0.67161378517411562949 -1.3873494369644832069 -0.2246836278393714581 0.09572593947795449143;0.24157345977987607077 0.12408787295070650147 1.1644867272271484104 -0.43810935687955265827 -0.68842047433033914938 1.5478267656026940546 0.94546935188201275402 -0.066709916065451929645;1.2585779359668953781 1.8045959906805923989 0.29544265071827197078 -0.18052909824675908435 -0.12687660684681253964 -0.23241346489403005848 -0.08584283471178802416 -0.84181457249048574631;-0.82652393735698037602 0.057014269294915362329 -0.80923150664463738835 -0.98521744208067507653 -0.84864432296005509748 1.3145875789339547257 -0.1218496264108544086 -1.3067674083279454678;0.755376039932765031 1.1441802804862721832 0.40752916012596479556 0.15362941463540230602 0.13820861653854427109 -0.26441863328414327849 -0.090092363260960039995 0.67574181671173016994];

% Layer 2
b2 = [1.3923791579902036197;-1.3879495829840076127;-4.597645417492921105;-2.7165908956410960151;-2.3553538850406208205;1.7911950657945634457;0.18111518263587858546;-0.35171281940413356981;0.24363331159436274476;-0.098207475154032691789;-3.3332329657076269491;-2.8424139402723489489;1.1096798873658375584;-1.3614930142183019335;2.9371177687785388066];
LW2_1 = [0.57547485151637878253 0.3217992075425781251 1.1729189053785871177 -0.80905301996199807668 0.48281696402581469885 -0.88946092782453445302 -1.2513655726545789282 -1.9345154326644420362 -0.11099808300737475097 -1.0465325472824615893 0.083288136128160675242 -0.63837045429288086229 1.7630732241680324712 -1.0002789283047401891 -0.11359439206015142787 -1.3274671034573874984 -1.3430060329135702357 -0.19392906843913826354 -0.82550986855829133759 0.39907581684858323134 0.10171368389692385381 0.40019085805813925383 -0.003445878468625493754 0.18641527914021099877 -0.67707596791365276534;-2.1442686357874403136 -1.0105511068963639598 -0.53473632788970226759 -0.91046547114215292673 -1.0965247668331030706 0.61200889147054415851 -0.43850388876031692886 0.26768654482682435747 -0.1151903179064676791 0.57498600169872826626 0.27182164264221220229 -1.5969951989666180125 -1.4930705646131188225 1.3796030343065490698 0.93211095853779635245 -0.32835128616664971402 0.26303222836834672105 0.063906172702124508511 0.76239277516172132554 -0.52290665464410079899 -0.16086792316618039878 -0.9719877373821463884 -0.89942544433595239539 0.54542858590231457772 -0.28626638293268524826;2.2068845019789931072 0.3486003994782087001 1.1276759738605328653 -0.24256005663099239333 -0.18887199510827079596 0.93737340479508057811 -0.99965022645003842605 2.8226840833910213924 -0.93305156401744770989 -0.684710932817544915 -0.18473847648712715674 -0.97216167548197562898 1.5488825257049891881 -0.70316241681695168531 -0.12812853559197193154 -1.4417563757910780708 -3.4788461598841875499 -1.5250530293138555571 1.09472746421949374 -2.8303185603583029817 1.3513855531440535085 0.15070174460345753253 0.46574393709874206682 0.35489840539056655144 1.079486826583711423;-2.1676927286843530318 0.22440824846732751463 -0.021699423062522149586 -1.7854955226184789918 -0.12324223872333825436 0.46571891999139070073 -0.92863344004746462002 -1.3521245094456062485 0.26523602896259951267 -0.98960394631334747118 1.2188995868723835336 0.2151569984857916451 1.7104853260819758987 0.17854338319914497824 0.3294364212219688115 0.44669732534225436993 2.1287306811420103614 0.056321555139435562631 -0.084252343976279620019 1.4476152266103206045 0.026945978648526248966 -1.1560217489550548287 -0.14881050416687227611 -0.54751116202696648028 0.59815489113309827651;-0.37818059258383407428 1.0424319039775356011 -0.061486989159808498606 0.6551182141309495055 0.35363815050411329111 1.0646412846760289206 -0.32318270636860169098 1.7594442599981499065 1.3209584198872155181 -0.47053579186761024244 1.549803980977227047 -2.5962971740095777129 1.6964033242021849812 -0.29984497326567310926 0.26982476858639514727 0.35744659624382268959 -0.36910307096700012419 -0.69047063645311312907 -0.56530517200778718223 -1.1955052424411845013 0.15459974210464130784 -0.49050460491821101794 -0.030293289282819602093 0.22440990554848513905 1.089186899427731392;-0.90830356656752775368 -0.11932288286907966057 0.060304416813076439086 -0.15895401678412193758 0.36528223993646236067 -0.4856437606840030563 -0.48652520209204391799 -1.2870785353523761518 -1.1882539605316919662 -1.0800510263611997619 0.058352007820292484119 1.6366890691354971654 0.71314600703282993166 0.8867441902228474504 -0.3620667168663763702 0.31311054789035985824 1.3823857574362687917 0.6644135173189208432 -0.705843706307641372 2.1234399842096882693 -0.47128510534202211479 -0.64304354108156103198 0.57603142520205552479 -0.16489410742370252949 -3.0966973679199614189;0.12343258876352382336 0.016604879733683759052 1.2490913320564420719 0.17527704146716605038 -0.61289194403554969615 0.187508140074810703 0.41307793401106102094 -0.45430424756757398219 1.2340188919694010661 0.84976373000093452514 0.10808732458036865443 -2.2367749364902760334 -1.7291238619960791478 -1.2979598097168287829 1.1614223156652674351 -1.2916898627186004234 -0.60146910775323392251 -1.6313553290888751501 0.36420341662877786515 -1.0841886511436500395 1.2920071135647863692 0.69949601092882041087 -0.74886555137731591891 0.4235014777102509842 1.3023589043098311446;-1.0749541585257600396 0.6175197871611035616 0.42740238511430389545 -1.0205549703854042143 0.28064892676927821746 -0.58489985438285885255 -1.250476975173645533 -0.19014111265212441104 1.2011461609475460044 -0.81459128972378092026 -0.063889682595363844753 -3.7484056197203541316 0.63221609292969571658 -0.46673201522991408785 0.46483429731713876709 -0.32271798304983861927 0.92618526817194557488 -0.71649354383420282844 -0.34489511649837062413 -2.6178172751171313593 -0.38067001758728058469 1.8904982459868038269 0.1314963941359060251 0.30445773610124393649 1.1863771670312552864;0.33855824408574941575 0.11943217228993002754 0.82777605775030949076 -0.38077294704790060065 -0.22410527969115343061 -0.7584423256311171313 -1.7142435515121177136 -2.2854660431900151174 -0.42371070248339959274 -0.30643047321386723825 -0.016988435963596882489 -0.43399721671883462193 0.69113287623376418178 -0.80773402981452613325 0.092998475125235971483 -1.9034960730989924915 -1.638737774709352335 0.0058599117099125287017 -0.74166611750357480659 0.48610470435920238863 0.75625851596689419587 -0.36794421223224654804 -0.014925445490393903294 0.12344157491607804511 0.0789622391970014742;0.24829690408744259367 0.2615003436515778712 1.3792380900240290842 0.043394677862163173743 0.91696808596595336116 0.65859217468338315804 0.33098533475333191678 -0.053302616826510239301 0.78106753674189621783 -1.1879867697686197214 0.20429702116820347424 -0.78046237745909297079 1.1242833588501213171 -1.1550310962813046256 0.16813277307866361654 -0.318959126369819157 -1.0734568982533907455 -0.3565983161302054838 -0.60191040913984217298 -0.036731526938350599065 0.10303331993876285733 -0.099933319901013159536 -0.66003826850106139812 0.1763013952599340739 0.40367629137779442727;0.46916289205188954625 -1.3043088560203288662 -1.8135216073984947727 -1.2022069451226540249 -0.14841967523076524116 1.0216148664956099523 -2.1024624305781052236 1.8775600077355174555 -0.48090743507006222801 -1.3865632916042212042 1.829946199636124371 -1.3882445489903567637 3.36120979701660616 -1.0771217945498303692 0.058235309375346537564 0.41553965544383036912 -1.1513848147446670289 -1.0527678942149165664 0.02713017630914135031 -1.5779432271905136087 -0.073113465890737680541 -0.92931416156989432586 0.12842964978005644294 -0.30660309368727051149 0.59263013430211297639;-0.62764545930724002787 -0.22062264183358765157 2.1026999006539712767 -1.1435302227559684951 -0.38787008300615977152 -0.045229330842947876312 -1.3560619694666820045 0.50991196808891781078 1.3567468098488739692 -1.5590171149303646558 -0.25240326954072883625 -4.4627328945277264793 -2.1879272615535256818 1.2104147668572140528 0.67802995055809534009 -0.43939642504152376379 0.99911029753123437036 -0.30990544481550308387 -0.24067804810379411773 -1.7215876081735042025 -0.2184553122909730638 0.55757059698036670614 -0.66721561129730555528 0.15552559739605903433 0.93142123160093415191;-0.83711180952230657137 -0.18977690764400498291 0.55810563956180858636 -0.057930795903403685398 0.37357115577881905288 -0.013993247205232423394 -0.32228466238299285918 -0.64691635021940863126 -0.99304248221808411046 -0.90811790391856961335 0.01801157090431679611 0.74895579622963026623 0.3141919762910810876 0.65507357739256633078 -0.2959165902246655433 -0.2298059817061129384 0.85696552263786429648 0.34722583495060749259 -0.54954794650367699838 1.1761561207603705714 -0.092394353476931909142 -0.50122684775656123346 0.16334410817852162978 -0.099387955384986806373 -2.1300728055614293055;-0.65159886418972823297 0.4341927898610052905 0.55969873212235232707 -0.19371222208647675966 0.10512024538235634208 0.20737884878711460024 0.31229719280563028594 0.74991380446708832963 0.83362270957938855442 -0.30058042268460183344 0.13212317065587017573 -1.6069375116890591837 -0.37222451209433604635 -0.075366217519454589691 0.49421808219503737813 -0.28847758845026272301 0.44161917040682119318 -0.28758085621758883166 -0.062200377920910050744 -0.21765173024770570476 0.34214596410589925846 0.63198805463080875988 -0.21352218222806462755 0.17126618649585426124 1.6816822084941482895;0.45016436085285493895 0.28380608542762986168 -2.0355168537149448404 0.95244411853933164114 0.4477307306925260022 -0.0061020379324032104559 1.4276247551199572339 -0.67585286283068535251 -1.2994294925831948451 1.3922768470548410047 0.1773364501540536653 4.3883617931561129666 1.9563236826210570474 -0.97479785637045579527 -0.58735328447757617631 0.60263278873101699951 -0.3137215704613002254 0.45300561938108058202 0.18002328306648421319 1.9074991143800064108 0.075062672394603929815 -0.37685450861083008567 0.6009975626596685272 -0.012729136320950174693 -0.56336451645412277678];

% Layer 3
b3 = 0.52362073888662397092;
LW3_2 = [-1.7771422867085999719 -0.26887939219589984674 0.60316230749510157949 -0.51458138401119912775 1.1523219347522450917 4.382899944950930049 0.36887184635819897371 0.7993440165196546987 1.030844026637275368 0.90481174123083241145 -0.2972636296291970015 2.9546500606533467881 -5.6435345111157033671 -2.2016887912257279147 3.3160956052647780901];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 1.90239167587468;
y1_step1.xoffset = 7.36300175998845;

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
