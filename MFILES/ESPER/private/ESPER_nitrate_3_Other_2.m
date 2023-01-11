function [Y,Xf,Af] = ESPER_nitrate_3_Other_2(X,~,~)
%ESPER_NITRATE_3_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:33.
% 
% [Y] = ESPER_nitrate_3_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-138.904784718693;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0042772329208737;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.3265092354726095714;-2.3161552604681818757;-1.8846638401297148846;2.389157258832082853;9.1588808100973295723;1.5720080387890329465;1.3545173481335204269;1.0545781139402550242;3.4700272387874622204;-0.044448009189516690665;-1.1993523623571398229;1.2512732089982494532;-1.1194756359996271833;-2.0350086452204965148;-3.4158554914553165283;1.2222335229654768263;1.7183238600827654974;-2.3549383531169865158;0.76160382679401383843;10.002417603286804848];
IW1_1 = [-0.057173504412800210672 -0.20827912978231508112 0.89784455362227333364 -0.84149338572890253207 -1.0467157715759867198 -0.65891582988912866803 -0.67001967740313639332 -2.0332320015387064238;1.2549769705060749292 0.087180088849849790567 -1.3931043860871719087 0.018611464056582378845 -4.7871214372588157815 2.4763659130554884769 -0.14820230399111097852 -0.99779626998884951572;0.38498978686952589623 1.0111986568416939658 1.3813667822828072929 -0.31938474043381442913 0.66382853710697398153 -1.1317171222353852311 -1.076473505926029528 0.20055708904799129932;0.15127858224629087136 0.023762700452306412258 0.98829744546915510028 0.62569560299760895727 -3.6918146244875060802 -0.092723325282549062809 -0.16311151562284725847 1.0959432357572886918;-2.0278713445689060535 1.5759556841229191981 4.8030923141886434635 2.472829567915165061 1.1499708589995938368 0.22390852640573208343 1.7986493023746308317 -0.42557707201660877061;0.1549646927699052934 -0.20009740153928115558 0.50487776070801226869 0.19468749171001736853 -4.2632910974922459957 -0.014911982323249259666 -0.24803547925544117136 0.64626876905946273322;-2.1098946122504802503 -0.8549952083393261848 5.3629427474917017449 -0.15780045570518008713 4.7394784239849947483 -4.0406901700943240741 -2.4333746439854944121 1.4999789428463801144;-0.1960325791943176299 -0.050865175169759947738 -0.94182266296805972061 -0.028556878169805545442 1.0962472888083631961 -0.097125601679991610538 0.20983121836453377806 -0.48622276697813027191;0.70822569504987065159 0.72978630858378124024 0.04590212244233188299 2.7601187880598980762 -0.71987775865903624872 -2.0263348190595378284 -0.32196574406607553209 -0.1089947319622828481;0.13459699871998148657 -0.88482097482803234723 -1.0649752809003170118 -0.76941701674551388201 0.95542086557600014451 1.9737158340235527731 -0.30352590513601029887 -0.88015147753322797719;0.23584838950535699786 0.20880410728186263314 -0.60007909794935165859 0.25122876391981441291 -0.86920781902494947957 0.74338993511760120825 1.6754741479502295309 -1.1540823213442681361;-0.085967521243144676135 -0.0093567904853623738021 0.45662615672272499623 0.050867029577650144101 -0.18083747293337545936 0.32943230660776245955 0.22600256084743536311 0.23857754449045695533;-0.52348821126873279486 -0.55262927393424732259 0.54656812167491553023 -1.1143591520381479487 0.97477514532763531729 0.027187974568471363901 -1.0778297930174696173 -0.19122761929842127526;0.13072900919626465055 -0.018730609262621517042 -0.8742087511263049171 -0.22542774069483595412 -5.0163684535000649944 -1.4632509221718943415 0.21244174833459528173 -0.091419968021515321532;0.62435714316938606405 0.29569289890067856907 -3.8046779060132220884 0.16663644265782309728 1.9562547309042062071 0.49937335492670853476 -1.6011183129109716194 -0.91199055239640525183;0.010065818650272196147 0.042894814803639935707 -1.8776277258163713846 -0.14772612890161046506 2.8763185947321363578 0.3137109220137589416 -0.13783132836939029908 0.90881231213420332615;0.068898904374461375499 0.027861979400537140172 0.091989770089399833242 0.35004178067058705492 4.5969245314888027565 1.4988923990067892777 0.11670158153954318658 -0.88094924282531295034;-0.13343562198339059544 -0.092116160311671746452 0.97960712068883815107 -0.25210287803351560942 -1.1779298520766341429 0.47813451461855255298 -0.17913201141738860356 -0.81735340633572073354;-0.091266414311782589897 -0.0073484565435930175592 0.30489547515115089649 0.032564555175389572894 -0.23355063450937996183 0.42893867850725198032 -0.22935265117223579967 0.437915305598200022;0.64800614559223856936 0.68810041718438597069 0.7704293330956858421 4.0306728646859371068 -0.40273724236731750503 -1.7681270372503106181 1.7302134976846814673 4.4230175019779744972];

% Layer 2
b2 = [-1.9724946709287067836;1.6042863098882114148;-10.302598081104811456;0.99817191486972378645;-2.7062135297935721034;-6.428296016685352221;-0.46735466465957198334;5.4231422647213918253;-2.9991809615367590602;-9.9726627863943750896;1.3897239632278537425;-1.9259386493517192118;0.4765660285273290997;2.4162582792726485792;-10.905505749953707806;-0.20339649026200362636;11.538878470076847194;-4.3433378569827292282;-10.410359289290184392;-0.78127757583517865037];
LW2_1 = [0.52690639562839325372 0.92088485112481743844 0.64874566658297339039 1.8489487387784919736 0.92260831988597358233 -1.8366211555489500018 -0.17206382869989739204 4.3893302752447116433 2.3459856214705356869 0.3052569176202734047 -0.75987161055158969347 5.7142345062288129043 0.27923125585464486642 2.1587622566013875414 1.0398662219300682796 -4.948111925278587897 1.338892442878547806 -1.4625032004173712519 1.4618035754135523252 -2.368118517955037472;0.092411523750975069214 0.42613957871594554394 -0.28894310592470356935 -0.31824093753379267824 0.59883372553416913675 -0.031885071677526119482 -0.15000247495420118704 1.0115674881871046242 0.077743763597123066722 0.16463827976812422937 -0.55553324067282738863 -0.28092346041873361706 -0.34948491548715804322 0.42900209609952272283 -0.30111398404182038169 -0.51736514957032941897 0.5605983171471617954 4.2247835161293325612 1.2579610201861257668 0.47436592442313374463;-1.2683013043312814094 -0.3029888167803386767 0.50094433505450963739 -3.3600255007123736029 -7.564809417371842315 2.2664332806045304203 0.0578618111111508307 -1.2418781545934554078 1.6648281326130685009 -0.98988626078570440825 1.5329436856336489114 20.726283379405519014 1.5042902686412855662 0.12922556254542591492 -2.9223004606775133674 0.41237543479283056413 -0.3418521714278064394 4.0083742586823731813 -4.3454817705962431873 2.4769169221072262488;0.19808982333751931626 0.3606086593362590742 0.34707104111496461796 -0.10314827029838251338 0.30957399486464370897 -0.23096093916065343898 0.049466894599600651217 1.2889600956741213977 -0.7190022823058114243 0.23464701324441342156 -0.63927413974426094967 0.44691179510703038513 -0.38216922120973112476 -0.64840646296043458197 -0.34458599496010783314 0.44606849729278796879 -0.21845677129716259413 3.1834905738841845491 -0.48976883230773449185 1.5528323402880817028;0.98378497568530476336 -1.5046317920428737835 2.1264839890096438602 4.8382220819446253302 0.83055827718421082473 -5.2762361486809545852 0.71712544143179290934 -1.0501179977462065107 -3.1201256497238114562 0.59426022666093114566 -1.5799797458980577858 10.997844722587869271 -2.7204462151043262708 -0.49829233537072431526 3.3721285151467812469 -0.13863176506536040167 -0.18896723383547869846 5.1166767319833530436 -10.009336998037003497 3.3777916806250853909;-0.702712427813438123 -0.32876874342867767886 0.7509476859419853545 1.218540175826481553 0.39764243464694132291 -1.2032422905967312676 0.033955923564319250019 6.3357502092384425652 -0.87348530841455718221 0.58545066548695112729 -1.1213570553175429279 6.5810627335229199275 -0.92364638141104338231 0.51757536487496058619 -0.18694653092159962959 -0.54071296081356878815 1.2024594518364741802 3.8913089715893653953 -4.8066524661863967438 -0.29957440134960344746;-0.93920732523646377476 -0.39918238984490594934 -0.60076040803723074202 1.7815412411865270137 0.010713208294218790062 -1.3160435263822172214 -0.029233095006276018507 0.20084677529213823544 -1.072746511016929416 -0.30385177588458889275 -2.5757755125313095412 1.9020744958807822567 0.18700338222290738277 -2.853798505548800879 -0.17820364399764401786 -0.075257106811145987879 -1.4107285083901623945 1.3339509744923971901 -10.906297382530626194 0.53071439357310579155;-0.50565135060633592534 0.16040743517619915348 -0.29657772701677964644 -1.3077499116108284305 1.0024446454441522381 0.51370012064411085806 -0.12582638719259428561 -2.2290321276251718885 0.54162629882550139548 -0.07538229092539941445 -0.34065606987088892232 -4.6220551868707913101 -0.020868191427438623903 1.6705120261261672265 0.27840989142868688422 -0.53246233237862716781 1.1386082336829248529 -1.562777625017726546 -0.0022987037346422246886 0.027342901838630931377;1.610401439801133705 -0.14187716164922792306 0.12330896572839185332 2.2084966801483698617 0.92222568542542182612 -0.8768525664486347404 -0.20448585634156324353 -0.36769574663972182105 1.7326840141214032709 1.9431016841902333248 -0.78791787721100658981 -1.0429879021350587642 -0.42966733655440608741 -1.7793977998469538981 0.9387165615940028518 1.3211406401866034876 -2.6280198444161846361 2.2310287938988455636 3.5875592464934888959 -0.77128236563518326463;-0.73744090643248838557 0.38838382552693090366 0.90217986392595139744 0.18093435817837696944 0.39701729171929639461 0.61170194557534363078 0.086092341556841922579 7.0038444982459999721 1.2838285744152855283 0.51677873271568075531 -0.90194774981141923043 6.3068065603596172863 -0.65177523872682330275 -0.25378834794614357762 0.24723845245595996167 -0.52267102029911993721 0.40185418818636425842 2.7778813914113187344 -6.2117186698448980664 0.81281578393391740178;1.4389151657387897476 0.67319939137931106288 -0.050750253285196798847 -3.4754815747474245846 0.88261499208175731024 1.5927560816164876289 -0.080317543409342564975 1.7007113592693028359 0.60152448035157801698 0.25098415740742180402 -0.97858896756092716984 -5.3394018847331485844 -0.87206533368804883821 -1.1271568718067455084 -2.528200542052024602 1.060416572078128894 -0.68969848186694404646 1.5014218328525228419 0.69854699748104331913 0.37790846745189976552;-0.2138623311034536234 -0.44856759620337888217 -0.58876836933632881177 -0.99790337689696162471 -0.27190155735408028104 1.0142620586840032093 -0.10486095144276574853 -1.0885837258577801823 1.2868556548287444219 -0.3182581005694400389 0.70501734981416486203 0.62026398696225126805 0.21126792499588553653 1.132594084135808421 0.37484883992981005152 -0.57385718286005227284 0.52880621777870639999 -2.9459096996894071907 0.99504578847526448104 -2.077741236426289273;-0.07521132427334759063 0.29025916971845205916 -0.56268274541916740716 4.9361819190007345526 1.1324907644928439154 -3.8113478718360642716 -0.099073576079664979477 -4.3330002942641279162 0.45964444964755862832 0.60605129781129996047 -0.95518827261661920502 1.3027972949821311488 0.35400617572613035655 0.38157719081667185179 0.045313615827980859019 0.33810309193400156502 0.75479286693690095245 0.57929520754431040874 -2.5176855157700273047 -0.21903532718615031927;-0.1745541428569603315 -0.18216393363955596385 0.3743391200460165491 1.0192585051708464849 -0.62020336497464423342 -0.54996046920263563784 0.089769895210907280081 1.4583509612056368976 -0.95273219308073420208 0.1890184204883197383 0.18578487680043537233 -4.3031556711108640556 -0.00077523904261126244143 -1.6289410419611218828 0.032119992225300146582 -0.37457476960201163685 -1.2025286734299041047 -0.40743432658843015126 -1.0440643777101763501 0.10547478985853478428;0.14114784194449758736 0.27765570464058403344 -0.37926701605616219837 0.14131724234951897889 -4.3938486859371792548 0.51316719183133119309 0.59776104335117818955 -3.2386551956032452182 -0.49662709904239471514 0.011208186683132776232 1.2661770554344404793 -4.0752319110720129913 0.56413403664026928741 -0.67004974527505445359 -16.17589631460033317 -1.0341632897310746309 -0.43999281960690067583 -5.5386287305663817193 -0.14147900690805723101 -1.3743715661499342762;0.12704012231703196156 0.27304070100447580538 0.014196017286552771106 -0.98331172241847875082 0.12756298831473183997 0.34787746925800588915 -0.032915276382122236887 -0.2021860559167599769 0.59388047923082587509 0.028908338242591410938 -0.53940088720190626503 6.8795195947755711785 -0.13118010249426584335 0.42504931836927362099 -0.010991154330839705833 -0.016476566147596388018 0.32023509620613693771 3.6519830197615483769 -1.8066456054553821087 0.72971790995959984194;17.720766141443164088 -0.075740769550542555444 -9.8222420070808666992 -15.951027699558521178 -2.1473934249713115818 8.6448963730710204345 2.2492353212424349707 -17.401017974897165175 0.85429419945397577774 -5.5885883541694996168 -6.0506283188662024486 5.0426711208409713549 -13.770261725058372093 6.0840698132820998367 -5.2578616280910281944 11.515624130424338745 6.5653292866076462531 10.577659232932516531 13.570441987378117688 2.5541733615348083397;1.7563858135865870747 -0.13929104516493490129 0.66400743140357254557 -1.5603525519793879361 1.3288305278075507321 1.0316194258503681436 0.22467946280335659992 8.7553743432159798488 1.0340287479823355721 1.3289531583788196656 0.10821387867892390267 -7.8819418040588873353 -0.8432973932277118001 -3.2363263265543071867 0.10151080529869346836 1.7713099331330128816 -5.0228129598968829228 4.5424684459834967498 10.763093975121607926 -0.24211877313363763808;0.15688337454533807236 0.35961510840509930409 -0.71506705127426084534 5.36504491828494956 -0.18701150439133329417 -4.2629742777999952708 -0.32393760617515765521 -5.947020697508289544 1.0750938775480092247 0.2545882151403182081 0.79146430047055060797 7.3050545873453573975 1.794532960240071251 3.4932363864174797996 -0.19136231040235823775 1.006098923297801262 1.5939933423138032964 -9.911135041336455842 2.4088599629584792616 -0.50072599387545313476;0.086827504501413646509 0.081244829011864166946 0.25702879512889992286 -3.8724955674408501061 0.32686132081167817987 3.0080675699151346869 0.064888928659070241323 2.268609819955729634 0.5664235101263557226 -0.3919454497142454974 0.13499508881087271539 1.0867739866492684975 -0.33870942965787548262 0.19998091924859437429 0.076785885483479032754 0.88749265230260565573 0.23848161208501900132 3.2198356204346434417 2.2093431022462639568 0.30952914639585349388];

% Layer 3
b3 = 1.0350926158353297435;
LW3_2 = [-0.3362538952412494897 -0.81120346536999810017 0.22773845201432632579 4.3526194287119110982 -0.22357870343964916793 -0.33102254246431006557 0.16539756695974761058 -0.94480889789136501822 0.29752703942270841697 0.42905809313335235844 -0.15966933532133542606 2.3427062489742485063 -0.67336031921817485202 -1.4941270393931824767 0.59847290481173587295 -1.6795810251068130636 0.53211887284484915561 -0.22676496191681971837 0.15896629303554107859 -1.1137185424566935499];

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