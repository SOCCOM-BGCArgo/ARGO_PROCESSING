function [Y,Xf,Af] = ESPER_silicate_6_Other_4(X,~,~)
%ESPER_SILICATE_6_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:40.
% 
% [Y] = ESPER_silicate_6_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-0.04979];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.568038194888224];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.96380810584696530174;-0.19853291303918677024;11.121881859634941137;0.046214286331925154283;-0.54145353421287922835;-6.7117943976300962206;1.4613619939491293565;-0.12994083064519679205;3.102932770459637446;3.1254211763489179177;-3.3138704492303756588;1.1040099525581876794;-1.0576996099191673384;2.4154262244343862598;-0.0092274876644413035554;-7.4381778378271627972;1.6857705071991919699;2.0114507965600045836;6.8748065870616175133;-0.93291599292236659924;1.294760973489802014;-1.4726479946232562401;-6.8068313158498616744;0.014437213507863868103;2.9906386043826822529;4.4335233186532114757;3.5947257827625831439;-9.9625581363144419811;3.1333648373360243333;1.778864819240630224];
IW1_1 = [0.43910581005138998867 0.67844347368514434926 -0.35387679752307393821 0.087400130863352021304 -0.040834856149644362111 0.3559261953641529308 -0.1403434449336841694;0.055502570200373552156 0.046914646710870734625 -0.46105081339036135812 -0.17213488171125770898 0.075863054474819266448 -1.0321693704328633867 0.536214878673074824;-2.9352941432500307961 -3.1038493854731501465 -3.0999913065391950795 2.1534588218543446736 -3.8597945219015179497 8.2321802243627217877 -1.6958663586079107954;-0.083945748764272049014 -0.076366722546591159215 0.071283550006022877721 0.88742435216719262048 1.1641193997868719556 -0.95754311375291722985 -0.11587822589864268941;-0.11028003479003760379 0.16984151287794815866 0.91855829753698969853 0.13655811577132029311 1.0715059889899158918 0.59217373195075317938 0.77251115621962229696;-5.4687339458760284572 -11.476334590274200664 2.517877831255088239 -1.5228635447869636188 6.1630732569345747507 1.6843510174350067032 3.2277932814294607589;-0.074204783174662544698 -0.10662143891444081367 0.26835864175120888175 1.5261019977094547606 -1.4148642719847372984 0.99442865950482073334 0.59189507070858748961;2.0444671419456543759 -0.61305566949725964765 0.20310117759236948065 -0.44232708981169716056 -2.4019353034227091115 4.7210852234336080357 3.001907806557970293;-0.32164811543606414146 -0.087735730593771502073 1.1512753464195384367 -0.0050273928523490336862 0.35492704636936867058 2.2108089564469883825 0.93688243830039752424;3.1658725284397153388 -0.57295286495743069644 0.81574044855259442155 -0.56137391714690443312 0.0052890991215935106648 -0.35273482762965002735 1.3147654212336163937;-1.3084226320777112562 1.1271792350152378681 -0.66882915823542232747 -0.09674749477130160924 -2.9647166852612394372 -4.2716392314609441527 -2.8813640237073490091;1.8647774236311636109 -1.8988801706280045156 -4.9329354242274749964 0.94238234280605215343 -9.6972062628096598047 5.0265629685532289272 1.288879143269413996;-1.8094728421748187408 0.07968046729259695149 -0.053601753277910529905 -0.15191325137135242018 -0.98499978261487886044 -0.31889744741563463837 -0.063105044873644919323;-0.048034066182128187794 0.17692033173303131632 0.8660751366530662354 -2.8755898138942241893 -19.881229722748507527 3.7946538500620330403 -2.4315301039007630735;0.016357453504404012579 -0.048921460730847105003 0.35140074270183691008 -0.0099641579268084541632 0.32996280528115778941 0.50064576516363723524 -0.49913752568249625252;2.7673009933359686841 -2.4444618080798896642 2.3078164855192593308 1.0229596303224146592 8.4539374345874946926 -5.1334753806657413833 6.252471045255016513;-0.093533755672968094119 -0.098759291829196349521 0.42290212752680017561 1.5579569391510024179 -1.4963102747606706888 1.2500998876845952346 0.74592822560476323979;-0.069937418589659308688 0.063001240559564436472 0.76155599356881220618 0.41946172126900399535 -1.9988634262969722233 1.5847454649422592077 0.24670965317001303241;0.11399215355936326288 0.00096957830194692516505 0.24191802486158298091 0.38679845955398589608 13.954582552379296345 7.0400989777401603575 0.7910520898059497652;-0.18751230897291107369 -0.019044137077983525064 1.0274171419820101736 0.018767788382728238517 -0.078046970939614379592 0.9582910144842001765 0.72517697572888140112;-0.13669167144227692412 0.012970270843832612181 -0.66306513421063717484 0.66288822516124912276 -1.5006304309424882693 1.0124626487066210601 -0.31022854749478262493;-0.58494531039950214968 1.3141814430585603279 0.28016756792460500991 0.13769137017416058222 -0.57938553490007083102 0.19205330054799144013 -0.03854305471837085606;-0.12902164974346930215 -0.0066824074837265635962 -0.22640958657627255346 -0.2821669252943684314 -13.722749913788110376 -7.1445370346583443322 -0.81587494913108538075;-3.992381509023148034 0.98599722389993116334 -5.490223453287848443 1.0285361214149453257 -1.6844499972231476992 3.8838823557576693801 -1.0745669101878496665;0.37491168444647393621 -0.086809032377034647809 -2.8143785936450202989 -0.10682538704800798035 2.8999751579646937394 -0.40411476250817618405 -0.75780814306606170128;2.153887938611830144 -1.8569682088479981807 -0.12877218495864189984 2.6361927129070230613 -5.5692562147992461519 -1.7183509798489535747 -1.8070763878134736835;0.93000326840964664399 0.95223120602076360974 0.56788047728481338972 0.29697070108340134498 4.204882076465547236 1.9053506083392033332 0.88858055139584846138;-2.0770327661399159425 -1.291745533282205427 -0.86755406346961450037 3.1260821500691906571 -1.9430478649969045168 -10.024517601542047629 -4.5532248422756840256;1.6316492315483412945 -1.9543023420849614524 0.13819878458518591868 2.0143059644233960448 -0.36425648703394886141 -2.3390535963732297375 -1.5971900742473221779;-0.098830362060216947473 -0.0099573203150175610548 0.38255591081378603757 0.026198095873490388763 -0.76905729197781624862 1.2132989135989564744 0.15704633008126589466];

% Layer 2
b2 = [187.84914924367168965;-9.9172415835338263435;-115.43316668633595157;-56.92905511707988353;30.95635798849645326;-163.31815199205024669;-57.089345468849678866;-40.481389179166392012;-94.63551642986631407;187.93978587806981295];
LW2_1 = [108.65663595740640801 -13.674055523042852656 8.0111143280604615313 39.427942980455043198 -7.2189413911917155176 22.310229251518887139 -15.004992489714238602 42.198165624159422293 -11.285754640325540876 16.573776803425545268 -13.242340124066210194 7.3547276593567891823 53.327946125968622937 7.0435269443128136757 17.697387117945670099 226.2220302614894365 12.24110856176475437 74.946006252783604396 -25.743286956792921671 32.01472649133417292 -14.768308595408228712 -20.18958997428903146 -9.5961131041326499513 -3.358957733511392707 42.521790160079937948 -4.2838480701057068956 12.001566254989059601 22.036604317904380679 5.1618991271528340548 -53.710341511101361789;1.5331340825858825738 -5.6217245263027075808 0.039407626551283353455 3.5932312130304726594 -3.5761853489890369495 -0.016289157680805490397 -5.1636888233645210633 -0.062561467675486565065 -2.3122518528084419209 0.28362594400125679117 -0.085108394693415620469 -0.0055422633620988140044 -0.28803509169021107583 -0.10736941716774274502 -8.4896567445320503253 -0.025938531068831476023 4.0238939279221712653 -3.7499504649188621919 4.723904899196464946 3.6410768971190541343 -1.8262268922257447468 -0.35320109897404849963 4.7144735779290671829 -0.030772879441196682865 0.46742532598504205055 0.25085324291059962754 -0.38165083836690960606 0.031406201222818126706 -0.18735376793027266884 16.315343582380780418;8.5768482970361521467 27.55998404022750492 3.2889395064426643067 -20.607773061508414258 -8.6447206258737807616 -6.0999909750516883022 -41.550187663864768695 70.514509905534140444 5.6377911384209591361 29.951523469649494302 -29.423697497559871294 -8.0862426170993604302 -47.187644588349002106 -55.289980217169983234 8.7469062108914368281 -171.66261126192398478 22.156771252795596894 18.053898339885407154 6.8250284578430484572 -4.805795955483779025 23.592693707336668041 -31.700606061361806809 12.863202875957984617 27.535410814812063762 20.706378484274527096 -4.5614496457191284762 -37.789931401384194487 40.979458614409971062 -3.6675300251648215166 54.02250798231205664;18.19811717375684168 -10.79656479632127386 1.3841810581364040189 -33.073504459620245655 -4.2572250361941401664 -1.3319128251776080063 15.769409464665544718 0.88559864070317728668 7.7748281905278684434 -6.7338899608175122324 17.064911003340675677 7.2895726151252784319 -16.179303526275258918 10.782333000306739024 3.8942743948134959098 -27.028389354975647763 -2.9740388586985808672 14.684313326784312892 -16.304353327310924016 -2.9120149719317796055 -5.4313555562783726316 -3.8282705292030501099 -20.785961268741754537 7.7566867690918197553 -11.899181114534471249 -1.5520636953561739446 16.094863256539248653 -9.7365239139430546089 12.203686983493623686 4.1085036640021250776;0.29983650528022198767 17.271355127941621532 28.998142422952550845 60.631008831328401243 10.617197785690791534 20.823019851697345928 -24.137986750781692535 -28.202646800361414137 -29.718683946802720186 -34.633813623327910136 -3.8229520419665807829 -2.0237585262518069662 -81.930382970476486548 0.78248224136750810231 -52.027748341235643181 118.2227765209070327 -43.846527414517744603 -19.676208811358009143 46.376501377825945838 -0.97652737785042820207 91.510176334297582912 8.4154424639614209269 51.040358309192576769 -10.284886547268586199 58.217911936613639057 8.639820532640859696 -16.835660276556108528 -37.378310147869761693 22.619107077302121667 34.874607520062895105;-13.836533479567791005 -11.266495912319372863 4.4785229869302458638 -8.7350103967379677528 -5.1701792909802239961 -4.8387713400260841112 15.742603822473697051 1.5568141552710850206 -9.1131275805618869867 -8.5659227251777139145 -4.5009809597181789087 17.470775438654015943 14.025389381149283352 35.897125706041236981 -6.3270336782499629891 -128.29783664793302478 -16.173265812573205835 -0.24689506823317911022 18.363874242507641554 -26.444598992322408293 -20.606238008303900955 -1.9963965452949286039 17.950651753912964637 -12.626442653862797627 -2.9609470099354022743 -1.5401074342698228836 1.0189324147624696515 -34.795272350383392279 8.2246207037266962914 -16.815589683666580356;-33.613218602487307862 -0.75039879198990266929 -37.531019900788564314 -19.280549934282991131 -7.5951405978791788343 -4.5649520236140688212 1.8562582575333668178 20.216093310213132384 20.357465552463086311 -25.493310799446437187 13.501448345552518049 2.6901323063284210413 -6.74769846965735276 -12.637596311845925356 12.181219346646408042 -177.87083778897991237 4.6789519428677026625 23.963742967135278406 -7.369046847550110968 -5.439825044659158948 62.972681983834505104 0.36445431158192592713 17.069169097544854452 6.0654471326652510399 21.760413069042304102 -28.15680836210274407 -26.145449913969425637 37.47135573163994593 28.602353093804627093 -9.3277095823102307293;6.7466792368955879056 -3.1610522749738598769 5.9828466395298427472 22.033798495617283209 1.2098726467376572646 1.7478575400251397109 3.4435723049203437185 4.2131264887879593672 -11.873128208918478421 -2.9672161725784924613 2.3960908506226674675 -3.0133793080659012631 -5.6122643544317591946 2.8312971959591917148 4.5462056014870144693 -10.548039572625661719 -3.6173076208735013282 -2.6267167737382659531 -2.769448046131727903 26.629237759532522034 -10.26803308215779964 -17.802058353025742576 -2.048928903291587833 -0.10781963776493209972 -5.9112750893707044497 7.909635171913124374 -4.5029717699250815244 -57.833414901105733463 -16.403202076427124467 23.570286802367142798;-6.454212825778458118 -11.812596281964731659 -18.077395631936823861 -3.9201575517017221983 -2.3564833720317808208 -4.2139248425915392104 -4.6713365831659832139 -7.2315417809059239218 -19.511316952112363055 1.6002184719251144607 1.3741047060996360862 11.328236279398373298 7.5408606872832519841 1.7134190938380282621 2.3525279482469048098 -112.74711281615873304 -5.1427106178555952454 26.175818183628464908 -1.2546351203563770316 1.4136815560310718176 -44.522373806873837054 -12.325172431366494763 1.3850112304583344081 -4.1447650572737773089 -15.57535695574226331 -0.45974432119193747726 3.2067878188686846919 -14.322473794307994055 15.003959441209731906 8.4972819419560128296;-0.82927125548527114152 3.3755549902158983322 14.743448417212222878 31.386656160850183284 -6.7296498296306426212 6.8665351214330065943 -0.6606281598807328681 23.40299100919658315 2.1062643637697027188 9.6043333064649338837 -10.480412851857307643 0.5928522853912846502 33.267505186173629284 -22.575031054151391885 -22.268832039200482598 200.55136034654120181 5.0801784752343674967 1.1393087948063387316 -10.402573898678523889 19.59804799040288259 24.068721843466349242 -23.906318512045817926 -11.561762205574202866 9.1797151125765719115 19.758637129764068163 -2.9649740839518727853 3.0233923775066422479 5.2661738800788571169 -26.667398665466230057 -8.7448246606303392525];

% Layer 3
b3 = 1.0248905435799531283;
LW3_2 = [0.007068229711723026526 2.0038476060933865242 -0.0089699538323165968479 -0.0091718292611840113171 0.0053331217461229088114 -0.0049463418488638879675 0.0071960647741908505964 0.019274401024576438141 0.0091022871585013983425 0.0034380695475735999385];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00819772922900357;
y1_step1.xoffset = -2.02;

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
