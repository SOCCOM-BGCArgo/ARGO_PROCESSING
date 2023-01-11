function [Y,Xf,Af] = ESPER_phosphate_8_Atl_1(X,~,~)
%ESPER_PHOSPHATE_8_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:28.
% 
% [Y] = ESPER_phosphate_8_Atl_1(X,~,~) takes these arguments:
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
b1 = [10.236032832650911573;-3.1213185434018333453;-3.3160071161601636369;3.6305137388862900316;1.1123668940254993753;1.8568565026301633303;-0.088407327806362581701;-0.97596208126849548492;3.7318990407280767663;3.2417136158077828334;-5.0501793944074488962;-3.7510889170707741869;4.8884860790863475799;2.3655013843337449053;-5.7891758323582669377;-0.75239489138212733987;-3.0804964717324936352;1.331666501815416126;-9.6834148676724076665;-7.7231227968425013586;-2.657127867261938281;-2.6531362242698204312;0.62198106346249315557;0.50250767096057713523;-2.2387745850554403582;11.381587091692610514;-0.69112488427216733911;3.3235341007497019028;7.6695186930452603491;2.0324160844912633195;0.01128478471246349181;-0.10042481050330696879;1.0624338111583786493;-8.2606316832723720722;2.1056774706408427633;1.1265287017213434506;-2.8666476405011382944;2.8585064029453048917;-0.010087264900922337033;-4.128954304392626895];
IW1_1 = [0.017122120544290634486 0.89482516495831865289 0.33185763650068283992 10.287350125916601584 3.3209053395435179645 -0.38437945796480765814;0.33706295071907355076 0.52866803714867527386 0.41603320177933733603 -0.87276123881319689346 3.2072123646871775193 0.17760833343679632845;0.27275368036922514881 1.5845302900483033692 2.1859058266241397916 0.55278417378578803554 -0.72008848865029295716 0.45927213327412547139;0.066815444087090905545 -2.1167412652554142838 1.0429635489478013 1.0741959641989153873 -2.8526902646009166808 -0.21546171775236475088;-0.35833072298881996076 -0.71955287771274789588 1.1307193493854763577 1.3505944608795941164 -0.091891398525120751573 -0.48012885668907545611;-0.55501789327504313931 0.33844258181966330179 -1.6354852449196715103 -1.2083445060693251971 -1.465537642551967723 -1.8726963112981800386;0.43069674460802664928 0.12113640165448399943 -1.0078841745122228968 0.72242270167361233302 -0.22556589484440633941 1.0525244240916118432;-0.2214903474055307242 -0.50990115194262064513 1.3911018061348621178 1.5467027732209508439 2.1091430342827930389 -1.7129503578513478601;-0.046104149068035588099 0.33447500032821669835 -1.115651475464036313 0.3987546624123721184 -1.4890122240033820677 1.5237474014255909527;-3.1740655816151925883 -0.74501265858799303743 1.7327393452865540446 -0.72652399075936158201 -3.8345563982222632404 0.70106321026164952759;1.2824389942961704758 -0.3696536228726921558 -4.9132588302001130032 -0.11275430342151021579 4.5821989981966346761 -1.8304428661228087094;0.95725981484988487136 3.447406338264781045 2.4607018787870522658 1.4696801023493957228 2.8286137969757110255 -3.0285728471532191719;-0.83064991091188211669 -1.315014915430627207 -2.2231513699492406744 -1.0808445830093071915 -4.451831241573334097 1.4706093300968940962;-0.29990354248575118445 -4.9524532118414317239 -1.3252809440566535315 2.1089953891270942776 3.4420715903854604001 6.9412983523705840128;-3.6067467138130475668 6.8672371917517995854 -7.4659235668386170204 1.2992800479338346697 16.754118279103700218 -1.9636967400786757842;-0.28244991352384107852 0.81929409133048369807 0.99911748886119411406 -0.16893294640032677223 2.7049711878081645899 -0.13858292608736733076;0.14526716257680427624 -0.14421640589003553723 0.11936283157833790791 -0.49451936702587900863 1.1010127366348723132 -1.7913722208764271127;-2.2804712473832369213 -1.3130867278012490207 0.45162041272039776807 0.68154906620310728993 -0.33842850633981275132 1.1050870236146161929;-0.37911192447678121464 1.2364961635204281265 3.1150845969347873421 0.51187310541777530926 8.5080189910239916173 -1.0322020207848869866;0.35358022602865341044 0.21254349175600101285 -0.58306087880066859164 -20.163935290346120865 -16.050454774487825205 2.3315424111812474983;-0.24864749513144049353 0.081233465376188765394 -0.34560183804047001077 3.2153605421787925245 7.1451874495460776515 -2.2606624608683278943;-4.2554067529834105343 -0.84768493973451308765 -4.9163497885628233419 -0.010678341648426362132 4.3225438678591032726 0.6433458578071736822;-0.42288187904194440847 -0.65077123144510995445 2.8565988395241244113 0.25273030591068584849 -4.5756364656889871867 -0.85901545547716062767;0.37007280341765413922 0.83267117420113945236 -0.52939582537143259344 -0.54226793679314144736 -2.412806844296182085 -1.9213306572571799968;-0.2762993395769856031 -0.43432321931758161027 -0.31459719594759966155 -2.2099107318183821569 -1.5562934426741408256 -0.19751643692485906789;0.38210350741497461913 1.5042988075289920324 1.7231458143662499172 2.7667165916689326899 -20.672999068364511288 -3.9459590415262377761;0.019200104870832249149 0.61204918767999572093 2.2759445024501339105 3.3254641824939001893 7.7806307642761085575 -0.92210913164303565104;-0.25086985235220271573 -0.31170826557851416538 -1.3256187694561913926 0.73467114983616232937 -3.1149779313953138704 -0.44965623307822449872;-0.19966151672611642809 0.39838337266177153984 -1.5862272918763742346 -3.1891525780264835532 -12.215192873838260468 2.3586344543633495086;0.57929923110031067424 0.73728599800129457353 0.2877096587604046074 3.2544866252228961123 3.8120324137122336694 1.2286267520500748773;-3.4946984783323209456 0.9480916496800138038 5.1104781670825847328 -0.40759030626848735679 -4.1197472019534000509 1.1662767124118773587;0.29863635947427619177 0.11335876462865691894 -0.041175437945495038661 -1.869263311784288506 0.3790079134331756916 1.8249452613606882423;-0.14588521321674841058 -0.79668776825850529111 0.9663798170583045577 -0.56223240020522113891 2.7747017073436066603 3.3225826752238867279;-0.41073247864987733058 0.99527625376676331914 2.2784578499504712923 0.28647333470200697647 7.2796180843658628845 -1.2411332472977514296;1.6457217451020842436 -0.62976262744604749599 -0.24443432087723079849 -0.18134399305628706722 0.00048281127262862896452 0.012432804130289947658;0.20294975954188368417 0.029149553584167646469 0.68758190855346756898 0.037111416433174268525 -1.3314458240509077225 1.3391754131685078555;-0.75280492848317459575 -0.53998716252557443784 1.0677295402535249558 -0.35564343694040478905 -1.2511286151895755214 -1.0672380091641127109;1.542572610625132068 1.3622384835465013442 -0.45605536996196943678 0.17114210803322568721 0.22432727713128761127 0.32028384375850188492;-0.84535478016080745078 -1.2011169345653018858 1.251937698699544077 -2.7728050631320724229 -8.0295193549192145355 -1.1675346338623391329;0.11867377619695747282 -0.099946605374503977304 3.2700331680598990758 0.8036105672401129274 2.2956615125365154029 1.202356604128281603];

% Layer 2
b2 = 5.9907017595374156826;
LW2_1 = [3.2908030655045275026 2.9260253022340405593 3.2241036760871466527 -0.40592832930546113301 -1.6022484062982045305 -1.6353093885474538904 -2.772047006518286949 1.648041104246813271 -7.4908306297275588648 0.1999847974298337161 0.34409362330685916431 0.26179907612614616852 0.65075880792926155127 0.11635059041670808844 -0.084916985975271280784 -2.8713766557560496118 -5.8922632294427144117 -0.21545947472674983891 2.2188583541383901654 -0.13874914202813842801 1.1203508241254473976 -0.10563650326156526216 -1.2288769222002151604 0.70236487148026705007 6.6196142730022113909 0.15377005513692085237 0.59540843647271013417 2.7103608537799055433 0.49515648557661684492 1.9592462604244016422 -0.28477693982444807208 -0.78167096282000814167 -1.0599767615916315666 -3.3663068156162840872 -4.7541599771042051259 2.3456084946416724257 -6.6115875998800932578 -3.7226862947685450678 0.73517151024168869711 -1.91144537638158174];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.604229607250755;
y1_step1.xoffset = -0.03;

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
    a2 = repmat(b2,1,Q) + LW2_1*a1;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a2,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

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