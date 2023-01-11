function [Y,Xf,Af] = ESPER_silicate_8_Other_1(X,~,~)
%ESPER_SILICATE_8_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:40.
% 
% [Y] = ESPER_silicate_8_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987];
x1_step1.ymin = -1;

% Layer 1
b1 = [-10.025925356682330758;-3.2899766267568950617;-17.601291545882720158;-9.2447301691717420624;-11.295209981399109367;-10.293268656673156158;-19.478832657876512968;-6.5788285399616315274;3.7364663958775481412;5.9203426510073731137;4.5181420475192233255;-15.935840548211311685;0.73902640923877149781;2.6170267807692613182;1.7275895428826331113;3.3252184531367299236;-2.8365814786096734679;-2.3872549061526733283;4.3944313523360420959;2.5617079731871745807;-0.75255380366728941155;-1.5933759583720632502;-1.5164110319581551867;-0.54146828310660399808;-0.73661366104584558023;1.1153926383925414356;4.2775988447034842821;-2.2417857256845152847;-2.6388072791399204142;-2.5758902432563091089;4.8345696171929173346;1.7571493812309855187;-4.714819345757452318;-2.0772098047284410782;0.81923947522339646277;-2.7887099273664706089;9.8697804341947499296;-2.1508453485169023445;9.9017293302266313759;8.7678248201199906475];
IW1_1 = [0.34160192620509227668 0.17305243272113574049 4.8874378403669975768 0.20470048025105291867 1.6049985358138607605 -8.1157275748205375265;0.25116868716887202906 0.35617754762994580098 -1.3320850488951498924 0.35269618765951776806 6.4500762871224175399 -3.0873688014940383617;4.2956970352572838578 2.7666896828730456193 6.0934021003315930187 -5.7352578240872995607 -3.2527310057769365237 -8.3194944074377605858;0.085643053087324630401 -0.04282857372589402023 -1.6658716590586866158 0.45517280507903373765 -6.0069090209609896647 -9.8080836344863477905;0.084208580455323894798 0.0047368631982752677431 -1.9186259572927875094 -0.84951108869104641208 -9.1870049370260975508 -9.9688187450728360517;0.4008462340549949765 5.6155270082892219108 1.4375528749937436235 -2.7001379434608474561 1.586466103489721613 -1.178300025875634649;3.0068418987661513775 8.4740876098765447466 3.5935754621146549503 -2.3276716905243479694 9.7971431578579224464 -14.505569828631275087;1.8851358137934342007 0.79619511319383562409 1.0001419269704954207 -4.9566643248084174544 5.0231726666994331865 -0.74299009401315341261;-0.50160538928177345053 -0.52049416813915316826 0.63795364243505825463 -0.80429649193642238636 -6.7790019829161023424 5.8030187726901711898;0.91061516662434616176 2.7182592589965648422 0.55311971710712870554 0.50094315849364512427 -14.847162049211895862 5.5476943503358278065;-1.191045865956710692 2.0384200203250655292 0.8546805141026146968 -0.28935994385982855004 -4.9937905547928993499 3.6841744141769048326;1.1883808317214668016 7.8545783220065228392 3.5787238232980707053 -5.8351062198983720819 1.9159500132879432321 -2.5971928932946108937;-0.18730963526583113699 -0.040665235168197838567 1.9085272119521379253 2.9778207730062722369 12.05607165082946608 -4.1035450928599965792;-2.6244328350477199585 -0.17355241849475433469 4.7217330354318134411 -0.68045317709182862043 2.7925210113015288549 -0.65302147595541437664;-0.19184796464352010803 -0.52039235437297148845 1.5741545551954017057 -0.1126764219021916813 0.11608543864703116388 0.71328538645889283387;-0.24702925498704636498 -0.41976433563759596446 1.6118724219323294911 -0.24671911041010388677 -6.4708278099517952953 2.791952571109409309;0.69390744319845754084 0.058324264157489834337 -0.12370650966710844409 -0.059591430145328631862 4.5561787648317748634 -3.3096077984120295667;1.9620351120752543839 0.030847886012436744052 -3.4986428569782677123 0.0027848871886272360057 -1.1733444144381079965 -0.14427172173624744356;-0.89704753108017409868 0.0036667900308679363119 3.5116477416771156861 1.9194595400394280915 0.34107017762592650723 -2.2871929409596343241;0.16736969116296596116 -0.30415085034220157301 0.022598345732924407514 -0.33761054691000791372 -4.9070984456176942601 3.3775689551918044451;-0.60704041238314465634 -0.48976819095933893378 5.6092235819182230117 2.1407019380980769618 -11.374758780931486513 5.2820989792362817994;0.20743499200186168996 0.46673602586895124311 -1.3287029624718016407 0.22655703868267906831 0.22141538465878171649 -0.97629259418734037901;-0.59190292241930508865 1.1222517724394038829 5.2698689417361936904 0.70148813987981795126 -2.4257271033268135163 4.2249329727939981538;0.1405373211555084878 0.027096000284773005817 -1.8344673525957650995 -2.7003340033325748948 -11.488431447240749605 4.0714536114861665794;0.091519775639570505321 0.029325978194533149079 -0.84603079245952850229 -0.044626834062895567867 -1.0069513287615965424 -1.0413540421920945978;-1.5235372927013997035 -1.0687024739194186385 -1.5854842007950606941 0.72498758583875511086 -3.4278418769262404275 2.6859968058428691506;4.0857294035490241413 1.5663599303443045407 -6.295631966542130975 7.0635783663754727613 5.8339939447634296954 1.6204664891746680677;0.2036044781807738735 0.19683364150377338775 -2.7221860754790876413 -0.58251532797539118391 -1.1001543734065579727 0.42899939846074847116;1.0894095256331308708 4.7553778707872531584 9.3470379568145869342 2.0398988958959263229 9.6120855720537932143 0.65818413099697214985;-2.0105931606842846193 4.9798384778475615065 2.8154788940426138666 4.8093605199739934974 -6.521164149629972151 -5.010529837775097306;3.1215072053925752016 -2.5059390202403286629 -1.1887005258163767518 0.2488962590374622641 0.26726549733112187068 0.36685110459101949454;-0.22349261358319855164 0.17703826822634696758 11.420426759394178617 3.856283260351749842 8.0132219163698881204 -3.1737307199233795174;-3.0041755920441359784 2.4021821283607134312 1.0717005479639452048 -0.31879651304673795043 -0.34745575512934734341 -0.37903472568444268598;0.55057300372426964952 0.61826366802276910839 -0.56047047964284824051 -0.42680240990149681535 6.4357503719496982342 -3.9311311665420176631;-0.033327341053662862258 -0.017436919920978129206 0.94980587221512768892 -0.010459756615172129124 1.1015832500079065781 0.91741271374839983288;-0.54773261997342759244 0.39964186264537310267 2.8664924501905342247 -2.1601576006504568106 -8.9316486522636147072 1.9364315923386199447;-8.0477359315857857069 -0.30799018844868053524 0.50078866986383163695 2.8412858316234430589 10.21628820757387146 10.361629149950157824;-0.54969577549704506936 1.0193169111579452046 4.967590794599002102 -0.11221746592274638166 -1.7151883500766178514 3.8461050393702640449;-0.15297896753880332876 0.0044296525324329354828 -1.441214505005746771 -0.40747814410025201548 -3.0840467854394688629 9.6410592476599408229;-0.36619209421841436569 -0.13800468549902600301 -4.4402204007097259009 -0.37814892813160078111 -1.4608046695586851538 6.9441529171433646894];

% Layer 2
b2 = -2.6122755224644498284;
LW2_1 = [-5.2378036700841077433 -7.9587455713132779778 0.79690396512884653646 0.56726300964886455791 -0.52224373088508824203 5.1619189196798913599 -0.10080351144545265929 0.19780742338275181247 0.64524618440109349216 -0.12462624128785479427 -0.26506494138748870615 -1.5083210569907068432 -1.5810508751907130787 -0.27696682909491560087 5.2840460203884864754 -6.055560762929570906 0.97808305466509237736 -0.43024565554224969421 0.18535366760487567706 -1.6171293913148765409 -0.067467083841287370238 6.6079352115497433928 -0.54301365967648640432 -1.8535976254659642937 3.73357798945248609 0.1333382255902991198 -0.031652561047937073646 0.67682209754070177254 -0.02661086349714944857 0.024849223358582733096 -2.7144544361528124554 -0.032819883053974249765 -2.7658724700512662409 -0.55082165158805529437 4.3501378541754256446 0.12093417247763739442 0.017611113499921626085 0.54640373576766165797 8.8087391873902607387 -6.8604962015256782948];

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
