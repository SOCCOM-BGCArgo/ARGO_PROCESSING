function [Y,Xf,Af] = ESPER_tco2_13_Atl_3(X,~,~)
%ESPER_TCO2_13_ATL_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:58.
% 
% [Y] = ESPER_tco2_13_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-0.2178;-177.232103709389];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0461352500991908;0.00425400504498441];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.0543854095167942;5.5585406725944191;-5.0493028388118217;1.0404301955912927;1.9121852955213654;-1.0079913688445237;0.87229130124680154;7.5226795551552001;3.5008602896051704;6.3088192702335695;-0.99458795029924085;-1.442575785143994;-3.2458849529162781;0.32767366975219026;1.2326350147242897;-4.6828187843418814;1.0477684852626865;4.2282368405630653;-11.437981654352596;5.6382788307959082;2.5743889521586061;-1.2141081650570138;2.0479993184597758;5.6804345627758694;-2.7188652674045759];
IW1_1 = [-0.059309864179448064 -0.34895313813897821 0.40256424568510335 -0.24816670715025177 1.5414127060150333 -1.8016999943838021 0.85263012663034476;4.0685084513848642 -10.103996961780693 -4.4170782997517195 2.5057076732639301 -23.912305270874626 -2.0630892778458292 1.5187100295042337;0.047045234670579902 0.40102925595682359 0.24864834444595241 0.2032708912512845 5.6308726942382465 -0.14460732070196938 1.0581929170266118;0.15861357289925121 -0.26688317275603091 0.10982543216765783 0.19202083902933631 -2.9097810319559221 -1.2562654978697749 -0.30164216537879318;-0.34676011404687507 1.5551589390753688 0.59759703384199181 1.2539084507547533 0.15734021862037373 1.2534264584890669 -0.72642526473341562;0.47297796401222408 1.5573031522660106 0.57950982107950644 0.41352350738238963 2.4941471986311954 0.84796635508488261 2.1199586605848273;-0.34126687003185052 0.67341542907882979 -1.3190920928651888 -0.10068641327487306 0.80885905334980746 0.35699461880182815 0.76525073168411217;3.0596733117225194 -1.2028998851679691 -0.52969189018709684 0.2474477020159154 -14.726327616105703 2.6591647455755232 -0.95802017074645451;-1.8457322874766191 3.0508238321803818 -1.369009288760116 -0.31744325021985459 -1.0624318408787878 -2.5360160788615622 -1.6849163580969819;0.432862492333766 -0.94229401191631612 0.5975389935405746 0.10354239320906973 -9.0615931374392851 -0.78650200031188922 -1.1277300842796656;-0.22247602717899248 -0.051467370838826208 -0.096066610954548923 -0.16720568209981204 1.6680741541769271 0.45732543280857629 0.3379409324522098;-1.3756296772249685 0.012649568781626387 0.004127257914395826 0.74848704409425981 4.0050610504143291 -0.072889003132199093 -1.2139769126386404;0.4628755250620975 -0.68719053689308651 0.029200453401852219 0.54945345733395146 4.0048064708641231 -0.014789776293774863 -0.44042348720082231;0.19712286556993519 -0.050180692318360323 0.23856492659622897 0.83412168964332567 0.32419360604908393 1.0981245014969021 1.3406973628556413;0.37295945555149018 0.23181463593114401 0.01256924426353194 0.28054364432556816 -1.8326663757648169 -0.31264244629853871 -0.59222501057150057;0.14616312201552981 -0.49306652666452305 -0.053279261661430072 -1.0195182269980023 4.9913891281935605 -0.15495769206885437 0.68964640022501078;-0.73471621063115289 -0.3150594539696871 0.26374252369962931 0.13028441316695902 -2.0552142781594114 -0.26179488184283572 1.2517311286813027;-0.8722973094155555 1.0216159927402539 -0.57826432174464382 1.2676111980129994 -3.3835130458066001 -0.20828146754472462 -3.4431556679974511;0.76113991569916517 -3.5146075156526964 -3.4552867085036456 -4.237565318725232 4.8944432289890383 -1.7468197565664891 2.7641520714252947;0.072661678239665281 -0.065590786310794402 -0.80190790564588765 -0.25592812707574492 -6.8957828175456637 0.0050903592580527892 -1.1979073208907818;-3.4132461859353493 0.86544912705332366 -0.26897824681107857 2.1108724744275404 1.9568435169129932 1.1036775388501072 -1.3343490285536193;-0.93459749036678896 6.1837565542324331 -1.2055972957300338 0.70056554102227675 12.930557397149171 -2.2018867540215945 1.2131742790715636;-1.7270428925075361 2.3544237812793116 -1.221181352684152 0.53602188029877662 -0.5811931392422881 -0.48072074308589835 -0.88022523714111967;-2.7004521259950733 -3.121485781903703 -2.9822144246172932 3.9396727595327325 -0.44395026151662331 -1.8440911488415173 0.36366173549878633;5.7758787700763028 3.4513303485469145 4.5446748851526788 1.3888786001904219 6.6415964058313994 -3.7472919191315901 2.6302656969679266];

% Layer 2
b2 = [47.082033620491714;-6.3317999489854495;-240.79066453452802;-5.1712471040419414;37.026219315039548;-50.798245226913629;-13.848811151442654;-166.61564844199836;33.577536501480957;27.907937547962856;-196.65491109013294;35.573224205222608;83.036461196007622;-31.564437115161439;21.811631708011234];
LW2_1 = [14.69529291687523 0.94110434715319646 22.703779809972808 97.627416694840804 15.064180127136474 10.198807660695588 -23.919120220358714 10.913823857288884 34.880076365911165 81.430321028135438 74.516225771019165 49.647382742200101 -106.73232173809969 2.8483174735779975 4.9854274749828988 -13.49337305401386 -64.887009445966456 30.256019923797101 11.416526753498905 -13.869929173983985 65.411997058407053 -9.3708000593095253 7.5613611087668424 -3.6540945401532317 -20.562573373242532;6.0916691138732686 1.0013003927071811 -2.5623054546029183 -15.938627513449079 -35.137672895195664 29.462514526190432 -16.07236911552685 -17.538459904697241 -8.6032122109386417 -3.9638358793401975 -48.480301879544498 26.420102313673929 -46.341235497811446 -11.754257689325964 -14.984383628948285 -15.249126016167979 14.398337600665586 23.593504954674291 -12.083556334777626 23.515698335765695 -1.9644166085674026 6.0803405057493283 -4.1089712200678647 12.11499771243645 -0.055731261742387898;150.70422832083116 -106.64639638195007 -20.600904607364203 -51.390034292268737 -160.61792616178249 163.77334472755655 3.7027880027840641 -43.992329431400762 -169.0025103979805 -3.9957250094117325 -131.45016114331042 -74.550423055067441 11.849854934376433 -15.356855465094883 50.850163400633214 114.01636065078907 -41.475027791434037 114.92952175885056 -26.330479266401984 231.17563003573252 33.773313062045702 -40.19628340739699 131.54789761053655 9.2147965290504033 1.7677242367187216;-0.4771354751088141 -0.01623987944138677 -3.7789420012468957 -1.414334625912169 0.1460566746129251 -0.18759344165198513 0.55266833735724319 -0.047694121407913681 0.04430538012121582 -0.67593345638835911 -8.4157551273196436 0.11966169814548828 -0.47506029843584235 0.38285895562240102 -3.7923461806181384 0.57796621842627027 -0.53237312555069272 0.17486225034619421 0.044099922558713811 -1.2250812901502031 0.042526611832034801 -0.086536027435509033 -0.200315523418054 0.059499023547177718 -0.054938245514578224;-26.274750798276976 -1.8946474318038986 -26.485153039104858 18.810967241352891 -7.42174516033249 -13.837104397441456 -30.091414561749538 -12.559667869620926 -0.11859622458605851 -5.257721654765672 18.400054327018893 -26.202490605064565 55.792886599641626 -17.493433908466258 -1.7593866999624967 -9.12772046135483 31.611545489293395 -15.655745353666337 -12.79257900122534 -20.24204287161399 -10.916477510754239 -17.994694641100224 15.805352289456243 -25.147799244723529 -5.8683789450247641;40.268210466743973 -1.1936185764985607 -20.28182794995352 -42.582689466137843 -0.6598396993451634 -3.8598480310676551 -11.168436748756259 10.572041355932852 -20.764578420826698 23.589927125555278 97.192019894235969 49.297011900609021 11.710448717366555 -24.791208192037335 4.0112368516519794 -68.589198490494155 8.2665493979935754 -23.678145985510806 -10.1475383274873 -11.112474231892163 -28.247385421762168 -5.0968012100840907 -32.894675431390048 -3.9571196528803774 -2.3681993893894697;5.7879841667948106 -14.084821089064627 -2.5296383163003613 22.914231411923705 -0.64071856120660931 -13.761003862726138 4.0299580492712543 0.057791035824585016 -17.175331380143437 1.8932217188259011 -13.828329082481929 8.6286132454198174 -9.8764876141804052 -4.7910062303326955 13.573453333340989 -9.2882644620398551 15.189580140769941 1.0759245693959258 6.0966481072667813 6.8468697424084031 -11.676806864361 1.0707038513767264 -34.156030785757153 20.104434614752925 2.513626867275879;45.414681568542129 -29.211759276977471 -90.611377568276836 -21.063601584138166 141.42980355087738 -56.305891306971311 20.143136286641944 48.92804591316105 48.692198889193996 -15.608227611490667 -48.951445803242002 91.389220352353888 31.831026575429245 -19.812172206155811 -1.0020489087215103 -33.707989462148092 -44.743990087519286 -38.27701077086288 -17.34843948616237 65.62838101854642 -83.56861956871046 68.628456336144197 10.591360829680367 -87.63365329871985 40.549652807635326;70.226369525950034 125.95745786783014 -62.523021386272987 -122.98180919973467 -150.58837994737269 16.773739472844714 35.004818799718265 123.73190215453361 21.407348053492203 -14.733108802549571 15.913315670898214 -90.367191753848914 -45.079157545907805 -45.432167133896179 -11.337456897268936 -4.1112634818421876 -130.14275829425787 -51.992889554082439 -75.038402555667147 52.044496277611358 -89.611057948432872 -57.469742300590397 96.285177828900501 53.294207693507325 86.515971838940573;-10.044211110603928 3.3874669695524955 7.2299782140419078 27.46479696584186 6.1091789709057949 -3.988273490296292 -2.8849573753637623 -2.8996277015879244 10.518212780418057 9.1971166693852862 15.291540822546072 -9.5274980679796073 7.2631002195663106 19.208567376498664 -8.3613719841966514 12.076039171470155 19.693282921568432 -3.0569720157423022 0.33781733579001316 4.0409523472221514 0.9609609354304276 3.975999037858267 13.882230269309186 -0.64260960470731932 9.9254400497654043;82.69879923172185 20.247100628325057 -106.81840526499228 12.266657048695905 16.941340612268206 80.957371179601736 19.661346209944298 -93.705918335593211 -38.274818072457819 -108.62699254699812 -59.382569345838689 -51.44113468100845 -39.175352097190299 11.043824517078619 -32.359147983788475 66.74737762641368 37.451700129324955 12.965019411893589 -19.96197568121794 0.72741288734696263 125.34131000185486 59.076106068750335 -62.397484029941367 -4.7418645726010489 0.50427695438481379;-2.6195589587304831 1.8512913618158169 30.880147585893678 14.07638289734472 -4.3993912520403589 -2.4505051800191944 14.139072811639339 -1.9379309495170673 -7.6366685306522504 -3.2664363096958402 36.928351461672207 -3.4309468473193667 25.683218783277319 -11.043063324138041 8.5759426332584745 10.416836587185822 33.250099960896172 23.037588437475598 -0.23638485793318495 -4.150605916012946 -6.2992994226320755 4.4436661033975993 -30.215443226517767 -17.389450273152722 1.9390657878387512;-25.155584949274822 3.8294000406100457 -55.256072871512394 44.206463262532147 111.45491184954423 41.376853005792206 -42.053613791886171 4.6587664413677947 26.923214216047047 20.70252878927511 -4.4031126374515175 71.637325700521515 38.901523046595464 -45.964774453645141 16.339939647076122 4.8548014112588556 58.05112070564762 34.642659227795662 57.712169084717274 -5.7055828642923512 -25.606111063898609 -3.7511543038213899 -100.24876607951879 27.663966946576821 -0.84461808405344363;-3.8302721890553899 -25.723491508735375 7.4579709467902422 -30.783457995040187 43.114389772933364 7.0596554771192386 -2.4197305501286457 -3.4814720811241928 -11.477929890435759 25.455588017038167 -44.601289345726514 2.9882711947299794 -70.562469412515512 22.663129115545672 -5.9811393480777326 36.007741733519673 -49.350916282061874 -47.427886217307261 16.75462311503831 14.936900885960823 -19.269448441526034 -3.2206133178914547 24.166595147413773 15.220821433618553 -2.6898122351537208;-34.145417548082015 -8.6890332785299247 41.143523221099805 12.313507599724455 23.079221701838936 1.3232453669115043 -19.126600137878309 -13.609696780225194 19.866368703551835 -11.795365926249474 64.838050737908873 -24.998849862157041 40.915858698876491 -28.232568366593025 4.4908675910913933 -33.590779839390748 2.8298068709569368 -24.340612264851472 16.287624318508364 35.089764894826317 19.153712730460136 16.598956377822233 0.82159268533093599 2.1007702495264562 -1.4710882108985943];

% Layer 3
b3 = -10.173709569948468;
LW3_2 = [0.0015981084144487527 -0.0051811271364958847 -0.0020435681511736816 -11.029100461995379 -0.0039925513803319651 0.0048779807661315736 -0.0045062396315832262 -0.0011268203081921981 -0.0010808195528572166 -0.0047898369480900574 0.0017513559452885308 0.0028509892950635946 0.0026520455700049651 0.0019024081703822131 0.0049486928182047784];

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
