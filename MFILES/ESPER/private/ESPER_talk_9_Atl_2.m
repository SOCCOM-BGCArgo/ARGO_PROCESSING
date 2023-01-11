function [Y,Xf,Af] = ESPER_talk_9_Atl_2(X,~,~)
%ESPER_TALK_9_ATL_2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:47.
% 
% [Y] = ESPER_talk_9_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-0.2178;-143.241523403244;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0479365703301392;0.00458552991639557;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.102168915554917;1.0241829297572467;-0.25695374755251454;1.1001581615980396;-2.5099828493856791;0.19628908847223264;1.0790372815325171;1.2920336726825663;-0.91836978118742008;-0.91825075449756444;-1.2133048739631036;-0.27090578184162878;2.7368509087344792;0.89563800638378188;0.56958576563547425;-2.4297835519867546;-2.3629007271535341;-2.9390731456018586;1.6771921349768417;2.3669988393569197];
IW1_1 = [-0.66681678811847367 0.0082687443759658572 -0.51425437137264851 -0.6019046417181696 0.64026322199568131 -0.84602119028985523 0.29247012836927089 0.5876330023387113;-1.3074008063596885 -0.65329959494352929 0.19586026733258385 -0.32717667611378815 2.3860666085753195 0.94705554145653326 -0.82404506326968563 1.2265775711318452;-0.7466300026690097 4.4071004667702001 -0.37432187664987976 -1.1513802797220294 -1.5619951497685109 -1.1683708072608743 -1.6628362850486675 2.1546195897002942;-3.4358174165348072 2.3045540701232938 -0.34676646184211263 -0.28675288480719213 1.8792452180801225 0.79734594761719491 -0.44537638728172452 -2.4633482031982057;3.8665200716602688 -1.4359731946374958 -2.91404216701754 0.23472400337202362 2.569138355223612 0.6817216992400833 0.46866013493513681 -1.5549450455193266;0.17215870399255265 -0.33124351643573208 -2.5188010843950432 0.66842191693336261 1.7564673047723169 0.83450898282556485 0.91261875094534139 -1.5379799355006414;-1.2694758185971742 -0.046130325267446134 -0.11766535461310511 -0.64139210510540068 -0.27808143069688923 -0.061481274507360004 -0.27044775262974602 -1.6342787730382522;-0.81756311932997416 -0.51458202630974204 0.49630547486011506 -0.15143345460955968 -3.922051610676728 -0.40764148009909446 0.46926451094078631 -1.7363263538498883;-0.46802880775445876 -3.0737443729647804 0.47814167316258771 0.044500479971561445 1.7485923251740136 0.62083987590367451 0.35885427566367017 -2.2410130551177936;0.33952006471850221 1.0472736703728565 2.9285573217584844 1.5518730506827727 1.0718438685322569 1.1514048418797238 1.1737195133485983 -0.81357863387286178;0.9950748322023405 -0.62490145673393294 0.59856192931414132 0.25614581737401426 1.9958976925225007 -0.1132378164022475 -0.34874109531075637 1.5543950201945451;-0.25245747090202031 0.45641721869973811 -0.34263675598683913 0.49609514719183911 -0.24392376103846491 -0.31065830313083032 -1.9892399780587933 -1.4007619643612901;0.45173230200477538 0.50493359981171893 -0.25641832661340908 -0.067210608640949368 -3.3307460641847837 0.095813354852626903 -0.047426741751208375 0.58994572559515535;0.75063705085835353 -1.0769934568900952 1.5617005954029841 0.33879410085152845 1.3327431671985854 -0.21676546859717546 1.4434141660848887 2.949742883888923;1.7516839777456428 2.0761244043698239 -0.08324462689602552 0.19025425661055434 1.9441780728118856 0.30417184680310244 0.40929892110412291 1.934937699864508;-2.2234140338109558 0.27874997586860201 -1.5902205024187999 -3.0394445748649432 -0.3863520362838172 -2.0440486444241346 0.023322860967082484 3.3091410390235385;-0.32445929813240987 0.77031851134410068 1.3552395034411848 -0.81280815968631048 1.7765232137831888 -0.66875141785912862 0.18428272641117635 0.77851533564164255;-1.7528243419274705 -0.36972355193963846 -5.5467363991235787 -1.0817024213543021 1.6160840322186552 -0.57758026322688438 0.0010535146715402698 -2.6926662906394658;0.55254432593183944 -0.61450638299023419 2.021658920989629 0.29888105470016996 -2.0910945170045681 -0.28720377537556308 -0.98525822054257284 1.0231461742696815;1.2040014750619288 -0.062394179361396906 1.0232683153057587 1.9631272745210315 2.1231385359538391 0.33864473196532152 -0.043553675330851603 0.048109766169979799];

% Layer 2
b2 = [0.87333122760635362;-0.67357126083211216;2.4574618835641244;1.6459652017726081;2.4054399838043343;0.28858725795802326;0.15227136089160798;1.2437443294101573;-0.81020699430629295;-0.7216750024382903;-0.35891682240110467;0.86546544416369198;-1.1119255200771976;0.49978765382510559;1.2683268725947483;1.2328809645896104;-1.4235629470615678;1.6485112749404942;-2.53752760538078;1.225945301043883];
LW2_1 = [-1.1159337949381833 -0.015964856726286912 0.42437403891901782 0.1042197511743783 1.1650599032854745 -1.2858733271234597 -0.49208877397778666 -0.70191894072619065 -0.56159924459403032 -0.30047192737808093 -0.77729287922519552 0.44446537333877478 -0.48723075326602849 0.78260093081222826 -0.11497520659141677 -0.043582953114625907 -0.50650692283059517 -0.55760327134307708 -0.19487545627584277 1.998457869403268;0.91997531412686551 -0.91472229026690866 1.2619741257672166 -0.51316093355171077 0.61307640430193144 1.4185358428513679 -0.99812766023793931 0.085932321634486275 0.47846275711876929 1.1464407241059158 -2.0838269717890165 0.78145181833055122 1.5625965578029819 0.34125776253724016 -1.0707683380935689 1.0027012944444553 -0.15165801758703312 0.39942733096084981 0.3735291049954495 0.96058296480572991;1.416865548172547 0.24996729995759265 0.42498749951975262 -1.3573951886408842 1.0978058962322466 0.33411532418093931 1.2538864086378589 2.3379009911036031 -1.4820419761587835 0.42851401411527001 0.21678652810425755 -1.4706784380765003 1.4136874146743927 0.063217044426187197 -1.8343863198514887 -0.11901423209719503 -3.0246585798568395 2.10282393354065 1.4540699384652773 0.86636790809793718;0.60190040758376262 -1.3560887938944239 1.1179682266541269 -0.78439077071903784 0.093099696465725384 -0.29253033817601759 0.15592455303756578 -0.54397631331957053 -0.97034370370231593 0.35679992990534903 -0.42490113843529431 -0.49378901644979811 0.076218006792650783 -0.62544846417630373 -0.57906499571827175 0.85879588031012621 0.065145671451214429 0.037957934641727996 -0.87838050129342349 -1.0949243727186395;0.29798668885897289 -0.22596987133104529 0.6906699755386565 0.036773728604436999 -0.96127240386631463 -0.79620128628655362 2.2379778007432889 1.2150451166722169 -1.3469347060179946 0.28996485106859171 1.2278794825800838 -0.065565790335456384 0.9907180058517937 -0.37211427075103004 0.5625891613362668 1.4823718167876474 -0.24715162873208729 0.3899489756127329 -0.87987533226281478 1.0297423908222922;0.0041107421763250523 -0.70588404192076881 -0.75043348844638802 0.25996900137801943 0.25302454975416155 0.21712413432421554 -0.46964795541303384 0.086586338065406179 0.45754321313902835 0.5292223340070531 -1.0042344077696677 -0.57078264951939439 -0.34357326705414576 -1.4844764906256103 0.61670353373984677 -0.26563637783160804 -1.1967812360975236 0.0098133481979415862 0.29949116671662995 0.21723573702358953;0.82315506506325775 -1.55415284272529 -0.92524438107239337 -0.61909585311334625 -0.76737605561642153 0.23755993202286937 -0.27980444701560309 -0.3030667166281461 -0.0020734306160397897 0.42710871898368646 -0.81306435850247349 -0.4231673488505226 1.3677772334029465 -0.52993411353442965 -0.87226146901013757 -0.94346648483604756 -0.38483986378747792 0.68775181569255295 0.55503665000746183 -1.1880086868244211;0.24791656086567992 2.0281221130546569 0.65934208231340641 -0.436508066632323 -0.96928384042454796 0.79244969846500413 -0.074571157355752787 -0.17117327722610304 0.77959308274059147 -0.36181753522963422 1.6543878775089842 0.44196327998505375 0.032409760871631568 0.42620282043659519 1.2006932361334108 -0.80656222606608474 0.61980541523251831 0.60796322155257065 -0.71350712340615363 0.18860131534005273;0.1151334894213337 0.46283332075038647 0.79810632055614084 -0.43540416692519507 -0.068699285487179337 -0.75243974589310758 0.91466197318239395 -0.3059332952362841 -0.84770537684510361 -0.45496564131197964 0.73693224624568654 0.27103740171197516 0.080145980054550953 1.3149536688549359 -0.80621804668457386 0.46614412452412851 0.5511920254413285 -0.021971960067662395 -0.28535372546455706 0.61701332439975598;-0.77766103242826012 0.64519027741198764 -1.5267035227046639 -0.35500544302316439 1.1809347629525764 -1.0190171782830719 -1.5548426388711862 -0.25272076819493178 -0.78677239246867736 -0.68060587945038875 -1.4121807521015581 -0.29106946658933874 0.29719848745916594 -1.0552347231990762 -0.036893337877536597 -2.3893365804394247 0.56101551077693979 -0.56758141096734627 1.00384315215435 -0.22333560802726019;-0.19598428580728855 0.14548294970359166 0.20712557854219718 -0.95800110358144108 0.83466708090754083 2.2512400813766384 0.54168257868558412 -0.36512076770285051 1.0023422268717765 -0.5498612631203561 -0.29847682116266344 1.1920364712732103 1.5529494141396176 1.0383279664794982 -0.39795324483601202 -0.30893235911391559 1.7154456870934278 0.19640815817752222 0.14387854099544634 0.34053820195939222;1.1131190916998592 -0.46324129957200516 -0.061493747503055585 -0.17521823862354546 1.1029082820916458 0.34349579116920359 -0.87592590555185768 0.045539686097594165 -0.79389328528232994 0.82928604242991366 -1.469363918789315 0.46167218808718474 0.10441985054602819 -0.019892575930447853 -0.89893674016336933 0.52124849989004829 0.50627706029649955 0.21629077307863245 0.13179995056894245 -0.70241830789606707;-1.7196476536332774 -0.93683921663175507 -0.51306828505822899 0.82259543192551299 0.13458208635188146 1.3256163897966977 -0.60990086330530868 -1.7353594970738624 -0.52853638611171683 -1.260003599408978 -1.9596631029499525 -1.0248199472999806 -0.23388296231621825 -1.6824978077554174 0.72551324674252216 -0.88018274013619546 1.7111616098362719 -3.835773505651225 0.99400657551295679 -1.5589914008322865;0.8329431642690881 0.74067074165460978 1.7300623966684514 0.81967014570647667 -0.27794913189304354 1.3193580235402451 0.2280998790994464 0.78681611420743125 -1.2469048793284612 -0.7858378759826411 -0.57087819808957407 -0.28895917578545077 0.41819922597106385 0.26694442467280044 0.28949828022733504 -2.4115133068991645 0.55801096274198947 -0.36602601157199782 0.154142189954696 1.1788847774677984;0.65914417602764552 0.93018860738545261 -0.90347526107670273 0.81931860900841258 0.92048345682548149 0.55292634484517234 0.25056768463314266 -0.07265185697703419 -0.32134985375910102 -0.37237401328170561 -1.3812653065218434 -0.20729487890521728 0.75492579063000975 -1.8310095071237671 -0.25048296051638091 0.52626679260545883 -0.010760450202851973 -0.44637230678519962 -0.12004767387082602 -0.27606420881899019;0.63463689594613071 0.28865683328893732 0.049194270117913177 0.2500558810473027 -0.67512317260647314 -0.029442954687419538 -0.22117362081481814 0.15903055404562641 -0.01431116866845945 -0.086120829448492087 0.55976911270101615 0.0094955995335043002 -0.76914570787638448 -0.14196517389439242 0.53894021494598998 0.52496023888877741 0.42930969564617871 0.019439571355106024 -0.22806463299242696 -1.1081529767631322;-0.47069507417836615 -0.56081930726440155 0.099841934966136414 -0.66169764950445364 0.69899263536066225 0.67597417501447887 1.2070145282099551 -0.32717280135604443 -1.2563459368314629 0.26179704556618411 0.21954942055480853 -0.29918326235335513 0.88552818975978831 -0.32785127344542919 -1.3464198007916504 -1.7693812281019561 -0.543367371216394 -0.48497692518569968 0.11185324354096886 1.3274104227650363;-0.16048214350515361 1.0318200982874401 -0.84445851318784715 -0.6549819235900628 -2.3175616580997511 -1.1321594771075565 -0.082040361578014753 -0.70666481355816035 0.12866465578753702 0.11369487893855872 -1.0656448355994477 0.029198931585470857 0.33993540329477417 1.0603394408904885 -1.7018234115626218 1.5633019063747502 -0.57947751566978567 0.46145398793554765 0.38275568720239422 -0.27787895363197102;-0.98886643050880829 -1.043582289106336 -1.8739544177587444 0.051790340400263757 0.52805300625608143 1.0667448386380518 -0.055496576885758653 -0.6982516764604294 -0.19003650257044538 -0.65428479537259887 -0.91171359623027348 -1.1653021788335125 0.13740775760132987 -2.1275450092202499 0.59037064925614691 0.49758282376076818 1.3540381115091651 -2.9310010856661242 0.28345218365228075 -1.0654343413049903;0.21878845974868763 0.19307780514985207 1.2524458397509004 0.61253892610064087 -0.63673877230941722 2.0230066777844722 -0.33319003830340344 0.41719698961370311 0.92362365501837551 0.57706293902786698 0.24835054727329969 0.30995380674581519 -0.24109328598283894 -0.23533913557044572 0.17689140092492248 -0.35960416267118206 -0.94196123147048383 -0.89483557348528275 -0.29334595647297024 2.0883571000565406];

% Layer 3
b3 = -0.32070072201360128;
LW3_2 = [0.10615094727818279 -0.043727640265031904 -0.029907953152895968 -0.060170496623470543 -0.061976401226798837 0.28436274529487365 0.067693914025582308 0.1120495184638518 0.26902989170802749 -0.030711454290900355 0.11550097077141468 0.052132280596475941 -0.080070168769117517 -0.93543960679902716 0.095314187260156918 0.30985525614269227 0.059856804907006352 -0.037362699753497 0.46496685742910687 1.7961951914894083];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00111018595614765;
y1_step1.xoffset = 1025.5;

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
