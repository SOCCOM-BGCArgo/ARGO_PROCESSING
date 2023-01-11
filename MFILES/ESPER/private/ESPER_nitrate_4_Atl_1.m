function [Y,Xf,Af] = ESPER_nitrate_4_Atl_1(X,~,~)
%ESPER_NITRATE_4_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:33.
% 
% [Y] = ESPER_nitrate_4_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.94465199867687765778;3.4366736161463076193;11.222734890380687389;22.394780797881669088;-9.51538697473547046;-5.847884962858957536;7.4221421678491692475;2.7983644159184355438;2.1198024066394367537;1.8175232537638990049;-2.5067889502461930817;7.3356653850747104428;5.8482626906330219541;9.7663332473217874252;-1.225799015119548363;0.1901173615379617432;-4.3751769506534845533;1.8120186702593767958;22.489536379349043926;-4.5442451114092392572;-6.707799515701116988;-1.0584128895425388706;6.7286314890424412383;-1.0826082108401067217;-4.9330541175577939939;1.0683297448890394232;1.6090914621741458923;-5.0883913816062751323;-0.62659544272166378676;5.3206628750836673092;-5.7356095848854788954;-6.8368650275224229773;3.709869638551431148;-1.4980515611849765101;2.82562651501068407;-1.9108553886609218608;5.3649010984295193438;-3.7423032301082796991;-1.5690697610197268386;13.78912618372926957];
IW1_1 = [-0.89580490745452490131 -0.44256236739789267531 -1.4592725008644402696 0.50394561628219269878 -3.0347848381731816758 1.0609218238886288788 -2.7681296537615236453;0.48244592481557063701 -1.1088794659571155066 4.0369933083123479278 3.866002302341780883 -0.92291208710116079406 0.52253334329826384153 1.7803293873047802975;-0.4540665391139214746 -0.12856029555767303774 0.25754837596632207397 7.4177745808054993404 -0.96202423355177879571 -0.14832885796232370135 3.0167351915341638069;-2.9503410172071844464 8.6090253721503238182 -0.14565125546466251016 -4.2052483372361502134 -12.820479187082700179 1.4791176657004245776 7.0957972310279542327;0.27992470364474303368 0.19482440659759947543 -0.23253042682811150832 -6.2476090928805945168 0.64355573653486375374 0.33515763800481479517 -2.6750483784332637072;0.028849159070205417088 0.35245761351079518686 0.87701251305817928827 -4.6643608159338976193 -1.6574948714349673207 -0.32328728998986078214 -0.98688985361205172619;-1.2316847040840597227 -0.18495728216321549731 1.5312640819236751355 -0.6833628014543432494 -1.4020785902040702187 1.2169689047194522935 2.7544356680685369021;1.1452993475828001646 -0.25113452424085530579 0.13581572829944948966 -0.083255101127600261712 -4.9183888838473057703 0.94176512540042067734 0.44300078070463205782;0.59480902772879773455 -0.31784828305440215868 -1.9877374864967924939 -0.60219074023592578104 -1.1986089636170920159 2.5304557363905417766 1.1983260514260507534;1.3410765065502059379 0.70543184538751912971 -3.4083947444789344239 0.31822341902260953272 -4.5953534673779730113 -0.40394919896456493724 -3.5527560620653200196;-2.1656844569536524148 -4.4965836191408241618 0.30624907557278913117 3.9061120498957304648 -1.538392099728094875 2.5859781361839839775 -6.3981673238918137514;0.33871748550417374313 -0.31883660825707083042 -0.092380775848086293944 -1.6120824580043024987 -12.899032020038170288 3.9914793892166109401 -4.0660004485967293775;-1.3251895770957591303 0.53728018841248892379 -0.37913112597183895591 0.24139007388080269134 -4.4019028648336213649 3.3067420744885809292 -1.399302560627116998;-0.31711005935947494949 -0.15220317636013067086 0.18231606076035319908 4.9120178250083936788 -11.349405739522230974 -0.9067678234727144515 0.64702695800788445624;0.22002255499895123636 0.67103125309157041212 0.17838014714563529828 0.43818772477305079072 -2.286829096921457527 1.242049591813861209 -3.9612352587550749838;-0.071764343649842438078 0.35240147477305161505 -2.2159302259811743951 -1.2470306241105757472 -0.71525762435384321325 2.4347853372771255387 -1.1985949807019151692;1.2867256752760447558 -0.86302366553574727526 -0.64202879490595043954 -0.25897424143845676481 7.2905084397123394169 1.7133248343799467595 -0.11211244614305411471;-1.8792771243816626114 1.4134464547865639439 -3.4640225738759582441 -0.48501400734213218069 0.7859052924880003399 -4.7897266190293787957 -1.4115818845797918346;-0.088321365942530774551 0.14584849783136544232 -6.5220422909194057937 8.9682647762351450638 -2.4437385652242080702 -2.7794413061895499872 9.7833365583549269928;-5.299181424995929035 -2.7807131503051056853 2.6546878287749953174 1.1577703391070632843 1.8772401065603026016 -2.4648057932377382961 -0.19463201883811737103;0.41304267047905085519 0.031942981861992368608 -1.534155214094516495 -8.1716924052165271775 0.12049709946527013826 -1.4213011627613252319 -0.72837328960565383973;0.21606949262125321742 1.0394022015972628203 2.4763529306041784928 -0.35789668259827783769 5.5015735229742093182 -0.41405917687299859109 1.7494692896804611149;0.67797071891556359713 -0.53375006408530245583 -0.2400567430830588822 -2.0491049228834841145 -12.430858291125860049 0.63490431909452549331 -0.14731843886496845397;-0.5906827013257348602 0.94984262649929973588 0.50054472892213475177 1.312810731591815383 1.5828498873742598008 -1.9634449998289253969 -0.2362988719246548186;-2.1638485358900099165 -0.52813076736974717829 0.24058034670550645928 0.4304252355047530898 5.6780973492650002044 -3.1290417905157421607 -0.17888368971706627986;-0.1455337167839792778 -0.23210006933472165924 -1.4995096827835263653 -0.22637470820648031289 -0.7854329081306204996 0.87997606799945726763 1.2654415199062034603;0.41690650389422956312 0.06804718352099940748 1.2007750636662604116 -0.19734473555237957987 0.61251421497716385112 -1.0945509337254348825 1.7111317937618386598;-7.1074389100748973647 0.24943033346126541883 1.7080111530890356164 -0.63924333822197909605 2.5191458100512216234 -2.9193293342667612222 -3.7501239941612145579;-0.43186735627640249202 0.01642882753745297103 -1.6856033212680503475 -0.18164692446924773916 -1.0448271799056114695 -1.8303232989984632439 1.111070564870337396;1.1363379158861242058 -2.2759017166141757293 8.4470721292214889786 6.9809367423478949632 8.3684602893008221969 3.9984323024069405861 5.3028801031640337982;-0.3973944522229787002 0.33522759146505831662 0.013914954194315563851 1.2216687953531326638 11.809937403798482336 -0.90818469555236791546 1.6320801039005172761;0.072024846734714723095 0.29815179294193849824 -0.21075423134916701495 -4.7282731891225404652 0.21198563725822369697 0.6006226770571482243 -1.7081199263084689832;0.28735825157613859027 1.359331243686948909 2.3479167362391151208 0.018233110913409220077 2.2114837376133289482 2.6979506805182862372 1.337542064101467032;-0.053991226635889788543 0.67490106709131936213 -0.57899972604902028017 0.58764073214461443051 -4.7339333589996703822 -0.84805032488673781632 -3.8655852427889612066;0.24969616011502543396 -0.68473025305564372278 0.28809439030211353439 4.6996149604779962061 0.83129473373204343378 -0.81127582392936692379 -1.3386507075553408797;1.8274198691949288254 -1.3771485224565740424 3.3196401605424696868 0.12951128687842120013 -0.75010433511664664774 4.5882608053741558152 1.5443056650926185736;-3.3230347951847432064 -7.382885570614241999 -0.039873728218083395736 -9.7865926064471064905 -14.971145603525686241 3.6284295264761672328 6.0355155583363266913;-0.22582992844089114226 0.67692822519282558336 -0.20665900804943304614 -4.3454660820238970587 -0.61368729979090264059 -0.1648840218881245101 1.0526312778464683717;0.86619726506964800361 0.1104198441487197202 -0.87640805574873126282 -0.27878544189040316281 3.2529585195216501603 0.9365144829366840451 -0.33458565490780606266;0.76477581920198844934 3.8875506686278429314 -10.712897287954458392 15.851915425023992867 10.729605447383498529 -1.5117897090017378137 -3.30203416373459957];

% Layer 2
b2 = 7.9130869985149772106;
LW2_1 = [1.2966936421116701528 0.13882933677628678248 -5.2249169553847423231 -0.040651526038373346494 -12.401233183642512259 -1.3714842844056436277 -6.6479206325719939485 0.5565575400149224139 -0.40875043197120575211 -0.10800573730690130048 -0.045008758434839250817 0.19459003439996466756 0.19615051621306800467 0.073220892806234949557 -0.43919419688982352934 0.37201346684709324375 -0.37191754302766028006 -5.2949610827065090035 -0.14214871269140016974 0.026485824691673980069 -0.063002416347621811865 0.54412268583854361825 0.66561132539234957317 0.2428311519599031798 0.13457134428194267306 1.8032535797668309829 2.6577162715801452286 -0.029357225816319609812 -1.1567542735737181836 0.094718084243613109852 0.815868500744352243 10.415513659773273147 -0.58786468009823422687 0.52046105248302010349 2.5234499338890756803 -5.6256431757928178072 0.019975109286886137683 3.0928644545461616566 0.60813738048134347025 0.051114766096800400008];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0470499670650231;
y1_step1.xoffset = -0.28;

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
