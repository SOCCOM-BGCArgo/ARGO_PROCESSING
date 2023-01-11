function [Y,Xf,Af] = ESPER_TA_8_Other_1(X,~,~)
%ESPER_TA_8_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:02.
% 
% [Y] = ESPER_TA_8_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228];
x1_step1.ymin = -1;

% Layer 1
b1 = [-8.3768729386792468716;22.926353300147109593;0.86697983712024306868;-1.2910212226010935321;-2.3851497286257288266;4.269350931086647094;-3.2486615552687667652;7.2503052820639775078;-3.0305578486451874021;-4.3635709478294657515;3.4413302861684700318;-6.6304934858510744533;-5.0300325058107437926;-5.9519305415699594874;7.4093586238153203993;-15.668327749131909954;1.1254992429970753776;-3.5599775625082608066;-1.3250881254155535416;-12.691979015088499239;1.729625335226537608;-0.14750856676734394601;-8.9072434101770081583;0.98946792053255405719;-0.7248200667613049486;2.0000364338398850528;-14.478486404742190174;9.4799777770788065112;6.9213444231005816221;4.2332466123821603787;2.8161881147584537644;-0.77844878558502905719;4.0370894701576656161;-5.3468575206158845958;4.5820495803931358125;-3.3519404942931161351;-7.1124300813923513331;3.2640876129766502878;-0.42263549434949715389;1.3795624985228516213];
IW1_1 = [-9.6988632616001471121 -0.29019642162073311376 -10.941819074360383368 -0.53025813280881406442 0.99611972943605675557 1.2008190606030180181;-1.3349680080317638531 -2.3339585952690158344 -3.5818839417323076191 -0.084456079258105304364 -28.896580415788331919 13.751378480225378453;-1.2720668798225494012 -1.3699010023520534762 -0.026575673296323883144 0.29016448358919072481 0.17611696168105653304 0.16301766108385518095;-2.9981065672529250143 -1.9245106353553023126 -2.0497597524951056691 -2.6562819805994108968 -4.0027413209009745643 -2.4234614877792983734;0.43869776445971653267 0.14332210300386263868 -4.9260773944522036061 0.33012716753926307778 -1.2884127919670151918 0.67262157342364525814;-0.44309500152672703921 0.0071550294470530736296 2.8049287864197793851 0.51048338030922757458 -3.8454167208900496 2.3020111149085300539;-0.20973188103436207763 0.058438603465101865431 0.42811209847976317011 -0.99601008023882042774 0.50271072229845148449 -2.3851622962548177931;3.0029072088100940974 0.0055997957581827408416 -2.809257226103424987 2.5556376000077758448 -13.434666416534948752 3.8623258723144093807;0.75807069180820341181 1.8433487578916145644 -0.83652788983701853454 -1.7832201061476817827 0.017609888193129519945 -2.9770862786519884935;11.587544259607248875 -5.9825406114124701773 -2.6829488714825635576 -13.56441234572385568 -6.5951913780142197297 -8.0942937128017824477;-2.5285007038640525678 1.5000663685376605105 3.8167165526950008747 2.802122177573181272 -6.9854024561004939287 -1.9470764446583568308;-3.5873521045435858845 0.017141958566020513988 -2.1905358294929784257 2.9186542447251264676 9.7588536722944674295 0.28207666783811174982;5.0498021166932298698 -10.755992799802564619 -2.4672457630787065774 4.0123330035729845022 -6.5417423057330426062 -14.022578159586826629;6.4375941849193809574 -4.9722907198523298433 -10.814628591613384856 -4.7137280311704445879 18.00264833589880098 -1.1803334473608595445;-2.0752346657813625086 24.995554073442420417 1.9271152541519209755 1.101246380546059056 -7.0545047582668241049 13.48542190270404717;-5.4556053271945526006 -1.5199462138427513036 -25.800815196943823082 -27.323558480226171952 -28.516805822949955029 -8.6803129126604741117;-1.2640338261917267726 0.18288568940934710194 -0.42825248173182495215 1.8687203034257870993 -0.17523179393898949407 1.3016562017918678063;-0.81996759961921039661 0.054301259880861034446 -1.1130063534897012278 -6.5349895705738063612 -11.895278353798651949 -5.4706966915337602586;-0.65526388373328514625 0.67497885952375746843 0.55966101749076058169 -1.1326774360679994214 1.1941753216822565964 -0.87328536506341958301;-1.2734901638800870671 -0.56088417429298365047 -6.5721219361959617444 -1.7479945230412881063 39.542715593411571717 -7.4503562230641611208;-6.1756025013735333928 -12.755360267237810135 -0.59429028807637518383 1.7267215444919625789 18.516096132773839145 7.347425642133015522;-0.10670716759353050807 -0.096178995184756796433 -0.15514949765873037313 1.2357430695549105337 1.1305222286399325693 -0.065071138795793154186;3.1318558480657681109 3.2787697162219222591 -18.024424904695521832 -0.31592394888868557956 10.99226180071643455 -3.5764890961245394685;0.89965408303318217254 -1.3964380528763575828 -0.58936546258527255926 1.5813091666154484027 -4.4917417259844665267 -0.46242478225632943634;-0.42447181034434161573 0.21978344330561638831 1.4240662411722848635 -1.1889463065614307435 -2.6982345233081841407 -0.024124586307482698971;7.3352447964714730944 -1.2366961456461629254 2.6732367485329189272 2.1418715748154544265 0.94854024389037772114 5.2210730449207609638;-12.512604611968205504 -6.5224785969665086682 -4.2513212604304326803 -2.7566851591428975787 -3.5613159692431208825 4.2795335007786610149;-2.2682541073115820573 -1.995463598864432031 -2.5817029987095612142 1.0966651623368648405 -6.7234133745880670929 9.7629251038857596257;-0.0012443358251978224593 0.27600792263108903102 0.17762622460554694515 -1.9172752172351934519 -12.7485001100325821 4.1439419562715666956;4.0244780009062877113 -1.5384135426172287353 0.78684918498054701175 1.7772900659975239268 -5.4520205994922443793 8.6356494910443988289;0.37856565083568338714 -0.46029166086168604721 -0.60265830750805049476 2.5129996229499393934 -1.0552326653153381031 0.38387443292352846358;2.034457579696013596 0.5911500330185094354 -1.7326987004399647674 0.090027956193819799746 4.6283068218891214585 -3.0224899588609877021;-0.13514586076422163696 0.089373242923655532577 0.22832014806489997483 -2.0163311958157472681 3.5457695343873942839 6.9555943849911017196;-1.8887296925333230124 -5.1454728917230321272 -3.372913115738265688 1.5897539826635616045 7.7849361594383168494 -2.0694103740498484001;0.048449206816557198463 -0.017549662511327411563 -0.3798170205918384168 0.25274380250821310012 -0.23845931300563469324 3.6852619059432107917;-3.1267661001760420625 1.9391245125134803917 -1.000852831075216276 -5.0024246755064698533 -5.2964559834485767098 -1.0932640307186807238;0.012212819490794887128 0.028291735479286903715 0.48796850148073545039 -0.31660852137419992491 -0.083982512538077119801 -6.1689712464962092042;0.10002647864019413582 -0.016020702202087695049 -0.3273986200947590075 0.38566682477498476223 -0.37927525995683369908 2.5294123538678419116;-0.11571832566999107339 0.24038555222079460716 0.59121311191015901709 -2.4579636492381293777 1.6005173435812141136 0.97266470313667130476;0.21389695451692306016 -2.7909658687960203771 -12.349809222074062731 -0.76722351642064678945 -9.1720560131925239489 0.67941987479216814361];

% Layer 2
b2 = -3.1383589792050434752;
LW2_1 = [-0.024386114566448360841 -0.36525416620447453964 0.117894925599126843 -0.02783250646670354031 0.1281445630968434668 0.14951539428525309483 -4.4274283156557805441 0.025577095341698580233 0.055200557638103692204 -0.012043168034940118435 0.028411828931692494071 -0.0728523298423688348 0.0069771618428416575311 -0.01270199351884527747 -0.0087345669518477126597 0.0022622740917601349628 -0.18402361478840231768 0.021810968989954420544 -0.37087265080925285732 -0.02508085464569586101 0.007997497400101060902 1.7536915093316802761 0.01552330847632420012 -0.098909566705951371124 -0.23857627672053796641 0.025035870745586122355 0.025396942668526507197 -0.046508622774052953142 0.61598919064208990104 0.013210278367684057457 -0.70853501245886174598 -0.067847559592971257159 -0.25376861203258799238 -0.021613605276844395642 27.158472764901386398 0.024060232985619146168 8.7512420426529438089 -19.176371481381430328 0.76701599046436563789 -0.019959162945237619857];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00492610837438424;
y1_step1.xoffset = 2075.6;

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
