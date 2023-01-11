function [Y,Xf,Af] = ESPER_oxygen_13_Other_4(X,~,~)
%ESPER_OXYGEN_13_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:48.
% 
% [Y] = ESPER_oxygen_13_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-0.04979;-0.9];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.568038194888224;0.0415843642790311];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.83083102017856114152;0.061441270795254590564;-6.0055804720210410608;-2.4242522724618456031;1.5622700825614586062;-0.38359530279017622911;-0.45273743297303953304;2.149097934466615456;2.6092989196412932351;0.91874399933700778842;3.4479233837559615239;0.94788444478494882084;1.1179498192611512231;7.7579555871343277573;-1.0374830290088004059;2.8497334454050502472;-0.98203741732533866848;0.78230829260225742683;0.35634733235279997254;3.3750331802563082739;1.4572685804698808365;1.895081466026877548;9.4027319598917369348;-7.9807684267875451667;-0.81750997424071325526;0.88197930053970363673;-10.458717150775283855;1.7048659863635986511;-0.30880171447222981573;-9.2616188487801753837];
IW1_1 = [0.05004298371995909156 0.15370911686523880402 0.39043110189272328681 -0.34597218570665660975 0.30585325581704675413 0.59377906560469184694 0.12196454133235666573;0.35921917781697987326 -0.345034846415124552 -0.8962057794480585704 0.84187619303606520127 3.6646468144233401709 -0.18423800843836776742 0.72738471597508669308;0.1639486294864145044 0.53241183127648594731 -3.8807054684913468989 -3.3603948066793782878 -1.8257248843738094202 0.3107354764966123728 -1.5415858565964437599;0.025965277453100881616 0.3325066139666843501 -1.1520027966203980618 -1.6397951977249809463 5.0610521716985052976 0.1734574294392191296 -0.97209247496269335365;-0.3372743818789268877 -0.67980318967220898507 -1.4970620257494917649 0.44929355274322330116 0.87592212625743448573 -0.86990436953926120189 -0.10313811658310069463;0.30250331264838659084 -0.47011290794799570625 -0.096412913590119275753 0.58690697553744697323 4.2254983366366074549 0.72197244757894507572 1.9949532636298881183;-0.18600960858392190378 -0.36473826235505985016 -0.83896889158699627931 0.7017723220212858104 -4.8864185806547819624 -0.13123287208172906837 0.013314969883573200132;-0.27461144214872956759 -0.52725706932962512763 -0.23338714400058446774 1.8325352668010879142 -3.6581816585388109431 0.38179091868348741245 -0.38759278852211631605;0.51225159382536589181 -0.040738339238017950983 -0.85277451333544695977 0.65604099359063305918 0.47857305965822899418 -0.40278641466072062638 -0.54030670802198865754;0.44193403329002756141 -0.60762503483473606636 -0.5226169402027596167 0.76312642628793703192 -1.3225979110421879348 -0.51683262880852665067 0.15924199457429422111;-0.28775604888869626974 -0.69549955922953687892 -0.33853716746488410339 4.1621243095727864514 0.099447709029091058142 0.15899745929932324051 0.60367110332530538042;-0.65933145760806333158 -1.561972524030419418 -0.12622965377804051812 0.15547797467483182032 0.29595189558302237298 -0.55622443781657360251 0.55486075062603079644;-11.546329424546444287 3.8828623687738508252 2.6585733962866173741 -0.98277658441978199555 -10.32642304525275101 1.8377485827718205158 13.346117675240011025;0.44725583082439479554 0.55261747437557973139 -1.741814096949335644 8.8568302169683477132 -2.7172217764535373519 -2.2029771932101516185 2.5930726893175082282;-0.30129840712368943789 0.86605187694524188302 -0.38399273577544995595 -0.31747233820754544587 2.3961138992331387776 -0.028624157026378729829 -1.1367389842493282881;2.6520669045686844356 2.6691081539474188489 0.86542141699607300076 0.46164666380229529841 5.3129313488642031515 1.7830712392737877448 0.62467170710895059127;0.76420643062138016166 0.79675739127132039208 0.68117814916384056279 0.17186440234953043182 -1.5027565656957442286 0.10311578459810674668 -0.93671336256838699796;0.97072471184264708111 0.25481954190052030773 0.39646438717291748732 0.52151964409874396189 1.2753881472323991275 0.38663218343047556225 1.1595423461560618872;-0.12082079193326841771 -0.4047225302843751038 0.065233943003832506324 0.27575910791813895706 -2.0370746858031041349 0.035224718812053475159 0.22945126126488121776;-0.32015952355909721438 -0.017571627465297591508 0.65893970148333413928 2.9208677687580952842 6.6532275551007789005 1.7336333895833708763 1.1485515573580502569;1.6477446204271883889 0.48037474085708786742 -1.8631645587541072473 0.68991358783987843939 1.6710120586024637124 -1.9059424940849001295 0.26628304740459479705;0.38148613501389311287 -0.038962325855646552653 -0.87425834483111108941 0.72671334219022709622 1.4830026905109703339 -0.39042949020238681967 0.28820130768941304433;-0.64176006782994632349 -1.1299951384404458565 -1.1002390221763678202 6.7733622483789925894 0.44982355288348091493 0.52090675590280355323 1.2725063734645623015;0.5676090706208302139 1.5295613202990323565 -2.9394361863346962771 -5.7469685697469179431 -5.0478607522132259433 -0.6362612801486596803 -0.53411980833687511261;-1.097966896698487016 -0.4200163318557345038 -2.7869563457521246974 0.1078760310524727456 -1.3635759048782079006 -0.76966969116746730695 -0.27259998250177258639;0.21485454463357317456 -0.73779310306945100617 1.2110125761509507036 0.33471853121797168207 -0.56473736854396860974 0.60688126664538100119 -0.13279513577609280972;0.68989224178517594144 -0.20492218440522586009 -12.049991578912946366 -0.34863725010369889246 -6.054779961374409325 -0.71021780982566817819 1.1394326986263436652;0.20631311017274639452 -0.021195767719149286035 -1.0454631394994895466 1.0250728815477345712 1.9914751050443273428 -0.40785181496305278337 0.66717569233120699312;-0.015129194300445411928 0.23313128961375023662 -0.50724611492039051175 -0.1870853277031930828 -1.2771454183940027516 -0.2126823321493452601 -0.41256267571578286724;-0.22095986288997551683 0.21158352081940648426 3.9953875495047572741 -8.501626406092237076 5.9732756357728176155 -0.60839912452392708175 0.32594649363516292606];

% Layer 2
b2 = [-15.129126903200521781;4.1684584661290919883;-2.5613707111408960948;-7.0483018544144124817;8.9818637867510062733;-3.3378421627781595937;-3.4795180148647082419;-3.0025852438927240762;-9.4241347523217395832;4.086338413982026907];
LW2_1 = [-18.669494373267369269 -9.6044472939434708536 0.16122611153794719208 8.3875496254858497025 -26.834831583166625535 -17.394134072542811964 0.26761602783813787099 2.3995195073180433099 29.749867985627354017 25.282062226807113348 5.8164929017482744911 12.448669583180672404 -0.11411546166452286866 3.7551093115825717739 6.7095411468677008671 1.2807743396350439458 6.1568023093377428623 4.8108246184797378575 -19.4384218050906874 7.0449904297881937865 -0.87341909462742073345 9.2222343797078494276 8.8049846874088384396 -7.3337630002332536705 -3.6805832233286248645 -9.3273708321167312363 16.783463209022823293 -29.981897706611015053 6.3737929278150335932 6.0628561005403946993;0.013445712902673455369 -2.3501442626757564902 -0.20708765771964443569 0.13809399712268261018 0.14348180919105366082 0.21975464884575418578 1.7447773829290729175 3.6314457753803841555 0.55590805286473565783 4.6664052296270437381 -2.0729008585036625689 -0.84250793795196132763 -0.018354494287752771975 1.3404896841198734592 -0.93564620902593065299 0.99718552627204359862 -3.431986206271436135 0.19404350296463485082 -12.983474207640643527 -0.68035513108278000693 -0.78172531557843527583 -10.95600374417527334 -3.6145623259992096798 -0.6415257526159819923 1.2654655824648179951 2.5726386467207853492 0.65702559322256104402 6.1416987776286946499 -1.6641252230661698341 0.48930716284963809981;2.9090442928609703266 0.53046530737377428633 0.54110674098020883349 -1.0573335260772023236 -0.7844767383938439087 0.45376855128324283806 -0.36212281247983152399 -0.8664741302558306435 7.1895875429794289246 -0.38696361273983465034 0.80504126626032213032 2.1073729905117617456 0.029924765789090784773 -0.9314292983260777925 -2.5663655161644052249 -0.34923004937436541439 2.2116457355162322251 1.491592016196913173 -3.04718461168704291 -0.78138358774900373849 0.40725518657983889037 -5.769604529771501511 0.36736160198244177666 0.16243305921463901953 1.387015702524370786 -2.4532181599683631212 -1.5688159678905044814 -1.7000042514130992988 -1.8330996526116438705 -0.2896764844184767651;6.8409582474622157022 3.4536149751940641295 -3.6483948265655272714 1.6693401042363993358 -3.86770946170353902 0.4636149916345398414 -3.9503245026750857249 -6.6017072316138252219 11.10681723983426572 -1.5660250517921681102 1.6626480447604297108 -0.17713070853570514473 -0.13059732547992941587 -1.2180112493619428271 1.7084393203851879761 -0.59946422966941137478 -0.22407489423857768207 0.71466005025441226017 17.209071430006613213 -1.4039077461473059572 -0.14206481849320792832 -11.75781714480320872 1.3280135664469892021 1.306680382616239644 -0.027782800620353487187 0.29753166265503683263 -0.61228153544131935782 4.4119539465565296155 10.360809744261478826 0.3007323781566834664;-4.8428038190714088884 0.30686450713732976725 0.26126748995812343335 -1.7497338975221872204 -1.2196065313599810942 -0.20493892931342380659 -3.058320435081767652 -2.3436919714414563565 3.7728356950538133319 0.83157138037702960087 1.0013036161294117665 0.25537989818279088761 -0.12519296453353265242 -0.037914102671618342411 2.9238942997984773697 -0.043959978357705639662 1.1285828711541177061 0.18352879248307191196 6.8282361783002523126 -0.19041492364355425138 0.14008407618538831763 -7.7425563151372847059 -4.0190173428890201279 0.1289939642328425673 0.061673959617211726592 0.36472718988725150169 -0.84791108442085627583 4.7840861789217381528 4.7950431880917063765 0.41654602375744687048;4.9385635959462836908 -0.9038257428416037742 -0.14227552339477606003 -1.7078961471530069716 -1.8788409163189974116 0.39011250158911714436 1.8586270652920251489 -0.24845208111363537018 -2.4427643549908157894 2.6334462275723815594 -0.012175621834886940151 0.59806992650307888582 -0.026525899999291346204 -0.69967402684888813091 1.1460446658145428955 0.11296254231304644899 -1.6963139910961961832 -1.1395916001534942552 -6.2129491177832170123 0.40504269654723512595 -0.27111235515767023729 3.4581872536494180359 0.7597982888843570759 0.64559844658762144398 -0.42324248946021447715 0.062490382919743861501 -0.26731646076627285646 -2.0866545090155215192 -1.9860093635737881002 -0.78423534543234285632;1.6824847568736429615 -1.6706102459069693822 -1.9706744748318254512 1.8333806464270974157 -0.064888624429243840774 1.7053002027517050188 0.17374441464640427712 3.6951723400708207379 1.3121951527342656796 2.1101321164966182486 -1.9302204052073612672 1.5974691924241513075 0.18725920507009646832 -0.44250173781250756733 -4.9306474620343720972 0.15827995317510873408 0.56184104086554809143 0.057905881150731225437 -8.3869661606593357561 0.3667618276944578648 0.69175016116281029532 1.7641504253685691417 -1.1449232879385713435 0.46545282824533062982 1.3663556089581965214 -3.3079098856506488069 -0.99792826333225193025 -3.8040780893855532518 -3.7890547986819749227 0.35466115036714601638;2.6379399995625276709 1.7033956337363664879 -0.53614160744093886102 -0.38416385112153805803 -0.92048884893416205166 -0.16551728431530876584 0.30527461366519842434 0.66761529042716405424 -1.8804107986568596722 0.21480514631621938637 -0.69495781052307137315 -0.71299288613796174552 0.0069986480449861603048 -0.87412432606405032764 1.1149760187651180932 -0.036260322004279560426 -2.0522171900059729133 -0.8878287308811011469 0.85404825818037399809 0.21883148185087997195 -1.0073135286405339173 4.7463962367442338319 0.41610992545831265721 0.54728175553949198751 -1.006022219558523112 0.021842227787928686356 0.60414956713651835063 -0.83937360225218304866 1.9259358826764318984 0.062291749669482553153;10.043043131448262884 6.0573206617669388407 3.7897797132925066244 -7.1953605195705145903 -0.58594703620113697262 -0.36684759125082577746 -5.3996743701650204983 -3.5160526392812077034 1.6609772335934918086 -1.6763021509584208335 -0.68407055920427417561 -4.9711537083042616203 -0.41038359816775710209 -2.3033022302141770155 12.143339251537192069 -0.037155683719844331825 -3.337305933218867704 -3.5924877121214633391 16.193064486333607732 -1.5636214418560190342 -0.65607920583732004527 -12.050840961436053433 3.5579836156517594681 -0.72957081281000513151 -0.99528120573060996357 6.2676754369300313741 -1.6186972330957334343 14.296154699077812111 7.6933956065660753509 0.51960318931454385627;-0.75265747622959866536 -4.6595857805924500994 -0.67246482635554860607 1.2940074929520337665 -1.4997034518814851278 1.8434235485973975432 4.6436518175421968024 3.857929221328118885 1.139581377176685395 2.9391545045270390801 -0.91863435545635052648 0.87011892352265729489 0.011526325071644622985 0.43957685039299493335 -0.59181223457968745727 0.11421169438224922088 -0.19327193039049492262 -1.8865906915709864666 -9.9739663788567227698 0.040068365023171875761 1.131189376571382299 -4.2369544814862738491 -0.12765372606047717396 0.24817287346758862276 -0.88363600901032413049 -0.63300777786430550975 1.4491521081680425187 1.2732993975564070865 -5.6066783079347048968 -0.25465208665642036223];

% Layer 3
b3 = -0.40505412756475533964;
LW3_2 = [0.15013608145115436843 0.17463814138868294368 0.23198167254158372219 -0.16087370609980006231 0.34017173624995300507 -0.33804204945800614412 -0.19378939816080484326 0.31987759545686272888 -0.051755361852240867881 0.20117992915551130562];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00454442172233583;
y1_step1.xoffset = 0.3;

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
