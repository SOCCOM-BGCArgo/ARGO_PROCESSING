function [Y,Xf,Af] = ESPER_silicate_12_Other_4(X,~,~)
%ESPER_SILICATE_12_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:42.
% 
% [Y] = ESPER_silicate_12_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-0.9];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0415843642790311];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.9757728591476184654;5.1627029941814477354;-0.62950624034240887372;3.5833163428901575465;-5.3849553126727354169;-1.9177743524318928703;0.19498567482482279112;1.8159410238923574532;1.0378005681154769402;0.79816679453274907985;13.31606169566476261;0.40483993939295948028;-0.13065372569267433289;0.15150649242020750029;-0.53261985694950964021;3.4560183295807971326;0.22678639610407755334;-0.37896105029101767592;0.29226819665599718112;1.1706098450699875535;-0.30695507512357628199;1.5416468274943864181;-2.1580252550777903053;-1.2136763978530240582;1.062666170601402893;2.4610137411310142319;-1.6636785923366230744;4.2980728782590782444;3.458034504554363231;-0.47507717807671823174];
IW1_1 = [3.2888755394408755173 -1.0385342944036641821 -2.134744614677981378 0.46002563877139535098 -0.7580010281448426035 0.029970178995301612057;-0.23259908571045262438 -0.15405705912762440724 2.5899786455746234282 2.0125393427945819091 2.4639114291707802096 -0.22738771941246249675;0.4557512775157870899 -0.090537459451230678509 0.49809198597804532627 0.018293817078328451775 1.6704491451959548876 -0.94333616852916735152;-2.489364919466554138 0.97769064505744995852 2.5527932304001268449 -0.75930456128036027952 1.3854165125146491455 0.10817469357494549631;0.52669696417109945763 0.66427411368258504343 -6.5896788571944293622 0.20719679721047762344 3.2704162135797223421 0.11440993961981824056;-0.77638625209319545029 -0.68615793735784280649 1.4937345777046724393 -1.2504622090367167875 3.9282199716899017794 -0.20906620181546969883;-0.49481764147092954786 -0.27551304644055479853 -2.1615430687467918069 0.83038167614818836348 1.6921398783761720797 -0.5055010850257027899;-0.20859163916187220544 1.212814539437525152 1.691346483459757799 1.2063284812476544605 3.5051969080482505703 0.72692127774001558116;-3.0619524058292468816 -2.0123012614027477873 -1.7486878279914168655 -0.39180392293596677522 -3.8783558925600494227 3.0712292956469320337;1.6771920484669300411 1.1737647121608301148 0.68689126643388276339 -0.58824864365739693017 2.7806618639365137113 0.93269649006068322539;-9.639925851752389363 10.768097721758490337 -31.58756735018991435 0.10101467873072617065 55.765574005026792292 1.787773470786868435;-0.97209344568629729721 0.5596655978038996837 -0.29227981277888592793 0.12353033038509264596 -1.7879102215875670723 0.78472655338176489082;-0.26008402956409865103 -0.54764221857887418476 0.74368955631637023718 -0.74501578541999224914 -0.35656878651679924097 0.9083517713835256302;9.0601213678720675659 6.0678649633009369424 -11.992699607730207489 3.5740876944888082534 3.3943475349971121702 -0.98101034228686834204;-0.20988092551444445522 -0.42587161899335668203 -0.46997916402773937117 -0.58936706083724843541 -1.6957514671741906742 0.12666521660289675744;1.7842345071621679153 0.91037459142132204359 1.0228501187953453755 1.3142851829131936281 4.9920921633039583298 1.4475712151817501372;0.22541659244966158471 0.43559757367826690277 -0.72985875870847716218 0.93477652416879530861 -0.087292810949198301751 -0.94785608754395989095;-0.072477070982655766862 -0.046144086382609750263 0.042696019908988726288 0.82558849817163237894 4.5923924868283831913 1.1459243943207475525;-1.7662261375547567699 -0.18100125344203979449 1.3526804373374110568 1.8358422607723619002 0.0075781060537447528236 -0.079476766342086696548;0.37630395180183201509 -4.2890637464297860504 2.8827661776881097211 0.59443278850979663908 -1.1867253684204359665 0.34797641301310938688;3.4651150819735003772 -0.0056227507445575017342 0.94010237700569609931 -1.3971064683363281222 2.0754475377179293361 -1.6353188238466593862;0.16766551447859673329 -0.048002227172051319115 -3.2013255090808905301 1.492916394320680018 -1.1489507365515005954 0.2285905150004384212;-1.1008820586897500871 -0.43559943930208167018 -2.1911907064800382194 0.40733706698617039299 -0.53417620470750937223 -0.19084445280733994266;-0.82881496324479075888 -0.058657765400162827085 0.2762963889079189661 0.8236683829521683764 0.16813076213431724093 -0.66710859369660602258;0.17506659330545720765 0.18611990300679698285 0.23097600059503770042 0.35356187028288449392 -0.31917345780756917195 -0.093397036744946898557;1.5180774137655446143 -1.6288042121154755471 1.956353168550852395 -4.4540600969349135596 14.692007316677358375 -3.4684315119883932788;-2.0976976348459923649 -1.3660568328466862198 0.46844398449280449936 0.67837838906208292045 1.9412792242766134798 0.1883500461808966242;2.5349714540821657316 -1.0129399140190851192 -3.3479483083067598237 0.77161881414931532408 -2.970541600108345115 0.66490297129296438428;-0.40332100220380362465 0.16182237822116299353 6.884565388546480591 -1.7628153520715337255 -35.66710786618983775 0.59337419969004245868;-0.066279794563596741708 -0.082316429071477820179 2.3783302613286965155 0.04615362279689651237 -2.648028176662526878 -0.2912184149752933271];

% Layer 2
b2 = [3.3101485196261775101;29.146589815693040748;1.7619937469075903635;-9.1285177064811637848;-7.0030769721145009754;-1.7358190702672433581;-8.2541236243116351545;18.907342985630787524;-7.0425203598534755756;24.595570512573406319];
LW2_1 = [9.7485585697100614055 2.9717731710674004653 15.61636515268511971 10.964572075978525234 -1.0289555782477251 4.5720339842640340677 0.58745387472823873409 2.6828309291366814193 -0.69360516554183782034 -2.1346369120888959436 2.3210454472936081771 6.694957173133121664 1.8212718859145784478 -0.80246717318110805106 11.174315100968495784 -0.11363429512374056962 1.2795389229077300808 3.0653419934435719441 0.1447526708939651241 1.0460221465276582364 -0.1586325251266701708 -2.6500815830676742735 -7.0874555145204372053 7.7618165674833781154 7.8346534047319815741 0.15930344392626807903 -4.1865640844699454348 6.5378317256404931967 -2.0931422228029141763 1.2446845374879624835;-15.718711621853262983 -1.4777700795815493784 5.6404415674408392434 -12.15907359419006184 -0.46156890203597733535 -5.0075341364775001196 -2.5024268874574220334 1.4914114910904843203 2.6112238868709782835 -1.7535831817696176937 -0.73868834777935754587 -5.3117065042700497557 8.2373323772342281757 -0.04504804461578815894 -0.51501864080038073368 -0.32140870109966213075 11.706467557579719241 4.0656878430014442571 0.23713789583056588128 0.43016992226943950861 2.5490059706147634166 -3.2845871804274695549 -6.1726260104848957866 -0.1398011254895855926 -35.965769823239405412 -0.070454408289538708599 -1.9418130582153325658 -1.1561731681953504314 -0.14134824648608024233 -3.1897748274682107272;-1.635424625316318048 2.6182821643904028441 1.2310884747122650129 -1.6871932232320878686 -0.12015194569054656049 -0.86741242256146677825 -0.18911270083076131554 0.085261950448362017196 0.079435278110778398397 0.52261189923943995872 0.0010213159748343552356 0.16341798856774591986 -3.4329365891683414347 -0.19170306362251224308 3.012440234967305841 -0.7643356491108693751 -2.1167490986953692023 -0.33382716751472502059 -0.02396879415984108988 -0.12100242131226861264 0.47490767032831238126 -2.0410244532036085907 -0.653102747117262461 3.8451723614905346516 0.19530581271096769669 -0.33073450710306073752 0.43231602368497840105 -0.056130676958806018095 -0.054187514728309453538 -1.4322335734227218396;-1.5439541302187329563 -1.6397805234995170842 -11.079788551237024663 -1.9144884372201227407 -3.2998807805378804758 0.47703963245888869604 1.4138408249016725104 0.39288401029109909501 1.1499719422269856217 0.65333427669905386193 0.086848075253725995348 -1.9413937090783026651 -21.903377360522668482 -0.64433360622538360385 5.358009748077646961 0.7820729369269401321 -22.592375191359479913 3.6922895996128546336 0.16882397686811287385 0.039883548312067210184 0.69692779282141625874 -0.42751466244643221426 1.3843284499017609601 5.975619832220620431 9.1166266128170825311 0.8618689255575430197 -1.1107199076317249009 0.53900732359334257993 -2.0096576830755488174 1.4325155457125131608;3.9640673162055906076 -3.5557431130049215184 -0.56871476252568120557 5.6841729888788981029 1.9514289047948192302 0.87617853713784532577 0.82413900325051259177 0.31905619100187643244 0.37236780120280738027 0.3179856676491131795 0.074799826182081696757 0.10077356629339367933 0.43739099963850225361 0.066800293180869288467 -0.55486210299199645313 0.067368456728272366196 -1.6093280269535308324 1.5752617954264700462 0.46991416688787746025 0.46939105607687892086 0.27140556578028834922 0.39958509600680852669 -1.0309231485702949893 -2.1506475245100395988 12.349878605618766159 0.38795442143518743539 -0.26086912728784256155 0.62158640348627325345 0.1268650300229556549 -0.22767510884574543306;1.6191064320611230354 -2.6095403318659284864 -1.3138614345009285156 1.691920353720704151 0.11609065257195651244 0.85433914163261426289 0.22241771370481075443 -0.090495273679374599318 -0.067628308464232025909 -0.51915350686212569098 -0.002057893215651635313 -0.22077141436834393784 3.1676547440834426084 0.19214338758601037793 -3.0018811124672510893 0.75380568784689849782 1.862378874158378661 0.3343282635791976487 0.039440397534962860771 0.11488370814214424187 -0.45633321887573530251 2.0131292554746069712 0.68570044090449078222 -3.7734890749056262393 -0.15680879262698776966 0.31967989665134521893 -0.41330987053088269478 0.068742411317516835534 0.049566681640497649297 1.4656060823440266727;0.92288463374638396441 -3.4166089001084656296 0.47406544312800180663 3.9446342069693933752 -5.4894372201459686877 0.10452717006084247253 2.7573580884399699897 0.48432948704449435562 1.1323792952927098643 1.2546357353381800692 0.14214173481291039747 -1.3028609895294491761 1.6604472584198106677 -0.057493707973429231117 -0.99331352666826544517 -0.64357466321405576437 2.4765559761760647106 -4.4580220587839578528 1.4991577587092461687 0.25397589955769583892 0.90446284862702652774 -0.12423149866688000531 0.63225443248196999679 1.6783617593250534927 5.135382123882513028 0.29532785555762292029 0.40627564499942281362 1.2633960780541144686 0.43617570514547776561 1.3980715399093244233;-20.260913381086993468 4.0512839568017033898 -1.9264544429280689997 -23.023381403712047444 -9.7922335802312758801 -1.1130312606600960557 1.1237172940069850569 -0.76049994753543215253 2.4945500552087742641 5.5801085428497732011 -1.2291941835003659111 4.6120368813911003514 7.4857779571425488996 2.353771573231322467 -13.224746774278376193 2.2036277700460642492 8.981285551695611602 -3.693515283057419829 -2.0490952963723043645 1.2493948328266464287 3.1836824205884672878 -3.9700523643663530038 17.368220409002020688 7.3572142164150466925 -39.563064298179781986 -0.7562396354565458223 -2.0979991780723854511 -1.4779557045424167594 1.4191194603927028695 3.0136949167405133743;-0.13195636842856123261 -2.1923692894244837248 -3.7402730891415107983 0.5164341364480923291 -0.19209028242308098799 1.4226973141254013466 0.90723206973802617004 -0.12925056233078721646 -0.12572341968759587205 1.4282566828340474974 -0.006893178922019175224 -1.6834895793576616363 -7.292091929650962534 -0.10274025292006798427 0.29358177993939726935 -0.04233290702520005333 -7.7046311368189206448 -1.4253286948777430521 1.4540526723028686096 0.14429550482647518495 0.29914925985617540904 -1.2400147191633481558 2.8359326565938065556 -3.998776316679113485 3.9821456310659031352 0.32283350619568201845 -0.12765123845584649964 1.0240466855152057235 -0.11564041214215327324 -1.0509326072157967324;28.983838680466941895 7.9385777090580065263 8.3745291143078741669 16.134952840689720688 7.9384621915362183486 -4.2716890967094176546 -4.725391918469645347 0.87660567577204306922 0.59783873469227877262 -5.1760851990002638701 -0.071970491548434892271 -0.84133713504986618048 -7.5626094005686050892 0.35448631208451536478 4.4892342130417679513 0.69765033804015663677 -1.4341261029287053486 3.8526313168475945226 -2.2086790266440208264 -0.73830877076407674942 -2.2128189122560120161 -2.1808540434379248829 -3.5841417028682207579 -1.9608711566057974718 -12.447373603482201787 0.24924317813703150692 -0.27499032920423593929 -2.6690692520356411599 1.7598712019538440554 -1.6066819069238693629];

% Layer 3
b3 = -0.99604670547752149901;
LW3_2 = [-0.26274133656473808163 0.48573955931348061776 9.2032867840929544201 0.15190196907569025742 0.11216465428650342073 9.5245299972498607133 -0.1101171934125904639 -0.1021622468278396284 -0.47075536496670389308 -0.084830552831265920721];

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
