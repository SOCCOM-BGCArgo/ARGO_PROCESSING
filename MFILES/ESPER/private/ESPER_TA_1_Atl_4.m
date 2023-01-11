function [Y,Xf,Af] = ESPER_TA_1_Atl_4(X,~,~)
%ESPER_TA_1_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:00.
% 
% [Y] = ESPER_TA_1_Atl_4(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 9xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-0.2178;-143.241523403244;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0479365703301392;0.00458552991639557;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.6825227848365664762;-0.36551837627289568422;1.9876825877657786634;-3.5388844220234374482;0.56110458808021823618;0.94969464561010730197;1.7276114502911841253;3.6730286942986163368;-1.1497726143630480067;-1.1796454263009128471;0.31025762403400525269;0.94311804862028647545;-0.75953106822299798484;0.80967891312841466789;-0.10479736505790754331;-0.22403450996450136046;-1.0418926016836536075;-1.2772818442484887846;-3.9185853208618492971;1.9109774103581333993;0.71097373210761616313;0.51085331223730279238;-1.2884390340921292317;6.3102297480891067494;-0.68667873668738943493;0.92988832625953088407;5.3058922646723924998;-0.20744305960465878047;-4.8956191498703969245;-0.44248791632002659746];
IW1_1 = [-0.4453933211043763496 1.7931308137490091781 0.23247939679775925237 1.8942274366279894249 -3.4393481108073253871 -3.1821908526517632687 0.55414019261429769969 -2.7724154922324357742 -1.4702855659012785683;-0.02640755307416955619 -2.680110526791660952 -0.7276796237335700912 1.0071885259252739697 -4.5752579489491047582 1.1208886126385033677 2.0555300044396851433 -1.4187207922899764956 0.66648215225143458973;0.541030050192465084 0.13792250153827498305 -0.46581083922056237823 0.86562769021059571273 1.9133283663108437977 0.63825108756810400301 0.52877696482075486895 0.86893095360091765489 -0.56030579082600084906;0.82846304980087515091 -0.99714472587595248676 0.13642411891106515087 0.15147213084036084685 0.8762894459669717806 -0.61614664335202939949 0.31103088162310210674 0.22888483843978033039 -0.65451173348679492214;-0.11362327549227738133 0.52058853403719029806 -0.035966452393207656857 -0.17473024762507349816 -1.0261562356361824122 0.10974634585672907172 0.10086956373407161835 -0.15845946773668895258 -0.063132845920218658153;1.5063480258915351584 0.066430112448582295803 0.14314664324275336305 0.49875504859694314019 -0.53378953440217880022 0.10048917767773432597 0.60097025936118730893 0.79408925102692673459 0.51324988449190378503;-0.97506623534014646992 0.31285274268211527016 -0.094080183361033492551 0.20585322284981127927 -2.5742319915250750917 0.36479033287335072 0.10448709977592325615 -0.098973597219608855369 0.19589081251637058267;-0.5548044428830268826 -0.77727971143263274723 -0.46452827664183465872 -1.3805506886598668537 0.3835001952213466514 -1.5574922904143149172 -0.25391951263333878419 -1.3804221430466507847 -0.2206688584225786065;1.2155456111090956473 -0.043924392374055071242 0.060884241760068208427 -0.6496998611736464202 1.3961023644410326394 -0.23534911922843768095 0.13414955419891763655 0.22927417304143016263 -1.03107557349425627;-0.90182399244534661165 -1.1821384784222694897 -0.91021538166208815213 -0.98497215004286087225 1.2341980045969751512 -1.3881991225074732377 0.079321898434737400652 0.51959159992381886983 -1.1351859252955331936;-3.2119403692216419444 1.2117245452821694407 0.67443537668886033298 -0.55402680448002783464 3.2781194549499037549 -0.2623912327123938959 0.074181198598157749302 -0.083449215215705752069 -0.20283009080702191218;-0.82478778437425737113 -0.14977795275104660577 -0.44660603196544568227 0.26430112253501514274 -0.62312268810965298904 0.81085726690558068697 -0.24414506244975220373 0.09749354318952790921 0.75150201543074945931;0.30243971456849466639 -1.7221954048145295957 0.66237716887865405457 -0.016632815040053528649 2.3912275740871935525 0.11430538702734058254 -1.0772096212061645915 -0.698814387633004519 0.87323227245432544219;0.98704042412490367031 0.24795286515695072471 -0.59467787677663119794 0.65035416985979122018 -1.7281301789773957367 0.26009886606331561465 0.81245191386372817011 0.046860556293339288014 0.35426228411158955378;-1.1177929550709877482 1.9882599736692934123 2.6833816826579552917 -0.82843895883896934151 -2.6914337682122502216 -0.33204867412898941836 0.20501029768487472049 -0.57019116749857134963 -0.71481572137679127987;-0.45154032765422380535 -0.12665710584687667195 -0.2275823338337122359 0.5245505718449376964 0.92310098493122094698 0.14160550187534431954 -0.31692920245503158894 -0.1210731890744639605 0.80742485991868984119;1.203388272163418371 -2.2701132186650849754 -1.880428933090915633 1.0167207807553166887 3.0914244939398534129 0.30269805238478031395 0.023048346381973290609 0.62676980040672414241 0.25728429115831008556;1.3825594565966399241 -0.64589565482443145328 0.26470591598135523315 -0.4440169188244906251 0.17078222911170160003 -0.71782920348570111724 -0.087668019425942655953 0.064855633306612411881 -0.63405048517287443399;-0.97908494599447704054 -0.7231261526857034605 1.7292706752909463574 -1.1143842867174988953 2.5207741765756654928 -3.0718011190336427241 -0.057111202357680378172 -0.45176540459887598589 -0.21995603703364857706;-0.71643754726600261673 -0.35329555435955656284 0.2392685744609757692 -0.37148278150918206952 -1.8230356356901833426 -0.73208581543402839742 0.12516760765080250928 0.034166689406896641468 -0.527936114579500404;1.2006258185616558443 -0.79404440222051431153 0.47277708572203086534 0.49430653239522082121 -1.8363285271207023808 -0.20832522088087193635 0.13296846666399561587 0.66982623652218375998 0.31458430282263066813;-0.9832654264357505447 0.81042357860536751168 -2.414241782849349871 -0.30219974814167799293 4.5303015373036474855 -0.025966452665100857311 0.39492854096346868831 -1.1422736157487318742 -1.2596588655171161442;-0.12946098627569613271 0.36928365626059717064 -0.65554075233890507501 0.34433322982524222633 1.0958920183768299239 0.80982459501107539879 0.045109604735090569338 0.23167060623867255509 0.2469517800985591427;1.6324950249343932462 -0.56149433702160911341 -2.8100693911300651529 1.0320975974569257527 -3.1616558994750518785 -3.0976238060393241902 -1.6056842960854105229 1.0758685467170374661 -1.45983678507066994;-0.43793673490589712927 -0.73712603220843075924 -0.32443284519921722353 -0.028039311915493844718 -0.1802638708193406436 0.32429697854046951599 -0.36427517411643128842 0.11524713134086211896 0.2192427151345823777;0.15450231584228266146 0.63530709934106655012 -0.42650037958418801631 0.18712331356621894618 -0.72276142306201829335 0.31514962562082154962 0.25828022190488436616 0.21946636724639104199 0.17020940654362232136;-0.68776308913616568397 -0.080955782944515025656 -0.82177290080594733723 0.050698386865448225036 2.7680115865328440705 0.66539527744848625535 1.0623104099368418929 0.2071163009323256643 2.9114012518630856263;0.39127764481728871671 0.05018604444026097261 -0.94527753837062700182 0.40779339700170080807 3.8382339475505542481 0.40186009418356410938 0.47005961795561257421 0.5592952190766954379 -0.13152488853600086127;-6.1184608226467451786 -0.36269179153541813987 -0.8317637108179540828 1.9537558065527473872 -11.315377716261613372 -4.0957043098480152565 13.830732916701464319 5.2858554841611500308 -1.5280409673835291251;-0.1248231611018385534 1.1865356652089655842 -0.48783798170743286615 0.29231558537413587784 1.687679770495032372 0.25856314359636917155 1.1346058026656646867 0.21799438771958765115 0.06094456410240121913];

% Layer 2
b2 = [1.5263386701173440319;7.3756334437127888393;9.9208656989010073346;-7.0485737030000548131;-10.73494252457859055;-0.51908253026695216903;6.9450826139116408342;-4.4889819382757156774;5.0674788830234120951;5.0678803614849519121];
LW2_1 = [0.94101387400700653441 -2.9106907432715436457 -0.4679709690999115379 1.3855277727792834064 -3.2491940422365388486 1.9468507710455540405 -4.5321703212071433597 2.0407431465362262379 2.6793539346945078883 1.7510977794152347009 2.3543019494671297842 -2.6694561990167979815 -0.24347342310385605701 1.1673314743474740851 2.7893761670474450654 4.4140659884591526563 -0.89784234754026515102 -0.73828358729649790071 -2.8917099772866952101 -2.0587977529825045586 -5.4491151996876983077 -3.0623095182925381685 -2.2728437024357237561 2.3045544525100685185 -2.3505635556428328137 -1.9201044254022885394 -2.794629677571289772 -1.0256235376516054103 2.8675934773622606677 6.4742846213037523739;3.7355766098720732948 -2.293469024865717909 13.064065427964909105 -3.1913246464779878586 -29.859520096440011372 1.8972939544784537258 6.8665223657546627933 -6.1251882048112493706 12.942478950303319962 -2.1359317277415788006 0.12420310017660426971 17.57092054112001378 -7.1193674564208757616 4.8799036667303070658 22.354447148468384654 1.8124606269395220348 13.272471797039026242 7.5282405232856781652 -4.5520303881002712743 -5.6165967069940654355 -1.3215480895884610391 -0.94708617663290395328 1.6715888208118447888 11.724148011259066138 -11.810819170915339171 -23.030157447117250769 -7.3221373649053260735 -1.7739019652250809234 32.995899358390119005 2.6882290227805936134;-0.56846494282459780756 1.6358855770935296636 -5.030118811458978989 0.36076949959987414385 -5.688129526205750075 -1.4807732470620189602 -3.8006343203103503114 -2.4007764923942707114 -2.4628814065488779939 2.1685771303491256923 -2.7724626702847459114 0.27296872466911226152 -1.4325361440234389843 0.9164964415897922656 2.5450789108839395603 -5.3828387327229529902 -2.0695635488779662303 -3.7586123306189640481 -1.823969245399891248 -6.4878286851905393462 3.3671651592334717051 2.9848152074786020727 5.0519813315712669421 -2.4986527149879012022 -0.47748998819365023882 -10.281702677149150205 10.379925201445141525 1.1169716162044320829 4.7734363832750892342 2.5403025769483624785;0.17740775578995102157 -0.95528358342309638296 -0.25897934740784522845 1.06623180331842482 0.81418191802583572692 -0.21425436059040289805 -1.5094851382596206069 -1.7193005535405512241 -1.4063959191013641892 0.33164956919051225004 0.6983752303723894661 -2.0403442464027632397 1.0282071571045541258 0.96451582432944305001 -3.5800198922788166733 -0.72229366828200702866 -3.2134018264449664137 -0.90492730831365875144 0.10436103296837145882 -1.0593243455658862828 -0.3761999363082711989 -0.69715420339152711637 -1.0521280812139179872 -3.7073196562368071483 2.9119119682928324977 2.0677078961027604009 1.0473654383053518924 -0.21320622784878703881 -12.543471797550695968 0.054035031780870303175;0.13304246959289972585 -1.6049622748188709931 5.7205680322953522676 -0.830917036504614126 7.047411495840782969 1.6517877826875499991 -2.0246902486204039207 1.037589467842503721 -0.41315998250920410451 -2.2522311653053925795 1.8317408567639648087 0.48559535283062910116 0.92432173613539447832 2.0781306492938962371 -0.74937990282549604704 -0.57759624063710490738 -0.17272152771790955095 -1.3532784294663984248 0.77791253403185867121 3.772618437515369294 -0.78016848255928861189 0.89812282959904021773 -1.6493128371025227441 1.4391399776813267319 0.051631296867952682084 0.4050050833807437578 0.16667653518846844984 0.34850657700749343437 2.8186161672831300251 -0.25537043276339710385;0.16454804785940380096 -0.90686331084902604971 0.13098471213487472564 1.7510162004940943437 0.21404947170171337301 -0.58814288283021953774 -0.23002709925539363556 1.0201001214469336276 -1.6029559550452714412 0.14375089611369842868 0.67261396912821969085 -0.17413183027317608031 0.91505048386972509533 0.72865557983664919206 -1.8687303822248388929 -0.87933480128546448462 -0.21386544830408041329 0.68198587270052457576 0.11441250871158403779 -0.19414277231995669126 -0.52658238829624304245 -0.67300557285752005132 -0.81174192373095699349 1.2804747273374024186 0.64693017102684491348 2.3152469786523361961 -1.1325232028003260343 -1.2667821194215651559 -0.055862965911261740215 0.4340246069923955341;1.3056930221729210473 -0.84683980781583845232 3.2525085427263422666 -0.89721558881038043776 2.339445357886407173 -0.53785208807112505625 1.0081376209301389224 -9.6333391489908368754 0.69645625465531313392 0.39808107599619119155 -1.9274289758749556523 -2.4526843235551614164 -1.1467499175200432937 -0.40776189874641644018 -0.64161621164573112708 0.97868349979550206541 1.6234579855362798462 -0.35160403220327646201 -0.077265469830244681848 0.72577779225252470052 -0.76626153156461751514 0.10170402997777301668 2.9307845542341377509 -1.4215224898829195155 0.30271210974345280631 -2.8707371369730410038 -0.86467136565594859832 2.5463425023985650775 -0.81875121336899525648 2.4323437748417573445;0.21640686902705597028 -0.52968094847530289293 -8.3040867937110238728 2.1736085027398797109 7.2134852625049399322 -0.54108325375633714582 0.028141399432165998401 1.5486372766365927856 -0.98625992161748443632 0.22046362730290169907 0.27547508616057786446 -4.2388218292059187675 1.1996869180556362089 0.1001212203229827763 -0.53050389023195676508 0.93844876523744935604 -0.39909858596221786664 -1.4243996679446717657 -0.47971393429421810284 -0.18958596361152271759 1.2806611634557161228 -1.1255089389558883273 0.48786612936226209092 -0.70680814596899765867 5.1603735647190038804 1.5952743173661758647 -1.5713019569568242861 2.9173151428212511505 -19.007948618813156116 -0.82580610947649590869;-0.15475085854910414218 -1.2523712331373568762 -4.2540006840036932445 -1.7221593393350322909 1.2364435643947642163 2.1450771776603456686 1.4574479053672795281 1.759911405709233545 2.6418674541290254609 -1.3623830156163252436 -0.912572273489074548 -0.28812948753061040641 -0.83623636300592707471 -0.017357622249729445091 0.20274319943219837326 -0.27879408759089424619 -1.4520477866398622258 -0.82449027599716684378 0.86269055038753628128 2.0841200529652295081 -0.073169704266449864249 0.066670817027785711506 5.695518024344887742 0.61748303116537905932 2.2790220043102884517 -6.7096361816188787586 -1.2224350859444248485 0.64039166670249414448 -0.0010306540776911779153 1.9811159400616218385;-1.0177077643257874673 2.1770589673208626813 -0.62805555561487924621 -7.8894038954073266723 -3.8714141832040747104 -2.0839441299580436784 2.3229474887719878673 -3.2835929103156198394 0.16850471285861715054 2.0875309756822399976 -1.2602223702109391912 -0.73609371134946399717 0.1882506610208051212 -0.4225795492196383063 0.86513482106241046399 -1.9496586943180358453 1.8174562283577511135 2.6565481672887543141 1.7561880878069107137 -1.9971922680191427357 6.6015574409532424838 -0.64418935345556493921 1.53725459017760957 -0.57385935470124938007 5.1785934993683850891 4.2741031513969423017 1.594836005249369304 -3.630843522896894715 2.981725169422071442 -5.3610916722742025797];

% Layer 3
b3 = 0.54673147852955628068;
LW3_2 = [0.0074902466504978698625 -0.010029684365669101417 0.014250220263146585326 0.37422575593669010763 0.04161188461621648238 -0.33372369527132933475 0.0099905738595786053313 -0.16156889639892949018 -0.084354924216659368796 -0.015466593091587476153];

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
