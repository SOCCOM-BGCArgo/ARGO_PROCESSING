function [Y,Xf,Af] = ESPER_nitrate_6_Atl_2(X,~,~)
%ESPER_NITRATE_6_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:34.
% 
% [Y] = ESPER_nitrate_6_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.95648195453551887457;-2.2364651658396383915;-0.5267279978561913989;2.8522729267819242338;-2.9601670732782965345;1.017711582043450802;-0.29311900810016877639;3.0576749317892675251;-2.9592386089345725786;-3.2021442063210296958;-1.9169803361967685529;0.46727939867985340827;2.8711099413997391849;-6.5806600495450435773;-5.0053638399489646105;2.3049638987992255856;-2.7685963665695161318;-0.023378186487250783576;-2.5570837033340345634;-6.8012665528537370463];
IW1_1 = [1.2514023085133236801 -0.67660784617671532182 0.5606812448679897587 0.39460247011900079217 -0.88294569250686127937 -0.20087644095670093014 0.40613535084142765585;0.035239861814027551024 0.64021390948078438221 1.2433350301492205325 -1.3810010199672115849 2.7930083835103243217 0.30977500410186908919 0.16451818144672741839;2.5359209679068386301 -2.1718759661397766259 0.79947448676458643391 -1.618848591346393917 -4.4065981423860520749 -0.72281223410257267048 -0.56405280342030417184;2.1130699523048099842 -2.0662965388897704955 1.6691694927248801328 1.1615478879454401717 -2.1074657945528514347 1.359000586481639905 2.6358722884447329804;-1.567361472420450097 2.7397939830148674289 -0.35049802276233121834 -0.46599179670636659001 3.6951602981722628805 -3.0731875856851633699 -2.5522403823926449462;-0.76807871712355901828 0.058986387914963939183 -0.16498098749699505183 -0.57534075289091357597 -0.25244092729082912019 -1.5977485815437200145 0.63711680902871514753;0.1685111047089148939 -0.43556934582489237329 -0.61797366328082270481 -0.69426787262578493198 -1.659609860795985492 -1.2442481787873591603 -0.24895532862560817988;-0.86916729812930149102 0.52815864668168299723 -1.0701059126883198935 1.1851380568146487882 -1.3967501466673044597 0.70808152593979190037 -1.0016842185185657055;0.24361797876887852921 0.26527831521934025494 0.40036561587195956813 -0.23319929209730821951 4.7827941263408035866 -0.48138283528796038668 -0.92916379374724977769;0.083071210578716969675 -3.0556534744463417219 0.49153443560664072942 2.481652997793076576 2.0425048285912739132 2.2576753313090121011 -3.8359065932812375088;4.0787944735863286283 0.4138888199410445945 2.9432485809726935955 0.052652488576166299461 -3.8333280960391231673 -0.029860958798663392311 2.8987422404678033772;0.091007536637455907491 0.29956334164079234661 0.13638364032215838373 -0.14099213790817327618 -1.6987245126242322968 0.65510513327311248322 -1.2833316045679583617;0.78236938619116125349 -0.35211910469880824959 1.1390588917556399018 2.9502414592522909764 0.17115995660877547646 -0.44910669948373405092 0.74660874437873170351;-0.45960109652772690447 -1.08711727576335071 -0.29260529123496370474 0.32030987205174721311 8.0295700439281532113 -0.089633173854507877665 1.4931581839989216398;-0.0029140987533855635186 0.25953970420865951629 -0.01071012989778371248 0.66879975273883640607 7.3182155426814698984 -0.71485284915961488927 0.66815200990616285814;-0.78841690808711639704 -0.071512396116152787684 -0.020433827889594153443 -0.29978843218932960113 -2.5307965831140424307 -0.6717719894090042887 0.85065376366353662174;0.22737550491685318277 0.033814743841575858885 -2.3201768116398442388 -1.0261979235865599236 -0.77077378932459650329 -0.68355677428860806 0.56631714394701226922;-0.063368126370308142103 0.65800674996247965343 -0.26097014621138109325 0.32828072161290350328 3.8243365211476190169 1.6522745771991398378 0.90986386389110507356;0.22808378629741765664 0.19244612550364575965 -2.2798102624554412721 0.24991192512988746888 1.3448277242240198603 -2.4903501605382847295 -0.84338048246702912447;0.65453518086975248647 -0.68899143037905485887 -0.12573200697123965153 0.35056296246206036793 7.2543331730043476568 -0.75764032252164204007 1.3893314166715011737];

% Layer 2
b2 = [0.69275101574457842446;2.4418043787593810556;-1.7934432071960519917;-6.5167147949715884536;-6.9125868033303818905;-1.2114520096066336841;7.4123001197923139571;1.3013931432915126774;-1.8952434050269792998;-2.4700538961993836473;-0.24001333509511815767;1.5065278011405502934;0.18438248143917890243;-0.66720695434184340566;-2.7460946970334418893;2.038051735770738393;-3.7242624474245791255;4.4411824831405795067;-1.3003903610403046986;1.0907357975145468565];
LW2_1 = [1.5204477421720226094 -0.72456665678749188064 0.68617884565797371277 -0.5961686429123115305 -0.65089580242457378834 -0.30208332396530340747 -0.64081383319965357082 1.104508804240203812 -1.6138404573588285285 -0.092593687016550685343 0.3840610129621647717 1.2110283148476461701 0.6051653283449178744 -1.1859541521271919073 -0.7246496399740212091 0.66629306852021519258 0.49076920013992869629 -0.84423385932359396033 -0.79973592643549462622 -1.0514429566836198671;3.6739130175115608523 0.41724682419289260693 2.7661667513494503368 -1.2990846551669792408 0.40010336088224307538 -2.5505315570234214917 -2.3173179466450273267 1.3143233111064980267 -4.8614473043799515395 0.40227278539625194176 -0.087699722032872781674 0.43124764018794997655 -0.13317640217349868692 2.2516470891695452572 0.95366425919918307308 2.7160743216400882183 -2.0190996146455075788 -1.685139207562595276 2.0781342207120117038 -6.5972796704380574084;-0.83444076644171405466 -0.42858162324765691764 -0.078311640542025648393 0.13042173982001500909 -0.40886253514234111028 1.7324174805570846658 -0.85478847078511654978 -0.46882904280970494337 2.1368809990908781415 0.0079035962437982842238 -0.33412601514532086489 0.34919445531793302795 0.46659012396781546217 0.87035620262944302095 -0.22299364890342854295 -1.4443719348828536475 0.61375490536218579862 -0.42394627144897401028 0.34226433911933351606 0.66021172883743906024;-2.5722334757360623136 -1.507588745137346864 -0.25444489701799916093 -1.0696360073038011507 2.837504271777327336 -4.0773136771518538879 8.6373346589335522339 1.6656377431854831883 -3.2052399689206607647 1.5902164830420186803 -3.5899276800508235219 -3.9652436101804067192 -2.6955603186917249658 -5.3840806361690720649 -0.079661845551187826509 1.3507595297484669139 -2.0726064235092676036 6.7455639316378812609 -1.0053046198642154874 3.8174591247945417471;-4.1404700718040086116 -0.081858490972692848908 0.31266121066703955744 -0.65159363462826525026 -0.14701712865639662842 -1.8623813466156036167 2.3762858547544776222 1.0205389368713864151 -0.77702210675986049804 -0.11950717067499282142 -2.5198355739514624041 -0.9736506249875703789 -1.0969330901924829735 -0.54455877940017327621 -0.069986131107483931513 1.5903552800886604501 -0.450951049058789466 1.0365044566116441427 -1.2562829342665897681 0.84411232345927988963;5.149477817534339863 -0.21150144721765823386 2.0374633284901766395 0.22916929844719352571 0.97944948705537493439 -1.1790786927436556208 -2.981491993236212501 -3.4257648581271231691 0.81924000088914206152 0.22229927902500692949 -1.403085945130073009 -3.6787513452672282455 -0.90711171620778041635 3.8321924654119556841 -1.5222360165142905242 0.82333147199357903379 -3.1121367393516696787 3.4433832310601442117 1.2455004571508303535 -3.8643022270417137776;4.1449156924207697372 -3.1606153291510246106 -0.41902747319325939213 -1.1300855649162786509 0.23801477626774419916 -3.9116857945362188964 -2.1659217256796847728 -1.4718378983643851665 -2.2828443867700713987 1.4492708028173513224 -1.5667914807659599763 1.2307387113993133454 -0.28709315833053206113 0.86680687313555293549 1.0242966885685507172 0.92824440271511621248 -1.1884268421345414168 -0.76124876967355470381 -0.0049810884855275364136 -5.2062637899789949714;-0.76575502056173272436 4.5402671614569527492 1.1142483218864365213 0.53498131052795905305 0.91752784878317783512 2.5289909364755689403 -1.0667270983997882983 1.3405864734541965699 -6.0109378057514959437 1.6834795429525930555 -1.3619112535177351209 -1.863011991455505223 3.7040911447231383846 -0.21361167990277957207 -1.5746866914213255573 -1.5410514749936350398 0.18143710644324728687 -5.4362994373730586517 -0.35995264648462899082 -2.0241189117630122851;-1.5512493326818079797 1.5039446939139291537 -1.4760791826587729147 0.65668381520512530969 -0.73745227864304763177 -2.7901446622310963619 4.89449910833633961 0.93027666910193951288 -0.53097990569744657474 -0.27516069024967282353 1.2745251231082230792 0.54121961935980522629 -0.18554226361995157513 -1.6399876682877079226 0.92974324426243004638 0.92962555780575972797 0.29427651697244144868 1.464453272819524976 -2.4948040909836040058 -0.89123285233689097939;-1.1704318440251477718 -1.1162916683080872371 1.2851223506092441795 -0.61374427027455669137 -1.430780080899481721 -0.81725970579074391775 -1.6430729218031647765 1.209692387538460423 5.0261979593401218125 -0.54090876410396027119 1.231433131037302342 0.83824685747871463537 2.2037035683243875184 0.72639188181105462938 -1.4472520074851658034 0.62531273367517681283 2.9886803249633406843 -1.6886256978428275666 -0.81362584770257961875 -0.36885313889132570475;-0.50938857403991910644 -0.44530232059738106276 -1.2707744482952878506 0.1796040222086513205 -0.60969761518411813572 0.99043428914175080102 0.70211227457899905158 0.66447594108554686176 3.1870692091151320291 -0.16682495577724271163 0.59335301575844945354 -1.3981037867985801615 1.022748922284266726 -1.8724154847335214846 -0.87747279453034610075 -1.9333915092968712024 0.38577817855377399825 0.37098215780879184589 -1.1089075186128500849 5.3358292976770576033;-2.3627784882742921546 -2.475226376334387357 -0.55532121794694278005 0.67529475196528321312 -0.45295891665317444819 -3.0206284291245304452 -3.2502968270044072518 -1.1485165466261546552 -3.9172918960951941258 -1.993562358675387447 0.53673734265020811574 1.826098120594671137 0.3481747233322845525 -3.7029421431263456199 0.5721629547048885378 3.5238464907776712032 -0.53687247400184068269 -0.12292910117501898437 -0.0095390580235230112288 3.5070005228368419381;-0.10405088409187347731 -0.10195036986874599083 -0.13801584561669077988 0.074856864884106005809 -0.23975593906434292113 1.5143037036982043286 -0.78810600593923285651 -0.077702611455345627078 -0.19894924991046791796 -0.058217982376465318706 -0.56590095435085352449 -1.2504460514314594821 0.50718054724903571806 0.73041593306404095642 -1.1858944911279123424 -2.1899384302443678507 0.69911412302043529632 -0.56699355258951555925 0.0029375100748383134452 0.9795831806507014905;4.7491292342757036593 1.0830788301308691413 0.48582108602286355747 0.25240821272152047516 0.9508409206466501562 -0.95329797760212287017 1.3057971045265333032 -1.1465267166334407456 -1.9430465932171225862 1.0690015554149558952 -1.846213681028987974 1.6810730096649646459 2.1482145837982700165 1.0934524015051545476 -0.021420136733283856589 -0.92471032486795534666 0.050294049496913420527 3.3948821517675189519 0.16853319066337726917 -5.911490095619017282;-3.7813031173271842711 -0.96290233034414685065 -0.80245338324059112267 0.26881486305037327611 -0.89266262145958841234 4.9833828575570500163 -3.2948811683574055209 -0.97634705938119015212 3.7937271104680725387 0.10031344518627463569 0.25416502563596110598 0.0013577845994309604982 0.95383376285675225414 -0.33455812343591301294 -1.2961098163501871827 -3.8612739284860029265 0.2497228954431715553 -1.6519219371828226528 1.4637210250725689598 5.3293321143378298288;-3.5914279738041900991 -0.1581313944266718452 -0.14824369669952466855 1.639811030273421677 0.53347845615251454632 -2.3233834761404885683 -1.7607815854448811699 -1.9859948545310193158 -0.50640892124850944711 0.50495331855450931169 -0.84505782463348355638 -1.9381324921836970976 0.24131501248335404664 2.1035016015475553885 -1.2753896696933322286 0.23009283244872916208 -0.27855858680516476511 -1.270776735005153979 0.20031907214008837226 0.59119699796025737548;-4.9447872521953373237 -3.0126678825917014848 2.0305789618846703526 1.9117853649557992046 0.4266223833542673316 3.5177146660089650965 -4.5156187072367286461 -2.2056531511121906242 -2.1279773057088631205 -1.7736714506678781245 -0.13872957787883247183 -1.6775254759811197847 -6.7930371835514913315 -2.4642748177340734195 4.5165362988685613388 -1.0092786623890219744 0.1710040507767898732 -4.3508049795051011444 -0.29888921216359254363 -5.0179024337394366739;2.3776737683463804451 1.1744702107136191227 0.015636546429505737699 -0.89472177555312415453 -1.0716808150593735149 -3.2968986665930963831 3.6124570310214716429 0.90479135931124665682 -2.4543700336903429005 1.3678426220768462596 0.67897068230625412522 2.236578962591059927 3.4925849893482867614 -2.2698268411382374055 -0.38864227398280604131 1.9759455491649622338 -1.8916391223073178995 0.98593545926586723116 6.2289644874035587918 -3.7537030959111885586;-3.8068213508774433329 -0.6487914229950462186 -1.052051909096895832 3.5725196235678109424 -0.10621806962414902242 2.2409992812969274922 1.3732286056525477846 -5.2881933976712218737 6.7775476450620155688 -0.35379738534560389862 -0.55674949273975748998 4.6210887455314173522 -0.40634325566118545403 -1.2974809448651010158 4.8916886790139546548 0.45330186731332006156 2.0558679299870878943 -0.55379741707143370633 1.7466517541353006049 -0.39503979408835659948;-2.7753493410732197333 1.2385411002472850317 -0.094863888651495448245 0.52511370379852961854 -2.5394082806024353971 3.6702600204940996065 -0.062338342366453705634 0.18297952757516175448 -4.7360176352169212066 0.023746967486553915982 -1.3798609202447222888 4.465012992560143168 1.2906500482438159061 4.6755907758736992008 3.2300544011741361139 2.6054925876304420029 -2.4368945846874963124 1.353377437914889958 5.917970392357490006 -5.1434720123167814165];

% Layer 3
b3 = -0.61207414684257888116;
LW3_2 = [-0.26624754479732470891 0.088872512781174334284 -1.0453399105335594577 -0.042390381139433601365 0.24232484742021703772 0.030116802422605917083 -0.14445455049844527218 0.045880600104816010587 -0.12625553634522340696 0.13673010387782127961 0.24562489303377402217 0.06938596557747669924 0.59362024997142515126 0.081140482317105633303 0.16386713834377847343 -0.073181930856821120202 -0.07456448599123376153 -0.058923497920355907831 0.047148424282391646067 0.054230092547534684833];

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