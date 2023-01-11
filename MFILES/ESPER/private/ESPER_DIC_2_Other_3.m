function [Y,Xf,Af] = ESPER_DIC_2_Other_3(X,~,~)
%ESPER_DIC_2_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:10.
% 
% [Y] = ESPER_DIC_2_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172;-0.9;-0.1386];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402;0.0419287211740042;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.5119639718160626796;3.2677821265466677758;0.791755245767149618;-1.4373146234448392899;1.3691796101900586269;-1.6488389118292259283;-0.97280279329367158603;-0.20545699660583846557;-1.2144422663750154001;-0.43862781874942513571;-1.2596999320682371604;1.152833978737149323;-0.081021482830191893232;0.26117241809682589793;-0.64064828347292779309;1.3357828727211851216;-0.28780695497462915888;1.113166755396309382;-0.43447997773755886453;-3.8055154040172136654;0.93580472232489453344;2.2153296218326712363;-1.9119318578331672143;5.1236810221539830223;-2.1267165424151066588];
IW1_1 = [-0.030512387479017539366 0.18797824711336702519 0.45414766077529677979 0.027278254098867381339 0.21466818960505934522 0.59191429742019796567 -0.050678395523011440793 0.15651360397612959496;0.22087514743646102788 0.51521210975828435341 1.0233447175460874679 0.97536238581642276113 0.19035943153389248139 -2.216970789491135907 0.13382691370509308681 0.58402448173349863492;-0.24231415813384343805 -0.36002956384151357794 0.8440059539970934388 -0.19211382322773432052 -1.7608408852531471922 -1.1715084591179962814 -0.40618656557716725786 0.78387251580656236438;0.16481581718325463704 0.54036878100454188711 -2.2897217580242243962 0.5057310784720870922 2.5267184734663015178 1.8262578427829581074 1.0350228826509799696 -0.74400576406461482737;-0.22704181240788839924 -0.66264862587610451339 3.177345942512330268 -1.0149522263002344058 -4.1828825241064917506 -2.1370970329871044058 -1.5656925737883455874 0.1844562852287889021;0.18578223534153312513 -0.40151068117031302185 1.5716504288882746643 0.050555395082877541135 0.12927186185925801376 0.083991443344634403312 -0.54078021363580752734 -0.32576315882915907185;-0.12973663033228385877 -0.61802620442936262446 0.20451534895303449813 -0.18407637567695936398 1.8061547683059686609 1.2411472517251731507 0.084828294906882872883 0.41077653307018763273;0.026097828062264724625 -0.052952989351721102551 0.098416631951486019392 -0.063676859423615075895 -0.5651505534447015533 0.23120423617215221657 -0.45147942685410252794 0.072207300030800058988;0.12148778406358914927 0.35868615213589061241 -1.6427706709049889611 0.33939202838642851257 2.0908234375612497757 1.5846962581206129883 0.84963606851619688509 -0.61236173651272540752;-0.07711081385700718771 0.091959294313007866606 -0.2872735244747674499 -0.13814796785783134414 0.47880148813675260655 1.263384210753691228 0.47224755769271786754 -0.070602491781870418963;0.058807854397075920438 0.28799455786075461861 0.36121860340947758994 -0.13069577702167844979 0.69024148820676012317 1.2475673693910953599 0.078901927218154657706 -0.39752805110806838096;0.054069492162342726316 0.30502087050380422006 -0.5578925902012780158 -0.059857672762694026647 -0.22918778136757031461 -0.52986721449660767469 0.21812122511974407768 -0.25082304709469832593;0.11271076178198347062 0.042398162478898451799 0.075457135785417808749 -0.15228676162230692093 0.43413448966280482333 -0.070212199576975595217 -0.00034869480799887454673 -0.28749549263567752977;0.0089320994103080902232 -0.34584704474126365392 0.48154181994364536612 -0.097290934701349446323 0.67656733069605834441 -1.3049153719611876756 0.69894024756372563623 0.26022429756359333197;-0.1145104750705504687 -0.065606976093803906491 0.57853011360511152894 0.0016236810248497408372 0.019040318635814436343 0.37119474990802442393 0.089885131300820098055 0.092884925056871417426;0.34635478499814131848 -0.10303906074214011501 -0.15845677011340977836 0.069070278233032397774 0.034377146208984160625 -1.0116578219321354926 -0.96493596490613786631 0.72334064561766842072;0.38193634305315510291 0.12675449543362374283 -0.21902534475211229981 0.35290501690827325021 1.3555107807103705664 -2.0629791623689190949 -0.38657929770034116235 -0.17174792860157164776;0.24909392919396811106 0.076320243607640975148 -0.9796203560278489153 -0.089600439597693815053 -1.1194818242264124031 -0.044560698937921028318 -0.40446439482262092069 0.26846966099114411008;0.037421259705469125945 0.19989253116937769073 0.055920092187641406989 -0.39391396921787674712 -0.66352431827298274936 -0.57432549362076945609 0.44391499419610813693 -1.2360677252037555185;0.0027683837532497162018 -0.14527384917567873202 -0.29163426565306466154 0.22555981493217433509 2.0740382771857710509 -3.2109383054224904797 0.12238343393538703152 0.064547170861940894593;-0.29678943769786225726 -0.40816247705158181791 0.90011220183798090755 -0.17599433666745314953 -1.7392550502066441709 -0.67869743436592389063 -0.1532793619719613476 0.62187187947122812837;0.15781166256112605373 -0.18895093158144835033 -1.1191891955808266346 0.034600535566385620667 -0.18879846318091994539 0.024100247499033181087 -0.0029730011953269935843 0.32398548145975242241;-0.032774477060561857389 -0.15372665145218142313 1.1233919061309731724 -0.17532664270569794929 0.024490266301455822229 0.14097527777051713227 -0.18014454693074097302 -0.4205893873529425453;-0.022170928886865375745 0.10961369107479752849 0.51981429748592078433 -0.21933302787644057275 -1.6302724051833115126 4.0415397712142970832 -0.54137849580937946747 0.36252884979422911327;-0.099137690739464726097 0.099123016733982119431 2.2958039493626114513 0.47554252185274314968 0.063639142150448108071 1.8558571356661353047 0.06642471780382672164 0.6646708629379859401];

% Layer 2
b2 = [10.07588801626588193;23.007737689032438766;-1.3748942143622779888;-6.8751424841166004498;4.9551768636558675141;-9.1565402608151256203;5.1627744803906514548;8.0634193114968830685;5.2349976997629390496;-13.305327525385738596;-5.7135604851623948974;9.5205890934842312134;1.0942585341143691302;16.295602224943031899;-12.628600536852321312];
LW2_1 = [-2.7533825318128135606 -0.76664682377263104751 6.7214162813139948227 -4.0932895654338281588 -0.3244133398268751467 4.8359736892841214129 -2.8439963482495764779 2.6852927864615301345 6.1901033776089926519 3.4998056880853591899 -1.4142818708504785175 -0.090570821848577232061 -0.6990917938527130282 -0.28015927921412964174 3.0767633229261446282 -0.5512721157410508388 0.26908946526555022061 0.63811564132061227017 -0.64908940585886054642 -3.2063950470172692953 -6.4472815494267212699 -1.0431684673920358364 -2.4953492056088788509 -3.3605598455681220393 4.389586322135776264;4.209768664511422287 -2.5792377454102854983 7.0017140604120315928 1.9865374775791044826 0.55297174798050796607 5.4745345935983031893 -2.2050595050418566601 -0.64993497132461797161 -5.5498123305340394396 6.2305173721799729236 -2.6099478751993903103 -8.8025331209592589943 -5.0092750861983086352 -2.4194541880750586316 -8.0399789801189545102 -1.9061289916723493931 2.6794637323488434966 -1.9206977110144334286 1.2190959859246532915 -2.0131005381768916074 -8.5561471330616978292 -5.0924206214445817409 4.6538421197031372145 -5.4481258439137274507 -4.1615850955513966625;-1.2888468480053438814 -4.7612075242481841642 -10.049039322536808427 -2.0148422500190954487 1.0781208977399607818 -1.3938821401258438826 -6.2637439532050409596 -0.7748076295030322802 1.3067590359231640029 11.576433050246976109 -2.6757892037109094474 7.7896750163672958678 -1.1088223577718294699 2.1153646094080564843 1.1239963429525903216 10.83208853454336662 -0.35981994722141663612 -3.7057104940051694975 0.46926991784039862043 -2.4637203713483617484 5.9814888021298555643 7.2737468295022331688 10.566618979880290752 -10.27803358577812709 -0.58181001152742839455;4.8240318433778961449 -1.3121420842389825978 0.56698318003281567012 -0.56078159182669828908 0.075493550462100295073 2.5845712293680693783 0.45131633754122102431 1.3999929224840568232 1.1928143554570502793 2.6631772231902903769 -2.1180349939341329701 9.505072879683863718 1.7051443896673150746 -0.43230824013228547065 3.3201349065239766922 0.88444743679893800259 0.10940412220432263468 -1.0740763015780272482 0.50347239273693156214 6.3577903784735454451 0.15045388259879807769 8.5361815715488074119 0.9843651016978663959 7.9911006989588351246 0.59152074488970740607;5.5556479953987025766 0.056482522719341878015 0.14507554359410498868 0.1134906065890050747 0.0008767536239165954438 -0.99783021276361760776 0.30533788307899778491 0.44661624379231851245 -0.81826522026132164633 -0.15554655645706691458 -0.55614083398752278686 -2.4061432122609507367 -0.76773415172198311573 0.16220601091776926106 -1.7752076922143251725 -0.61813232533317086226 0.1603919879915674529 0.13706962279994588938 0.51489525516933953853 0.70776818693859278042 -0.4308294516837496424 -2.5205554066041075956 0.77467896172330219517 0.82317669715143204634 -2.3903182542634957564;-7.3310767910162040195 11.374952869133451827 -6.965425593461329612 3.4914318412710114714 1.2132540778263340986 12.544555681264851543 1.7605662430822648368 -9.0071072964836655927 3.3215232309681703704 -4.4863213335811629179 9.6365192034595459347 -26.433126395609367165 -11.379928798495127396 4.6907458792924776603 1.1880630889284531371 6.1642376051048044516 -3.4171942941980035791 4.9292903213033048004 2.738035815990595534 -9.551719134555396451 8.3633627505049670248 -3.9289155511801263643 -31.847565716931935498 -12.364372081290635919 -2.8880056860397278484;15.092612242285129653 -11.474085704228418336 -13.861400484899803587 -7.1634029662216072154 -1.7357552745990292298 4.9019375998785292836 -10.122688101716525111 4.1405926531474683827 2.4208368835191831359 8.6144718969527662722 2.3497172769087208977 -6.6498434066870890646 -18.566063496602360772 17.354833358246960273 -0.51098237010472713493 17.256236451433064616 -8.2710599309157597503 -18.695350853913353717 2.6093344535950491192 4.0233653433827258894 2.1408262120499319536 10.437594589251972366 -9.2577629464427229067 8.2251580403098305538 -2.213703573881275144;6.7139540784054396738 -0.45829682969048279872 -1.5559505801443715978 4.7257587762345938387 -0.53042031425259106303 -1.5304478874674787292 -0.043280850546631194353 -0.011821940380616239447 -6.8933447989361242847 0.064108945230794078807 -3.1378581728967946951 -8.3647765697566249798 5.236176085733382557 -1.1367665249687795015 -11.760985962037414865 -1.6909259066642863267 -1.784776088477757483 -5.9000107821347347326 -1.4619489828439087287 -3.3339185910534419044 3.1254345371338230208 6.6290660232931211127 3.4726461907349084335 -1.5967694557769751551 4.1691390826325722685;2.0616703597702530359 2.5693164055894706088 -2.436824885331638324 4.7221879122794154782 0.18278772408876042821 -6.8722568588123218447 0.35112369221828260146 0.056019198622725058234 -7.6127961110687110136 1.4927454247071620941 1.7464938701353378558 -5.9712463891433165131 -0.91672873148339095728 1.1298458166432689964 -1.8553292713827067573 -2.2286526532586456995 0.97079635016543752712 1.4117653410457688956 -2.9309634549532144199 -11.997908369600803979 1.630035866611926787 2.9527234374662452154 14.51085284235159456 -8.686724818613068777 -2.6857764497106044743;6.2171685320893423921 24.983261904907916318 5.438439388572876787 8.043217149942277544 6.8266958954092427092 -4.748087373643841147 5.5268600021615581497 6.8823073292062622031 -2.0725737218048423394 -15.114440021363185096 5.878346857451953511 2.7137727567218319358 13.768749938021946022 4.9207501738835111027 10.254596778919387035 -14.306810004597116048 -1.346162786571404979 7.3268465572934360353 -6.6272245160347065251 -2.0493357685187789663 -6.5271421788041203982 -6.2185423651176829907 -2.3027544817298091218 -2.9326297217867858258 -4.6115058666365600359;-5.9047539497564622835 -1.4828088130110315124 -3.3857111315544989871 -0.063769695658276018269 0.1958893389347954872 -8.6325302050365220197 -1.1028731196059264885 -8.9097393709132077788 1.0261541724388159036 -12.68985325469520653 3.4886896295754366015 -19.053056093482084066 2.4928237165883309068 -0.31189268546871784515 -4.9424538322626272446 -3.5454460972386692141 -6.5314411617736736204 8.3328397084635916769 -1.7148000014341293618 -3.0348711194311528416 2.2619075270581840975 1.2202245235878255158 6.984159068903620593 0.56509111499482467256 -2.2483667451616735633;-9.1639610799718855816 1.7208625913717863209 2.9617731963362050607 1.1699643847579792588 -0.18682005887761429785 -0.86639874135070726879 0.3408586441027150471 2.7966248982418302482 -2.4679644519841108519 2.0466021441602229025 0.094751318819690663608 -8.4219168826030994524 -0.30560595494165759822 0.41083780793254121022 7.0545943870758769378 -1.1442157875085328111 0.96406033164521420797 1.1703688732884709456 0.90700594081059171536 0.36043036449379262098 -2.7783136828943777452 -16.033251873200001114 -11.249349727403240351 0.77711053683425379557 1.4347429591928708703;-10.015046794391263774 -0.053330336924846608038 -0.071187930231150026161 -0.52375467162580846558 0.16554021654621089987 0.28063921044262823967 -1.3259633758915398971 -0.97136660042767108969 2.0388898081001625151 -4.7689406642672365066 3.3238415534565581844 -7.4631345495207996876 -2.5812779663722795753 -0.020871182412306812232 -0.04437067530455353348 -0.63834819786876706438 -0.42760225316157063347 0.63589370055965688255 0.10394817626580871162 -4.3797179530627650124 -0.27297543763532056138 -1.6419213872196680715 -5.4391202556012698821 -5.5198252036247250629 4.7857646145265544391;10.031211943509406481 -1.886387838641351733 -1.1728637801753627468 -2.4319660403340566113 1.2842874656675591361 -3.5090217844468680752 -2.4949957622125813472 -2.6224128788109579702 8.3146329202983881856 -0.048082668408596936083 0.45152838059481464894 -15.979909414653652888 2.7719292238059347611 -0.41408686882905854088 -15.029255145619501377 6.3841975016929239928 -2.3983654197046515577 -5.7234260588336791642 2.3343935950857730433 10.970251334403533505 3.1931160001968570405 4.2984621981180728412 3.2564389992878872704 9.3544038814875634102 5.8235887530926051525;19.829365466167299559 -1.0638443862352777725 -3.4493085786917112578 -5.0071735223715796437 0.15398321655018987197 12.824594121896053522 -1.5206137286453478286 1.3963893253479624512 10.916290936364354991 -3.9667350763059205931 -2.5729283291878188855 5.8222906033219024025 0.48705959451388997072 2.4662077676588296526 -3.7138302332409032225 4.1296166511454091008 -3.8396485313157837638 -5.3857672211301634846 -0.80511778892789021 -1.7098368026020789312 3.8926640467117108457 -8.2342525856963444397 -31.180073168660243255 0.19217230202967491848 -13.242994017930895367];

% Layer 3
b3 = -1.0074995939938131695;
LW3_2 = [-0.20077805371233717335 0.10492945843038987974 -0.041793717137972204168 -0.75103008550533623122 -2.4475490802599328966 0.013206459545162698224 -0.011798507498311839453 -0.32395958028297872078 -0.32681937143298089854 0.083154896850972895295 0.064845529538880300469 0.80502371009354700693 -0.59982809644452450559 -0.25668507363131320975 0.12547981932584370557];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00314743155299998;
y1_step1.xoffset = 1774.26765577007;

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
