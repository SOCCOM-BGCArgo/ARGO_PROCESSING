function [Y,Xf,Af] = ESPER_nitrate_7_Other_1(X,~,~)
%ESPER_NITRATE_7_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:34.
% 
% [Y] = ESPER_nitrate_7_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-138.904784718693];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0042772329208737];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.5824722406977742573;5.7192084967986698274;11.445513092821563816;5.5583652148419178118;2.7184843753575265168;4.6076935356338504235;1.3191073739807599452;1.1318598766880241246;2.4746855509378935523;-5.7894960923390845764;-3.6652629988532230421;1.1633311766244576191;-1.1505742015517175236;0.67106736529336175057;-2.3787475485808355913;-1.2291656695959758405;4.1522584053565356754;-0.091785622463348079392;2.7369216084813081658;-1.1199488480889256969;2.2186965124655766424;0.7862562915812583153;-2.2865558856332821591;1.947663870490969229;0.16960716282803159349;4.6748488557177454439;-6.2373486178399017987;-2.5747577583343601582;-1.4339979718475004411;2.0270114795181850553;-1.3342474075162169544;-2.3082026300945654995;1.540447042305381764;-1.3407754584805571252;-0.93231905793248381897;0.90550722912193637448;-2.762949724564163656;2.3206192999849686309;-1.9347569337910373033;2.5884690965878278668];
IW1_1 = [-0.10272094900104655757 -0.22656771757219373864 -1.4521254320921945791 0.27358589851321230002 -2.7346098458891110461 -1.4970502774342113739 2.371283188615659121;-3.3846750562652250238 -1.2596772597518164361 -7.0632760301551007487 0.93629243384845972642 2.0420347799670883759 -1.165912010744525551 -0.41653464560168268482;2.2032554848817209603 5.5653953250406056341 -2.4737993205801238616 3.3692821972516631135 2.1875533750462348159 3.5799790054374374115 0.85505422783871276593;-0.39475119699337557488 1.4673541586163310058 -2.2381693665305775198 4.3363723388866013408 -2.9466903906741346297 -2.6600013498962900016 0.14954811938279771732;-0.1775777881998401031 0.04586126771064855806 -0.85854576004150007229 0.78800568563210848794 3.9313172183973623675 1.4300573464666437307 0.37047776162250067866;-1.9375400748035023213 -2.1404602979291951215 -0.87431997031887753735 -0.17813399744118749157 -2.1408295897380700268 -2.8565184679204067564 -2.2140098657494462486;0.27149204018123351823 0.17284076610516288519 -1.4365896373596596103 0.45402260588805726238 -4.0049719829094518886 1.1499203375549842399 -0.98605581138008424613;-0.03603334043435181544 0.065114922443052017376 -0.55129147738564554526 0.67426923107132330504 0.18068500148888896883 -0.12381632986249732631 0.95730149740495773525;-2.5711026933948950557 -2.7372537293950092163 -0.51915141099859007934 0.14732056523467468301 -1.7892677723015335545 0.075863875094327237747 0.52174463796432546392;3.4726340636373667081 1.1767029647916449342 7.2976279576931943893 -1.0114073005997314958 -2.085730339738046446 1.0990910105902389482 0.30487461027927464929;0.26146141047737514462 0.14527067769689128163 -0.67650812667310611381 -0.90209600865574723905 -0.22422073467277736647 2.0819248260912344151 -2.7160216325292130257;0.13757623250838915685 0.26412205935867694162 -1.2807912918927641943 0.40894636116848265273 -3.4068172414646369361 1.0183792819014121367 -0.91066683458841435339;1.0849514484946101156 1.4598319201531799738 -0.12073035819602488672 0.75800957043059780105 -8.5892494263261838228 1.0165529427774870896 -0.0021455183838437232036;-0.15638968658529439626 -0.52745708605438801353 0.0021711137746106556889 -0.081674474331214386158 -0.45839082991700086378 -0.18243658426359984914 0.22714786389580593484;-0.99281154631444812342 -1.3703575470026585492 1.8771943197950893278 -1.0059917042196715453 1.9859192023758620582 -3.2143797274730365743 0.66099290030351087744;0.046592140720863123304 -0.094669815668830176181 -0.35439221779117302535 -0.035328299146913808015 1.8216632030555275357 -2.8920662415704896375 -0.26134309828567953149;-0.30627975584597172398 -0.1420122108259693483 0.72328361724498146934 1.2319523415259587829 0.29717176334325601239 -2.3683527638360346579 2.6942086022265887557;-0.61144680628766723274 1.9082057409245851698 -2.7276255786869558762 0.42610497473021335946 0.17592960440120991583 -2.0198200866200908266 0.62044511913218902688;-1.1874791925635499457 -1.713793470079300052 -1.7730938798011406377 0.82447270156728902357 -0.97643880866670607066 1.4647410240542939608 0.74743231792199404939;-1.4711938411382703418 0.71919453586123127664 1.2191654684548194343 0.076436674943651666814 -1.6736855828038748051 0.71042159616025735147 1.0227413922853436645;0.80478968154942776003 1.0460936463139656194 -1.6330095515451736876 0.88493019666261429901 -1.8650942257520546352 2.8224457468605832489 -0.68454706309757329397;1.0742180311241908974 -0.14302197715835693326 -1.9031815046264477509 0.4965654411443050642 -2.4562813854234524236 -1.418280866228511794 -0.70300769634238058359;-0.90195442860169039445 -1.2065972527131503256 1.7447600319912273115 -0.94561757665123791217 1.9987881695979132779 -2.9997937083282306681 0.67470005263667187645;0.22116920601001105462 0.30126056704558840105 -1.5087071198478043677 0.30248600605136044539 5.543141843882405162 0.74584597402230112806 -0.80564590063145380405;0.37256602695501728206 0.8006765599447527304 0.83403444971257045815 -0.36042500595937632113 -0.06936939195999668184 -0.67481323028535300512 -0.0088031204033802330777;-2.3503314984626961781 0.33053383902410560236 0.61617131478584985516 0.83581965115443712655 9.0550950768679729919 7.3520080329101515915 0.36755785695527631862;-0.004554897384579694071 -0.088825499238040256378 -0.56435278045026682214 -1.1296431005706730755 -4.4549328380105590597 -4.3521840718875610676 1.0231833792309124043;-19.383444934655894798 1.4214112821602373771 -7.2285447673719840367 2.0347121011154105652 -10.119274918863347068 2.5102561934011786704 -1.6906955732818516047;0.13272525943659499448 0.11842565601643353923 -0.87387795729073025619 -0.45421192135708410298 0.93276620645418317768 1.5966957856283652273 1.7305629966364115546;-0.050092340295204401679 0.56247113804236215096 0.44689510516026148546 -0.21791234622331193127 5.0390523163706681942 0.67056162938464125567 -0.15545818658203738871;0.039021431877336150429 -0.083377210999690365423 -0.34027237285775779885 -0.088374745441754373743 1.6706847156269288845 -2.8945892662526642169 -0.31398843996590669603;0.18961595052387847371 -0.29003910850917202779 0.2907629929341530306 -0.46156816877868878102 0.38352893212958527602 -1.6899928372770829377 0.16385139799061154697;-0.049232681733349374642 -0.68917552882028099148 0.57927988023844989574 -0.56634836571917901438 -3.7695923051035595819 3.2724576880398150713 -0.092688647894070261279;1.0554549670192370225 1.7839806205807593464 1.647820357747551645 -0.75570056445059707873 0.2613931967704439252 -1.3360592535911155476 -0.77953295537971800666;0.33780588585998583273 -0.19404822709522451252 0.54449499909002185838 -0.26142853883312494601 -0.27269357914031894241 -0.36000889447744011118 0.31923105101233240077;-0.10610601031378066872 0.37877007651513444619 -1.1767504793437173483 0.31653313848544450959 -2.3325049882367108189 0.63489173420670952375 -0.88066506622124873793;-0.53289151011777702305 -1.2499730739316614603 -3.9622566635162090876 0.54975482974673983705 3.864880092051508953 3.8732746366938992644 2.8165280374776862082;-0.126547119374488648 0.11145080687587333168 -0.63389309984654995667 0.70984399764233052643 3.4057917945523725756 0.9315732900782690562 0.20937130454513658284;0.032113695837217862139 -0.62913497141510643651 -0.59588796875984773038 0.27490144236176672354 -4.8994677318906276753 -0.53625085000144878666 0.13705122715264400624;0.11708965017489790517 0.23925872644710993731 1.5233759421408792356 -0.22906000560316217007 2.3680394026427955545 1.4532512805795279043 -2.4583946279129964552];

% Layer 2
b2 = -1.7386199345475219502;
LW2_1 = [4.4044692105489025025 -2.6817353964168706426 -0.16965038312948821364 -0.22623561101989492217 2.3810597851987265194 0.60033087507244109471 -3.0449551682639146577 1.4557441260542627326 0.10550263011450525252 -2.5727578133280712969 4.4128315833377040889 5.9392848165352392087 0.085237605892122231266 2.4886040628978638445 2.7821450108514880206 5.6316303040526483059 3.8273212740434887635 -0.12303301322724016775 -0.47221950525144684718 -0.29139747873228172281 -4.6606840217321581221 -0.27551460916505954302 -7.0220192848818117071 0.51449554138569797423 1.0279614759838511251 0.054988957609931610193 -0.66811132280829610419 -0.013626153344059879768 0.65520051030400972003 4.4323262558627156693 -5.485075842970096538 1.6348656664458109944 0.16181574504438547835 0.2306052209042238621 -2.9388337898855456665 -3.0970525936419628366 -0.36484498437387258285 -4.3201804563062298215 3.792251273994728944 3.9947473490039207];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0415843642790311;
y1_step1.xoffset = -0.9;

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
