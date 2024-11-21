%Retrieve directory that the file will go into
userDir = getenv('USERPROFILE');
userDir = userDir+"\Documents\";

%Identify the directory of the desired file on the GDAC
webtarget = "https://data-argo.ifremer.fr";
web_dir  = "/";
fname = "argo_synthetic-profile_index.txt";

%Full name of the file url
dataUrl = webtarget+web_dir+fname;

%Name of the file to be put into your local directory
sProfFile = userDir+"argo_synthetic-profile_index.txt";

%Save the file into userDir as sProfFile
fullPath = websave(sProfFile,dataUrl);







