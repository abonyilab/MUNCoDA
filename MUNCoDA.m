%% Processing of Voluntary National Reviews for network analysis
% Janos Abonyi, Viktor Sebestyen, 11/02/2020

close all
clear all



K=0.5;
tr_nmin=5
tr_nmax=75; %do not use freq. words
tr_tfidf=0.4;

%% Read the knime CSV file

if 1
    M = readtable('data_term_csv.csv');
else 
  load M %saved as table to speed up
end  

%% filter the words - based on characters / remove numbers, etc. 
ascfind=[65 90;%asccapit 48 57; %ascnumber
        97 122]; %vascsmall
delindex=[];
for i=2:3
    t = table2array(M(:,i));
    %coltext = regexp(t, '\w','match','once');
    dum=double(char(t));
    [id1,idd]=find(dum==45); %remove "-" 
    [id2,idd]=find(dum==39); %remove "-" 
    [id3,idd]=find(dum==44); %remove "-" 
    [id4,idd]=find(dum==96); %remove "-" 
    ok=zeros(size(dum,1),1);
    for k=1:size(ascfind,1)
        dok=((dum(:,1)>=ascfind(k,1)) & (dum(:,1)<=ascfind(k,2)));
        ok=max(ok,dok);
    end 
    dumind=find(ok==0);
    delindex=[delindex; dumind; id1; id2; id3; id4] ; %find(ismember(t,words_common)); find(ismember(t,words_unique))]
end
M(delindex,:)=[];

%% Generate list of country names (xxx.pdf -> xxx)
Countriestext = table2array(M(:,1));
Countries = cellfun(@(Countriestext) Countriestext(1:end-3), Countriestext, 'Uniform', 0);
%listcountries={'aut','bel','bgr','cyp','cze','dnk','gbr','est','fin','fra','grc','nld','hrv','irl','pol','lva','ltu','lux','hun','mlt','deu','ita','prt','rou','esp','swe','svk','svn'}
Countries = upper(Countries);
listcountries=unique(Countries);
listcountries=upper(listcountries);
%xlswrite('listcountries.xlsx',listcountries)
N=length(listcountries);


%% Calculate the word frquencies, select the improtant words
%TF=readtable('results_from_knime.xlsx','Sheet','term_frequency','ReadVariableNames',false);
load TF 
TF_n=table2array(TF(:,3));
TF_countries=table2array(TF(:,2));
TF_countries = cellfun(@(TF_countries) TF_countries(1:end-3), TF_countries, 'Uniform', 0);
TF_countries=upper(TF_countries);
TF_words=table2array(TF(:,1));
ascfind=[        65 90;%asccapit 48 57; %ascnumber
          97 122]; %vascsmall
delindex=[];
%coltext = regexp(t, '\w','match','once');
dum=double(char(TF_words));
ok=zeros(size(dum,1),1);
for k=1:size(ascfind,1)
    dok=((dum(:,1)>=ascfind(k,1)) & (dum(:,1)<=ascfind(k,2)));
    ok=max(ok,dok);
end
delindex=find(ok==0);
TF_n(delindex,:)=[];
TF_n=min(TF_n,N); %valami francia bekeveredett, 
TF_countries(delindex,:)=[];
TF_words(delindex,:)=[];
% calculate the num. of supporting documents
[words,index,i3]=unique(TF_words);
w_num=accumarray(i3,1); %Total number of occurence of words
sel=find(and(w_num>10,w_num<70));

words_sel=words(sel);

%% select the word-pairs based on N_t and filter the words 
[wp,dum,gi]=unique(M(:,2:3),'rows','stable');
Nt=accumarray(gi(:), 1);
index=(find(Nt>=tr_nmin));
wp=wp(index,:);
Nt=Nt(index);
index=(find(Nt<=tr_nmax));
wp=wp(index,:);
Nt=Nt(index);
wp(find(Nt==75),:)
word_list=table2array(wp);
word_list=unique(word_list(:));

%word_list=intersect(word_list,words_sel);

%% Deterime the informative word pairs and code the net
edges=[];
[im,countries_loc]=ismember(Countries,listcountries);
edges=countries_loc;
for i=2:3
    t = table2array(M(:,i));
    [im,words_loc]=ismember(t,word_list);
    edges=[edges words_loc];
end
edges=[edges table2array(M(:,4))];


remind=unique([find(edges(:,2)==0);find(edges(:,3)==0)]);
edges(remind,:)=[];
Countries(remind,:)=[];
listcountries(edges(end,1)) %OK
word_list(edges(end,2:3)) %OK
%% calc n_t after pruning 
dum=edges(:,2:3);
[word_pairs,iwp,iied]=unique(dum,'rows','stable');
Nt=accumarray(iied,1); %Total number of occurence of word-pairs


%% Data structure based on the list of countries 
%from, to, t, n_t, tf, idf, tfidf
X={};
XX=[];
for c=1:N;
    index=find(strcmp(cellstr(Countries(:,1)),listcountries(c)));
    X{c}=[edges(index,:) Nt(iied((index)))]; 
    t=X{c}(:,4);
    tn=X{c}(:,5);
    tf=K+(1-K)*t/max(t); %term frequency
    X{c}(:,6)=tf;
    idf = log10(N./tn); %OK.
    tfidf = tf .* idf;  
    X{c}(:,7)=tfidf;
    XX=[XX;X{c}];
end

%% Pruning the data 

%index=find(XX(:,5)>=5); %Min. number of the supporting countries 
%XX=XX(index,:);

index=find(XX(:,7)>=tr_tfidf); %remove the not informative word-pairs
XX=XX(index,:);

dum=[XX(:,2) XX(:,3)];
[oindex,nindex,pindex]=unique(dum);
rindex=reshape(pindex,size(dum));
XXnew=XX;
XXnew(:,[2:3])=rindex;

%%

figure(200)
subplot(3,1,1)
hist(XXnew(:,6),30)
ylabel('Number of {\it t} word pairs with a given {\it tf}')
xlabel('{\it tf} measure') 


subplot(3,1,2)
hist(XXnew(:,5),30)
ylabel('Number of {\it t} word pairs with a given {\it n_t}')
xlabel('{\it n_t} number of documents (VNRs) with a given {\it t} word pair') 

subplot(3,1,3)
hist(XXnew(:,7),30)
ylabel('Number of {\it t} word pairs with a given {\it tf-idf}')
xlabel('{\it tf-idf} measure') 


%% Generate sparse adjacency matrices 
DX=[];
Nnode=length(word_list); % max(max(XXnew(:,[2,3])));
network_layer={};
for c=1:N;
    index=find(XXnew(:,1)==c);
    network_layer{c}.A=sparse(XXnew(index,2),XXnew(index,3),XXnew(index,6),Nnode,Nnode); %4: number, 7; tf-idf 
    network_layer{c}.A=network_layer{c}.A+network_layer{c}.A'; %make it symetric
    network_layer{c}.G=graph(network_layer{c}.A,word_list); % (oindex)
    network_layer{c}.C = centrality(network_layer{c}.G, 'eigenvector');
end    


%% % projection based on the max. 
network_layer{N+1}.A=network_layer{1}.A; %N+1: TOTAL
for c=2:N
    network_layer{N+1}.A=max(network_layer{N+1}.A,network_layer{c}.A);
end



%% draw a network
trplot=0.65

%for ii=1:length(listcountries);
 %   ccode=listcountries{ii}
%il=find(ismember(listcountries,ccode))

%il=21% N+1 %21-DEU 
il=N+1


A=(full(network_layer{il}.A)>trplot).*(network_layer{il}.A);

selnode=find(sum(A)>0); %just the connected

G=graph(A(selnode,selnode),word_list(selnode));
figure(500)
clf
HG=plot(G,'Layout', 'force','NodeLabel',word_list(selnode))

writetable(G.Edges,'WholeG_v2.xlsx'); %'UseExcel',false
C=centrality(G, 'eigenvector');
[dum,ic]=sort(C,'descend')

names=word_list(selnode)

%graphtogml([ccode 'all.gml'], A(selnode,selnode),names)
%end
graphtogml('TELJES_065.gml', A(selnode,selnode),names)

%% NODE similarity ... 
tic
for c=1:N+1
c
S=network_layer{c}.A;
for k =1:N
    for i =1:N
        for j =1:N
            if S(i,j)< S(i,k)*S(k,j)
                S(i,j)=S(i,k)*S(k,j);
            end 
        end
    end
end  
S=S-diag(diag(S))+eye(length(S));
Sim{c}=S;
end
toc


%% Plot the  similarity of the projected / summarized layer
figure(240)
clf
indp=find(sum((network_layer{N+1}.A))>60);

clf
S=Sim{N+1};
Y = cmdscale(S,2); %+randn(length(S),2)*0.01
plot(Y(indp,1),Y(indp,2),'w.')
%plot(Y(:,1),Y(:,2),'w.')
hold on 
%text(Y(:,1),Y(:,2),word_list(oindex(indp)),'Fontsize',11);
text(Y(indp,1),Y(indp,2),word_list(oindex(indp)),'Fontsize',10);


%% Similarity of the layers - based on similarities 
L_sim=zeros(N,N);
for i=1:N
    for j=i+1:N
        maxsim=max(Sim{i}(:),Sim{j}(:));
        minsim=min(Sim{i}(:),Sim{j}(:));
        ind=find(maxsim>0);
        L_sim(i,j)=mean(minsim(ind)./maxsim(ind));
    end
end  
L_sim=L_sim+L_sim';
L_sim=L_sim+eye(N);


%%
figure(250)
clf
Y =  mdscale(L_sim,2); %,'Start','random' ,'criterion','sammon'
plot(Y(:,1),Y(:,2),'w.')
hold on 
text(Y(:,1),Y(:,2),listcountries,'Fontsize',8);



%% clustering based on degree correllation Country
%D = pdist(K','pearson');
CK=1-L_sim;
ZCK = linkage(CK,'complete');
figure(70)
[H,T,outperm_ck]= dendrogram(ZCK,0,'Orientation','left','Labels',listcountries); %
figure(71)
h=heatmap(listcountries(outperm_ck),listcountries(outperm_ck),CK(outperm_ck,outperm_ck));

%%


K=[];
for c=1:N;
    K=[K network_layer{c}.C];
end
CK=corr(K, 'type','Kendall');

