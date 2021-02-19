function [nodes,elements,Np,Nelem]=TriPartition2DSW

model=createpde(1);
geometryFromEdges(model,@SolveDomain);
M=generateMesh(model,'JiggleIter',10,'Hgrad',1.8,'Hmax',0.5);
nodes=M.Nodes;
elements=M.Elements;
Np=length(nodes);
Nelem=length(elements);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=SolveDomain(bs,s)
xiDt=1;nbs=15;S=25.5000;D1=12.6500;D2=0.1;H1=0.01;H2=0.1;
if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0% start parameter value
  S/2 S/2 S (S/2)-H1 (S/2)-H1 S D1 H2 D2 H2 S-(D1+D2) S-(D1+D2) H2 D2 H1% end parameter value
  1 2 2 2 1 1 1 1 1 1 1 2 2 2 2% left hand region
  0 0 0 0 0 0 2 0 0 0 0 0 0 0 0% right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error(message('pde:lshapeg:InvalidBs'))
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
  error(message('pde:lshapeg:SizeBs'));
end

if ~isempty(s),

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=interp1([d(1,1),d(2,1)],[0 0],s(ii));
y(ii)=interp1([d(1,1),d(2,1)],[S/2 0],s(ii));
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=interp1([d(1,2),d(2,2)],[0 0],s(ii));
y(ii)=interp1([d(1,2),d(2,2)],[0 -S/2],s(ii));
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=interp1([d(1,3),d(2,3)],[0 S],s(ii));
y(ii)=interp1([d(1,3),d(2,3)],[-S/2 -S/2],s(ii));
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=interp1([d(1,4),d(2,4)],[S S],s(ii));
y(ii)=interp1([d(1,4),d(2,4)],[-S/2 -H1],s(ii));
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=interp1([d(1,5),d(2,5)],[S S],s(ii));
y(ii)=interp1([d(1,5),d(2,5)],[H1 S/2],s(ii));
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=interp1([d(1,6),d(2,6)],[S 0],s(ii));
y(ii)=interp1([d(1,6),d(2,6)],[S/2 S/2],s(ii));
end

%boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=interp1([d(1,7),d(2,7)],[0 D1],s(ii));
y(ii)=interp1([d(1,7),d(2,7)],[0 0],s(ii));
end

%boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=interp1([d(1,8),d(2,8)],[D1 D1],s(ii));
y(ii)=interp1([d(1,8),d(2,8)],[0 H2],s(ii));
end

%boundary segment 9
ii=find(bs==9);
if length(ii)
x(ii)=interp1([d(1,9),d(2,9)],[D1 D1+D2],s(ii));
y(ii)=interp1([d(1,9),d(2,9)],[H2 H2],s(ii));
end

%boundary segment 10
ii=find(bs==10);
if length(ii)
x(ii)=interp1([d(1,10),d(2,10)],[D1+D2 D1+D2],s(ii));
y(ii)=interp1([d(1,10),d(2,10)],[H2 H1],s(ii));
end

%boundary segment 11
ii=find(bs==11);
if length(ii)
x(ii)=interp1([d(1,11),d(2,11)],[D1+D2 S],s(ii));
y(ii)=interp1([d(1,11),d(2,11)],[H1 H1],s(ii));
end

%boundary segment 12
ii=find(bs==12);
if length(ii)
x(ii)=interp1([d(1,12),d(2,12)],[S D1+D2],s(ii));
y(ii)=interp1([d(1,12),d(2,12)],[-H1 -H1],s(ii));
end

%boundary segment 13
ii=find(bs==13);
if length(ii)
x(ii)=interp1([d(1,13),d(2,13)],[D1+D2 D1+D2],s(ii));
y(ii)=interp1([d(1,13),d(2,13)],[-H1 -H2],s(ii));
end

%boundary segment 14
ii=find(bs==14);
if length(ii)
x(ii)=interp1([d(1,14),d(2,14)],[D1+D2 D1],s(ii));
y(ii)=interp1([d(1,14),d(2,14)],[-H2 -H2],s(ii));
end

%boundary segment 15
ii=find(bs==15);
if length(ii)
x(ii)=interp1([d(1,15),d(2,15)],[D1 D1],s(ii));
y(ii)=interp1([d(1,15),d(2,15)],[-H2 0],s(ii));
end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RbW,CbW,Rbu,Cbu,Rbd,Cbd,Rbiu,Cbiu,Rbid,Cbid,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes)
xiDt=1;nbs=15;S=25.5000;D1=12.6500;D2=0.1;H1=0.01;H2=0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Rbu,Cbu]=find(nodes(2,:)>0 & abs(nodes(2,:))>=S/2);

[Rbd,Cbd]=find(nodes(2,:)<0 & abs(nodes(2,:))>=S/2);

[Rbiu,Cbiu]=find(nodes(2,:)>0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[Rbid,Cbid]=find(nodes(2,:)<0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[RbW,CbW]=find(abs(nodes(2,:))==0 & abs(nodes(1,:))<D1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[RbIv,CbIv]=find( (nodes(1,:)==D1 & abs(nodes(2,:))<H2) |...
                  (abs(nodes(2,:))==H2 & nodes(1,:)>=D1 & nodes(1,:)<=D1+D2) |...
                  (nodes(1,:)==D1+D2 & abs(nodes(2,:))<H2 & abs(nodes(2,:))>H1));

[Rv,Cv]=find(nodes(2,:)~=0 & abs(nodes(2,:))<S/2 & ...
              ( (abs(nodes(2,:))>H2) |... 
                (nodes(1,:)<D1 & abs(nodes(2,:))<=H2) | ...
                (abs(nodes(2,:))>H1 & nodes(1,:)>D1+D2 & nodes(1,:)<=S & abs(nodes(2,:))<=H2) ));

for l=1:length(Cbid) %% align the Cbiu and Cbid
    a(:,l)=nodes(:,Cbid(l));
    b(:,l)=nodes(:,Cbiu(l));
end  
aa(1,:)=sort(a(1,:),'ascend');bb(1,:)=sort(b(1,:),'ascend');
for l=1:length(aa)
    aaa(l)=Cbid(find(a(1,:)==aa(l)));
    bbb(l)=Cbiu(find(b(1,:)==bb(l)));    
end   
Cbid=aaa;Cbiu=bbb;
end
