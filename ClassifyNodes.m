function [RbW,CbW,Rbu,Cbu,Rbd,Cbd,Rbiu,Cbiu,Rbid,Cbid,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes)
xiDt=1;nbs=15;S=20.0;D1=9.9;D2=0.1;H1=0.01;H2=0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Rbu,Cbu]=find(nodes(2,:)>0 & abs(nodes(2,:))>=S/2);

[Rbd,Cbd]=find(nodes(2,:)<0 & abs(nodes(2,:))>=S/2);

[Rbiu,Cbiu]=find(nodes(2,:)>0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[Rbid,Cbid]=find(nodes(2,:)<0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[RbW,CbW]=find(abs(nodes(2,:))==0 & abs(nodes(1,:))<D1);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B]=ABMatix2D_SpinWave_U(Nelem,elements,nodes,A,B,alpha,gamma1)

%Coef2=6*gamma2*gamma2+(gamma1*gamma1+1);
%Coef1=3*gamma1*gamma1+(2*gamma2*gamma2+1);


for k=1:Nelem
    clc
    Nelem
    length(nodes(1,:))
    k
    clear Ek n1 n2 n3 x1 y1 x2 y2 x3 y3 xi1 xi2 xi3 eta xi omega D Si a2xlij a2ylij L1 L2 L3 a2lij A2lij
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ek=elements(:,k);

    n1=nodes(:,Ek(1));x1=n1(1);y1=n1(2);
    n2=nodes(:,Ek(2));x2=n2(1);y2=n2(2);
    n3=nodes(:,Ek(3));x3=n3(1);y3=n3(2);
    
    xi=[x2-x3 x3-x1 x1-x2];
    eta=[y2-y3 y3-y1 y1-y2];
    omega=[x2*y3-x3*y2 x3*y1-x1*y3 x1*y2-x2*y1];
    
    Di=sum(omega);
    Si=abs(Di)/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% gamma2 for different domain %%%%
 
    gamma2=gamma2SD(y1,y2,y3,gamma1);

    Coef1=3*gamma1*gamma1+(2*gamma2*gamma2+1);
    Coef2=6*gamma2*gamma2+(gamma1*gamma1+1);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    akji=zeros(3,3);
    for o=1:3
        akji(o,1)=Coef1*(1/(1*(Di.^2)))*eta(o)*eta(1)*Si+Coef2*(1/(1*(Di.^2)))*xi(o)*xi(1)*Si+...
	    Integral2D_kvtu(o,1,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1)+...
                  Integral2D_kvUiu(o,1,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1,gamma2);

        akji(o,2)=Coef1*(1/(1*(Di.^2)))*eta(o)*eta(2)*Si+Coef2*(1/(1*(Di.^2)))*xi(o)*xi(2)*Si+...
	    Integral2D_kvtu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1)+...
                  Integral2D_kvUiu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1,gamma2);

        akji(o,3)=Coef1*(1/(1*(Di.^2)))*eta(o)*eta(3)*Si+Coef2*(1/(1*(Di.^2)))*xi(o)*xi(3)*Si+...
	    Integral2D_kvtu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1)+...
                  Integral2D_kvUiu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1,gamma2);

	
        % akji(o,2)=(1/(1*(Di.^2)))*eta(o)*eta(2)*Si+(1/(1*(Di.^2)))*xi(o)*xi(2)*Si+...
	%   Integral2D_kvtu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha)+...
        %           Integral2D_kvUiu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha);
              
        % akji(o,3)=(1/(1*(Di.^2)))*eta(o)*eta(3)*Si+(1/(1*(Di.^2)))*xi(o)*xi(3)*Si+...
	%   Integral2D_kvtu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha)+...
        %           Integral2D_kvUiu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha);
            
    end  
    
    bkji=zeros(3,3);
    for o=1:3
        
        for f=1:3
            bkji(o,f)=Integral2D_kvu(o,f,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di);
        end            
    end
    
   %%%%%%%%%%%%matrix symmetrized%%%%%%%%%
   Ak=(1/2)*(akji+akji');
   Bk=(1/2)*(bkji+bkji');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
      for u=1:3 %% The lth element locaizes inside Domain, us the original index 
  
          A(Ek(u),Ek(1))=Ak(u,1)+A(Ek(u),Ek(1));
          A(Ek(u),Ek(2))=Ak(u,2)+A(Ek(u),Ek(2));
          A(Ek(u),Ek(3))=Ak(u,3)+A(Ek(u),Ek(3));
          
          B(Ek(u),Ek(1))=Bk(u,1)+B(Ek(u),Ek(1));
          B(Ek(u),Ek(2))=Bk(u,2)+B(Ek(u),Ek(2));
          B(Ek(u),Ek(3))=Bk(u,3)+B(Ek(u),Ek(3));
      end

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ivu=Integral2D_kvu(o,f,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       integrant=@(x,y) ((eta(o).*x-xi(o).*y+omega(o))./Di).*((eta(f).*x-xi(f).*y+omega(f))./Di);
       Ivu=TriEleIntegral2D(integrant,x1,x2,x3,y1,y2,y3);
     
end