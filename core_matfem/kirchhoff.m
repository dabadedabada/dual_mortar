function [C,S] = kirchhoff(rCG,lambda,mu)
% rCG ::= [rCG(1) rCG(4) rCG(6); rCG(4) rCG(2) rCG(5); rCG(6) rCG(5) rCG(3)]

S = zeros(6 ,size(rCG,2));
C = zeros(21,size(rCG,2));
J2 = rCG(1,:).*rCG(2,:).*rCG(3,:)+2*rCG(4,:).*rCG(5,:).*rCG(6,:)-rCG(1,:).*rCG(5,:).^2-rCG(3,:).*rCG(4,:).^2-rCG(2,:).*rCG(6,:).^2;
J = sqrt(J2); 
I = [rCG(1,:)+rCG(2,:)+rCG(3,:);   rCG(1,:).*rCG(2,:)+rCG(1,:).*rCG(3,:)+rCG(2,:).*rCG(3,:)-rCG(4,:).^2-rCG(5,:).^2-rCG(6,:).^2; J2];         
rCGinv = (1./J2).*[rCG(5,:).^2-rCG(2,:).*rCG(3,:); rCG(6,:).^2-rCG(1,:).*rCG(3,:); rCG(4,:).^2-rCG(1,:).*rCG(2,:); rCG(4,:).*rCG(3,:)-rCG(6,:).*rCG(5,:); rCG(1,:).*rCG(5,:)-rCG(4,:).*rCG(6,:); rCG(6,:).*rCG(2,:)-rCG(4,:).*rCG(5,:)];
al1 = (I(1,:).*lambda)./4 - mu./2 - (3.*lambda)./4;
al2 = mu./2;

S(1,:) = 2.*al1 + 2.*al2.*rCG(1,:);
S(2,:) = 2.*al1 + 2.*al2.*rCG(2,:);
S(3,:) = 2.*al1 + 2.*al2.*rCG(3,:);
S(4,:) = 2.*al2.*rCG(4,:);
S(5,:) = 2.*al2.*rCG(5,:);
S(6,:) = 2.*al2.*rCG(6,:);

ze1 = lambda./4;
ze7 = mu./2;

C( 1,:) = 4.*ze1 + 4.*ze7;
C( 7,:) = 4.*ze1;
C(12,:) = 4.*ze1;
C(16,:) = 0;
C(19,:) = 0;
C(21,:) = 0;
C( 2,:) = 4.*ze1 + 4.*ze7;
C( 8,:) = 4.*ze1;
C(13,:) = 0;
C(17,:) = 0;
C(20,:) = 0;
C( 3,:) = 4.*ze1 + 4.*ze7;
C( 9,:) = 0;
C(14,:) = 0;
C(18,:) = 0;
C( 4,:) = 2.*ze7;
C(10,:) = 0;
C(15,:) = 0;
C( 5,:) = 2.*ze7;
C(11,:) = 0;
C( 6,:) = 2.*ze7;
