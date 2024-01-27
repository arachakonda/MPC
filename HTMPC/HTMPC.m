clc;
clear all;
close all;
%% System data
n=2;m=1;

A0=[1 1; 0.3 1.12];
B0=[0.65;1];
A=[1 1;0 1.15];
B=[0.5;1];
Q=eye(n); R=0.1;
W_ric=Q;
Qsqrt=sqrtm(Q);Rsqrt=sqrtm(R);Wsqrt=sqrtm(W_ric);

K=[-0.5 -1.2999];

Phi=A0+B0*K;
Rfor=11;
%% Constraints
xminf=15;uminf=20;

F=[0 1/xminf;1/xminf 0;0 -1/xminf;-1/xminf 0;0 0;0 0];
G=[0;0;0;0;1/uminf;-1/uminf];
X=Polyhedron(F(1:2*n,:),ones(2*n,1));
U=Polyhedron(G(2*n+1:end,:),ones(2*m,1));

FGK=F+G*K;
nsFGK=size(FGK);
nr=nsFGK(1);
OneM=ones(nr,1);

%% constants
denom=1+xminf^2+uminf^2;
adaptl=0.1;

%% Terminal set
% % XtF=[0 1/1;1/1 0;0 -1/1;-1/1 0];
% % XtG=[1;1;1;1];
% % PhXt=Polyhedron(XtF,XtG);
% % Xtermhalf=PhXt.minHRep;Xtermhalf=Xtermhalf.H;
% % figure;
% % plot(PhXt);

%% W_theta and W_xtilde
Wt1=0.01;
Wx1=0.1;

%% Tube cross section
nDJ=8;
TubeD_V=[0.5 1;1 0.5;1 -0.5;0.5 -1;-0.5 -1;-1 -0.5;-1 0.5;-0.5 1]; %vertices of D
TubeD=Polyhedron(TubeD_V); % polyhedron D

D_V=TubeD_V'; %vertices of D in column [x(1) x(2) x(3)]

disp('  ');
x0k=input('Enter x0|k as a column matrix = ');
x0kt=x0k';
XXT2=[];

N_again=1;
% % % if contains(PhXt,x0k)==1
% % %     N=0;J=x0k'*W_ric*x0k;disp(' is not the correct cost');
% % % else
% % %     N=input('Enter prediction horizon length, N(>0) = ');
% % % end
N=input('Enter prediction horizon length, N(>0) = ');
xtilda=[];Thetat=norm([A B]-[A0 B0],2);
while N_again==1
    
    ustore=[];xstore1=x0k(1);xstore2=x0k(2);
    
    for ki=1:1:Rfor
        if ki>1
            Y0=[x0k_Y;uka];
            [Anew,Bnew]=ABcall(A0,B0,adaptl,xpt,Y0,denom);
            A0=Anew;B0=Bnew;
            Thetat=[Thetat norm([A B]-[A0 B0],2)];
        end
        disp('k= ');disp(ki);
        disp('A0= ');disp(A0);
        disp('B0= ');disp(B0);
        
        
        cvx_begin
        variable Tube_c((N+1)*n)
        variable Tube_alpha(N+1)
        variable Utube_v(nDJ*N*m)
        XTube_v=[];
        for Tube_i1=1:1:N+1
            for Tube_i2=1:1:nDJ
                XTube_vt=Tube_c(2*Tube_i1-1:2*Tube_i1)+...
                    Tube_alpha(Tube_i1)*D_V(:,Tube_i2);
                XTube_v=[XTube_v XTube_vt]; %contains x11 x12... x14 x21 x22...xN4
            end
        end
        
        %             utemp=Utube_v(nDJ+1:end);
        %needed later while implementing
        % % %     UThalfA=[];UThalfb=[];
        % % %     for Tube_i4=1:1:N-1
        % % %         vtemp=utemp(:,nDJ*(Tube_i4-1)+1:nDJ*Tube_i4);vtemp=vtemp';
        % % %        [Atemp,btemp]=vert2con(vtemp);
        % % %        UThalfA=[UThalfA;Atemp];
        % % %        UThalfb=[UThalfb;btemp];
        % % %     end
        %             TubeTheta=[Tube_alpha;Tube_c;Utube_v];
        
        % state with disturbance
        xnext1=[];
        xnext2=[];
        xnext3=[];
        xnext4=[];%A0*x0k+B0*Utube_v(1)+[-1;1]*Wx1;
        Tube_i51=1;i52=0;
        for Tube_i5=1:1:N*nDJ
            if Tube_i5>=1 || Tube_i5<=nDJ
                xn1t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[1;1]*(Wx1);
                xn2t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[-1;-1]*(Wx1);
                xn3t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[1;-1]*(Wx1);
                xn4t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[-1;1]*(Wx1);
                xnext1=[xnext1 xn1t];xnext2=[xnext2 xn2t];
                xnext3=[xnext3 xn3t];xnext4=[xnext4 xn4t];
            else
                xn1t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[1;1]*((Tube_i51-1)*Wt1+Wx1);
                xn2t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[-1;-1]*((Tube_i51-1)*Wt1+Wx1);
                xn3t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[1;-1]*((Tube_i51-1)*Wt1+Wx1);
                xn4t=A0*XTube_v(:,Tube_i5)+B0*Utube_v(Tube_i5)+[-1;1]*((Tube_i51-1)*Wt1+Wx1);
                xnext1=[xnext1 xn1t];xnext2=[xnext2 xn2t];
                xnext3=[xnext3 xn3t];xnext4=[xnext4 xn4t];
                i52=i52+1;
                if i52==nDJ
                    Tube_i51=Tube_i51+1;
                    i52=0;
                end
            end
        end
        Jv=0;
        if N>=1
            for iJ=1:1:N*nDJ
                Jv=Jv+norm(Qsqrt*XTube_v(iJ))+norm(Rsqrt*Utube_v(iJ));
            end
        end
        for iJ2=nDJ*(N-1)+1:1:nDJ*N
            Jv=Jv+norm(Wsqrt*XTube_v(iJ2));
        end
        J=Jv;
        %             J=Jv+norm(Qsqrt*x0k)+norm(Rsqrt*Utube_v(1));
        Tube_i8=1;Tube_i81=1;
        
        minimize J
        
        subject to
        -Tube_alpha<=0;
        Utube_v<=uminf*ones(nDJ*(N)*m,1);
        -Utube_v<=uminf*ones(nDJ*(N)*m,1);  %u in U
        %   XNbarH(:,1:n)*x0k<=XNbarH(:,end); % x0k is in ROA ??
        F(1:4,:)*XTube_v<=ones(4,(N+1)*nDJ);   %x is in X
        % % %             Xtermhalf(:,1:2)*XTube_v(:,(N-1)*nDJ+1:N*nDJ)<=[Xtermhalf(:,end)...
        % % %                 Xtermhalf(:,end) Xtermhalf(:,end) Xtermhalf(:,end)  Xtermhalf(:,end)...
        % % %                 Xtermhalf(:,end)  Xtermhalf(:,end)  Xtermhalf(:,end)];  % xN is in terminal set mRPI
       
        % x0k in X0
        norm(x0k-Tube_c(1:2),inf)<=Tube_alpha(1);
        norm(x0k-Tube_c(1:2),1)<=1.5*Tube_alpha(1);
        
        % inclusion
        %             norm(xnext1(:,1)-Tube_c(1:2),inf)<=Tube_alpha(1);
        %             norm(xnext2(:,1)-Tube_c(1:2),inf)<=Tube_alpha(1);
        %             norm(xnext3(:,1)-Tube_c(1:2),inf)<=Tube_alpha(1);
        %             norm(xnext4(:,1)-Tube_c(1:2),inf)<=Tube_alpha(1);
        %
        for Tube_i6=1:1:N
            for Tube_i7=1:1:nDJ
                norm(xnext1(:,Tube_i8)-Tube_c(2*Tube_i6-1:2*Tube_i6),inf)<=Tube_alpha(Tube_i6);
                norm(xnext2(:,Tube_i8)-Tube_c(2*Tube_i6-1:2*Tube_i6),inf)<=Tube_alpha(Tube_i6);
                norm(xnext3(:,Tube_i8)-Tube_c(2*Tube_i6-1:2*Tube_i6),inf)<=Tube_alpha(Tube_i6);
                norm(xnext4(:,Tube_i8)-Tube_c(2*Tube_i6-1:2*Tube_i6),inf)<=Tube_alpha(Tube_i6);
                Tube_i8=Tube_i8+1;
            end
        end
        
        %             norm(xnext1(:,1)-Tube_c(1:2),1)<=1.5*Tube_alpha(1);
        %             norm(xnext2(:,1)-Tube_c(1:2),1)<=1.5*Tube_alpha(1);
        %             norm(xnext3(:,1)-Tube_c(1:2),1)<=1.5*Tube_alpha(1);
        %             norm(xnext4(:,1)-Tube_c(1:2),1)<=1.5*Tube_alpha(1);
        %
        for Tube_i61=1:1:N
            for Tube_i71=1:1:nDJ
                norm(xnext1(:,Tube_i81)-Tube_c(2*Tube_i61-1:2*Tube_i61),1)<=1.5*Tube_alpha(Tube_i61);
                norm(xnext2(:,Tube_i81)-Tube_c(2*Tube_i61-1:2*Tube_i61),1)<=1.5*Tube_alpha(Tube_i61);
                norm(xnext3(:,Tube_i81)-Tube_c(2*Tube_i61-1:2*Tube_i61),1)<=1.5*Tube_alpha(Tube_i61);
                norm(xnext4(:,Tube_i81)-Tube_c(2*Tube_i61-1:2*Tube_i61),1)<=1.5*Tube_alpha(Tube_i61);
                Tube_i81=Tube_i81+1;
            end
        end
        
        cvx_end
        if isnan(cvx_optval)
            %   disp('*** Increase value of N ***');
            N_again=1;
            N=N+1;
            break
        elseif cvx_optval==Inf
            % disp('*** Increase value of N ***');
            N_again=1;
            N=N+1;
            break
            
        elseif cvx_optval==-Inf
            %  disp('*** Increase value of N ***');
            N_again=1;N=N+1;
            break
        else
            N_again=0;
            if ki==1
                TUBE1=XTube_v(:,1:nDJ);
                
                TUBE1=TUBE1';
                TUBE1=Polyhedron(TUBE1);
                XXT2=[XXT2;TUBE1];
            end
            TUBE1=XTube_v(:,nDJ+1:2*nDJ);
            TUBE1=TUBE1';
            TUBE1=Polyhedron(TUBE1);
            XXT2=[XXT2;TUBE1];
            
            % has to be changed. take u as a policy element
                cvx_begin
                variable Lambda0(nDJ)
                Jmin=[0;0];
                for Lambdac=1:1:nDJ
                    Jmin=Lambda0(Lambdac)*XTube_v(:,Lambdac)+Jmin;
                end

                minimize norm(x0k-Jmin,2)
                subject to
                Lambda0<=ones(nDJ,1);
                -Lambda0<=zeros(nDJ,1);
                sum(Lambda0)==1;

                cvx_end
            ukat=0;
            for Ulambda=1:1:nDJ
                ukat=ukat+Lambda0(Ulambda)*Utube_v(Ulambda);
            end
            uka=ukat;
            ustore=[ustore;uka];
            x0k_Y=x0k;
            x0khat=A0*x0k+B0*uka;
            x0k=A*x0k_Y+B*uka;
            xpt=x0k-x0khat;
            xstore1=[xstore1;x0k(1)];
            xstore2=[xstore2;x0k(2)];
            xtilda=[xtilda norm(xpt,2)];
        end
        
        
    end
end

disp('     ');
disp(' N = ');disp(N);
disp(' J* = ');disp(J);
t=0:1:Rfor-1;
xs1=xstore1(1:Rfor);xs2=xstore2(1:Rfor);
figure
plot(t,xs1,'-r');
hold on;
grid on;
plot(t,xs2,'-g');
hold on;
plot(t,ustore,'-*b');

%% figure for X, XT, trajectory, tubes
figure
plot3(xs1,t,xs2,'-k','LineWidth',1.5,'Marker','*');
hold on
grid on;
set(gca,'Ydir','reverse');



%         Xbarv=Xbar.V;Xbarv1=Xbarv(:,1);
%         Xbarv2=Xbarv(:,2);
%         Xbarc=convhull(Xbarv1,Xbarv2);
%         Xbarv1=Xbarv1(Xbarc);Xbarv2=Xbarv2(Xbarc);
%         lXbar=length(Xbarv1);

Xv=X.V;Xv1=Xv(:,1);Xv2=Xv(:,2);
Xc=convhull(Xv1,Xv2);
Xv1=Xv1(Xc);Xv2=Xv2(Xc);lX=length(Xv1);

for iplotc=1:1:length(t)
    
    XXT1=XXT2(iplotc);
    
    
    Xbarv=XXT1.V;Xbarv1=Xbarv(:,1);Xbarv2=Xbarv(:,2);
    
        
        Xbarc=convhull(Xbarv1,Xbarv2);
    
    Xbarv1=Xbarv1(Xbarc);Xbarv2=Xbarv2(Xbarc);
    lXbar=length(Xbarv1);
    
    patch(Xbarv1,(iplotc-1)*ones(1,lXbar),...
        Xbarv2,'y','EdgeColor','k','FaceAlpha',1);
    
    patch(Xv1,(iplotc-1)*ones(1,lX),...
        Xv2,'c','EdgeColor','b','FaceAlpha',.1);
    
    
end
% % iplotc=N+1;
% % XXT1=XXT(iplotc);Xbarv=XXT1.V;Xbarv1=Xbarv(:,1);Xbarv2=Xbarv(:,2);
% % Xbarc=convhull(Xbarv1,Xbarv2);
% % Xbarv1=Xbarv1(Xbarc);Xbarv2=Xbarv2(Xbarc);
% % lXbar=length(Xbarv1);
% % 
% % patch(Xbarv1,(iplotc-1)*ones(1,lXbar),...
% %     Xbarv2,'y','EdgeColor','k','FaceAlpha',1);
ylim([0 Rfor+1]);zlim([-xminf-0.5 xminf+0.5]);
xlim([-xminf-0.5 xminf+0.5]);
hold off
%% figure XT, trajectory
% % % figure
% % % plot3(xs1,t,xs2,'-k','LineWidth',2,'Marker','^');
% % % hold on
% % % grid on;
% % % set(gca,'Ydir','reverse');
% % %
% % % Zv=PhXt.V;Zv1=Zv(:,1);Zv2=Zv(:,2);
% % % Zc=convhull(Zv1,Zv2);
% % % Zv1=Zv1(Zc);Zv2=Zv2(Zc);lZ=length(Zv1);
% % %
% % % for iplotc=1:1:length(t)
% % %     patch(Zv1,(iplotc-1)*ones(1,lZ),...
% % %         Zv2,'b','EdgeColor','k','FaceAlpha',0.2);
% % %
% % % end
% % % ylim([0 15]);zlim([-11 11]);xlim([-11 11]);
% % % hold off


%% ki korbe
% xtilda norm and theta tilda norm plot
% plot U Tube
% separate code for (i-1)* W
%
%%

figure
plot(t,xtilda,'-r');
hold on;
grid on;
plot(t,Thetat,'-*b');