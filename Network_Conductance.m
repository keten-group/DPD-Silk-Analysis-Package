A=load('cluster_connectivity.txt');
B=load('cluster_coordinate.txt');
L=load('boundary_condition.txt');
Lx=L(1);

[m]=size(A,1);
NetList=zeros(1,3);
Vnod=zeros(1,2);

totconnections=0;
newnodes=m;
Vcount=0;

for i=1:m
    for j=i+1:m
        if (A(i,j)>0)
            test=abs(B(i,1)-B(j,1))< Lx/2;
            if test>0
                totconnections=totconnections+1;
                NetList(totconnections,:)=[i j 1/A(i,j)];
            else                
                left_id=(B(i,1)<0)*i+(B(j,1)<0)*j;
                right_id=(B(i,1)>0)*i+(B(j,1)>0)*j;

                totconnections=totconnections+1;
                Vcount=Vcount+1;
                newnodes=newnodes+1;                
                NetList(totconnections,:)=[newnodes left_id 0.5/A(i,j)];
                Vnod(Vcount,:)=[newnodes 1];
                
                totconnections=totconnections+1;
                Vcount=Vcount+1;
                newnodes=newnodes+1;
                NetList(totconnections,:)=[newnodes right_id 0.5/A(i,j)];
                Vnod(Vcount,:)=[newnodes 0];
            end
        end
    end
end

%break;

    %disp('test if any nanowire touches left and right electrodes');
    v_left_index=find(Vnod(:,2)==0); 
    if (length(v_left_index)==0)
        disp('Not Touching Left Electrode');
        %break;
    end
    for i=1:length(v_left_index)-1
        totconnections=totconnections+1;
        NetList(totconnections,:)=[Vnod(v_left_index(i),1) Vnod(v_left_index(i+1),1) 0];
    end

    v_right_index=find(Vnod(:,2)==1);
    if (length(v_right_index)==0)
        %disp('Not Touching Right Electrode');
        %break;
    end
    for i=1:length(v_right_index)-1
        totconnections=totconnections+1;
        NetList(totconnections,:)=[Vnod(v_right_index(i),1) Vnod(v_right_index(i+1),1) 0];
    end
    % sort NetList
    [B,IX]=sort(NetList(:,2));
    NetList=NetList(IX,:);
    [B,IX]=sort(NetList(:,1));
    NetList=NetList(IX,:);
    zero_index=find(NetList(:,3)==0);

    %NetList

    %disp('beautify resistance graph');
    % get rid of zero resistance connections
    for i=1:length(zero_index)
        indd=find(NetList(:,1)==NetList(zero_index(i),2));
        for j=1:length(Vnod(:,1))
            if Vnod(j,1)==NetList(zero_index(i),2)
                Vnod(j,1)=NetList(zero_index(i),1);
            end
        end
        NetList(indd,1)=NetList(zero_index(i),1)*ones(length(NetList(indd,1)),1);    

        indd=find(NetList(:,2)==NetList(zero_index(i),2));
        NetList(indd,2)=NetList(zero_index(i),1)*ones(length(NetList(indd,2)),1);    

        NetList(zero_index(i),2)=NetList(zero_index(i),1);
    end


    % sort again before compute
    for i=1:totconnections
        NetList(i,:)=...
            [min(NetList(i,1),NetList(i,2)) max(NetList(i,1),NetList(i,2)) NetList(i,3)];
    end
    [B,IX]=sort(NetList(:,1));
    NetList=NetList(IX,:);

% compute resistance
    l=size(NetList,1);
    N=max([NetList(:,1); NetList(:,2)]);
    A=zeros(N,N);
    B=zeros(N,1);

    for i=1:l
        n1=NetList(i,1);
        n2=NetList(i,2);
        if n1==n2
        else
            A(n1,n2)=A(n1,n2)-1/NetList(i,3);
            A(n2,n1)=A(n2,n1)-1/NetList(i,3);
            A(n1,n1)=A(n1,n1)+1/NetList(i,3);
            A(n2,n2)=A(n2,n2)+1/NetList(i,3);
        end
    end
    A0=A;

    if (N>0) %<-- add case: if N=0, kappa=0
        for i=1:size(Vnod,1)
            A(Vnod(i,1),:)=zeros(1,N);
            A(Vnod(i,1),Vnod(i,1))=1;
            B(Vnod(i,1),1)=Vnod(i,2);
        end 

        %%% Part VI. solve for the resistance %%%
        %disp('solving AX=B equation via singular value decomposition');
        V=pinv(A)*B;

        Right_index=find(B==1);
        indices=find(1:length(B)~=Right_index);
        Conductance=A0(Right_index,indices)*(V(indices)-V(Right_index));
        Resistance=1/Conductance;
    else
        Conductance=0.0;
    end

    fout=fopen('conductance.txt','w');
    fprintf(fout, '%12f', Conductance);
    fclose(fout);