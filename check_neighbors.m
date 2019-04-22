function [index_temp,temp,boundary,minID_matrix,list,TS_matrix,cross_matrix,tunnel_list,TS_list,tunnel_cluster,tunnel_cluster_dim]=...
    check_neighbors(i,j,k,in,jn,kn,level,iC,index_list,index_temp,temp,boundary,...
    minID_matrix,list,TS_matrix,level_matrix,E_matrix,cross_matrix,crossi,crossj,crossk,tunnel_list,TS_list,tunnel_cluster,tunnel_cluster_dim)
%check if neighbor has no minID
if minID_matrix.L(in,jn,kn)==0 %If min_ID has not been assigned to neighbor, assign min_ID, crossing and add info to temp_list
    if level_matrix(in,jn,kn)==level
        minID_matrix.L(in,jn,kn)=list.C(iC).ID(1);
        minID_matrix.C(in,jn,kn)=list.C(iC).ID(2);
        cross_matrix.i(in,jn,kn)=crossi;
        cross_matrix.j(in,jn,kn)=crossj;
        cross_matrix.k(in,jn,kn)=crossk;
        index_temp=index_temp+1;
        temp.info(index_temp,1:10)=[in jn kn level 1 0 E_matrix(in,jn,kn) crossi crossj crossk];
    else
        boundary=1;
    end
else %check if min_ID has already been assigned to neighbor and if it is in current list or temp, if not it is a TS
    in_list=0;
    in_temp=0;
    if minID_matrix.C(in,jn,kn)==iC
        in_list=1;
        idiff=crossi-cross_matrix.i(in,jn,kn);
        jdiff=crossj-cross_matrix.j(in,jn,kn);
        kdiff=crossk-cross_matrix.k(in,jn,kn);
        if idiff==0 && jdiff==0 && kdiff==0
            in_list=1;
        else
            tunnel_list=[tunnel_list;idiff jdiff kdiff];
            if isempty(find(tunnel_cluster==iC, 1))==1
                tunnel_cluster=[tunnel_cluster iC];
                tunnel_cluster_dim=[tunnel_cluster_dim;idiff jdiff kdiff];
            end
            in_list=0;
        end
    end
    
    if isempty(temp)==0 %checking if temp list exists
        if sum(temp.info(:,1)==in & temp.info(:,2)==jn & temp.info(:,3)==kn) %checking if neighbor is in temp-list
            idiff=crossi-cross_matrix.i(in,jn,kn);
            jdiff=crossj-cross_matrix.j(in,jn,kn);
            kdiff=crossk-cross_matrix.k(in,jn,kn);
            if idiff==0 && jdiff==0 && kdiff==0
                in_temp=1;
            else
                tunnel_list=[tunnel_list;idiff jdiff kdiff];
                if isempty(find(tunnel_cluster==iC, 1))==1
                    tunnel_cluster=[tunnel_cluster iC];
                    tunnel_cluster_dim=[tunnel_cluster_dim;idiff jdiff kdiff];
                end
                in_temp=0;
            end
            
        end
    end
    if in_list==1 %if neighbor is not in any list it is a part of another cluster
    elseif in_temp==1
    else
        C_connect=minID_matrix.C(in,jn,kn); %get min_ID cluster of connecting cluster
        for iM=1:length(list.M)
            if sum(list.M(iM).clusters==C_connect)
                connect_Mx=iM;
            end
            if sum(list.M(iM).clusters==iC)
                current_Mx=iM;
            end
        end
        if connect_Mx==current_Mx %check that connecting clusters are not already merged, if they are check for tunnels
            idiff=crossi-cross_matrix.i(in,jn,kn);
            jdiff=crossj-cross_matrix.j(in,jn,kn);
            kdiff=crossk-cross_matrix.k(in,jn,kn);
            if jdiff==0 && idiff==0 && kdiff==0
            else %we have a tunnel
                tunnel_list=[tunnel_list;idiff jdiff kdiff];
                if sum(tunnel_cluster==iC)==0
                    tunnel_cluster=[tunnel_cluster iC];
                    tunnel_cluster_dim=[tunnel_cluster_dim;idiff jdiff kdiff];
                end
            end
            list_connect=list.C(C_connect).info(:,1:3); %get list of coordinates for connected cluster
            cluster_pair=sort([minID_matrix.C(i,j,k) minID_matrix.C(in,jn,kn)]);
            if TS_matrix(in,jn,kn)==0 && TS_matrix(i,j,k)==0 %if neighbor or host has not already been identified as TS
                if E_matrix(in,jn,kn)>=E_matrix(i,j,k) %if neighbor has higher or equal energy compared to host check
                    TS_matrix(in,jn,kn)=level_matrix(in,jn,kn); %neighbor is TS, entry in matrix
                    find_line=(list_connect(:,1)==in & list_connect(:,2)==jn & list_connect(:,3)==kn); %find line of TS in list of merging cluster
                    line=find(find_line);
                    list.C(C_connect).info(line,6)=1; %Set to TS in list of merging cluster
                    TS_list=[TS_list;in jn kn TS_matrix(in,jn,kn) E_matrix(in,jn,kn) cluster_pair crossi crossj crossk];
                    
                else %if neighbor is not higher level than host
                    TS_matrix(i,j,k)=level_matrix(i,j,k); %host is TS
                    list.C(iC).info(index_list,6)=1; %update TS in list
                    TS_list=[TS_list;i j k TS_matrix(i,j,k) E_matrix(i,j,k) cluster_pair crossi crossj crossk];
                    
                end
            end
        else
            list_connect=list.C(C_connect).info(:,1:3); %get list of coordinates for connected cluster
            cluster_pair=sort([minID_matrix.C(i,j,k) minID_matrix.C(in,jn,kn)]);
            if TS_matrix(in,jn,kn)==0 && TS_matrix(i,j,k)==0 %if neighbor or host has not already been identified as TS
                
                if E_matrix(in,jn,kn)>=E_matrix(i,j,k) %if neighbor is higher level than host
                    TS_matrix(in,jn,kn)=level_matrix(in,jn,kn); %neighbor is TS, entry in matrix
                    find_line=(list_connect(:,1)==in & list_connect(:,2)==jn & list_connect(:,3)==kn); %find line of TS in list of merging cluster
                    line=find(find_line);
                    list.C(C_connect).info(line,6)=1; %Set to TS in list of merging cluster
                    TS_list=[TS_list;in jn kn TS_matrix(in,jn,kn) E_matrix(in,jn,kn) cluster_pair crossi crossj crossk];
                    
                else %if neighbor is not higher level than host
                    TS_matrix(i,j,k)=level_matrix(i,j,k); %host is TS
                    list.C(iC).info(index_list,6)=1; %update TS in list
                    TS_list=[TS_list;i j k TS_matrix(i,j,k) E_matrix(i,j,k) cluster_pair crossi crossj crossk];
                    
                end
            end
            
            merge=[connect_Mx current_Mx];
            nM=length(list.M); %Get number of lines in list of merged lists
            
            idiff=crossi-cross_matrix.i(in,jn,kn);
            jdiff=crossj-cross_matrix.j(in,jn,kn);
            kdiff=crossk-cross_matrix.k(in,jn,kn);
            
            list_all_connect=[];
            length_connect=length(list.M(merge(1)).clusters);
            
            for index=1:length_connect %Loop through all clusters previously merged to connected cluster
                iC_connect=list.M(merge(1)).clusters(index);
                list_all_connect=[list_all_connect;list.C(iC_connect).info(:,1:3)]; %get list of coordinates for ALL connected cluster
            end
            
            merge=sort(merge);
            merged_clusters=[list.M(merge(1)).clusters list.M(merge(2)).clusters]; %merge the two lists of merged clusters
            list.M(merge(1)).clusters=merged_clusters; %put the new list in the position of first merged cluster list
            
            if merge(2)==nM % if second merged cluster list is the last list
                list.M(merge(2))=[]; %remove second merged cluster
            else %if not
                
                list.M(merge(2)).clusters=list.M(nM).clusters; %move last list to position of second merged cluster list
                list.M(nM)=[]; %and remove last list
            end
            if jdiff==0 && idiff==0 && kdiff==0
            else %if current PI of neighbor differs from existing PI in cross matrix,
                %loop through all coordinates of all connecting cluster and add on the
                %difference in PI in eachdirection.
                for ncoords_switch=1:length(list_all_connect(:,1))
                    icoord=list_all_connect(ncoords_switch,1);
                    jcoord=list_all_connect(ncoords_switch,2);
                    kcoord=list_all_connect(ncoords_switch,3);
                    cross_matrix.i(icoord,jcoord,kcoord)=cross_matrix.i(icoord,jcoord,kcoord)+idiff;
                    cross_matrix.j(icoord,jcoord,kcoord)=cross_matrix.j(icoord,jcoord,kcoord)+jdiff;
                    cross_matrix.k(icoord,jcoord,kcoord)=cross_matrix.k(icoord,jcoord,kcoord)+kdiff;
                end
                
            end
        end
    end
end
end