close all
clear all

%%%read user defined input%%%%%
input=load('input.param');
E_unit=input(1);
T=input(2);
run_kmc=input(3);
plot_msd=input(4);
nsteps=input(5);
print_every=input(6);
nruns=input(7);
n_particle=input(8);
per_tunnel=input(9);
m_particle=input(10);
energy_step=input(11);
energy_cutoff=input(12);

%%%initiate timing%%%%
tstart = cputime;
twrite_start = cputime;

%%%read grid data%%%%%
[ngrid,grid_size,pot_data]=cube2xsfdat(E_unit);
save('grid_kJmol.dat','pot_data', '-ASCII');

R=8.3144621; 
RT=T*R;
beta=1/RT;
ave_grid_size=mean(grid_size);
msd_steps=[];
for log_order=0:log10(nsteps/10)
    msd_steps=[msd_steps 10^log_order*(1:9)];
end
tlevel=0;
kappa=1/2;
prefactor=kappa*sqrt(1/(beta*2*pi*m_particle));
level_stop=ceil(energy_cutoff/energy_step);
data_size=size(pot_data);
min_value=min(min(pot_data));
max_value=max(max(pot_data));
nlevel=(max_value-min_value)/energy_step;
level_scale=nlevel/(max_value-min_value);
pot_value_shifted=pot_data-min_value;
pot_value_scaled=pot_value_shifted*level_scale+1/(nlevel*1000);
pot_value_scaled_int=ceil(pot_value_scaled);

%%%%convert grid data into matrices%%%%
ind=0;
disp('Writing matrix...');
for z=1:ngrid(3)
    for y=1:ngrid(2)
        for x=1:ngrid(1)
            ind=ind+1;
            row_ind=ceil(ind/6);
            col_ind=ind-(row_ind-1)*6;
            level_matrix(x,y,z)=pot_value_scaled_int(row_ind,col_ind);
        end
    end
end
ind=0;
for z=1:ngrid(3)
    for y=1:ngrid(2)
        for x=1:ngrid(1)
            ind=ind+1;
            row_ind=ceil(ind/6);
            col_ind=ind-(row_ind-1)*6;
            E_matrix(x,y,z)=pot_data(row_ind,col_ind);
        end
    end
end
disp('Done');
disp('...finished writing potential grid.')
ttot=cputime-tstart;
twrite=cputime-twrite_start;
disp(['Time write: ',num2str(twrite/60)]);
disp(['Total time elapsed: ',num2str(ttot/60)]);
%------------------------

[xsize, ysize, zsize]=size(level_matrix);
minID_matrix.L=zeros(xsize,ysize,zsize);
minID_matrix.C=zeros(xsize,ysize,zsize);
TS_matrix=zeros(xsize,ysize,zsize);
cross_matrix.i=zeros(xsize,ysize,zsize);
cross_matrix.j=zeros(xsize,ysize,zsize);
cross_matrix.k=zeros(xsize,ysize,zsize);
level_max=max(level_matrix(:));
level_min=min(level_matrix(:));
tunnel_cluster=[];
tunnel_cluster_dim=[];
TS_list_all=[];
group_list=[];
list.M=[];
list.TSgroup=[];
list.tunnels=[];
list.tunnel_directions=[];
list.TS_all=[];
list.TS_all=[];
BT=[0 0 0];
tot_tunnel_min=0;
tunnel_out='tunnel_info.out';
tunnel_file=fopen(tunnel_out,'w');
%%%%%start growing clusters%%%%%%%%%%%%%
disp('Starting TS search.......')
for level=level_min:level_stop %level_max-1
    tlevel_start=cputime;
    TS_list=[];
    TS_tunnel_list=[];
    tunnel_list=[];
    disp(strcat('level: ',num2str(level),'(',num2str(level/level_scale),' kJ/mol)'));
    nCluster=0;
    if level>level_min %grow all existing clusters one layer at a time checking for current level points
        nCluster=length(list.C);
        temp=[];
        %go through boundary points of each cluster looking for neighbors at current level. Also update the boundary.
        filled=zeros(nCluster);
        while sum(filled)<nCluster
            for iC=1:nCluster
                if filled(iC)==0
                    length_list=length(list.C(iC).info(:,1));
                    index_temp=0;
                    for index_list=1:length_list
                        if list.C(iC).info(index_list,5)==1 %check if grid point is a boundary point
                            boundary=0;
                            i=list.C(iC).info(index_list,1);
                            j=list.C(iC).info(index_list,2);
                            k=list.C(iC).info(index_list,3);
                            [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(i,j,k,xsize,ysize,zsize);
                            neighbor_coords=[i ip im j jp jm k kp km];
                            for l=1:3
                                for m=4:6
                                    for n=7:9
                                        checkneighbor=1;
                                        if l==2 && m==4 && n==7
                                        elseif l==3 && m==4 && n==7
                                        elseif l==1 && m==5 && n==7
                                        elseif l==1 && m==6 && n==7
                                        elseif l==1 && m==4 && n==8
                                        elseif l==1 && m==4 && n==9
                                        else
                                            checkneighbor=0;
                                        end
                                        if checkneighbor==1 && TS_matrix(i,j,k)==0 %%%%fix to not grow TS
                                            crossi=cross_matrix.i(i,j,k)+cross_vec(l); %periodic image coordinate of neighbor
                                            crossj=cross_matrix.j(i,j,k)+cross_vec(m);
                                            crossk=cross_matrix.k(i,j,k)+cross_vec(n);
                                            in=neighbor_coords(l); %neighbor coordinate
                                            jn=neighbor_coords(m);
                                            kn=neighbor_coords(n);
                                            if i==29 && j==3 && k==43
                                                TS_matrix(i,j,k)
                                                TS_matrix(in,jn,kn)
                                            elseif in==29 && jn==3 && kn==43
                                                TS_matrix(i,j,k)
                                                TS_matrix(in,jn,kn)
                                            elseif i==29 && j==4 && k==43
                                                TS_matrix(i,j,k)
                                                TS_matrix(in,jn,kn)
                                            elseif in==29 && jn==3 && kn==43
                                                TS_matrix(i,j,k)
                                                TS_matrix(in,jn,kn)
                                            end
                                            
                                            [index_temp,temp,boundary,minID_matrix,list,TS_matrix,cross_matrix,tunnel_list,TS_list,tunnel_cluster,tunnel_cluster_dim]=...
                                                check_neighbors(i,j,k,in,jn,kn,level,iC,index_list,index_temp,temp,boundary,...
                                                minID_matrix,list,TS_matrix,level_matrix,E_matrix,cross_matrix,crossi,crossj,crossk,...
                                                tunnel_list,TS_list,tunnel_cluster,tunnel_cluster_dim);
                                        end
                                    end
                                end
                            end
                            list.C(iC).info(index_list,5)=boundary;
                        end
                    end
                    if index_temp==0
                        filled(iC)=1;
                    else
                        list.C(iC).info=[list.C(iC).info;temp.info];
                        temp=[]; %empty temp file
                    end
                end
            end
        end
    end
    %%%%Looking for new clusters!%%%
    [unexplored_X, unexplored_Y, unexplored_Z] = ind2sub(size(minID_matrix.L),find(minID_matrix.L == 0));
    %step through list of unexplored
    for unexplored_ind=1:length(unexplored_X)
        i=unexplored_X(unexplored_ind);
        j=unexplored_Y(unexplored_ind);
        k=unexplored_Z(unexplored_ind);
        %level_matrix(i,j,k)
        if level_matrix(i,j,k)==level
            [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(i,j,k,xsize,ysize,zsize);
            %check if level point is connected to an existing cluster
            neighbor_id=[minID_matrix.L(ip,j,k) minID_matrix.L(im,j,k) minID_matrix.L(i,jp,k)...
                minID_matrix.L(i,jm,k) minID_matrix.L(i,j,kp) minID_matrix.L(i,j,km)];
            if sum(neighbor_id)==0
                length_list=1;
                index_list=1;
                nCluster=nCluster+1;
                minID=[level nCluster];
                list.C(nCluster).ID=minID;
                minID_matrix.L(i,j,k,1)=level;
                minID_matrix.C(i,j,k,1)=nCluster;
                boundary=0;
                list.C(nCluster).info(length_list,1:10)=[i j k level 1 0 E_matrix(i,j,k) 0 0 0];
                first_empty=length(list.M)+1;
                list.M(first_empty).clusters=nCluster; %Add cluster to list of list of merged cluster
                
                while index_list<=length_list
                    if minID_matrix.L(ip,j,k)==0
                        if level_matrix(ip,j,k)==level
                            length_list=length_list+1;
                            minID_matrix.L(ip,j,k)=level;
                            minID_matrix.C(ip,j,k)=nCluster;
                            crossi=cross_matrix.i(i,j,k)+cross_vec(2);
                            crossj=cross_matrix.j(i,j,k)+cross_vec(4);
                            crossk=cross_matrix.k(i,j,k)+cross_vec(7);
                            cross_matrix.i(ip,j,k)=crossi;
                            cross_matrix.j(ip,j,k)=crossj;
                            cross_matrix.k(ip,j,k)=crossk;
                            list.C(nCluster).info(length_list,1:10)=[ip j k level 1 0 E_matrix(ip,j,k) crossi crossj crossk];
                        else
                            boundary=1;
                        end
                    end
                    if minID_matrix.L(im,j,k)==0
                        if level_matrix(im,j,k)==level
                            length_list=length_list+1;
                            minID_matrix.L(im,j,k)=level;
                            minID_matrix.C(im,j,k)=nCluster;
                            crossi=cross_matrix.i(i,j,k)+cross_vec(3);
                            crossj=cross_matrix.j(i,j,k)+cross_vec(4);
                            crossk=cross_matrix.k(i,j,k)+cross_vec(7);
                            cross_matrix.i(im,j,k)=crossi;
                            cross_matrix.j(im,j,k)=crossj;
                            cross_matrix.k(im,j,k)=crossk;
                            list.C(nCluster).info(length_list,1:10)=[im j k level 1 0 E_matrix(im,j,k) crossi crossj crossk];
                            
                        else
                            boundary=1;
                        end
                    end
                    if minID_matrix.L(i,jp,k)==0
                        if level_matrix(i,jp,k)==level
                            length_list=length_list+1;
                            minID_matrix.L(i,jp,k)=level;
                            minID_matrix.C(i,jp,k)=nCluster;
                            crossi=cross_matrix.i(i,j,k)+cross_vec(1);
                            crossj=cross_matrix.j(i,j,k)+cross_vec(5);
                            crossk=cross_matrix.k(i,j,k)+cross_vec(7);
                            cross_matrix.i(i,jp,k)=crossi;
                            cross_matrix.j(i,jp,k)=crossj;
                            cross_matrix.k(i,jp,k)=crossk;
                            list.C(nCluster).info(length_list,1:10)=[i jp k level 1 0 E_matrix(i,jp,k) crossi crossj crossk];
                            
                        else
                            boundary=1;
                        end
                    end
                    if minID_matrix.L(i,jm,k)==0
                        if level_matrix(i,jm,k)==level
                            length_list=length_list+1;
                            minID_matrix.L(i,jm,k)=level;
                            minID_matrix.C(i,jm,k)=nCluster;
                            crossi=cross_matrix.i(i,j,k)+cross_vec(1);
                            crossj=cross_matrix.j(i,j,k)+cross_vec(6);
                            crossk=cross_matrix.k(i,j,k)+cross_vec(7);
                            cross_matrix.i(i,jm,k)=crossi;
                            cross_matrix.j(i,jm,k)=crossj;
                            cross_matrix.k(i,jm,k)=crossk;
                            list.C(nCluster).info(length_list,1:10)=[i jm k level 1 0 E_matrix(i,jm,k) crossi crossj crossk];
                        else
                            boundary=1;
                        end
                    end
                    if minID_matrix.L(i,j,kp)==0
                        if level_matrix(i,j,kp)==level
                            length_list=length_list+1;
                            minID_matrix.L(i,j,kp)=level;
                            minID_matrix.C(i,j,kp)=nCluster;
                            crossi=cross_matrix.i(i,j,k)+cross_vec(1);
                            crossj=cross_matrix.j(i,j,k)+cross_vec(4);
                            crossk=cross_matrix.k(i,j,k)+cross_vec(8);
                            cross_matrix.i(i,j,kp)=crossi;
                            cross_matrix.j(i,j,kp)=crossj;
                            cross_matrix.k(i,j,kp)=crossk;
                            list.C(nCluster).info(length_list,1:10)=[i j kp level 1 0 E_matrix(i,j,kp) crossi crossj crossk];
                        else
                            boundary=1;
                        end
                    end
                    if minID_matrix.L(i,j,km)==0
                        if level_matrix(i,j,km)==level
                            length_list=length_list+1;
                            minID_matrix.L(i,j,km)=level;
                            minID_matrix.C(i,j,km)=nCluster;
                            crossi=cross_matrix.i(i,j,k)+cross_vec(1);
                            crossj=cross_matrix.j(i,j,k)+cross_vec(4);
                            crossk=cross_matrix.k(i,j,k)+cross_vec(9);
                            cross_matrix.i(i,j,km)=crossi;
                            cross_matrix.j(i,j,km)=crossj;
                            cross_matrix.k(i,j,km)=crossk;
                            list.C(nCluster).info(length_list,1:10)=[i j km level 1 0 E_matrix(i,j,km) crossi crossj crossk];
                        else
                            boundary=1;
                        end
                    end
                    list.C(nCluster).info(index_list,5)=boundary;
                    if index_list<length_list
                        index_list=index_list+1;
                        i=list.C(nCluster).info(index_list,1);
                        j=list.C(nCluster).info(index_list,2);
                        k=list.C(nCluster).info(index_list,3);
                        [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(i,j,k,xsize,ysize,zsize);
                    else
                        break
                    end
                end
            end
        end
    end
    
    %list.TS=TS_list;
    TS_list_all=[TS_list_all; TS_list];
    echelon=rref(tunnel_list);
    if ~isempty(echelon)
        for q=1:length(echelon(:,1))
            if sum(echelon(q,:))==0
                break
            else
                list.tunnel_directions(1:q,:)=echelon(1:q,:);
            end
        end
        disp(strcat({'Level '},num2str(level),{'('},num2str(level/level_scale),{' kJ/mol)'},{', cluster '},{' breaks through in directions:'}));
        disp(echelon(1:q,:));
        list.tunnel_directions=echelon(1:q,:);
    end
    
    
    
    fprintf(tunnel_file,'%s\n','------------------------------------------------------------');
    fprintf(tunnel_file,'%s %d %s %f %s\n','Level is ',level,' (',level/level_scale,' kJ/mol)');
    if isempty(list.tunnel_directions)
        fprintf(tunnel_file,'%s\n','Does not break through');
    else
        fprintf(tunnel_file,'%s\n','Breaks through in directions: ');
        fprintf(tunnel_file,'%d %8d %8d\n', [sum(list.tunnel_directions(:,1)) sum(list.tunnel_directions(:,2)) sum(list.tunnel_directions(:,3))]);
        fprintf(tunnel_file,'%s %d %s\n', 'Overall min value: ',tot_tunnel_min,' [kJ/mol]');
        for abc=1:3
            if BT(abc)==0 && sum(list.tunnel_directions(:,abc))>0
                BT(abc)=level*energy_step;
            end
        end
    end
    
    ttot=cputime-tstart;
    tlevel=cputime-tlevel_start;
    disp(['Time level ',num2str(level),': ',num2str(tlevel/60)]);
    fprintf(tunnel_file,'%s %d %s %f\n','Time level ',level,': ',tlevel/60);
    disp(['Total time elapsed: ',num2str(ttot/60)]);
    fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
end

disp('Breakthrough in:');
disp(strcat('A direction: ',num2str(BT(1)),'kJ/mol'));
disp(strcat('B direction: ',num2str(BT(2)),'kJ/mol'));
disp(strcat('C direction: ',num2str(BT(3)),'kJ/mol'));
BT_out = 'BT.dat';
fid_BT = fopen(BT_out,'w');
fprintf(fid_BT,'%s',num2str(BT));
fclose(fid_BT);

if sum(BT)>0
    %tmp=list; When debugging list->tmp so it does not overwrite list
    %Print cluster list   %%%%are cluster mins defined somewhere else?%%%%
    for iC=1:length(list.C)
        cluster_min=min(list.C(iC).info(:,7));
        index_min_cluster=find(list.C(iC).info(:,7)==min(list.C(iC).info(:,7)),1);
        coord_min_cluster=list.C(iC).info(index_min_cluster,1:3);
        list.C(iC).min=[cluster_min coord_min_cluster];
    end
    %group TS of cluster pairs into connected surfaces
    if isempty(TS_list_all)==0
        TS_list_all(:,11)=0;
        [list,TS_list_all]=organize_TS(TS_list_all,xsize,ysize,zsize,level_scale,beta,list,level,cross_matrix);
    end
    
    %Find clusters and TS which are in tunnels
    iT=0;
    idTSG=0;
    listed_clusters=[];
    process_data=[];
    for iTC=1:length(tunnel_cluster) %Go through coordinates in tunnel_cluster list
        TS_list_tunnel=[]; %initiate list of TS in tunnels
        tunnel_dim=tunnel_cluster_dim(iTC,:); % list of merging tunnel dimensions
        if iT>0 && ismember(tunnel_cluster(iTC),listed_clusters)%sum(list.tunnels(:).clusters==tunnel_cluster(iTC)) %check if tunnel cluster has already been checked
        else
            iT=iT+1;
            for iM=1:length(list.M) %loop through list of merged clusters
                if sum(list.M(iM).clusters==tunnel_cluster(iTC)) %find list of merged clusters containing tunnel cluster
                    list.tunnels(iT).clusters=list.M(iM).clusters; %Add list of merged clusters to tunnel cluster list
                    listed_clusters=[listed_clusters list.M(iM).clusters];
                    tunnel_energies=[];
                    for iiTC=1:length(list.tunnels(iT).clusters)
                        tunnel_energies=[tunnel_energies; list.C(list.tunnels(iT).clusters(iiTC)).info(:,7)];
                        list.tunnels(iT).min=min(tunnel_energies);
                        list.tunnels(iT).vol_proc=length(tunnel_energies)/(ngrid(1)*ngrid(2)*ngrid(3));
                        if list.tunnels(iT).min<tot_tunnel_min
                            tot_tunnel_min=list.tunnels(iT).min;
                        end
                        if ismember(list.tunnels(iT).clusters(iiTC),tunnel_cluster)
                            iiiTC=find(list.tunnels(iT).clusters(iiTC)==tunnel_cluster);
                            tunnel_dim=[tunnel_dim; tunnel_cluster_dim(iiiTC,:)];
                        end
                        echelon_TC=rref(tunnel_list);
                        list.tunnels(iT).dimensions=[sum(echelon_TC(:,1)) sum(echelon_TC(:,2)) sum(echelon_TC(:,3))];
                    end
                    break
                end
            end
            %Make process list of cluster pairs and TSgroup in tunnels
            iTSG=0;
            iP=0;
            iCM=0;
            for iG=1:length(list.TSgroup) %loop throup all TS groups
                if isempty(list.TSgroup(iG))==0
                    if ismember(list.TSgroup(iG).data(:,6:7),list.tunnels(iT).clusters) %check if in tunnel
                        iTSG=iTSG+1; %increment TSgroup index in tunnel
                        idTSG=list.TSgroup(iG).info.group(1);
                        TS_list_tunnel=[TS_list_tunnel; list.TSgroup(iG).data]; %add TS group list to list of tunnel TS
                        list.tunnels(iT).TSgroup(iTSG)=list.TSgroup(iG);
                        id_C1=list.tunnels(iT).TSgroup(iTSG).info.cluster1(1);
                        id_C2=list.tunnels(iT).TSgroup(iTSG).info.cluster2(1);
                        dE_1=list.tunnels(iT).TSgroup(iTSG).info.cluster1(6);
                        dE_2=list.tunnels(iT).TSgroup(iTSG).info.cluster2(6);
                        process1_info=[id_C1 id_C2 dE_1 iT iTSG];
                        process2_info=[id_C2 id_C1 dE_2 iT iTSG];
                        process_data=[process_data;process1_info;process2_info];
                    end
                end
            end
            list.tunnels(iT).TS=TS_list_tunnel; %Bookkeep whole tunnel TS list
        end
    end
    
    
    %write final clusters
    for iC=1:length(list.C)
        cluster_points=(list.C(iC).info(:,6)==0);
        final.C(iC).info=list.C(iC).info(cluster_points,[1 2 3 7 8 9 10]);
        final.C(iC).min=list.C(iC).min;
        final.C(iC).merged=0;
    end
    
    
    disp('Merging clusters with barriers < cutoff....')
    tmerge_start=cputime;
    TS2keep=[];
    processes_start=process_data; %%%%change order when done debugging%%%%
    clus2merge=find(process_data(:,3)<energy_step);
    while isempty(clus2merge)==0
        i=clus2merge(1);
        C1=process_data(i,1);
        C2=process_data(i,2);
        iT=process_data(i,4);
        iTSG=process_data(i,5);
        
        min_cluster_C1=final.C(C1).min(1);
        min_cluster_C2=final.C(C2).min(1);
        if min_cluster_C1<min_cluster_C2
            C_keep=C1;
            C_change=C2;
            min_cluster_pair=min_cluster_C1;
        else
            C_keep=C2;
            C_change=C1;
            min_cluster_pair=min_cluster_C2;
        end
        disp(strcat('merging cluster: ',num2str(C_change),' into cluster:',num2str(C_keep)))
        minID_matrix.C(minID_matrix.C==C_change)=C_keep; %change cluster number in minID_matrix.C
        
        %Go through TS coordinates to see if it remains TS after merging,
        %make list of these and remove from current TSgroup which will be
        %merged into cluster
        ind_TSG=0;
        while ind_TSG<length(list.tunnels(iT).TSgroup(iTSG).data(:,1))
            ind_TSG=ind_TSG+1;
            coords=list.tunnels(iT).TSgroup(iTSG).data(ind_TSG,1:5); %TS coordinates + levels and energies
            [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(coords(1),coords(2),coords(3),xsize,ysize,zsize);
            neighbor_clusters=[minID_matrix.C(coords(1),coords(2),kp) minID_matrix.C(coords(1),jp,coords(3)) minID_matrix.C(ip,coords(2),coords(3))...
                minID_matrix.C(coords(1),coords(2),km) minID_matrix.C(coords(1),jm,coords(3)) minID_matrix.C(im,coords(2),coords(3))];
            unique_clusters=unique(neighbor_clusters(neighbor_clusters>0));
            if length(unique_clusters)<2 %if not in interface of two clusters, remove from TS_matrix. Will be merged into cluster
                TS_matrix(coords(1),coords(2),coords(3))=0;
            elseif length(unique_clusters)==2 %if TS is in interface of two cluster-pairs (most often, maybe always?)
                TS2keep=[TS2keep;[coords unique_clusters iT 0 0]];
                list.tunnels(iT).TSgroup(iTSG).data(ind_TSG,:)=[];
                ind_TSG=ind_TSG-1;
            else %if TS is in interface of several cluster-pairs
                for iUC=1:length(unique_clusters)-1
                    for jUC=2:length(unique_clusters)
                        TS2keep=[TS2keep;[coords unique_clusters(iUC) unique_clusters(jUC) iT 0 0]];
                    end
                end
                list.tunnels(iT).TSgroup(iTSG).data(ind_TSG,:)=[];
                ind_TSG=ind_TSG-1;
            end
            if ind_TSG==length(list.tunnels(iT).TSgroup(iTSG).data(:,1))
                break
            end
        end
        %check if TS2keep clusters have merged into another
        if isempty(TS2keep)==0
            TS2keep((TS2keep(:,6)==C_change),6)=C_keep; 
            TS2keep((TS2keep(:,7)==C_change),7)=C_keep;
            TS2keep((TS2keep(:,6)==C_keep & TS2keep(:,7)==C_keep),:)=[]; %!!!Add this point to C_keep cluster?
        end
        %merge C1, C2 and TS coords to one cluster
        cluster_pair=[final.C(C1).info;list.tunnels(iT).TSgroup(iTSG).data(:,[1 2 3 5 8 9 10]);final.C(C2).info]; %Merge cluster pair and TS
        n_cluster_pair=length(cluster_pair(:,1)); %compute number of grids in merged cluster
        final.C(C_keep).info=cluster_pair;
        final.C(C_change).merged=C_keep;
        
        %find and make list of which processes to be removed
        remove1 = find(process_data(:,1)==C1 & process_data(:,2)==C2);
        remove2 = find(process_data(:,1)==C2 & process_data(:,2)==C1);
        remove=[remove1;remove2];
        
        %find processes which end in cluster which has been merged and update cluster number and cross vector
        clus2update1=find(process_data(:,1)==C_change);
        for j=1:length(clus2update1)
            j_cluster=clus2update1(j);
            process_data(j_cluster,1)=C_keep;
        end
        clus2update2=find(process_data(:,2)==C_change); %find clusters which end in cluster which has been merged and update cluster number and cross vector
        for k=1:length(clus2update2)
            k_cluster=clus2update2(k);
            process_data(k_cluster,2)=C_keep;
            
        end
        proc2update=find(process_data(:,1)==C_keep);  %find processes which start from cluster with updated volume and Bsum
        for l=1:length(proc2update)
            l_process=proc2update(l);
            process_data(l_process,3)=list.tunnels(process_data(l_process,4)).TSgroup(process_data(l_process,5)).info.group(5)-min_cluster_pair; %update dE
        end
        for iT=1:length(list.tunnels) %%remove merged clusters from tunnel cluster bookkeeping
            if sum(list.tunnels(iT).clusters==C_change)
                list.tunnels(iT).clusters=list.tunnels(iT).clusters(list.tunnels(iT).clusters~=C_change);
            end
        end
        %remove processes
        process_temp=[];
        i_temp=0;
        for i_process=1:length(process_data(:,1))
            if sum(remove==i_process)==0
                i_temp=i_temp+1;
                process_temp(i_temp,:)=process_data(i_process,:);
            end
        end
        process_data=process_temp;
        %check if there are any barriers <cutoff left. If not, break.
        clus2merge=find(process_data(:,3)<energy_step);
        if isempty(clus2merge)
            break
        end
    end
    
    TS_tunnel_new=TS2keep;
    for iT=1:length(list.tunnels) %step through each tunnel
        for iTSG=1:length(list.tunnels(iT).TSgroup) %step throup each TS group
            NP=find(process_data(:,4)==iT & process_data(:,5)==iTSG); %find row in list of processes corresponding to tunnel and TS group
            if isempty(NP)==0 %If it is in the process list, keep TS group
                disp('Keep TS group')
                lTSG=length(list.tunnels(iT).TSgroup(iTSG).data(:,1));
                C1=process_data(NP(1),1); % can be removed
                C2=process_data(NP(1),2);
                C=sort([C1 C2]);
                TS_tunnel_new=[TS_tunnel_new;[list.tunnels(iT).TSgroup(iTSG).data(:,1:5) C.*ones(lTSG,1)  iT*ones(lTSG,1) zeros(lTSG,1) zeros(lTSG,1)]];
            else
                disp('Remove TS group')
            end
        end
    end
    %%%%Recheck groups after merging%%%%
    new=[];
    if isempty(TS_tunnel_new)==0
        [new,TS_tunnel_new]=organize_TS_tunnel(TS_tunnel_new,xsize,ysize,zsize,list);
    end
    
    %%%%Update process list, TSgroup info and Tunnel info%%%
    iP=0;
    for iT=1:length(new.tunnels)
        for iTSG=1:length(new.tunnels(iT).TSgroup)
            new.tunnels(iT).TS=(TS_tunnel_new(:,8)==iT); %Update tunnel TS list
            new.tunnels(iT).min=list.tunnels(iT).min;
            new.tunnels(iT).dimensions=list.tunnels.dimensions;
            new.tunnels(iT).vol_proc=list.tunnels.vol_proc;
            C1=new.tunnels(iT).TSgroup(iTSG).info.cluster1(1);
            C2=new.tunnels(iT).TSgroup(iTSG).info.cluster2(1);
            idTSG=new.tunnels(iT).TSgroup(iTSG).info.group(1);
            TS_min=new.tunnels(iT).TSgroup(iTSG).info.group(5);
            dE_1=new.tunnels(iT).TSgroup(iTSG).info.cluster1(6);
            dE_2=new.tunnels(iT).TSgroup(iTSG).info.cluster2(6);
            iP=iP+1;
            process_data_new(iP,:)=[C1 C2 dE_1 iT iTSG idTSG];
            iP=iP+1;
            process_data_new(iP,:)=[C2 C1 dE_2 iT iTSG idTSG];
        end
    end
    
    disp('...finished merging clusters.')
    ttot=cputime-tstart;
    tmerge=cputime-tmerge_start;
    disp(['Time merge: ',num2str(tmerge/60)]);
    fprintf(tunnel_file,'%s %d %s %f\n','Time level ',level,': ',tmerge/60);
    disp(['Total time elapsed: ',num2str(ttot/60)]);
    fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
   
    %%%Check process boundaries%%%%%
    disp('Checking process boundaries....')
    tboundary_start=cputime;
    clear i
    %checking which clusters are at boundary
    for iC=1:length(final.C)
        if sum(final.C(iC).info(:,1)==1)>1
            final.C(iC).boundary=1;
        elseif sum(final.C(iC).info(:,1)==ngrid(1))>1
            final.C(iC).boundary=1;
        elseif sum(final.C(iC).info(:,2)==1)>1
            final.C(iC).boundary=1;
        elseif sum(final.C(iC).info(:,2)==ngrid(2))>1
            final.C(iC).boundary=1;
        elseif sum(final.C(iC).info(:,3)==1)>1
            final.C(iC).boundary=1;
        elseif sum(final.C(iC).info(:,3)==ngrid(3))>1
            final.C(iC).boundary=1;
        else
            final.C(iC).boundary=0;
        end
    end
    
    %%%Check boundary crossing and process direction%%%%
    TS_tunnel_list=[]; 
    ind_group=0;
    processes=[];
    for i=1:2:length(process_data_new(:,1)) %processes come in forward/backward pairs
        C1=process_data_new(i,1);
        C2=process_data_new(i,2);
        iT=process_data_new(i,4);
        iTSG=process_data_new(i,5); %row number in TSgroup list of respective tunnel
        idTSG=process_data_new(i,6); %TSgroup ID
        if final.C(C1).boundary==1
            coord_list=[final.C(C1).info(:,1:3);new.tunnels(iT).TSgroup(iTSG).data(:,1:3);final.C(C2).info(:,1:3)];
            start_min=final.C(C1).min(2:4);
            end_min=final.C(C2).min(2:4);
            [TS_cross_vector,coord_list]=get_TS_cross_vector(start_min,end_min,coord_list,ngrid,C1,C2,iT,idTSG);
        elseif final.C(C2).boundary==1
            coord_list=[final.C(C1).info(:,1:3);new.tunnels(iT).TSgroup(iTSG).data(:,1:3);final.C(C2).info(:,1:3)];
            start_min=final.C(C2).min(2:4);
            end_min=final.C(C1).min(2:4);
            [TS_cross_vector,coord_list]=get_TS_cross_vector(start_min,end_min,coord_list,ngrid,C1,C2,iT,idTSG);
            TS_cross_vector=-1*TS_cross_vector;
        else
            TS_cross_vector=[0 0 0];
        end
        
        %%%Compute boltzmann sums and rates for final processes%%%%
        Bsum_TS1=sum(exp(-1000*beta*(new.tunnels(iT).TSgroup(iTSG).data(:,5)-final.C(C1).min(1))));
        Bsum_TS2=sum(exp(-1000*beta*(new.tunnels(iT).TSgroup(iTSG).data(:,5)-final.C(C2).min(1))));
        Bsum_cluster1=sum(exp(-1000*beta*(final.C(C1).info(:,4)-final.C(C1).min(1))))+Bsum_TS1;
        Bsum_cluster2=sum(exp(-1000*beta*(final.C(C2).info(:,4)-final.C(C2).min(1))))+Bsum_TS2;
        k1=prefactor*Bsum_TS1/(Bsum_cluster1*ave_grid_size*1e-10);
        k2=prefactor*Bsum_TS2/(Bsum_cluster2*ave_grid_size*1e-10);
        
        %%%Write processes%%%%
        processes(i,:)=[C1 C2 k1 TS_cross_vector iT idTSG];
        processes(i+1,:)=[C2 C1 k2 -TS_cross_vector iT idTSG];
        
        %%%write list for printing%%%%
        TS_coords=new.tunnels(iT).TSgroup(iTSG).data(:,1:3);
        TS_levels=new.tunnels(iT).TSgroup(iTSG).data(:,4);
        TS_groupID=new.tunnels(iT).TSgroup(iTSG).data(:,10);
        TS_zeroes=zeros(length(new.tunnels(iT).TSgroup(iTSG).data(:,1)),1);
        TS_tunnel_list=[TS_tunnel_list;[TS_groupID ceil(TS_levels*energy_step) TS_zeroes TS_coords-0.5]];
    end

    disp('...finished checking process boundaries.')
    ttot=cputime-tstart;
    tboundary=cputime-tboundary_start;
    disp(['Time merge: ',num2str(tmerge/60)]);
    
    %%%Print tunnel info%%%%
    fprintf(tunnel_file,'%s %d %s %f\n','Time level ',level,': ',tmerge/60);
    disp(['Total time elapsed: ',num2str(ttot/60)]);
    fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
    tot_vol_proc=0;
    tot_n_TS=0;
    for iT=1:length(new.tunnels)
        if isempty(new.tunnels(iT).TS)
        else
            fprintf(tunnel_file,'%s\n','------------------------------------------------------------');
            fprintf(tunnel_file,'%s %d %s\n','Tunnel info ',iT,':');
            fprintf(tunnel_file,'%s %d %d %d\n','Breaks through in directions: ',new.tunnels(iT).dimensions);
            fprintf(tunnel_file,'%s %d %s\n','Min value:',new.tunnels(iT).min,' kJ/mol');
            fprintf(tunnel_file,'%s %d %s %d %s\n','This tunnel occupies ',new.tunnels(iT).vol_proc*100,'% of the total volume and has ',length(new.tunnels(iT).TSgroup),' transition states');
            fprintf(tunnel_file,'%s\n','------------------------------------------------------------');
            tot_vol_proc=tot_vol_proc+new.tunnels(iT).vol_proc*100;
            tot_n_TS=tot_n_TS+length(new.tunnels(iT).TSgroup);
            for iTSG=1:length(new.tunnels(iT).TSgroup)
                fprintf(tunnel_file,'%s \n','group ID, min x, min y, min z, TS_min, B_min_abs, Bsum_TS');
                fprintf(tunnel_file,'%s \n',num2str(new.tunnels(iT).TSgroup(iTSG).info.group));
                fprintf(tunnel_file,'%s \n','cluster_1, min x, min y, min x, cluster1_min, dE_1');
                fprintf(tunnel_file,'%s \n',num2str(new.tunnels(iT).TSgroup(iTSG).info.cluster1));
                fprintf(tunnel_file,'%s \n','cluster_2, min x, min y, min x, cluster2_min, dE_2');
                fprintf(tunnel_file,'%s \n',num2str(new.tunnels(iT).TSgroup(iTSG).info.cluster2));
                fprintf(tunnel_file,'%s\n','------------------------------------------------------------');
            end
            fprintf(tunnel_file,'%s\n','------------------------------------------------------------');
        end
    end
    TS_out='TS_data.out';
    TS_file=fopen(TS_out,'w');
    for iTS=1:length(TS_tunnel_list(:,1))
        fprintf(TS_file,'%s \n',num2str(TS_tunnel_list(iTS,:)));
    end
    fclose(TS_file);

    %%%Compute and organize basis sites%%%%
    tunnel_cluster_list=unique(processes(:,1:2));
    for cluster_row=1:length(tunnel_cluster_list)
        for iT=1:length(new.tunnels)
            if sum(tunnel_cluster_list(cluster_row)==list.tunnels(iT).clusters)
                basis_tunnels(cluster_row)=iT;
            end
        end
        iC=tunnel_cluster_list(cluster_row);
        %%%find cluster basis site
        if final.C(iC).boundary==0
            cluster_coords=final.C(iC).info(:,1:3);
            cluster_weights=exp(-final.C(iC).info(:,4)*1000/(R*T));
            cluster_basis=[round(wmean(cluster_coords(:,1),cluster_weights)) round(wmean(cluster_coords(:,2),cluster_weights)) round(wmean(cluster_coords(:,3),cluster_weights))];
            final.C(iC).basis_site=cluster_basis;
       else
            final=update_crossing(final,ngrid,iC);
            cluster_coords=final.C(iC).info(:,1:3)+final.C(iC).info(:,5:7).*ngrid;
            cluster_weights=exp(-final.C(iC).info(:,4)*1000/(R*T));
            cluster_basis=[round(wmean(cluster_coords(:,1),cluster_weights)) round(wmean(cluster_coords(:,2),cluster_weights)) round(wmean(cluster_coords(:,3),cluster_weights))];
            for dim=1:3
                if cluster_basis(dim)>ngrid(dim)
                    multiple=floor(cluster_basis(dim)/ngrid(dim));
                    cluster_basis(dim)=cluster_basis(dim)-multiple*ngrid(dim);
                elseif cluster_basis(dim)<0
                    multiple=ceil(-cluster_basis(dim)/ngrid(dim));
                    cluster_basis(dim)=cluster_basis(dim)+multiple*ngrid(dim);
                end
            end
            final.C(iC).basis_site=cluster_basis;
        end
                
        basis_sites(cluster_row,:)=[cluster_basis basis_tunnels(cluster_row)];
        update_cluster1_index=find(processes(:,1)==tunnel_cluster_list(cluster_row));
        processes(update_cluster1_index,1)=cluster_row;
        update_cluster2_index=find(processes(:,2)==tunnel_cluster_list(cluster_row));
        processes(update_cluster2_index,2)=cluster_row;
    end
    %%%Print basis%%%%
    basis_out = 'basis.dat';
    fid_basis = fopen(basis_out,'w');
    for i=1:length(basis_sites(:,1))
        fprintf(fid_basis,'%s\n',num2str(basis_sites(i,:)));
    end
    fclose(fid_basis);
    
    %%%%Print processes%%%%
    proc_out = 'processes.dat';
    fid_proc = fopen(proc_out,'w');
    for j=1:length(processes(:,1))
        fprintf(fid_proc,'%s\n',num2str(processes(j,:)));
    end
    fclose(fid_proc);
    
    %%%Run KMC%%%%%
    if run_kmc==1
        disp('running KMC')
        tstart_kmc = cputime;
        if plot_msd==1
            [msd,process_exec,distances,D_ave]=kmc_plot(basis_sites,processes,ngrid,n_particle,nsteps,msd_steps,print_every,grid_size,nruns,BT,per_tunnel);
        else
            [msd,process_exec,distances,D_ave]=kmc_noplot(basis_sites,processes,ngrid,n_particle,nsteps,msd_steps,print_every,grid_size,nruns,BT,per_tunnel);
        end
        disp('Diffusion coefficients in')
        disp(strcat('A direction: ',num2str(D_ave(1)),'+/-',num2str(D_ave(2)),'cm^2/s'))
        disp(strcat('B direction: ',num2str(D_ave(3)),'+/-',num2str(D_ave(4)),'cm^2/s'))
        disp(strcat('C direction: ',num2str(D_ave(5)),'+/-',num2str(D_ave(6)),'cm^2/s'))
        D_out = 'D_ave.dat';
        fid_D = fopen(D_out,'w');
        fprintf(fid_D,'%s',num2str(D_ave));
        fclose(fid_D);
        t_kmc=cputime-tstart_kmc;
        ttot=cputime-tstart;
        disp(['Time kmc : ',num2str(t_kmc/60)]);
        fprintf(tunnel_file,'%s %d %s %f\n','Time kmc: ',t_kmc/60);
        disp(['Total time elapsed: ',num2str(ttot/60)]);
        fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
    end
else
    disp('Diffusion coefficients in')
    disp(strcat('A direction: ',num2str(0),'cm^2/s'))
    disp(strcat('B direction: ',num2str(0),'cm^2/s'))
    disp(strcat('C direction: ',num2str(0),'cm^2/s'))
    
    D_out = 'D_ave.dat';
    fid_D = fopen(D_out,'w');
    fprintf(fid_D,'%s',num2str([0 0 0 0 0 0 0]));
    fclose(fid_D);
    ttot=cputime-tstart;
    disp(['Total time elapsed: ',num2str(ttot/60)]);
    fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
end

fclose(tunnel_file);




