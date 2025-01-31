close all
clear all


    %%%read user defined input%%%%%
    input=load('input.param');
    E_unit=input(1);
    nTemp=input(2);
    T=[];
    for iTemp=1:nTemp
        T=[T input(2+iTemp)];
    end
    run_kmc=input(3+nTemp);
    plot_msd=input(4+nTemp);
    nsteps=input(5+nTemp);
    print_every=input(6+nTemp);
    nruns=input(7+nTemp);
    n_particle=input(8+nTemp);
    per_tunnel=input(9+nTemp);
    m_particle=input(10+nTemp);
    energy_step=input(11+nTemp);
    energy_cutoff=input(12+nTemp);
    
    %%%initiate timing%%%%
    tstart = cputime;
    twrite_start = cputime;
    
    %%%read grid data%%%%%
    
    [ngrid,grid_size,pot_data]=cube2xsfdat(E_unit);
    
    R=8.3144621;
    RT=T*R;
    beta=1./RT;
    ave_grid_size=mean(grid_size);
    msd_steps=[];
    for log_order=0:log10(nsteps/10)
        msd_steps=[msd_steps 10^log_order*(1:9)];
    end
    
    tlevel=0;
    kappa=1/2;
    prefactor=kappa*sqrt(1./(beta*2*pi*m_particle));
    level_stop=ceil(energy_cutoff/energy_step);
    data_size=size(pot_data);
    pot_data(find(pot_data==max(pot_data)))=10*energy_cutoff;
    pot_data(find(pot_data==inf))=10*energy_cutoff;

    %%%%convert grid data into matrices%%%%

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

    min_value=min(min(pot_data));
    max_value=max(max(pot_data));
    nlevel=(max_value-min_value)/energy_step;
    level_scale=nlevel/(max_value-min_value);
    pot_value_shifted=E_matrix-min_value;
    pot_value_scaled=pot_value_shifted*level_scale+1/(level_stop*1000);
    level_matrix=ceil(pot_value_scaled);
    
    disp('Done');
    disp('...finished writing potential grid.')
    ttot=cputime-tstart;
    twrite=cputime-twrite_start;
    disp(['Time write: ',num2str(twrite/60)]);
    disp(['Total time elapsed: ',num2str(ttot/60)]);
    %------------------------
    
    [ngrid(1), ngrid(2), ngrid(3)]=size(level_matrix);
    minID_matrix.L=zeros(ngrid(1),ngrid(2),ngrid(3));
    minID_matrix.C=zeros(ngrid(1),ngrid(2),ngrid(3));
    TS_matrix=zeros(ngrid(1),ngrid(2),ngrid(3));
    cross_matrix.i=zeros(ngrid(1),ngrid(2),ngrid(3));
    cross_matrix.j=zeros(ngrid(1),ngrid(2),ngrid(3));
    cross_matrix.k=zeros(ngrid(1),ngrid(2),ngrid(3));
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
save('restart1.mat', '-v7.3')
    
    %load('restart1.mat')
    %energy_step=input(11+nTemp);
    %energy_cutoff=input(12+nTemp);
    tunnel_out='tunnel_info.out';
    tunnel_file=fopen(tunnel_out,'w');
    disp('Starting TS search.......')
    removed_clusters=[];
    for level=level_min:level_stop %level_max-1
        E_volume(level)=sum(sum(sum(level_matrix<=level)));
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
            for iC=1:nCluster
                if list.C(iC).ID==0
                    filled(iC)=1;
                end
            end
            while sum(filled)<nCluster
                for iC=1:nCluster
                    if filled(iC)==0 && sum(list.C(iC).ID)>0 %check that cluster has not been removed
                        %length_list=length(list.C(iC).info(:,1));
                        index_boundary=find(list.C(iC).info(:,5)==1); %list index for boundary points in cluster
                        index_temp=0;
                        for boundary_list=1:length(index_boundary)
                            %if list.C(iC).info(index_list,5)==1 %check if grid point is a boundary point
                            index_list=index_boundary(boundary_list);
                            boundary=0;
                            i=list.C(iC).info(index_list,1);
                            j=list.C(iC).info(index_list,2);
                            k=list.C(iC).info(index_list,3);
                            [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(i,j,k,ngrid(1),ngrid(2),ngrid(3));
                            neighbor_coords=[i ip im j jp jm k kp km];
                            for l=1:3
                                for m=4:6
                                    for n=7:9
                                        %if C2remove==0
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
                                            
                                            [index_temp,temp,boundary,minID_matrix,list,TS_matrix,cross_matrix,tunnel_list,...
                                                TS_list,tunnel_cluster,tunnel_cluster_dim,C2remove]=...
                                                check_neighbors(i,j,k,in,jn,kn,level,iC,index_list,index_temp,temp,boundary,...
                                                minID_matrix,list,TS_matrix,level_matrix,E_matrix,cross_matrix,crossi,crossj,crossk,...
                                                tunnel_list,TS_list,tunnel_cluster,tunnel_cluster_dim,energy_step);
                                            if C2remove>0
                                                removed_clusters=[removed_clusters; [C2remove iC]];
                                                filled(C2remove)=1;
                                                if isempty(TS_list_all)==0
                                                    TS_list_all(TS_list_all(:,6)==C2remove,6)=iC;
                                                    TS_list_all(TS_list_all(:,7)==C2remove,7)=iC;
                                                    
                                                    TS2remove=find((TS_list_all(:,6)-TS_list_all(:,7))==0);
                                                    for iTS2remove=1:length(TS2remove)
                                                        iTS=TS2remove(iTS2remove);
                                                        TS_matrix(TS_list_all(iTS,1),TS_list_all(iTS,2),TS_list_all(iTS,3))=0;
                                                        %minID_matrix.C(TS_list_all(iTS,1),TS_list_all(iTS,2),TS_list_all(iTS,3))=iC;
                                                    end
                                                    TS_list_all(TS2remove,:)=[];
                                                end
                                            end
                                            
                                        end
                                    end
                                end
                            end
                            list.C(iC).info(index_list,5)=boundary; %%KOlla detta
                            %end
                        end
                        if index_temp==0% && C2remove==0
                            filled(iC)=1;
                        else
                            list.C(iC).info=[list.C(iC).info;temp.info];
                            temp=[]; %empty temp file
                        end
                    end
                    %end
                end
            end
        end
        %%%%Looking for new clusters!%%%
        [list,minID_matrix,cross_matrix,nCluster] =initiate_cluster(level,minID_matrix,list,E_matrix,level_matrix,cross_matrix,ngrid(1), ngrid(2), ngrid(3),nCluster);
        
        %list.TS=TS_list;
        TS_list_all=[TS_list_all; TS_list];
        echelon=rref(abs(tunnel_list));
        if ~isempty(echelon)
            for q=1:length(echelon(:,1))
		    if sum(abs(echelon(q,:)))==0
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
    BT_out = strcat('BT.dat');
    fid_BT = fopen(BT_out,'w');
    fprintf(fid_BT,'%s',num2str(BT));
    fclose(fid_BT);
save('restart2.mat', '-v7.3')
    
    %load('restart2.mat')
    %Print cluster list   %%%%remove cluster with single point??
    ifinC=0;
    list.finM=list.M;
    for iC=1:length(list.C)
        list.C(iC).tunnel=0;
        if sum(list.C(iC).ID)>0
            ifinC=ifinC+1;
            cluster_min=min(list.C(iC).info(:,7));
            index_min_cluster=find(list.C(iC).info(:,7)==min(list.C(iC).info(:,7)),1);
            coord_min_cluster=list.C(iC).info(index_min_cluster,1:3);
            list.C(iC).min=[cluster_min coord_min_cluster];
            list.finC(ifinC)=list.C(iC);
            list.finC(ifinC).info((list.finC(ifinC).info(:,6)==1),:)=[];
            if isempty(TS_list_all)==0
                TS_list_all(TS_list_all(:,6)==iC,6)=ifinC;
                TS_list_all(TS_list_all(:,7)==iC,7)=ifinC; %Update cluster index in TS_list
                tunnel_cluster(tunnel_cluster==iC)=ifinC; %Update cluster index in tunnel cluster list
            end
            for iM=1:length(list.M) %Update cluster index in Merge list
                list.finM(iM).clusters(list.finM(iM).clusters==iC)=ifinC;
            end
        else
            for iM=1:length(list.M) %remove cluster from merge list
                list.finM(iM).clusters(list.finM(iM).clusters==iC)=[];
            end
        end
    end
    if sum(BT)>0
        
        %%%%Find and organize tunnels%%%
        iT=0;
        for iM=1:length(list.finM)
            if sum(ismember(list.finM(iM).clusters(:),tunnel_cluster(:)))>0
                iT=iT+1;
                iT_list(iT)=iM;
                tunnel_dim_list=tunnel_cluster_dim(ismember(tunnel_cluster(:),list.finM(iM).clusters(:)),:);
                tunnel_dim=rref(tunnel_dim_list);
                list.tunnels(iT).dimensions=[sum(tunnel_dim(:,1)) sum(tunnel_dim(:,2)) sum(tunnel_dim(:,3))];
                list.tunnels(iT).clusters=list.finM(iM).clusters;
                n_tunnel_grids=0;
                tunnel_min=1000;
                for iCM_list=1:length(list.finM(iM).clusters)
                    iCM=list.finM(iM).clusters(iCM_list);
                    n_tunnel_grids=n_tunnel_grids+length(list.finC(iCM).info(:,1));
                    list.finC(iCM).tunnel=iT;
                    if list.finC(iCM).min(1)<tunnel_min
                        tunnel_min=list.finC(iCM).min(1);
                    end
                end
                list.tunnels(iT).TSgroup=[];
                list.tunnels(iT).min=tunnel_min;
                list.tunnels(iT).vol_proc=n_tunnel_grids/(ngrid(1)*ngrid(2)*ngrid(3));
            end
        end
        
        %group TS of cluster pairs into connected surfaces
        if isempty(TS_list_all)==0
            TS_list_all(:,[11 12])=0;
            [list,TS_list_all]=organize_TS(TS_list_all,ngrid(1),ngrid(2),ngrid(3),level_scale,beta,list,level,cross_matrix);
        end
        
        %%%%Update process list, TSgroup info and Tunnel info%%%
        iP=0;
        for iT=1:length(list.tunnels)
            list.tunnels(iT).TS=(TS_list_all(:,12)==iT); %Update tunnel TS list
            for iTSG_list=1:length(list.tunnels(iT).TSgroup)
                iTSG=list.tunnels(iT).TSgroup(iTSG_list);
                C1=list.TSgroup(iTSG).info.cluster1(1);
                C2=list.TSgroup(iTSG).info.cluster2(1);
                TS_min=list.TSgroup(iTSG).info.group(5);
                dE_1=list.TSgroup(iTSG).info.cluster1(6);
                dE_2=list.TSgroup(iTSG).info.cluster2(6);
                iP=iP+1;
                process_data(iP,:)=[C1 C2 dE_1 iT iTSG];
                iP=iP+1;
                process_data(iP,:)=[C2 C1 dE_2 iT iTSG];
            end
        end
        
        
        %%%Check process boundaries%%%%%
        disp('Checking process boundaries....')
        tboundary_start=cputime;
        clear i
        
        %checking which clusters are at boundary
        for iC=1:length(list.finC)
            if sum(list.finC(iC).info(:,1)==1)>1
                list.finC(iC).boundary=1;
            elseif sum(list.finC(iC).info(:,1)==ngrid(1))>1
                list.finC(iC).boundary=1;
            elseif sum(list.finC(iC).info(:,2)==1)>1
                list.finC(iC).boundary=1;
            elseif sum(list.finC(iC).info(:,2)==ngrid(2))>1
                list.finC(iC).boundary=1;
            elseif sum(list.finC(iC).info(:,3)==1)>1
                list.finC(iC).boundary=1;
            elseif sum(list.finC(iC).info(:,3)==ngrid(3))>1
                list.finC(iC).boundary=1;
            else
                list.finC(iC).boundary=0;
            end
        end
        
        %%%Check boundary crossing and process direction%%%%
        TS_tunnel_list=[];
        ind_group=0;
        processes=[];
        for i=1:2:length(process_data(:,1)) %processes come in forward/backward pairs
            C1=process_data(i,1);
            C2=process_data(i,2);
            iT=process_data(i,4);
            iTSG=process_data(i,5); %TSgroup ID
            coord_list=[list.finC(C1).info(:,1:3);list.TSgroup(iTSG).data(:,1:3);list.finC(C2).info(:,1:3)];
            TS_coord_list=[list.finC(C1).info(:,1:3);list.TSgroup(iTSG).data(:,1:3)];
            if list.finC(C1).boundary==1 || list.finC(C2).boundary==1
                start_min=list.finC(C1).min(2:4);
                end_min=list.finC(C2).min(2:4);
                TS_coord_min=list.TSgroup(iTSG).info.group(2:4);
                [process_cross_vector,coord_list]=get_TS_cross_vector(start_min,end_min,coord_list,ngrid,C1,C2);
                [TS_cross_vector,TS_coord_list]=get_TS_cross_vector(start_min,TS_coord_min,TS_coord_list,ngrid,C1,C2);
                
            else
                TS_cross_vector=[0 0 0];
                process_cross_vector=[0 0 0];
            end
            
            %%%Write processes%%%%
            processes(i,:)=[C1 C2 0 process_cross_vector TS_cross_vector iT iTSG C1 C2];%*
            processes(i+1,:)=[C2 C1 0 -process_cross_vector TS_cross_vector-process_cross_vector iT iTSG C2 C1];%*
            
            %%%write list for printing%%%%
            TS_coords=list.TSgroup(iTSG).data(:,1:3);
            TS_levels=list.TSgroup(iTSG).data(:,4);
            TS_groupID=ones(length(list.TSgroup(iTSG).data(:,1)),1)*iTSG;
            TS_zeroes=zeros(length(list.TSgroup(iTSG).data(:,1)),1);
            TS_tunnel_list=[TS_tunnel_list;[TS_groupID ceil(TS_levels*energy_step) TS_zeroes TS_coords-0.5]];
        end
        
        disp('...finished checking process boundaries.')
        ttot=cputime-tstart;
        
        %%%Print tunnel info%%%%
        
        disp(['Total time elapsed: ',num2str(ttot/60)]);
        fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
        tot_vol_proc=0;
        tot_n_TS=0;
        for iT=1:length(list.tunnels)
            if isempty(list.tunnels(iT))%.TS)
            else
                fprintf(tunnel_file,'%s\n','------------------------------------------------------------');
                fprintf(tunnel_file,'%s %d %s\n','Tunnel info ',iT,':');
                fprintf(tunnel_file,'%s %d %d %d\n','Breaks through in directions: ',list.tunnels(iT).dimensions);
                fprintf(tunnel_file,'%s %d %s\n','Min value:',list.tunnels(iT).min,' kJ/mol');
                fprintf(tunnel_file,'%s %d %s %d %s\n','This tunnel occupies ',list.tunnels(iT).vol_proc*100,'% of the total volume and has ',length(list.tunnels(iT).TSgroup),' transition states');
                fprintf(tunnel_file,'%s\n','------------------------------------------------------------');
                tot_vol_proc=tot_vol_proc+list.tunnels(iT).vol_proc*100;
                tot_n_TS=tot_n_TS+length(list.tunnels(iT).TSgroup);
                for iTSG=1:length(list.tunnels(iT).TSgroup)
                    fprintf(tunnel_file,'%s \n','group ID, min x, min y, min z, TS_min, B_min_abs, Bsum_TS');
                    fprintf(tunnel_file,'%s \n',num2str(list.TSgroup(iTSG).info.group));
                    fprintf(tunnel_file,'%s \n','cluster_1, min x, min y, min x, cluster1_min, dE_1');
                    fprintf(tunnel_file,'%s \n',num2str(list.TSgroup(iTSG).info.cluster1));
                    fprintf(tunnel_file,'%s \n','cluster_2, min x, min y, min x, cluster2_min, dE_2');
                    fprintf(tunnel_file,'%s \n',num2str(list.TSgroup(iTSG).info.cluster2));
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
        tunnel_cluster_list=unique(processes(:,1:2)); %Find all clusters involved in processes
        for cluster_row=1:length(tunnel_cluster_list)
            iC=tunnel_cluster_list(cluster_row);
            %%% removed weighted basis sites. Should be fixed and added%%%
            basis_sites(cluster_row,:)=[list.finC(iC).min(2:4) list.finC(iC).tunnel iC];
            update_cluster1_index=find(processes(:,1)==iC);
            processes(update_cluster1_index,1)=cluster_row;
            update_cluster2_index=find(processes(:,2)==iC);
            processes(update_cluster2_index,2)=cluster_row;
        end
        
        %%%Print basis%%%%
        basis_out = 'basis.dat';
        fid_basis = fopen(basis_out,'w');
        for i=1:length(basis_sites(:,1))
            fprintf(fid_basis,'%s\n',num2str(basis_sites(i,:)));
        end
        fclose(fid_basis);
        
        %%%%Print processes%%%
        
        for iTemp=1:nTemp
            temp=T(iTemp)
	    D_out = strcat('D_ave_',num2str(temp),'.dat');
            proc_out = strcat('processes_',num2str(temp),'.dat');
            fid_proc = fopen(proc_out,'w');
            list.temp(iTemp).processes=processes;
            for j=1:length(processes(:,1))
                %%%Compute boltzmann sums and rates for final processes%%%%
                %processes=[C1 C2 0 process_cross_vector TS_cross_vector iT iTSG C1 C2];%*
                iTSG=processes(j,11);
                C1=processes(j,12);
                C2=processes(j,13);
                Bsum_TS1=sum(exp(-1000*beta(iTemp)*(list.TSgroup(iTSG).data(:,5)-list.finC(C1).min(1))));
                Bsum_TS2=sum(exp(-1000*beta(iTemp)*(list.TSgroup(iTSG).data(:,5)-list.finC(C2).min(1))));
                Bsum_cluster1=sum(exp(-1000*beta(iTemp)*(list.finC(C1).info(:,7)-list.finC(C1).min(1))))+Bsum_TS1;
                Bsum_cluster2=sum(exp(-1000*beta(iTemp)*(list.finC(C2).info(:,7)-list.finC(C2).min(1))))+Bsum_TS2;
                k=prefactor(iTemp)*Bsum_TS1/(Bsum_cluster1*ave_grid_size*1e-10);
                %k2=prefactor(itemp)*Bsum_TS2/(Bsum_cluster2*ave_grid_size*1e-10);
                list.temp(iTemp).processes(j,3)=k;
                fprintf(fid_proc,'%s\n',num2str(list.temp(iTemp).processes(j,:)));
            end
            fclose(fid_proc);
            %figure(10)
            %hold on
            %plot(E_volume)
            Evol_out = strcat('Evol_',num2str(temp),'.dat');
            fid_Evol = fopen(Evol_out,'w');
            fprintf(fid_Evol,'%s',num2str(E_volume));
            fclose(fid_Evol);
            %%%Run KMC%%%%%
            if run_kmc==1
                disp('running KMC')
                tstart_kmc = cputime;
                
                if plot_msd==1
                    [msd,process_exec,distances,D_ave]=kmc_plot(basis_sites,temp,ngrid,n_particle,nsteps,msd_steps,print_every,grid_size,nruns,BT,per_tunnel);
                else
                    [msd,process_exec,distances,D_ave]=kmc_noplot(basis_sites,temp,ngrid,n_particle,nsteps,msd_steps,print_every,grid_size,nruns,BT,per_tunnel);
                end
                disp('Diffusion coefficients in')
                disp(strcat('A direction: ',num2str(D_ave(1)),'+/-',num2str(D_ave(2)),'cm^2/s'))
                disp(strcat('B direction: ',num2str(D_ave(3)),'+/-',num2str(D_ave(4)),'cm^2/s'))
                disp(strcat('C direction: ',num2str(D_ave(5)),'+/-',num2str(D_ave(6)),'cm^2/s'))
                fid_D = fopen(D_out,'w');
                fprintf(fid_D,'%s',num2str(D_ave));
                fclose(fid_D);
            end
        end
        t_kmc=cputime-tstart_kmc;
        ttot=cputime-tstart;
        disp(['Time kmc : ',num2str(t_kmc/60)]);
        fprintf(tunnel_file,'%s %d %s %f\n','Time kmc: ',t_kmc/60);
        disp(['Total time elapsed: ',num2str(ttot/60)]);
        fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
    else
        disp('Diffusion coefficients in')
        disp(strcat('A direction: ',num2str(0),'cm^2/s'))
        disp(strcat('B direction: ',num2str(0),'cm^2/s'))
        disp(strcat('C direction: ',num2str(0),'cm^2/s'))
	  for iTemp=1:nTemp
		      temp=T(iTemp)
		      D_out = strcat('D_ave_',num2str(temp),'.dat');        
        fid_D = fopen(D_out,'w');
        fprintf(fid_D,'%s',num2str([0 0 0 0 0 0]));
        fclose(fid_D);
end
        ttot=cputime-tstart;
        disp(['Total time elapsed: ',num2str(ttot/60)]);
        fprintf(tunnel_file,'%s %f\n','Total time elapsed: ',ttot/60);
    end
    
    for iC=1:length(list.finC)
        sites_min(iC,:)=[list.finC(iC).min sum(exp(-1000*beta(iTemp)*(list.finC(iC).info(:,7)))) list.finC(iC).tunnel iC];
    end
    %%%Print basis min info%%%%
    sites_min_out = 'sites_min.dat';
    fid_sites_min = fopen(sites_min_out,'w');
    for i=1:length(sites_min(:,1))
        fprintf(fid_sites_min,'%s\n',num2str([sites_min(i,:) (sites_min(i,2:4)-0.5)./ngrid]));
    end
    fclose(fid_sites_min);
    fclose(tunnel_file);
save('E_matrix.mat','E_matrix', '-v7.3')
save('list.mat','list', '-v7.3')


