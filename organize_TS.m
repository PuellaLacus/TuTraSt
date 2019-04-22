function [list, TS_list_all]=organize_TS(TS_list_all,xsize,ysize,zsize,level_scale,beta,list,level,cross_matrix)
first_cluster=min([TS_list_all(:,6);TS_list_all(:,7)]);
last_cluster=max([TS_list_all(:,6);TS_list_all(:,7)]);
group_index=0;
TS2delete=[];
for cluster1=first_cluster:last_cluster
    for cluster2=cluster1:last_cluster
        pair_index=[find(TS_list_all(:,6)==cluster1 & TS_list_all(:,7)==cluster2);...
            find(TS_list_all(:,6)==cluster2 & TS_list_all(:,7)==cluster1)];
        pair_index=sort(pair_index);

        for N=1:length(pair_index)
            TS_index=pair_index(N);
            if TS_list_all(TS_index,11)==0 %Find first unchecked TS in cluster1 cluster2 interface
                group_index=group_index+1; %initiate new TS group
                TS_surface_list=TS_index; %Initiate list for attached TS points
                q=1;
                p=1;
                while q<=p
                    TS_list_all(TS_surface_list(q),11)=group_index; %Add cluster group_index
                    i=TS_list_all(TS_surface_list(q),1); %get coordinates
                    j=TS_list_all(TS_surface_list(q),2);
                    k=TS_list_all(TS_surface_list(q),3);
                    [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(i,j,k,xsize,ysize,zsize);
                    neighbor_coords=[i ip im j jp jm k kp km]; %get neighbor coords
                    for l=1:3
                        for m=4:6
                            for n=7:9
                                checkneighbor=1; %check this!!
                                if l==1 && m==4 && n==7
                                else
                                    in=neighbor_coords(l); %neighbor coordinate
                                    jn=neighbor_coords(m);
                                    kn=neighbor_coords(n);
                                    neighbor_match=(TS_list_all(:,1)==in & TS_list_all(:,2)==jn & TS_list_all(:,3)==kn);
                                    if sum(neighbor_match) %check if neighbor is a transition state
                                        neighbor_index=find(neighbor_match); %get the index in TS_list_all
                                        
                                        if length(neighbor_index)>1
                                            if TS_list_all(neighbor_index,10)==0
                                                TS2delete=[TS2delete;neighbor_index(2:end)]; %Find all redundant points
                                                TS_list_all(TS2delete,11)=-1; %Set group ID of redundant points to -1 for deletion                                            
                                            end
                                        end
                                        neighbor_index=neighbor_index(1);
                                       
                                        
                                        if sum(neighbor_index==pair_index) %check if the neighbor index is in pair index list
                                            if TS_list_all(neighbor_index,11)==0 %check that the neighbor has not already been added to group
                                                TS_list_all(neighbor_index,11)=group_index; %Add group ID
                                                TS_surface_list=[TS_surface_list neighbor_index]; %Add neighbor_index to group list
                                                p=length(TS_surface_list);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if q==p
                        
                        break
                    else
                        q=q+1;
                    end
                end
                TS_list_group=TS_list_all(TS_surface_list,:);
                TS_min=min(TS_list_group(:,5)); %get minimum TS value in cluster
                index_TS_min=find(TS_list_group(:,5)==min(TS_list_group(:,5)),1); %find listindex of TS min
                coord_TSgroup_min=TS_list_group(index_TS_min,1:3);
                cluster1_min=min(list.C(cluster1).min(1));
                coord_cluster1_min=list.C(cluster1).min(2:4);
                cluster2_min=min(list.C(cluster2).min(1));
                coord_cluster2_min=list.C(cluster2).min(2:4);
                dE_1=TS_min-cluster1_min;
                dE_2=TS_min-cluster2_min;
                list.TSgroup(group_index).data=TS_list_group;
                list.TSgroup(group_index).info.group(1:5)=[group_index coord_TSgroup_min(1) coord_TSgroup_min(2) coord_TSgroup_min(3) TS_min];
                list.TSgroup(group_index).info.cluster1=[cluster1 coord_cluster1_min(1) coord_cluster1_min(2) coord_cluster1_min(3) cluster1_min dE_1];
                list.TSgroup(group_index).info.cluster2=[cluster2 coord_cluster2_min(1) coord_cluster2_min(2) coord_cluster2_min(3) cluster2_min dE_2];
                list.TS_all=TS_list_all;
            end
        end
        pair_index=[];
    end
end
TS_list_all(TS2delete,:)=[];

end


