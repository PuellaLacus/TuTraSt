function [cross_vector,coord_list]=get_TS_cross_vector(start_min,end_min,coord_list,grid,C1,C2)
if C1==C2
    disp('Tunnel has only one cluster')
    cross_vector=[0 0 0]; %%%Should this not be direction of tunnel?
    length_list=length(coord_list);
    coord_list=[coord_list zeros(length_list,4)];
    start_pos=find(start_min(1)==coord_list(:,1) & start_min(2)==coord_list(:,2) & start_min(3)==coord_list(:,3),1);
    index_list=start_pos;
    coord_list(start_pos,4)=1;
    index_pos=1;
    endmin_found=0;
    while index_pos<=length(index_list)
        %if coord_list(index_list(index_pos),4)==0 %check if grid point is unchecked
        i=coord_list(index_list(index_pos),1);
        j=coord_list(index_list(index_pos),2);
        k=coord_list(index_list(index_pos),3);
        [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(i,j,k,grid(1),grid(2),grid(3));
        neighbor_coords=[i ip im j jp jm k kp km];
        for l=1:3
            for m=4:6
                for n=7:9
                    if endmin_found==0
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
                        if checkneighbor==1
                            in=neighbor_coords(l); %neighbor coordinate
                            jn=neighbor_coords(m);
                            kn=neighbor_coords(n);
                            if sum(in==coord_list(:,1) & jn==coord_list(:,2) & kn==coord_list(:,3))>0
                                neighbor_pos=find(in==coord_list(:,1) & jn==coord_list(:,2) & kn==coord_list(:,3),1);
                                crossi=coord_list(index_list(index_pos),5)+cross_vec(l); %periodic image coordinate of neighbor
                                crossj=coord_list(index_list(index_pos),6)+cross_vec(m);
                                crossk=coord_list(index_list(index_pos),7)+cross_vec(n);
                                if coord_list(neighbor_pos,4)==0
                                    coord_list(neighbor_pos,5:7)=[crossi crossj crossk];
                                    index_list=[index_list neighbor_pos];
                                    coord_list(neighbor_pos,4)=1;
                                else
                                    cross_diff=coord_list(neighbor_pos,5:7)-[crossi crossj crossk];
                                    if abs(sum(cross_diff))>0
                                        disp('crossing has been found')
                                        cross_vector=cross_diff
                                        endmin_found=1;
                                    end
                                end                                
                            end
                        end
                    end
                end
            end
        end
        if endmin_found==1
            disp('End point is found!')
            disp([C1 C2 cross_diff])
            break
        elseif index_pos==length(index_list)
            disp([C1 C2])
            disp('End point was not found!!!!!!!!!!!!')
            break
        else
            index_pos=index_pos+1;
        end
    end
    
else
    cross_vector=[0 0 0];
length_list=length(coord_list(:,1));
    coord_list=[coord_list zeros(length_list,4)];
    start_pos=find(start_min(1)==coord_list(:,1) & start_min(2)==coord_list(:,2) & start_min(3)==coord_list(:,3),1);
    index_list=start_pos;
    coord_list(start_pos,4)=1;
    index_pos=1;
    endmin_found=0;
    while index_pos<=length(index_list)
        %if coord_list(index_list(index_pos),4)==0 %check if grid point is unchecked
        i=coord_list(index_list(index_pos),1);
        j=coord_list(index_list(index_pos),2);
        k=coord_list(index_list(index_pos),3);
        [ip, im, jp, jm, kp, km, cross_vec]=pbc3d(i,j,k,grid(1),grid(2),grid(3));
        neighbor_coords=[i ip im j jp jm k kp km];
        for l=1:3
            for m=4:6
                for n=7:9
                    if endmin_found==0
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
                        if checkneighbor==1
                            in=neighbor_coords(l); %neighbor coordinate
                            jn=neighbor_coords(m);
                            kn=neighbor_coords(n);
                            if sum(in==coord_list(:,1) & jn==coord_list(:,2) & kn==coord_list(:,3))>0
                                neighbor_pos=find(in==coord_list(:,1) & jn==coord_list(:,2) & kn==coord_list(:,3),1);
                                if coord_list(neighbor_pos,4)==start_min
                                elseif coord_list(neighbor_pos,4)==0
                                    crossi=coord_list(index_list(index_pos),5)+cross_vec(l); %periodic image coordinate of neighbor
                                    crossj=coord_list(index_list(index_pos),6)+cross_vec(m);
                                    crossk=coord_list(index_list(index_pos),7)+cross_vec(n);
                                    coord_list(neighbor_pos,5:7)=[crossi crossj crossk];
                                    index_list=[index_list neighbor_pos];
                                    coord_list(neighbor_pos,4)=1;
                                    if [in jn kn]==end_min
                                        if endmin_found==0
                                            cross_vector=[crossi crossj crossk]
                                        end
                                        endmin_found=1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        if endmin_found==1
            disp('End point was found!')
            disp([C1 C2])
            break
        elseif index_pos==length(index_list)
            disp('End point was not found!!!!!!!!!!!!')
            disp([C1 C2])
            break
        else
            index_pos=index_pos+1;
        end
    end
end

end
    
