function [list,minID_matrix,cross_matrix,nCluster] =initiate_cluster(level,minID_matrix,list,E_matrix,level_matrix,cross_matrix,xsize, ysize, zsize,nCluster)    
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
                if nCluster==67
                    disp('check problem cluster')
                end
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
    end