function [E_matrix_coarse]=make_coarse(E_matrix,ngrid_coarse,ngrid,R,T)
apos_grid=round([1:ngrid_coarse(1)]*ngrid(1)/ngrid_coarse(1));
bpos_grid=round([1:ngrid_coarse(2)]*ngrid(2)/ngrid_coarse(2));
cpos_grid=round([1:ngrid_coarse(3)]*ngrid(3)/ngrid_coarse(3));
for a=1:ngrid_coarse(1)
    for b=1:ngrid_coarse(2)
        for c=1:ngrid_coarse(3)
            E1=E_matrix(apos_grid(a),bpos_grid(b),cpos_grid(c));
            E2=E_matrix(apos_grid(a)-1,bpos_grid(b)-1,cpos_grid(c)-1);
            E_matrix_coarse(a,b,c)=-R*T/1000*log((exp(-E1*1000/(R*T))+exp(-E2*1000/(R*T)))/2);

        end
    end
end

E_matrix_coarse=E_matrix_coarse-min(min(min(E_matrix_coarse)));

end