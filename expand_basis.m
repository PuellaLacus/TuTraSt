function [basis_exp n_image]=expand_basis(basis,ngrid,image)
basis_exp=basis
n_image=1
for a_image=-image(1):image(1)
    for b_image=-image(2):image(2)
        for c_image=-image(3):image(3)
            if a_image==0 && b_image==0 && c_image==0
            else
                n_image=n_image+1;
                basis_exp=[basis_exp;[basis(:,1)+a_image*ngrid(1) basis(:,2)+b_image*ngrid(2) basis(:,3)+c_image*ngrid(3) basis(:,4:9)]];
            end
        end
    end
end