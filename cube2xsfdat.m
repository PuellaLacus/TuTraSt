function [grid,grid_size,pot_data]=cube2xsfdat(E_unit)
disp('writing pot data...')
cube=importdata('grid.cube');
Natoms=cube.data(1,1);
if E_unit==1 %kJ/mol
    conv2kJmol=1;
elseif E_unit==2 %kcal/mol
    conv2kJmol=4.184;
elseif E_unit==3 %Ry
    conv2kJmol=1312.7497;
elseif E_unit==4 %eV
    conv2kJmol=96.4853;
elseif E_unit==5 %Hartree
    conv2kJmol=627.5096;
end
grid=[cube.data(2,1) cube.data(3,1) cube.data(4,1)];
a_grid=sqrt(cube.data(2,2)^2+cube.data(2,3)^2+cube.data(2,4)^2)*0.529177249;
b_grid=sqrt(cube.data(3,2)^2+cube.data(3,3)^2+cube.data(3,4)^2)*0.529177249;
c_grid=sqrt(cube.data(4,2)^2+cube.data(4,3)^2+cube.data(4,4)^2)*0.529177249;
grid_size=[a_grid b_grid c_grid]; %grid size in Angstrom
cube_data=cube.data(2*Natoms+5:end,:);
shift_cube_data=(cube_data-min(min(cube_data)))*conv2kJmol;
cube_size=size(cube_data);
line=0;
for row=1:cube_size(1)
    for column=1:cube_size(2)
        if shift_cube_data(row,column)>=0
            line=line+1;
            data(line)=shift_cube_data(row,column);
        end
    end
end
line=1;
for a=1:grid(1)
    for b=1:grid(2)
        for c=1:grid(3)
            data_mat(a,b,c)=data(line);
            line=line+1;
        end
    end
end
row=1;
column=1;
for c=1:grid(3)
    for b=1:grid(2)
        for a=1:grid(1)
            pot_data(row,column)=data_mat(a,b,c);
            if column==6
                column=1;
                row=row+1;
            else
                column=column+1;
            end
        end
    end
end
fid_xsf = fopen('xsf_data.dat','w');
for i=1:length(pot_data(:,1))
    fprintf(fid_xsf,'%f %30f %30f %30f %30f %30f\n',pot_data(i,:));
end
fclose(fid_xsf);
disp('writing pot data done')
end



