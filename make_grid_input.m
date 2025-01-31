close all
clear all

cutoff=12.5;
grid_size=0.4;
masses=load('masses.dat');
natoms=load('natoms.dat');
natomtypes=load('natomtypes.dat');
elements=importdata('elements.dat');
epsilon=load('epsilon.dat');
sigma=load('sigma.dat');
epsilon_CH4=0.294106;
sigma_CH4=3.73000;
atoms=load('atoms.dat');
cell_SC=load('cell_SC.dat');
cell_UC=load('cell_UC.dat');

A_UC=cell_UC(1);
B_UC=cell_UC(2);
C_UC=cell_UC(3);
alpha=cell_UC(4)*pi/180;
beta=cell_UC(5)*pi/180;
gamma=cell_UC(6)*pi/180;

f=sqrt(1-(cos(alpha))^2-(cos(beta))^2-(cos(gamma))^2+2*cos(alpha)*cos(beta)*cos(gamma)); %volume conversion factor for non-orthogonal cells 
V_UC=A_UC*B_UC*C_UC*f; %cell volume
cellvector_UC=[A_UC 0 0; B_UC*cos(gamma) B_UC*sin(gamma) 0; C_UC*cos(beta) C_UC*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma) V_UC/(A_UC*B_UC*sin(gamma))];

UC_mult=round([cell_SC(1)/cellvector_UC(1,1) cell_SC(2)/cellvector_UC(2,2) cell_SC(3)/cellvector_UC(3,3)]);
nUC=UC_mult(1)*UC_mult(2)*UC_mult(3);
natoms_UC=natoms/nUC;

A_SC=A_UC*UC_mult(1);
B_SC=B_UC*UC_mult(2);
C_SC=C_UC*UC_mult(3);

V_SC=A_SC*B_SC*C_SC*f; %cell volume
cellvector_SC=[A_SC 0 0; B_SC*cos(gamma) B_SC*sin(gamma) 0; C_SC*cos(beta) C_SC*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma) V_SC/(A_SC*B_SC*sin(gamma))];


for i_atomtype=1:natomtypes
    atomtype.mass{i_atomtype}=masses(i_atomtype);
    atomtype.element{i_atomtype}=elements{i_atomtype};
    atomtype.epsilon{i_atomtype}=epsilon(i_atomtype);
    atomtype.sigma{i_atomtype}=sigma(i_atomtype);
end

Li_qtot=0;
for i_atom=1:natoms
    atom.type{i_atom}=atoms(i_atom,3);
    atom.charge{i_atom}=atoms(i_atom,4);
    atom.x{i_atom}=atoms(i_atom,5);
    atom.y{i_atom}=atoms(i_atom,6);
    atom.z{i_atom}=atoms(i_atom,7);
    element=atomtype.element{atom.type{i_atom}}(1,1:2);
    mass=atomtype.mass{atom.type{i_atom}};

end


field_file = fopen('FIELD','w');
fprintf(field_file,'%s\n','MOF_CH4_out');
fprintf(field_file,'%s\n','UNITS   kcal');
fprintf(field_file,'%s\n','molecular types 2');
fprintf(field_file,'%s\n','&guest UFF CH4');
fprintf(field_file,'%s\n','UNITS   kcal');
fprintf(field_file,'%s\n','NUMMOLS 0');
fprintf(field_file,'%s\n','ATOMS 1');
fprintf(field_file,'%s %12f %12f %12f %12f %12f\n','CH4',16.000000,0.00000,0.000000,0.000000,0.00000);
fprintf(field_file,'%s\n','rigid 1');
fprintf(field_file,'%s\n','1 1');
fprintf(field_file,'%s\n','finish');
fprintf(field_file,'%s\n','Framework');
fprintf(field_file,'%s %d\n','NUMMOLS',nUC);
fprintf(field_file,'%s %d\n','ATOMS',natoms_UC);

config_file = fopen('CONFIG','w');
fprintf(config_file,'%s\n','COF_CH4_out');
fprintf(config_file,'%12s %12d %12d\n',num2str(0),3,natoms);
fprintf(config_file,'%16f %12f %12f\n',cellvector_SC(1,:));
fprintf(config_file,'%16f %12f %12f\n',cellvector_SC(2,:));
fprintf(config_file,'%16f %12f %12f\n',cellvector_SC(3,:));
index=0;

for i_atom=1:natoms
    element=atomtype.element{atom.type{i_atom}}(1,1:2);
    fprintf(config_file,'%s %12s\n',element,num2str(index));
fprintf(config_file,'%16f %20f %20f\n',atom.x{i_atom},atom.y{i_atom},atom.z{i_atom});
end
for i_atom=1:natoms_UC
    mass=atomtype.mass{atom.type{i_atom}};
    element=atomtype.element{atom.type{i_atom}}(1,1:2);
    fprintf(field_file,'%s %12s %12s %8s %8s\n',element,num2str(mass),num2str(atom.charge{i_atom}),num2str(1),num2str(1));
end
fprintf(field_file,'%s\n','finish');

fprintf(field_file,'%s %d\n','VDW',natomtypes);
for j_atomtype=1:natomtypes
    element1=atomtype.element{j_atomtype}(1,1:2);
    sigma1=atomtype.sigma{j_atomtype};
    epsilon1=atomtype.epsilon{j_atomtype};
    sigma_lb=(sigma1+sigma_CH4)/2;
    epsilon_lb=sqrt(epsilon1*epsilon_CH4);
    sigma_lb=(sigma1+sigma_CH4)/2;
    fprintf(field_file,'%s %8s %8s %8s %8s\n',element1,'CH4','lj',num2str(epsilon_lb),num2str(sigma_lb));
end


fprintf(field_file,'%s\n','close');
fclose(field_file);
fclose(config_file);
control_file = fopen('CONTROL','w');
fprintf(control_file,'%s\n','GCMC Run');
fprintf(control_file,'%s\n','temperature 500');
fprintf(control_file,'%s\n','steps 1');
fprintf(control_file,'%s\n','equilibration 1');
fprintf(control_file,'%s %s %s\n','cutoff ',num2str(cutoff),' angstrom');
fprintf(control_file,'%s\n','delr            1.0 angstrom');
fprintf(control_file,'%s\n','overlap         0.0 angstrom');
fprintf(control_file,'%s\n','ewald precision  1d-6');
fprintf(control_file,'%s\n','averaging window 1');
fprintf(control_file,'%s\n','numguests 0');
fprintf(control_file,'%s\n','history 0');
fprintf(control_file,'%s\n','henry''s coefficient');
fprintf(control_file,'%s\n','widom ins 1');
fprintf(control_file,'%s\n','surface tol 2.5 Angstroms');
fprintf(control_file,'%s\n','&guest 1');
fprintf(control_file,'%s\n','  pressure 0.1');
fprintf(control_file,'%s\n','  probability 1');
fprintf(control_file,'%s\n','  energy');
fprintf(control_file,'%s\n','&end');
fprintf(control_file,'%s %d %d %d\n','grid factors',UC_mult);
fprintf(control_file,'%s %s\n','grid spacing ',num2str(grid_size));
fprintf(control_file,'%s\n','finish');

fclose(control_file);

