function [ip, im, jp, jm, kp, km, cross]=pbc3d(i,j,k,x,y,z)
cross=[0 0 0 0 0 0 0 0 0];
if i==x
    im=i-1;
    ip=1;
    cross(2)=1;
elseif i==1
    im=x;
    ip=i+1;
    cross(3)=-1;
else
    ip=i+1;
    im=i-1;
    cross(2)=0;
    cross(3)=0;
end
if j==y
    jm=j-1;
    jp=1;
    cross(5)=1;
elseif j==1
    jm=y;
    jp=j+1;
    cross(6)=-1;
else
    jp=j+1;
    jm=j-1;
    cross(5)=0;
    cross(6)=0;
end
if k==z
    km=k-1;
    kp=1;
    cross(8)=1;
elseif k==1
    km=z;
    kp=k+1;
    cross(9)=-1;
else
    kp=k+1;
    km=k-1;
    cross(8)=0;
    cross(9)=0;
end

end
