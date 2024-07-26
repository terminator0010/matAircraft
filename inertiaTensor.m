function iT = inertiaTensor(Ixx,Ixz,Iyy,Izz,p,q,r)

 x = Ixx*p_dot - Ixz*r_dot - q*(Ixz*p - Izz*r) - Iyy*q*r;
 y = Iyy*q_dot + p*(Ixz*p - Izz*r) + r*(Ixx*p - Ixz*r);
 z = Izz*r_dot - Ixz*p_dot - q*(Ixx*p - Ixz*r) + Iyy*p*q;
iT = [x,y,z];
end
