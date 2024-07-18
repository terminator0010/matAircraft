function extForces
LEB = [x; y; z];

m_z = 0;

m_N = [0; 0; m_z]
W_N = m_N*.g;


LEB = [LBE]*(-1);

W_N_B = LEB*W_N;


GFM_N = [0;0];
AFM_Nm = [0;0];
TFM_Nm = [0;0];

GAT_FM_N = GFM_N + AFM_Nm + TFM_Nm ;

ExtForces_N = W_N_B + GAT_FM_N[1;2;3];

n = F(i)/W_N_B;

nx = Fx/W_N_B;
nx = (T-D)/W;

ny = y/W_N_B;

nz = (-1*Z)/W_N_B;

LF = [nx; ny; nz];

end



