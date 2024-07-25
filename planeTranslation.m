function pT = planeTranslation(plV, params, m);


u_dot = plV.Lpsi;
v_dot = plV.Ltheta;
w_dot = plV.Lphi;

pT = m*([u_dot; v_dot; w_dot] * cross([p;q;r],[u;v;w]));

end
