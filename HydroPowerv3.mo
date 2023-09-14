model HydroPowerv3
  constant Real g = 9.81 "gravity, m/s2";
  parameter Real rho = 997 "density, kg/m3";
  parameter Real mu = 0.89 * 10 ^ (-3) "dynamic viscosity of water, Pa.s";
  parameter Real eps = 1.5 * 10 ^ (-5) "pipe roughness height, m";
  parameter Real h_s0 = 82 "min water height in the surge tank, m";
  parameter Real H_r = 60 "resurvior height, m";
  parameter Real w_0 = 91.6;
  parameter Real C_v = 5.7 "minor loss coeff. of the gate";
  parameter Real V_rate = 24.6 "volum. flow rate through the turbine, m3/s";
  parameter Real W_rate = 103 * 10 ^ 6 "Rated electrical power generated, W";
  parameter Real H_net = H_r + H_i + H_d + H_p - H_t;
  parameter Real p_a = 1.013 * 10 ^ 5;
  //disturbe
  parameter Real t_start_u = 500.0;
  parameter Real t_end_u = 550.0;
  parameter Real u_t_start = 0.6;
  parameter Real u_t_end = 0.99;
  parameter Real W_ed_start = W_rate;
  parameter Real k_ed = 1.76;
  parameter Real W_ed_end = k_ed * W_rate;
  parameter Real t_start_W = 800;
  parameter Real t_end_W = 810;
  // surge tank
  parameter Real z_s = 120 "Vertical component of the length of the surge shaft, m";
  parameter Real x_s = 68 "Horizontal component of the length of the surge shaft, m";
  parameter Real cos_theta_s = z_s / (z_s ^ 2 + x_s ^ 2) ^ 0.5;
  parameter Real D_s = 3.4 "Diameter of the surge shaft, m";
  parameter Real A_s = D_s ^ 2 * Modelica.Constants.pi / 4;
  // intake
  parameter Real H_i = 22.5 "Height over which water fall in the intake race, m";
  parameter Real L_i = 6600 "length of the intake, m";
  parameter Real cos_theta_i = H_i / L_i;
  parameter Real p_n1 = p_a + rho * g * H_r;
  parameter Real D_i = 5.8 "diametr of the intake, m";
  parameter Real A_i = D_i ^ 2 * Modelica.Constants.pi / 4;
  // panstock
  parameter Real H_p = 429 "net height, m";
  parameter Real L_p = 600 "length of the penstock, m";
  parameter Real cos_theta_p = H_p / L_p;
  parameter Real D_p = 3.3 "diametr of the penstock, m";
  parameter Real A_p = D_p ^ 2 * Modelica.Constants.pi / 4;
  parameter Real V_p = L_p * A_p;
  // discharge
  parameter Real H_t = 10 "Height of water in the tail, m";
  parameter Real p_d2 = p_a + rho * g * H_t;
  parameter Real H_d = 4.5 "Height over which water fall in the discharge race";
  parameter Real L_d = 600 "Length of the discharge race, m";
  parameter Real cos_theta_d = H_d / L_d;
  parameter Real D_d = 5.8 "Diameter of the discharge race, m";
  parameter Real A_d = D_d ^ 2 * Modelica.Constants.pi / 4;
  parameter Real V_d = L_d * A_d;
  // agregate
  parameter Real J = 183000 "Moment of inertia of the aggregate, kg-m2";
  parameter Real k_b = 0.5 * 43.86 * w_0 "Friction factor in the aggregate bearing box, W-s3/rad3";
  parameter Real D_t1 = 1.85 "Turbine blade inlet diameter, m";
  parameter Real D_t2 = 1.48 "Turbine blade outlet diameter, m";
  parameter Real D_t_ = 0.5 * (D_t1 + D_t2);
  parameter Real B_t = 0.238 "Turbine blade width, m";
  parameter Real theta_h = 0.93 "Turbine hydraulic efficiency";
  parameter Real theta_e = 0.99 "Electricity generator efficiency";
  parameter Real alpha = 25 * Modelica.Constants.pi / 180;
  Real m_s, m_i, m_p;
  Real m_s_dot, m_p_dot, m_i_dot;
  Real R_e_s, R_e_i, R_e_p, R_e_d;
  Real f_s, f_i, f_p, f_d;
  Real F_f_s, F_f_i, F_f_p, F_f_d;
  Real M_s, M_i, M_p, M_d;
  Real p_n, dp_v, p_tr2, p_p2, p_d1, p_tr1, dp_tr;
  Real v_p, v_s, v_d;
  Real u_t, W_ed;
  Real K_a, K_tr1_dot, K_p_dot, K_d1_dot;
  Real W_fa, W_g, W_ts_dot;
  Real h_s, V_p_dot, V_i_dot, V_s_dot, V_d_dot, w;
  //Real h_s(start = h_s0), V_p_dot(start = V_rate), V_i_dot(start = V_rate), V_s_dot(start = 0), V_d_dot, w(start = w_0);
initial equation
  V_p_dot = V_rate;
  V_i_dot = V_rate;
  V_s_dot = 0;
  h_s = h_s0;
  w = w_0;
equation
  if time < t_start_u then
    u_t = u_t_start;
  elseif time < t_end_u then
    u_t = u_t_start + (time - t_start_u) / (t_end_u - t_start_u) * (u_t_end - u_t_start);
  else
    u_t = u_t_end;
  end if;
  if time < t_start_W then
    W_ed = W_ed_start;
  elseif time < t_end_W then
    W_ed = W_ed_start + (time - t_start_W) / (t_end_W - t_start_W) * (W_ed_end - W_ed_start);
  else
    W_ed = W_ed_end;
  end if;
// surge tank
  m_s = rho * A_s * h_s / cos_theta_s;
  m_s_dot = rho * V_s_dot;
  m_i_dot = rho * V_i_dot;
  m_p_dot = rho * V_p_dot;
  M_s = m_s * v_s;
  v_s = V_s_dot / A_s;
  R_e_s = rho * abs(v_s + 1 * 10 ^ (-10)) * D_s / mu;
  1 / f_s = (-2 * Modelica.Math.log10(eps / 3.7 / D_s + 5.74 / R_e_s ^ 0.9)) ^ 2;
  F_f_s = 0.5 * Modelica.Constants.pi * f_s * rho * v_s * abs(v_s) * h_s / cos_theta_s * D_s / 4;
  der(m_s) = m_s_dot;
  m_s_dot = m_i_dot - m_p_dot;
  der(M_s) = rho * V_s_dot ^ 2 / A_s + (p_n - p_a) * A_s - F_f_s - m_s * g * cos_theta_s;
// intake
  M_i = rho * L_i * V_i_dot;
  R_e_i = rho * abs(V_i_dot / A_i) * D_i / mu;
  1 / f_i = (-2 * Modelica.Math.log10(eps / 3.7 / D_i + 5.74 / R_e_i ^ 0.9)) ^ 2;
  F_f_i = 0.5 * Modelica.Constants.pi * f_i * rho * L_i * V_i_dot * abs(V_i_dot) / A_i ^ 2 * D_i / 4;
  m_i = rho * A_i * L_i;
  der(M_i) = (p_n1 - p_n) * A_i - F_f_i + m_i * g * cos_theta_i;
// penstock
  M_p = rho * V_p * v_p;
  v_p = V_p_dot / A_p;
  R_e_p = rho * abs(v_p) * D_p / mu;
  1 / f_p = (-2 * Modelica.Math.log10(eps / 3.7 / D_p + 5.74 / R_e_p ^ 0.9)) ^ 2;
  F_f_p = 0.5 * Modelica.Constants.pi * f_p * rho * L_p * v_p * abs(v_p) * D_p / 4;
  m_p = rho * V_p;
  der(M_p) = (p_n - p_p2) * A_p - F_f_p + m_p * g * cos_theta_p;
// discharge
  M_d = rho * V_d * v_d;
  v_d = V_d_dot / A_d;
  R_e_d = rho * abs(v_d) * D_d / mu;
  1 / f_d = (-2 * Modelica.Math.log10(eps / 3.7 / D_d + 5.74 / R_e_d ^ 0.9)) ^ 2;
  F_f_d = 0.5 * Modelica.Constants.pi * f_d * rho * L_d * v_d * abs(v_d) * D_d / 4;
  V_d_dot = V_p_dot;
  der(M_d) = (p_d1 - p_d2) * A_d - F_f_d + rho * V_d * cos_theta_d;
//guide vane
//dp_v = rho * C_v * (V_p_dot/A_p) * abs(V_p_dot/A_p) / 2 / u_t^2;
  V_p_dot = C_v * u_t * (dp_v / p_a) ^ 0.5;
  dp_v = p_p2 - p_tr1;
//turbine
  K_tr1_dot = K_p_dot + dp_v * V_p_dot;
  K_p_dot = rho / 2 / A_p ^ 2 * V_p_dot ^ 3;
  W_ts_dot = theta_h * K_tr1_dot;
  // dp_tr * V_p_dot = W_ts_dot;
  dp_tr * V_p_dot = 0;
  dp_tr = p_tr1 - p_tr2;
//draft tube
  p_d1 * V_p_dot = p_tr2 * V_p_dot + K_tr1_dot - K_d1_dot;
  K_d1_dot = rho / 2 / A_d ^ 2 * V_p_dot ^ 3;
//aggregate
  K_a = 0.5 * J * w ^ 2;
  W_fa = k_b * w ^ 2;
//W_fa = 0.5 * k_b * w^3;
  W_g = W_ed / theta_e;
  der(K_a) = W_ts_dot - W_fa - W_g;
  annotation(experiment(StartTime = 0, StopTime = 2000, Tolerance = 1e-06, Interval = 0.4));
end HydroPowerv3;