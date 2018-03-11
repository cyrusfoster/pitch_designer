
/*
 
All units base SI [m], [s], [kg], [rad], unless otherwise specified

Assumptions:
-3D spherical Earth with entire trajectory in XY plane (so effectively 2D)
-non-rotating Earth

*/

var r_Earth = 6378100.;
var GM_Earth = 3.986004418e14;
var Z_ = [0., 0., 1.];


// Define vehicle: 2-stage, roughly modeled after F9
// Ref: http://www.spaceflight101.net/falcon-9-v11.html

var m_payload_ = 13e3;

var s2_ = {};
s2_.id = "s2";
s2_.model_aero = false;
s2_.f_thrust = 801e3;
s2_.mdot = -s2_.f_thrust/(345.*g0);
s2_.m_dry = 3.9e3 + m_payload_;
s2_.m_prop0 = 92.67e3;
s2_.dt = s2_.m_prop0/abs(s2_.mdot);
s2_.dt_step = 20.;  // Propagator step size

var s1_ = {};
s1_.id = "s1";
s1_.model_aero = true;
s1_.f_thrust = 9.*700e3;
s1_.mdot = -s1_.f_thrust/(300.*g0);
s1_.m_dry = 23.1e3 + s2_.m_dry + s2_.m_prop0;
s1_.m_prop0 = 395.7e3;
s1_.ref_area = sq(3.66*0.5)*pi;
s1_.dt = s1_.m_prop0/abs(s1_.mdot);
s1_.qalpha_max = 240. *1e3*d2r;  // 240 kN-deg ~= 5000 psf-deg
s1_.dt_step = 2.;  // Propagator step size


// Insertion criteria
var r_target_ = r_Earth + 400e3;
var v_target_ = sqrt(GM_Earth/r_target_);
var FPA_target_ = 0.*d2r;
var r_tol_ = 1.;
var v_tol_ = 0.001;
var FPA_tol_ = 1e-5;


function get_pitch(ts, pitchs, t) {
  var n = ts.length;
  if (t < ts[0]) {
    alert("ERROR: requested pitchprog t before t[0]: "+t);
  } else if (t > ts[n-1]) {
    alert("ERROR: requested pitchprog t after t[end]: "+t);
  }

  var i_lo = binary_search_lowerbound(ts, t);
  var i_hi = i_lo+1;
  var x = (t-ts[i_lo])/(ts[i_hi]-ts[i_lo]);
  var pitch = interp_lin(pitchs[i_lo], pitchs[i_hi], x);
  return pitch;
}


function eval_frame(vehicle, pitchprog, t, y) {
  var f = {};
  f.t = t;
  f.y = y.slice();
  f.R = y.slice(0,3);
  f.V = y.slice(3,6);
  f.m_prop = y[6];

  f.r = norm(f.R)  
  f.v = norm(f.V)  
  f.h = f.r - r_Earth;
  f.FPA = half_pi - anglepi(f.R, f.V);
  
  // Mass update
  f.m = vehicle.m_dry + f.m_prop;

  // Gravity
  f.A_grav = s_x_v(-GM_Earth/sq(f.r), unit(f.R));
  
  // Aerodynamic
  var model_aero = (f.v > 0.1 && f.h < 100e3 && vehicle.model_aero);
  f.pitch_unconstrained = get_pitch(pitchprog.t, pitchprog.pitch, t);
  var axis_unconstrained = rot_z(unit(f.R), half_pi-f.pitch_unconstrained);
  if (model_aero) {
    var density = get_density(f.h);
    f.q = 0.5 * density * sq(f.v);
    
    // Get qalpha constraint given wind dispersions
    // Ref. "Ascent Guidance Comparisons", J.M. Hanson, NASA-TM-112493
    var dv_wind_dispersion = get_dv_wind_dispersion(f.h);
    var dV_wind_dispersion_unit = unit(cross(f.R, Z_));
    var dV_wind_dispersion = s_x_v(dv_wind_dispersion, dV_wind_dispersion_unit);
    var V_dispersion = v_p_v(f.V, dV_wind_dispersion);
    var dalpha = anglepi(f.V, V_dispersion);
    var qdalpha = f.q * dalpha;
    f.qalpha_constraint = max(0., vehicle.qalpha_max - qdalpha);
   
    // Get pitch and vehicle axis 
    var alpha_constraint = f.qalpha_constraint/f.q;
    var alpha_unconstrained = anglepipi(Z_, axis_unconstrained, f.V); 
    f.qalpha_violation = abs(alpha_unconstrained) > alpha_constraint;
    if (use_qalpha_constraint_ && f.qalpha_violation) {
      var alpha_constrained = alpha_constraint * sign(alpha_unconstrained);
      f.qalpha_violation = false;
      f.axis = rot_z(unit(f.V), -alpha_constrained);
      f.pitch = half_pi - anglepipi(Z_, unit(f.R), f.axis);
    } else {
      f.pitch = f.pitch_unconstrained;
      f.axis = axis_unconstrained;
    }

    // Now, eval aerodynamic forces
    f.alpha = anglepipi(Z_, f.axis, f.V);  // angle of attack
    f.qalpha = f.q * f.alpha;
    var f_drag = f.q * get_Cd(f.alpha) * vehicle.ref_area;
    var f_lift = f.q * get_Cl(f.alpha) * vehicle.ref_area;
    f.F_drag = s_x_v(f_drag, unit(s_x_v(-1., f.V)))
    f.F_lift = s_x_v(f_lift, unit(cross(Z_, f.F_drag)));  
    f.F_aero = v_p_v(f.F_drag, f.F_lift)  
    f.A_aero = s_x_v(1./f.m, f.F_aero);

  } else {
    f.pitch = f.pitch_unconstrained;
    f.axis = axis_unconstrained;
    f.A_aero = [0., 0., 0.];
  }

  // Thrust
  var a_thrust = vehicle.f_thrust/f.m;
  f.A_thrust = s_x_v(a_thrust, f.axis);

  // Compile output
  f.A_proper = v_p_v(f.A_thrust, f.A_aero);
  f.A = v_p_v(f.A_grav, f.A_proper);
  f.ydot = concat(f.V, f.A);
  f.ydot.push(vehicle.mdot);
  return f;
}


// RK4 propagation given vehice, tspan and pitch commands
function prop_traj(vehicle, pitchprog, y0, tspan) {
  var get_frame = function(t, y) {
    return eval_frame(vehicle, pitchprog, t, y);
  }
  var t = tspan[0];
  var y = y0.slice();
  var frame = get_frame(t, y);
  var frames = [];
  frames.push(frame);
  
  var h = vehicle.dt_step;
  var stop_prop = false;
  while (!stop_prop) {
    // Check if we're reaching end of simulation
    if (t+h > tspan[1]) {
      h = tspan[1] - t;
      stop_prop = true;
    }

    // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    var k1 = frame.ydot;
    var k2 = get_frame(t+0.5*h, v_p_v(y, s_x_v(0.5*h, k1))).ydot;
    var k3 = get_frame(t+0.5*h, v_p_v(y, s_x_v(0.5*h, k2))).ydot;
    var k4 = get_frame(t+h,     v_p_v(y, s_x_v(    h, k3))).ydot;
    
    t += h;
    for (var i = 0; i < y.length; ++i) {
      y[i] = y[i] + (h/6.)*(k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]);
    }
    frame = get_frame(t, y);
    frames.push(frame);
  }

  return frames;
}


function eval_terminal_residuals(pitchprog, frame_s1_handover) {
  // Prop through end of s1
  var frames_s1_rest = 
    prop_traj(s1_, pitchprog, frame_s1_handover.y, [frame_s1_handover.t, s1_.dt]);
  var frame_s1_end = frames_s1_rest[frames_s1_rest.length-1];
  
  // Prop s2 to end of pitch program
  var y0_s2 = frame_s1_end.y;
  y0_s2[6] = s2_.m_prop0; 
  var t_end = pitchprog.t[pitchprog.t.length-1];
  var frames_s2 = prop_traj(s2_, pitchprog, y0_s2, [frame_s1_end.t, t_end]);
  var frame_s2_end = frames_s2[frames_s2.length-1];

  // Evaluate residuals
  var residuals = [
    frame_s2_end.r - r_target_,
    frame_s2_end.v - v_target_,
    frame_s2_end.FPA - FPA_target_,
  ];
  
  return {
    residuals: residuals,
    frames: concat(frames_s1_rest, frames_s2)
  };
}
      

function perturb_pitchprog(pitchprog, x) {
  // Create perturbed pitch program
  var pitchprog_delta = {};
  pitchprog_delta.t = pitchprog.t.slice();
  pitchprog_delta.pitch = pitchprog.pitch.slice();
  var n_pitch = pitchprog.t.length;
  pitchprog_delta.pitch[n_pitch-2] += x[0];
  pitchprog_delta.pitch[n_pitch-1] += x[1];
  pitchprog_delta.t[n_pitch-1] += x[2];
  return pitchprog_delta;
}


function solve_terminal_guidance(pitchprog, frame_s1_handover) {
  var h = 1e-5;
  var dx_max = [10.*d2r, 10.*d2r, 50.];
  for (var k = 0; k < 15; ++k) {
    // Nominal solution
    var r = eval_terminal_residuals(pitchprog, frame_s1_handover);
    var c = r.residuals;    

    // Check for convergence
    if (abs(c[0]) < r_tol_ && abs(c[1]) < v_tol_ && abs(c[2]) < FPA_tol_) {
      return {
        solved: true,
        frames: r.frames
      };
    }

    // Partial derivatives. Use finite forward difference
    var c_x = [[0.,0.,0.], [0.,0.,0.], [0.,0.,0.]];
    for (var i = 0; i < 3; ++i) {
      // Create perturbed pitch program
      var dx = [0., 0., 0.];
      dx[i] = h;
      var pitchprog_delta = perturb_pitchprog(pitchprog, dx);
      
      // Eval residuals of this perturbation
      var rp = eval_terminal_residuals(pitchprog_delta, frame_s1_handover);
      var cp = rp.residuals;

      // Eval partial derivative of this perturbation
      for (var j = 0; j < 3; ++j) {
        c_x[j][i] = (cp[j]-c[j])/h;
      }
    }

    // Solve
    var dx = solve_linear_system(c_x, s_x_v(-1., c));

    // Control step size
    var alpha = 1.;
    for (var i = 0; i < 3; ++i) {
      var update_ratio = abs(dx[i])/dx_max[i];
      if (update_ratio > 1.) alpha = min(alpha, 1./update_ratio);
    }
    dx = s_x_v(alpha, dx);

    // Update
    pitchprog = perturb_pitchprog(pitchprog, dx);
  }
  return {
    solved: false,
    frames: r.frames
  };
}


function channels_from_frames(frames) {
  channels = {};
  channels.t = {label:"time", data:[]};
  channels.h = {label:"altitude [km]", data:[], option:true};
  channels.a_proper = {label:"acceleration (proper) [g]", data:[], option:true};
  channels.pitch = {label:"pitch [deg]", data:[]};
  channels.pitch_unconstrained = {label:"pitch unconstrained [deg]", data:[]};
  channels.q = {label:"q [kPa]", data:[], option:true};
  channels.alpha = {label:"angle of attack [deg]", data:[], option:true};
  channels.qalpha = {label:"qalpha [kPa-deg]", data:[]};
  channels.qalpha_constraint = 
    {label:"qalpha constraint [kPa-deg]", data:[]};
  channels.qalpha_constraint_negative = 
    {label:"qalpha constraint (negative) [kPa-deg]", data:[]};
  channels.qalpha_violation = {label:"qalpha violation [0 or 1]", data:[]};
  //channels.rz = {label:"rz [km]", data:[]};
  for (var i = 0; i < frames.length; ++i) {
    channels.t.data.push(frames[i].t);
    channels.h.data.push(frames[i].h*1e-3);
    channels.a_proper.data.push(norm(frames[i].A_proper)/g0);
    channels.pitch.data.push(frames[i].pitch*r2d);
    channels.pitch_unconstrained.data.push(frames[i].pitch_unconstrained*r2d);
    channels.q.data.push(frames[i].q*1e-3);
    channels.alpha.data.push(frames[i].alpha*r2d);
    channels.qalpha.data.push(frames[i].qalpha*r2d*1e-3);
    channels.qalpha_constraint.data.push(frames[i].qalpha_constraint*r2d*1e-3);
    channels.qalpha_constraint_negative.data.push(
      -frames[i].qalpha_constraint*r2d*1e-3);
    channels.qalpha_violation.data.push(frames[i].qalpha_violation ? 1 : 0);
    //channels.rz.data.push(frames[i].R[2]*1e-3);
  }
  var i_end = frames.length-1;
  channels.t_final = {label:"mission duration [s]", 
                           data:round(frames[i_end].t, 3)};
  channels.m_prop_final = {label:"propellant margin [kg]", 
                           data:round(frames[i_end].m_prop, 1)};
  channels.q_max = {label:"max q [kPa]", 
                    data:round(max(remove_NaNs(channels.q.data)), 2)};
  channels.qalpha_violated = {
    label:"qalpha constraint", 
    data:(max(remove_NaNs(channels.qalpha_violation.data)) == 1) ? 
          "violated" : "ok"
  };
  return channels;
}


function run_pitch_designer(pitchprog_toggles) {
  // Compile pitch program
  var pitchprog = {};
 
  // Pitch profile through vertical ascent 
  pitchprog.t = [0., 10.];
  pitchprog.pitch = [90.*d2r, 90.*d2r];
  
  // Pitch profile through User-controlled toggles
  pitchprog.t = pitchprog.t.concat(pitchprog_toggles.t);
  pitchprog.pitch = pitchprog.pitch.concat(s_x_v(d2r, pitchprog_toggles.pitch));
  
  // Pitch profile through (initial guess for) terminal guidance
  var i_last_toggle = pitchprog.t.length-1;
  pitchprog.t = 
    pitchprog.t.concat([pitchprog.t[i_last_toggle], s1_.dt+s2_.dt]);
  pitchprog.pitch = 
    pitchprog.pitch.concat([pitchprog.pitch[i_last_toggle], -10.*d2r]);

  // Prop through toggled pitch profile
  var y0 = [r_Earth, 0., 0., 0., 0., 0., s1_.m_prop0];
  var frames_toggled = 
    prop_traj(s1_, pitchprog, y0, [0., pitchprog.t[i_last_toggle]]);
  
  // Prop through rest via closed-loop terminal guidance
  var frame_s1_handover = frames_toggled[frames_toggled.length-1];
  var tg = solve_terminal_guidance(pitchprog, frame_s1_handover);
  
  // Compile channels and return
  var frames = concat(frames_toggled, tg.frames);
  var channels = channels_from_frames(frames);
  channels.terminal_guidance = {label:"terminal guidance", 
                        data: tg.solved ? "converged" : "unconverged"};
  return channels;
}


