
var hs_km_ = [  // base altitude (km)
  0,
  25,
  30,
  40,
  50,
  60,
  70,
  80,
  90,
  100,
  110,
  120,
  130,
  140,
  150,
  180,
  200,
  250,
  300,
  350,
  400,
  450,
  500,
  600,
  700,
  800,
  900,
  1000
];
var hs_ = s_x_v(1e3, hs_km_);

var rho0s_ = [  // base density (kg/m2)
  1.225E+00,
  3.899E-02,
  1.774E-02,
  3.972E-03,
  1.057E-03,
  3.206E-04,
  8.770E-05,
  1.905E-05,
  3.396E-06,
  5.297E-07,
  9.661E-08,
  2.438E-08,
  8.484E-09,
  3.845E-09,
  2.070E-09,
  5.464E-10,
  2.789E-10,
  7.248E-11,
  2.418E-11,
  9.518E-12,
  3.725E-12,
  1.585E-12,
  6.967E-13,
  1.454E-13,
  3.614E-14,
  1.170E-14,
  5.245E-15,
  3.019E-15
];

var Hs_km_ = [ // scale height (km)
  7.249,
  6.349,
  6.682,
  7.554,
  8.382,
  7.714,
  6.549,
  5.799,
  5.382,
  5.877,
  7.263,
  9.473,
  12.636,
  16.149,
  22.523,
  29.740,
  37.105,
  45.546,
  53.628,
  53.298,
  58.515,
  60.828,
  63.822,
  71.835,
  88.667,
  124.640,
  181.050,
  268.000
];
var Hs_ = s_x_v(1e3, Hs_km_);


function get_density(h) {
  // Check if outside atmosphere model
  if (h <= hs_[0]) {
    return rho0s_[0];
  } else if (h >= hs_[hs_.length-1]) {
    return 0.;
  }

  // Get base idx to use
  var idx = 0;
  while (hs_[idx+1] < h) idx++;

  // Evaluate rho
  var h0 = hs_[idx];
  var rho0 = rho0s_[idx];
  var H = Hs_[idx];
  var rho = rho0 * exp(-(h - h0)/H);
  return rho;
}


// Based on V2 Rocket diagrams in Rocket Propulsion Elements
// but assumed constant w.r.t Mach
function get_Cd(alpha) {
  return 0.15 + abs(alpha*r2d)*0.4/10.;
}
function get_Cl(alpha) {
  return (alpha*r2d)*1.5/10.;
}

// Assumed constant with altitude and azimuth
function get_dv_wind_dispersion(h) {
  return 70.;
}


