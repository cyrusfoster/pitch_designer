
var pi = Math.PI;
var two_pi = 2.*pi;
var half_pi = 0.5*pi;

var r2d = 57.295779513082323;  // radians to degrees
var d2r = 0.017453292519943;  // degrees to radians

var g0 = 9.81;


function sq(a) {
  return a*a;
}


function cube(a) {
  return a*a*a;
}


function abs(a) {
  return Math.abs(a);
}


function pow(base, exponent) {
  return Math.pow(base, exponent);
}


function sqrt(x) {
  return Math.sqrt(x);
}


function floor(x) {
  return Math.floor(x);
}


function sin(x) {
  return Math.sin(x);
}


function cos(x) {
  return Math.cos(x);
}


function asin(x) {
  return Math.asin(x);
}


function acos(x) {
  return Math.acos(x);
}


function acos_safe(x) {
	if (abs(x) > 1. && abs(x) < 1.0001) x = sign(x);
	return acos(x);
}


function exp(x) {
  return Math.exp(x);
}


function sign(x) {
  return x == 0. ? 0. : (x > 0. ? 1. : -1.);
}


function min(x, y) {
  if (typeof y == "undefined") {
    return Math.min.apply(null, x);
  }
  return Math.min(x, y);
}


function max(x, y) {
  if (typeof y == "undefined") {
    return Math.max.apply(null, x);
  }
  return Math.max(x, y);
}


function clamp(x, lo, hi) {
  return Math.max(Math.min(x, hi), lo);
}


function remove_NaNs(x) {
  var x_new = [];
  for (var i = 0; i < x.length; ++i) {
    if (!isNaN(x[i])) x_new.push(x[i]);
  }
  return x_new;
}


function norm(x) {
  var out = 0.;
  for (var i = 0; i < x.length; ++i) {
    out += x[i]*x[i];
  }
  return sqrt(out);
}


function unit(v) {
  return s_x_v(1./norm(v), v);
}


function dot(x, y) {
  var out = 0.;
  for (var i = 0; i < x.length; ++i) {
    out += x[i]*y[i];
  }
  return out;
}


function cross(R1, R2) {
	return [R1[1]*R2[2]-R1[2]*R2[1], 
          R1[2]*R2[0]-R1[0]*R2[2], 
          R1[0]*R2[1]-R1[1]*R2[0]];
}


// Return angle on [-pi/2..pi/2]
function anglepi(R1, R2) {
  var theta = dot(R1, R2)/(norm(R1)*norm(R2));
  //console.log("theta="+theta+", acos="+acos(theta));
  return acos_safe(theta);
}


// Return angle on [-pi..pi], given polarity vector N
function anglepipi(N, R1, R2) {
  return sign(dot(N, cross(R1, R2))) * anglepi(R1, R2);
}

  
// scalar times vector
function s_x_v(s, v) {  
  out = v.slice();
  for (var i = 0; i < v.length; ++i) {
    out[i] *= s;
  }
  return out;
}


// vector plus vector
function v_p_v(v1, v2) {  
  out = v1.slice();
  for (var i = 0; i < v1.length; ++i) {
    out[i] += v2[i];
  }
  return out;
}
  

// rotate about z axis by angle
function rot_z(v, angle) {
  var c = cos(angle);
  var s = sin(angle);
  return [c*v[0]-s*v[1], s*v[0]+c*v[1], v[2]];
}


// Returns idx of lowerbound on x in array
function binary_search_lowerbound(array, x) {
  var left = 0;
  var right = array.length-1;
  while (right-left > 1) {
    var mid = floor((left+right)/2);
    if (array[mid] < x) {
      left = mid;
    } else {
      right = mid;
    }
  }
  return left;
}


function interp_lin(a, b, x) {
  return a + (b-a)*x;
}


function concat(v1, v2) {
  return v1.concat(v2);
}


// Return determinant of a 3x3 matrix
// ref: https://www.mathsisfun.com/algebra/matrix-determinant.html
function det33(A) {
  var a = A[0][0];
  var b = A[0][1];
  var c = A[0][2];
  var d = A[1][0];
  var e = A[1][1];
  var f = A[1][2];
  var g = A[2][0];
  var h = A[2][1];
  var i = A[2][2];
  return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}


// Solve Ax=b via Cramer's rule
// ref: https://en.wikipedia.org/wiki/Cramer%27s_rule
function solve_linear_system(A, b) {
  var d = det33(A);
  var dd = [0., 0., 0.];
  for (var i = 0; i < 3; ++i) {
    var Ap = [[0.,0.,0.], [0.,0.,0.], [0.,0.,0.]];
    for (var j = 0; j < 3; ++j) {
      for (var k = 0; k < 3; ++k) {
        Ap[j][k] = (k == i) ? b[j] : A[j][k];
      }
    }
    dd[i] = det33(Ap);
  }
  return [dd[0]/d, dd[1]/d, dd[2]/d];
}


function round(number, digits) {
  var n = Math.pow(10, digits)
  return Math.round(number*n)/n;
}


