function [pt1 pt2] = line_from_rho_theta(rho,theta_rad,dist)
  %Return the two endpoints determined by the rho, and theta
  pt = [rho*cos(theta_rad) rho*sin(theta_rad)];
  vec_perp = [sin(theta_rad) -cos(theta_rad)];
  pt1 = pt + dist*vec_perp;
  pt2 = pt - dist*vec_perp;
end
