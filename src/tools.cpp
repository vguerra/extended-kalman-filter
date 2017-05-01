#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse = VectorXd::Zero(4);

  if (estimations.size() == 0
      || estimations.size() != ground_truth.size()) {
    return rmse;
  }

  for(int i=0; i < estimations.size(); ++i) {
    VectorXd c = estimations[i] - ground_truth[i];
    VectorXd c_2 = c.array()*c.array();
    rmse += c_2;
  }

  rmse /= estimations.size();

  return rmse.array().sqrt();

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  if (fabs(px) < 0.001) {
    px = 0.001;
  }

  if (fabs(py) < 0.001) {
    py = 0.001;
  }

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //compute the Jacobian matrix
  MatrixXd Hj(3,4);
  Hj << (px/c2), (py/c2), 0, 0,
  -(py/c1), (px/c1), 0, 0,
  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}

VectorXd Tools::ToCartesian(const VectorXd& polar) {
  assert(polar.rows() == 3 && polar.cols() == 1);

  VectorXd x = VectorXd(4);

  float ro = polar(0);
  float phi = polar(1);
  float ro_dot = polar(2);

  float px = ro*cos(phi);
  float py = ro*sin(phi);
  float vx = ro_dot*cos(phi);
  float vy = ro_dot*sin(phi);

  x << px, py, vx, vy;
  return x;
}

VectorXd Tools::H_x(const Eigen::VectorXd& x_state) {
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float ro = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  float ro_dot = (px*vx + py*vy)/ro;

  VectorXd h = VectorXd(3);
  h << ro, phi, ro_dot;

  return h;
}

float Tools::NormalizeAngle(float angle) {
  int factor = angle > M_PI ? -1 : 1;

  while (fabs(angle) > M_PI) {
    angle += factor*2*M_PI;
  }

  return angle;
}
