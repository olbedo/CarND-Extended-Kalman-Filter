#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() == 0) {
	  cout << "Error: Size of the estimation vector is zero" << endl;
	  return rmse;
	}

	if (estimations.size() != ground_truth.size()) {
	  cout << "Error: Sizes of the estimation vector and the ground_truth vector differ" << endl;
	  return rmse;
	}

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
    VectorXd diff = estimations.at(i) - ground_truth.at(i);
    rmse.array() += diff.array() * diff.array();
	}

	//calculate the mean
	rmse /= estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  if (px == 0 && py == 0) {
    cout << "Division by Zero" << endl;
    return Hj;
  }
  //compute the Jacobian matrix
  float px2_py2 = px * px + py * py;
  float px2_py2_sqrt = sqrt(px2_py2);
  float px2_py2_sqrt3 = px2_py2 * px2_py2_sqrt;
  Hj(0,0) = px/px2_py2_sqrt;
  Hj(0,1) = py/px2_py2_sqrt;
  Hj(0,2) = 0;
  Hj(0,3) = 0;
  Hj(1,0) = -py/px2_py2;
  Hj(1,1) = px/px2_py2;
  Hj(1,2) = 0;
  Hj(1,3) = 0;
  Hj(2,0) = py*(vx*py-vy*px)/px2_py2_sqrt3;
  Hj(2,1) = px*(vy*px-vx*py)/px2_py2_sqrt3;
  Hj(2,2) = px/px2_py2_sqrt;
  Hj(2,3) = py/px2_py2_sqrt;

  return Hj;
}
