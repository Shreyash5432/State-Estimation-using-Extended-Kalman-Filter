#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  VectorXd res(4);
 
  if((estimations.size() == 0) || (estimations.size()!=ground_truth.size()))
  {
      cout<<"RSME cannot be computed";
  }
  else
  {
    for (unsigned int i=0; i < estimations.size(); ++i) 
    {
    res = (estimations[i] - ground_truth[i]);
    res = res.array()*res.array();
    rmse += res;
    }
    // Calculating the mean
    rmse = rmse/estimations.size();
    // Calculating the squared root
    rmse = rmse.array().sqrt();
  }
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
  MatrixXd Hj(3,4);
  
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  float pxy = pow((px*px+py*py),(float(0.5)));
  float pxy2 = (px*px+py*py);
  float pxy32 = pow((px*px+py*py),(float(1.5)));
  
  if (px == 0 && py == 0)
  {
      cout<< "CalculateJacobian() - error - division by 0";
  }
  
  else
  {
      Hj << px/pxy, py/pxy, 0, 0,
            -py/pxy2, px/pxy2, 0, 0,
            (py*(vx*py-vy*px))/pxy32, (px*(vy*px-vx*py))/pxy32, px/pxy, py/pxy;
  }

  return Hj;
}
