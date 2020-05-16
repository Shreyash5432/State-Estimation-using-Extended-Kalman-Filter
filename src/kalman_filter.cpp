#include "kalman_filter.h"
#include<cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;

}

void KalmanFilter::Predict() 
{
  
  x_ = F_*x_;
  P_ = F_*P_*(F_.transpose()) + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) 
{
  
  VectorXd y_ = z - H_ * x_;
  MatrixXd Ht_ = H_.transpose();
  MatrixXd S_ = H_ * P_ * Ht_ + R_;
  MatrixXd Si_ = S_.inverse();
  MatrixXd K_ =  P_ * Ht_ * Si_;
  
  MatrixXd I_ = MatrixXd::Identity(4, 4);
  x_ = x_ + (K_ * y_);
  P_ = (I_ - K_ * H_) * P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
  
  VectorXd z_pred_(3); 
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  
  z_pred_ << sqrt(px*px+py*py), atan2(py,px), (px*vx+py*vy)/sqrt(px*px+py*py);
  
  VectorXd y_ = z - z_pred_;

  while(y_[1] < (-M_PI) || y_[1] > M_PI)
  {
    if(y_[1] < (-M_PI))
    {
      y_[1]+= (2*M_PI);
    }
    else
    {
      y_[1]-= (2*M_PI);
    }
  }
  MatrixXd Ht_ = H_.transpose();
  MatrixXd S_ = H_ * P_ * Ht_ + R_;
  MatrixXd Si_ = S_.inverse();
  MatrixXd K_ =  P_ * Ht_ * Si_;
  
  MatrixXd I_ = MatrixXd::Identity(4, 4);
  x_ = x_ + (K_ * y_);
  P_ = (I_ - K_ * H_) * P_;
  
}
