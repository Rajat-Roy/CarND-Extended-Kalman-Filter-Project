#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in; // Object state
  P_ = P_in; // Object covariance matrix
  F_ = F_in; // State transiction matrix
  H_ = H_in; // Measurement matrix
  R_ = R_in; // Measurement covariance matrix
  Q_ = Q_in; // Process covariance matrix
}

void KalmanFilter::Predict() 
{
  // predict the state
  x_ = F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
  // Update the state using Kalman Filter equations
  VectorXd y = z - H_*x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Ht*Si;

  // New state
  x_ = x_ + ( K*y );
  MatrixXd I_ = MatrixXd::Identity(4,4);
  P_ = ( I_ - K*H_ )*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  H_ = tools.CalculateJacobian( x_ );
  VectorXd h(3);
  
  float c1 = sqrt(px*px + py*py); 
  
  if (fabs(c1) < 0.0001) {
    cout << "Division by Zero" << endl;
    return;
  }

  h << c1, atan2( py, px ),  (px*vx + py*vy )/c1;

  // Update the state using Extended Kalman Filter equations
  VectorXd y = z - h;

  while ( y(1) > M_PI ) 
    y(1) -= 2 * M_PI;

  while ( y(1) < -M_PI )
    y(1) += 2 * M_PI;

  MatrixXd Hjt = H_.transpose();
  MatrixXd S = H_*P_*Hjt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Hjt*Si;

  // Compute new state
  x_ = x_ + ( K*y );
  MatrixXd I_ = MatrixXd::Identity(4,4);
  P_ = ( I_ - K*H_ )*P_;
}