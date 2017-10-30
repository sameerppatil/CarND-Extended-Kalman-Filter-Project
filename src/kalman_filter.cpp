#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>


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
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
  /**
  TODO:
    * predict the state
  */
}

void KalmanFilter::Update(const VectorXd &z)
{
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
}

VectorXd KalmanFilter::CalcZPred (const VectorXd &x)
{
    VectorXd calculated_z(3);
    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);

    float phi, rho, rho_dot;
    
    rho = sqrt(px*px + py*py);
    
    if (px == 0)
    {
        phi = 0.0;
    }
    else
    {
        phi = atan2(py, px);
    }
      
    if (rho == 0)
    {
        rho_dot = 0.0;
    }
    else
    {
       rho_dot = (px*vx + py*vy) / rho;
    }

    calculated_z << rho, phi, rho_dot;
    return calculated_z;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
   */
    MatrixXd Hj_ = tools.CalculateJacobian(x_);
    MatrixXd Hj_t = Hj_.transpose();
    MatrixXd S = Hj_ * P_ * Hj_t + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHj_t = P_ * Hj_t;
    MatrixXd K = PHj_t * Si;

    VectorXd y = z - CalcZPred(x_);

    // Noramalize phi here
    y(1) = atan2(sin(y(1)),cos(y(1)));

    // new estimate
    x_ = x_ + K * y;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj_) * P_;
}