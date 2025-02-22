#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // 1measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  
  // Finish initializing the FusionEKF.

  // Set the process and measurement noises
  noise_ax = 9.;
  noise_ay = 9.;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) 
{


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) 
  {
    // Initialize the state ekf_.x_ with the first measurement.
    // Create the covariance matrix.

    cout << "Initializing FusionEKF" << endl;

    MatrixXd P(4, 4);
    P << 100, 0, 0, 0,
	       0, 100, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;

    // Initialize transition matrix
    MatrixXd F(4, 4);
    F << 1, 0, 0, 0,
	       0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;

    // Initialize measurement matrix for laser measurements
    H_laser_ << 1, 0, 0, 0,
               0, 1, 0, 0;

    // Initialize ekf_ with the first state vector, 
    // estimated initial state covariance matrix,
    // and an empty matrix for Q
    MatrixXd Q(4,4);

    VectorXd x(4);
    if( measurement_pack.sensor_type_ == MeasurementPackage::RADAR ) 
    {
      // Convert radar from polar to cartesian coordinates and initialize state.
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      x << rho*cos(phi), rho*sin(phi), 0.f, 0.f;

      ekf_.Init(x, P, F, Hj_, R_radar_, Q);
      
    }
    else if( measurement_pack.sensor_type_ == MeasurementPackage::LASER ) 
    {
      // Initialize state.
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.f, 0.f;
      
      ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  // Update the state transition matrix F according to the new elapsed time.

  // Compute the time from the previous measurement in seconds.
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;
  
  if( dt > 0. ) 
  { 
    // Update the motion model matrix for a timestep dt.
    // We use a motion model with a constant velocity.
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;
  
    // Update the process noise covariance matrix for a timestep dt.
    // Our motion model uses Gaussian random accelerations in the x and y directions.
    float dt2 = dt*dt;
    float dt3 = dt2*dt;
    float dt4 = dt3*dt;
    float dt4over4 = dt4/4.;
    float dt3over2 = dt3/2.;
    ekf_.Q_ << dt4over4*noise_ax,                 0, dt3over2*noise_ax,                0,
			     0, dt4over4*noise_ay,                 0, dt3over2*noise_ay,
	     dt3over2*noise_ax,                 0,      dt2*noise_ax,                 0,
			     0, dt3over2*noise_ay,                 0,      dt2*noise_ay;
    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // Radar updates
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF( measurement_pack.raw_measurements_ );
  } 
  else 
  {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update( measurement_pack.raw_measurements_ );
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}