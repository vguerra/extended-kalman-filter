#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
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

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ = MatrixXd::Zero(2, 4);

  noise_ax_ = 9;
  noise_ay_ = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    previous_timestamp_ = measurement_pack.timestamp_;

    MatrixXd Q = MatrixXd::Zero(4, 4);
    MatrixXd P = MatrixXd(4, 4);

    P << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

    MatrixXd F = MatrixXd(4, 4);
    F << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;


    VectorXd x = VectorXd(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x = Tools::ToCartesian(measurement_pack.raw_measurements_);
      ekf_.Init(x, P, F, Hj_, R_radar_, Q);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
      ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // only if dt has an acceptable value ( > 0.001 ).
  if (dt > 0.001) {
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    ekf_.Q_(0, 0) = dt_4/4*noise_ax_;
    ekf_.Q_(0, 2) = dt_3/2*noise_ax_;
    ekf_.Q_(1, 1) = dt_4/4*noise_ay_;
    ekf_.Q_(1, 3) = dt_3/2*noise_ay_;
    ekf_.Q_(2, 0) = dt_3/2*noise_ax_;
    ekf_.Q_(2, 2) = dt_2*noise_ax_;
    ekf_.Q_(3, 1) = dt_3/2*noise_ay_;
    ekf_.Q_(3, 3) = dt_2*noise_ay_;
    
    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = Tools::CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
