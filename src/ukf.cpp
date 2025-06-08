#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;      // corrected due to what was learned in the lesson 

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;    // corrected due to what was learned in the lesso
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

   // first measurement?
  initialized_ = false;

  // definition of the state vectors dimensions: 2 positions, 3 velocitties
  // depends on the model chosen: CTRV
  n_x_ = 5; 

  // augmented state dimension: addining 2 more coming from noise vector
  n_aug_ = 7;

  // augmented sigma points 
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //predited sigma points 
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Define spreading parameter
  lambda_ = 3 - n_aug_;

    //init weights 
  weights_ = VectorXd(2*n_aug_+1);

  //Quality of the approximation: NIS calc Radar
  NIS_radar = 0.00;

  //Quality of the approximation: NIS calc Lidar
  NIS_lidar = 0.00;

  // Start time
  time_us_ = 0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  //first measument for initialisation of state and covariance matrices
  if(!initialized_)
  {
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      
      // Get the measured data for both dimensions
      double x = meas_package.raw_measurements_(0);
      double y = meas_package.raw_measurements_(1);

      // Init state vector
      x_ << x,y,0,0,0 ;

      // Init co-variance matrix based on sensor accuracy
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0,std_laspy_*std_laspx_,0,0,0,
            0,0,1,0,0,
            0,0,0,1,0,
            0,0,0,0,1;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Get the measured data for the three dimensions
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);

      // convert it into catresian coordinates
      double x = rho*sin(phi);
      double y = rho*cos(phi);

      // Init state vector
      x_ << x,y,rho_dot,phi,0;

      // Init co-variance matrix based on sensor accuracy
      P_ << std_radr_*std_radr_,0,0,0,0,
      0,std_radr_*std_radr_,0,0,0,
      0,0,std_radr_*std_radr_,0,0,
      0,0,0,std_radphi_*std_radphi_,0,
      0,0,0,0,1;
    }

    // Set Initialization to true
    initialized_ = true;

    // Get time 
    time_us_ = meas_package.timestamp_;

    return ;
  }

  // measurement of time delta
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;

  // new time
  time_us_ = meas_package.timestamp_;

  // do the prediction
  Prediction(delta_t);

  // do the update
  if(meas_package.sensor_type_==MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_==MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // START SIGMA POINTS AND ITS VALUES FOR POSTERIOR

   // 1. initalize the augmented state vector
    VectorXd x_aug = VectorXd(n_aug_); 
   
    // 2. init the augmented covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    // 3. fill the state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    // 4. fill the covariance
    // lazy implementation to set it all to default value of 0
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(5,6) = 0.0;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    P_aug(6,5) = 0.0;

    // 5. calculate square root matrix
    MatrixXd A = P_aug.llt().matrixL();

    // 6. create augmented sigma points add its values
    
    //1st is the x_aug
    Xsig_aug_.col(0) = x_aug ; 
    
    //next to be calculated by applying the above
    for(int i = 0 ; i < n_aug_ ; ++i)
    {
      Xsig_aug_.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*A.col(i) ;
      Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i) ;
    }

  // START PREDICTION
   for(int i = 0 ; i < 2*n_aug_+ 1 ; ++i)
    {
      //extract all values from the existing 
      double px = Xsig_aug_(0,i) ;
      double py = Xsig_aug_(1,i) ;
      double v = Xsig_aug_(2,i) ;
      double yaw = Xsig_aug_(3,i) ;
      double yawd = Xsig_aug_(4,i) ;
      double nu_a = Xsig_aug_(5,i) ;
      double nu_yawdd = Xsig_aug_(6,i) ;

      //predicted state values 
      double px_p , py_p ;
      
      if(fabs(yawd) > 0.001)
      {
          px_p = px + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw)) ; 
          py_p = py + v/yawd * (-cos(yaw + yawd * delta_t) + cos(yaw)) ; 
      }
      else
      {
          px_p = px + v*delta_t*cos(yaw) ; 
          py_p = py + v*delta_t*sin(yaw) ; 
      }
      
      double v_p = v ;
      double yaw_p = yaw + yawd*delta_t;
      double yawd_p = yawd ; 

      //add noise
      px_p = px_p + 0.5 *delta_t*delta_t * nu_a*cos(yaw) ; 
      py_p = py_p + 0.5*delta_t*delta_t * nu_a*sin(yaw) ;
      v_p = v_p + delta_t*nu_a;
      yaw_p = yaw_p + 0.5*delta_t*delta_t*nu_yawdd;
      yawd_p = yawd_p + delta_t*nu_yawdd ; 
  
      // prediction for the new state which has only physical properties (5 values) at the sigma points
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
    }
  
  // MEAN AND COVARIANCE for all results from the prediction step
    
     // set weights for the backward propagation into the mean and covariance
     // first position
     double weight_0 = lambda_/(lambda_+n_aug_);
     weights_(0) = weight_0;
     // rest of the weights
     for(int i=1 ; i<2*n_aug_+1;++i)
     {
       weights_(i) = 0.5/(lambda_+n_aug_) ;
     }
   
     // predict state mean
     x_.fill(0.0) ; 
     for(int i=0 ; i<2*n_aug_+1;++i)
     {
       x_ = x_ + weights_(i)*Xsig_pred_.col(i);
     }
   
     // predict state covariance matrix
     P_.fill(0.0);
     for(int i=0 ; i<2*n_aug_+1;++i)
     {
       //state difference
       VectorXd x_diff = Xsig_pred_.col(i) - x_ ;
       // angle normalization like mentioned in the solution of the lesson
       // M_PI is in C++ the number Pi
       while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
       while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
   
       P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ; 
     }

}

// Helper function to normalize an angle to be within [-PI, PI]
void NormalizeAngle(double& angle) {
    while (angle > M_PI) {
        angle -= 2. * M_PI;
    }
    while (angle < -M_PI) {
        angle += 2. * M_PI;
    }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
    // Measurement dimension for Lidar (px, py)
    const int n_z = 2;

    // 1. Transform sigma points into measurement space
    MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);
    // For Lidar, px and py are directly extracted from the predicted state sigma points.
    // If your Xsig_pred_ already has px and py in the first two rows,
    // this block operation directly gets them. No explicit loop is needed here
    // if Xsig_pred_ is structured such that its first two rows are px and py.
    // If not, the original loop was correct, but you should clarify what Xsig_pred_ contains.

    // 2. Calculate mean predicted measurement
    VectorXd z_pred = Zsig * weights_;

    // Calculate innovation covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z, n_z);

    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd z_pred_diff = Zsig.col(i) - z_pred;
        S += weights_(i) * z_pred_diff * z_pred_diff.transpose();
    }

    // Add measurement noise R (Lidar specific)
    MatrixXd R = MatrixXd::Zero(n_z, n_z);
    R(0, 0) = std_laspx_ * std_laspx_;
    R(1, 1) = std_laspy_ * std_laspy_;
    S += R;

    // 3. Update measurement
    // Create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // Residual for measurement space (Lidar doesn't have angles that need normalization here directly)
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // State difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        NormalizeAngle(x_diff(3)); // Normalize yaw angle (index 3)

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K
    MatrixXd K = Tc * S.inverse();

    // Incoming Lidar measurement
    VectorXd z = meas_package.raw_measurements_;

    // Residual between actual measurement and predicted measurement
    VectorXd z_diff = z - z_pred;
    // For Lidar (px, py), there's typically no angle normalization needed for z_diff.
    // If z_diff(1) in the original code was meant to be an angle from a different sensor,
    // that's a potential area for confusion. For Lidar, this is usually cartesian.
    // I've removed the angle normalization for z_diff here as it's not standard for Lidar (px, py).
    // If your 'py' measurement is somehow an angle, please clarify.

    // Update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    // Calculate NIS for Lidar
    NIS_lidar = z_diff.transpose() * S.inverse() * z_diff;

    std::cout << "NIS Lidar : " << NIS_lidar << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    // Measurement dimension for Radar (rho, phi, rho_dot)
    const int n_z = 3;

    // Create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

    // 1. Transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // Extract values from predicted sigma points
        double px = Xsig_pred_(0, i);
        double py = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        // Calculate common terms to avoid redundant computations
        double p_squared_sum = px * px + py * py;
        double rho = std::sqrt(p_squared_sum);

        // Avoid division by zero for rho_dot
        if (rho < 0.0001) { // A small threshold to prevent division by zero
            rho = 0.0001; // Set to a tiny value to proceed, though this case might need specific handling
        }

        // Measurement model: polar coordinates (rho, phi, rho_dot)
        Zsig(0, i) = rho;                                                // r (rho)
        Zsig(1, i) = std::atan2(py, px);                                 // phi (angle)
        Zsig(2, i) = (px * v * std::cos(yaw) + py * v * std::sin(yaw)) / rho; // r_dot (radial velocity)
    }

    // 2. Calculate mean predicted measurement
    VectorXd z_pred = Zsig * weights_; // More compact way to calculate weighted sum

    // Normalize predicted angle for z_pred (phi)
    NormalizeAngle(z_pred(1));

    // Calculate innovation covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z, n_z);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // Difference
        VectorXd z_pred_diff = Zsig.col(i) - z_pred;

        // Angle normalization for phi in the residual
        NormalizeAngle(z_pred_diff(1));

        S += weights_(i) * z_pred_diff * z_pred_diff.transpose();
    }

    // Add measurement noise R (Radar specific)
    MatrixXd R = MatrixXd::Zero(n_z, n_z);
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;
    S += R;

    // 3. Update measurement
    // Create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // Residual for measurement space
        VectorXd z_diff = Zsig.col(i) - z_pred;
        NormalizeAngle(z_diff(1)); // Angle normalization for phi

        // State difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        NormalizeAngle(x_diff(3)); // Angle normalization for yaw

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K
    MatrixXd K = Tc * S.inverse();

    // Incoming radar measurement
    VectorXd z = meas_package.raw_measurements_;

    // Residual between actual measurement and predicted measurement
    VectorXd z_diff = z - z_pred;
    NormalizeAngle(z_diff(1)); // Angle normalization for phi

    // Update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    // Calculate NIS for radar
    NIS_radar = z_diff.transpose() * S.inverse() * z_diff;

    std::cout << "NIS Radar : " << NIS_radar << std::endl;
}