#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
  * A helper to convert from polar to cartesian coordinates
  */
  static Eigen::VectorXd ToCartesian(const Eigen::VectorXd& polar);

  /**
   * A helper to compute h(x)
   */
  static Eigen::VectorXd H_x(const Eigen::VectorXd& x_state);

  /**
   * A helper to normalize angle between compute h(x)
   */
  static float NormalizeAngle(float angle);
};

#endif /* TOOLS_H_ */
