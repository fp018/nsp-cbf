#include <qpOASES.hpp>
#include <OsqpEigen/OsqpEigen.h>

struct QPParams{
    
    // ========== Core Solver Settings ==========
    bool verbosity = false;
    bool warmStart = true;
    
    
    // ========== Active-Set Parameters ==========
    int maxWorkingSetRecalculations = 1000;       ///< qpOASES nWSR
    double terminationTolerance = 1e-6;          ///< qpOASES epsilon
    double maxCPUTime = 0;                    ///< Maximum CPU time in seconds
    
    qpOASES::BooleanType enableEqualities = qpOASES::BT_TRUE;    ///< Enable equality constraints
    
    enum mode { Default, Reliable, MPC } solverMode = Reliable; ///< Preset solver modes for different performance profiles

    /**
     * @brief Reset parameters to qpOASES-style defaults
     */
    void resetToDefaults() {
        verbosity = false;
        warmStart = true;
        maxWorkingSetRecalculations = 100;
        terminationTolerance = 1e-6;
        maxCPUTime = 0;
        enableEqualities = qpOASES::BT_FALSE;
        solverMode = Default;
    }
};

struct cbf{
	double h;
  Eigen::VectorXd dhdq;
};

typedef Eigen::Matrix<double,7,1> Vector7d;

inline Eigen::SparseMatrix<c_float> dense2sparse(Eigen::MatrixXd& A ){
  
  int n = A.rows();
  int m = A.cols();
  Eigen::SparseMatrix<c_float> As(n,m);

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      if(A(i,j) != 0)
        As.insert(i,j) = A(i,j);
    }
  }
  return As;
}


inline Eigen::MatrixXd pinv(const Eigen::MatrixXd &vector, double tolerance=0){
    
  // Perform Singular Value Decomposition (SVD) on the vector 
  Eigen::BDCSVD<Eigen::MatrixXd> svd(vector, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Extract singular values, and the U and V matrices
  Eigen::VectorXd singularValues = svd.singularValues();
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::MatrixXd V = svd.matrixV();
  
  Eigen::MatrixXd singularValuesInv = Eigen::MatrixXd::Zero(singularValues.size(), singularValues.size());

  // Compute the inverse of singular values with a tolerance
  tolerance = tolerance == 0 ? std::max(vector.rows(), vector.cols()) * std::numeric_limits<double>::epsilon() : tolerance;

  for (int i = 0; i < singularValues.size(); ++i) {
      if (singularValues(i) > tolerance) {
          singularValuesInv(i, i) = 1.0 / singularValues(i);
      }
  }

  // Compute the pseudo-inverse: V * Sigma+ * U^T
  return V * singularValuesInv * U.transpose();

}

inline bool initQP(qpOASES::SQProblem*& solver, double lu, double ld, double gamma, const vector<cbf> &cbfs){
  
  constexpr double tolerance = 1e-6;

  // Hessian matrix cost function (dense row-major required by qpOASES)
  Eigen::Matrix<double, 10, 10, Eigen::RowMajor> H_dense;
  H_dense.setZero();
  for(int i=0;i<7;i++){
    H_dense(i,i) = lu;
  }
  H_dense(7,7) = ld;
  H_dense(8,8) = ld;
  H_dense(9,9) = ld;


  // Constraint matrix A

  Eigen::SparseMatrix<double> A_s(19,10);
  // Eigen::MatrixXd A;
  // A.resize(18,12);

  // A <<  MatrixXd::Identity(7,7), MatrixXd::Zero(7, 5),
  //       -cbfs[0].dhdq(0),                0,                0,                 0,                  0,                 0,                 0,  0,  0,  0,
  //                     0,  -cbfs[1].dhdq(1),                0,                 0,                  0,                 0,                 0,  0,  0,  0,
  //                     0,                 0,  -cbfs[2].dhdq(2),                0,                  0,                 0,                 0,  0,  0,  0,
  //                     0,                 0,                 0,  -cbfs[3].dhdq(3),                 0,                 0,                 0,  0,  0,  0,
  //                     0,                 0,                 0,                 0,  -cbfs[4].dhdq(4),                 0,                 0,  0,  0,  0,
  //                     0,                 0,                 0,                 0,                 0,  -cbfs[5].dhdq(5),                 0,  0,  0,  0,
  //                     0,                 0,                 0,                 0,                 0,                 0,  -cbfs[6].dhdq(6),  0,  0,  0,
  //      -cbfs[7].dhdq(0),  -cbfs[7].dhdq(1),  -cbfs[7].dhdq(2),  -cbfs[7].dhdq(3),  -cbfs[7].dhdq(4),  -cbfs[7].dhdq(5),  -cbfs[7].dhdq(6),  0,  0,  0,
  //      -cbfs[8].dhdq(0),  -cbfs[8].dhdq(1),  -cbfs[8].dhdq(2),  -cbfs[8].dhdq(3),  -cbfs[8].dhdq(4),  -cbfs[8].dhdq(5),  -cbfs[8].dhdq(6),  0,  0,  0,
  //      -cbfs[9].dhdq(0),  -cbfs[9].dhdq(1),  -cbfs[9].dhdq(2),  -cbfs[9].dhdq(3),  -cbfs[9].dhdq(4),  -cbfs[9].dhdq(5),  -cbfs[9].dhdq(6), -1,  0,  0,
  //     -cbfs[10].dhdq(0), -cbfs[10].dhdq(1), -cbfs[10].dhdq(2), -cbfs[10].dhdq(3), -cbfs[10].dhdq(4), -cbfs[10].dhdq(5), -cbfs[10].dhdq(6),  0, -1,  0,
  //     -cbfs[11].dhdq(0), -cbfs[11].dhdq(1), -cbfs[11].dhdq(2), -cbfs[11].dhdq(3), -cbfs[11].dhdq(4), -cbfs[11].dhdq(5), -cbfs[11].dhdq(6),  0,  0, -1,

  // cout<<A;
  // A_s = dense2sparse(A);			

  // Input bounds 
  for(int i=0; i<7; i++){
    A_s.insert(i,i) = 1;
  }

  // Joint limits constraints
  for(int i=0; i<7; i++){
    A_s.insert(7+i,i) = -cbfs[i].dhdq(i);
  }

  // Hard task constraint
  for(int j=0; j<7; j++){
    A_s.insert(14,j) =  -cbfs[7].dhdq(j);
    A_s.insert(15,j) =  -cbfs[8].dhdq(j);
  }

  // Soft task constraints
  for(int i=0; i < 3; i++){
    for(int j=0; j<7; j++){
      A_s.insert(i+16,j) =  -cbfs[i+9].dhdq(j);
    }
    A_s.insert(i+16,i+7) = -1;
  }

  
  // Gradient cost function

  Eigen::Matrix<double, 10, 1> gradient;
  gradient.setZero();

  // Constraints rhs terms / Lower and Upper Bounds;

  Eigen::Matrix<double, 7, 1> u_lowerBound; 
  Eigen::Matrix<double, 7, 1> u_upperBound;
  u_lowerBound <<  -1.2, -1.2, -1.5, -1.2, -1.5, -1.5, -1.5;
  u_upperBound <<   1.2,  1.2,  1.5,  1.2,  1.5,  1.5,  1.5;


  Eigen::Matrix<double, 19, 1> lowerBound;
  lowerBound <<  u_lowerBound,
                -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, 
                -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, 
                -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, 
                -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY;

                
  Eigen::Matrix<double, 19, 1> upperBound;
  upperBound <<     u_upperBound, 
                 gamma*cbfs[0].h, 
                 gamma*cbfs[1].h,
                 gamma*cbfs[2].h,
                 gamma*cbfs[3].h,
                 gamma*cbfs[4].h,
                 gamma*cbfs[5].h,
                 gamma*cbfs[6].h, 
                 gamma*cbfs[7].h,
                 gamma*cbfs[8].h,
                 gamma*cbfs[9].h,
                 gamma*cbfs[10].h,
                 gamma*cbfs[11].h;

  // Solver settings
  solver = new qpOASES::SQProblem(10, 19);

  QPParams params;
  qpOASES::Options options;

  switch(params.solverMode)
  {
    case QPParams::Default:
        options.setToDefault();
        break;
    case QPParams::Reliable:
        options.setToReliable();
        break;
    case QPParams::MPC:
        options.setToMPC();
        break;
  }
  options.printLevel = params.verbosity ? qpOASES::PL_HIGH : qpOASES::PL_NONE;
  options.terminationTolerance = params.terminationTolerance;

  solver->setOptions(options);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  A_row;
  A_row = Eigen::MatrixXd(A_s);

  qpOASES::returnValue ret = solver->init(
      H_dense.data(),
      gradient.data(),
      A_row.data(),
      nullptr,              // variable lower bounds
      nullptr,              // variable upper bounds
      lowerBound.data(),
      upperBound.data(),
      params.maxWorkingSetRecalculations, 0
  );

  return true;
}

inline Eigen::Vector2d computeBounds(const vector<cbf> &cbfs, double gamma, double ld){

  Eigen::MatrixXd W, A, B, dhdq;

  Eigen::MatrixXd N1(7,7), N2(7,7);

  Eigen::VectorXd task_gradient_1 = cbfs[9].dhdq;
  Eigen::VectorXd task_gradient_2 = cbfs[10].dhdq;

  if(task_gradient_1.norm() > 0)
    N1 = Eigen::MatrixXd::Identity(7,7) - task_gradient_1*(task_gradient_1.transpose()*task_gradient_1).inverse()*task_gradient_1.transpose();
  else
    N1.setIdentity();
  
  if(task_gradient_2.norm() > 0)
    N2 = N1 - N1*task_gradient_2*(task_gradient_2.transpose()*N1*task_gradient_2).inverse()*task_gradient_2.transpose()*N1;
  else
    N2.setIdentity();
  
  dhdq.resize(3,7);
  A.resize(3,7);

  dhdq << cbfs[9].dhdq.transpose(),
          cbfs[10].dhdq.transpose(),
          cbfs[11].dhdq.transpose();

  A << cbfs[9].dhdq.transpose(),
       cbfs[10].dhdq.transpose()*N1,
       cbfs[11].dhdq.transpose()*N2;
  
  W =  dhdq*A.transpose()*(Eigen::MatrixXd::Identity(3,3) + ld*A*A.transpose()).inverse();

  double c1, c2, c2_lb, c3_lb;

  c1 = gamma;
  c2 = gamma;

  c2_lb = std::max(0.0, c1*pow(W(1,0), 2)/(4*(W(0,0)*W(1,1))));

  c3_lb = std::max(0.0, (W(0,0)*pow(W(2,1),2)*pow(c2,2) + c1*c2*pow(W(2,0),2)*W(1,1)-c1*c2*W(1,0)*W(2,0)*W(2,1))/(W(2,2)*(4*c2*W(0,0)*W(1,1)-c1*W(1,0)*W(1,0))));

  return Eigen::Vector2d(c2_lb, c3_lb);

}

inline bool updateQP(qpOASES::SQProblem* &solver, const vector<cbf> &cbfs, double gamma, double lu, double ld, int taskset){
 
  Eigen::Matrix<double, 10, 10, Eigen::RowMajor> H_dense;
  H_dense.setZero();
  for(int i=0;i<7;i++){
    H_dense(i,i) = lu;
  }

  Eigen::SparseMatrix<double> A_s(19,10);

  // Input bounds
  
  for(int i=0; i<7; i++){
    A_s.insert(i,i) = 1;
  }

  // Update joint limits constraints

  for(int i=0; i<7; i++){
    A_s.insert(7+i,i) = -cbfs[i].dhdq(i);
  }


  // Update task constraints

  // Hard task constraint
  for(int j=0; j<7; j++){
    A_s.insert(14,j) =  -cbfs[7].dhdq(j);
    A_s.insert(15,j) =  -cbfs[8].dhdq(j);
  }

  Eigen::MatrixXd N(7,7);
  N.setIdentity();
  bool stop_proj = false;

  // Soft task constraints
  for(int i=0; i < 3; i++){
    Eigen::VectorXd task_gradient = cbfs[i+9].dhdq;
    Eigen::VectorXd task_gradient_proj = task_gradient.transpose()*N;

    H_dense(i+7,i+7) = 1.0/(1.0/ld+task_gradient_proj.squaredNorm()-task_gradient_proj.transpose()*N*task_gradient_proj);


    if(task_gradient.norm() > 0){
      N = N - N*task_gradient*(task_gradient.transpose()*N*task_gradient).inverse()*task_gradient.transpose()*N;
        
      //N = N - pinv(task_gradient.transpose()*N)*task_gradient.transpose()*N;
    }

    //H_s.insert(i+7,i+7) = ld;
    
    for(int j=0; j<7; j++){
      A_s.insert(i+16,j) = -task_gradient_proj(j);
    }
    A_s.insert(i+16,i+7) = -1;
  }

  
  Eigen::Matrix<double, 7, 1> u_lowerBound; 
  Eigen::Matrix<double, 7, 1> u_upperBound;
  u_lowerBound <<-1.2, -1.2, -1.5, -1.2, -1.5, -1.5, -1.5;
  u_upperBound << 1.2,  1.2,  1.5,  1.2,  1.5,  1.5,  1.5;


  Eigen::Matrix<double, 19, 1> lowerBound;
  lowerBound <<  u_lowerBound,-qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY,
                              -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY, -qpOASES::INFTY;

  // Eigen::Vector2d bounds = computeBounds(cbfs, gamma, ld);
  // double gain2 = std::max(gamma, bounds(0));

  Eigen::Matrix<double, 19, 1> upperBound;
  upperBound <<  u_upperBound, gamma*cbfs[0].h, gamma*cbfs[1].h,gamma*cbfs[2].h,gamma*cbfs[3].h,gamma*cbfs[4].h,gamma*cbfs[5].h,gamma*cbfs[6].h,
                    gamma*cbfs[7].h, gamma*cbfs[8].h, 
                    gamma*cbfs[9].h, gamma*cbfs[10].h, gamma*cbfs[11].h;
  
  Eigen::VectorXd gradient;
  gradient.resize(10);
  gradient.setZero();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  A_row;
  A_row = Eigen::MatrixXd(A_s);
 
  QPParams params;
  qpOASES::int_t nWSR = params.maxWorkingSetRecalculations;

  auto lastSolveStatus = solver->hotstart(
        H_dense.data(),        // updated H (dense row-major)
        gradient.data(),       // gradient
        A_row.data(),          // updated constraint matrix
        nullptr,               // variable lower bounds (unchanged)
        nullptr,               // variable upper bounds (unchanged)
        lowerBound.data(),     // constraint lower bounds (size 19)
        upperBound.data(),     // constraint upper bounds (size 19)
        nWSR
  );

  return (lastSolveStatus == qpOASES::SUCCESSFUL_RETURN);
}



Eigen::Affine3d kdlFrameToEigenAffine(const KDL::Frame& f)
{
    Eigen::Affine3d T = Eigen::Affine3d::Identity();

    T.translation() << f.p.x(), f.p.y(), f.p.z();

    T.linear() <<
        f.M(0,0), f.M(0,1), f.M(0,2),
        f.M(1,0), f.M(1,1), f.M(1,2),
        f.M(2,0), f.M(2,1), f.M(2,2);

    return T;
}