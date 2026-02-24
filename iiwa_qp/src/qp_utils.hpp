#include <OsqpEigen/OsqpEigen.h>

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


inline bool initQP(OsqpEigen::Solver& solver, double lu, double ld, double lk, double gamma, double k, const vector<cbf> &cbfs){
  
  constexpr double tolerance = 1e-6;

  // Hessian matrix cost function
  Eigen::SparseMatrix<c_float> H_s(12,12);
  for(int i=0;i<7;i++){
    H_s.insert(i,i) = lu;
  }
  H_s.insert(7,7) = ld;
  H_s.insert(8,8) = ld;
  H_s.insert(9,9) = ld;
  H_s.insert(10,10) = lk;
  H_s.insert(11,11) = lk;


  // Constraint matrix A

  Eigen::SparseMatrix<c_float> A_s(19,12);
  // Eigen::MatrixXd A;
  // A.resize(19,12);

  // A <<  MatrixXd::Identity(7,7), MatrixXd::Zero(7, 5),
  //       -cbfs[0].dhdq(0),               0,               0,                0,                 0,                0,                0,  0,  0,  0, 0, 0,
  //                     0, -cbfs[1].dhdq(1),               0,                0,                 0,                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0, -cbfs[2].dhdq(2),               0,                 0,                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0, -cbfs[3].dhdq(3),                0,                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0,                0, -cbfs[4].dhdq(4),                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0,                0,                0, -cbfs[5].dhdq(5),                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0,                0,                0,                0, -cbfs[6].dhdq(6),  0,  0,  0, 0, 0,
  //      -cbfs[7].dhdq(0), -cbfs[7].dhdq(1), -cbfs[7].dhdq(2), -cbfs[7].dhdq(3), -cbfs[7].dhdq(4), -cbfs[7].dhdq(5), -cbfs[7].dhdq(6), -1,  0,  0, 0, 0,
  //      -cbfs[8].dhdq(0), -cbfs[8].dhdq(1), -cbfs[8].dhdq(2), -cbfs[8].dhdq(3), -cbfs[8].dhdq(4), -cbfs[8].dhdq(5), -cbfs[8].dhdq(6),  0, -1,  0, 0, 0,
  //      -cbfs[9].dhdq(0), -cbfs[9].dhdq(1), -cbfs[9].dhdq(2), -cbfs[9].dhdq(3), -cbfs[9].dhdq(4), -cbfs[9].dhdq(5), -cbfs[9].dhdq(6),  0,  0, -1, 0, 0,
  //                     0,                0,                0,                0,                0,                0,                0,  k, -1,  0, 1, 0,
  //                     0,                0,                0,                0,                0,                0,                0,  0,  k, -1, 0, k;
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

  //  task constraints

  for(int i=0; i < 3; i++){
    for(int j=0; j<7; j++){
      A_s.insert(i+14,j) =  -cbfs[i+7].dhdq(j);
    }
    A_s.insert(i+14,i+7) = -1;
  }

  //  priority matrix

  A_s.insert(17,7) = k;
  A_s.insert(17,8) = -1;
  A_s.insert(17,10) = 1;

  A_s.insert(18,8) = k;
  A_s.insert(18,9) = -1;
  A_s.insert(18,11) = k;


  // Gradient cost function

  Eigen::Matrix<c_float, 12, 1> gradient;
  gradient.setZero();

  // Constraints rhs terms / Lower and Upper Bounds;

  Eigen::Matrix<c_float, 7, 1> u_lowerBound; 
  Eigen::Matrix<c_float, 7, 1> u_upperBound;
  u_lowerBound <<  -1.2, -1.2, -1.5, -1.2, -1.5, -1.5, -1.5;
  u_upperBound <<   1.2,  1.2,  1.5,  1.2,  1.5,  1.5,  1.5;


  Eigen::Matrix<c_float, 19, 1> lowerBound;
  lowerBound <<  u_lowerBound,
                -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, 
                -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY;

                
  Eigen::Matrix<c_float, 19, 1> upperBound;
  upperBound <<  u_upperBound, gamma*cbfs[0].h, gamma*cbfs[1].h,gamma*cbfs[2].h,gamma*cbfs[3].h,gamma*cbfs[4].h,gamma*cbfs[5].h,gamma*cbfs[6].h, 
                                                                                            gamma*cbfs[7].h,gamma*cbfs[8].h,gamma*cbfs[9].h, 0,0;

  // Solver settings

  solver.settings()->setVerbosity(false);
  solver.settings()->setAlpha(1.0);
  solver.settings()->setWarmStart(true);
  
  solver.data()->setNumberOfVariables(12);
  solver.data()->setNumberOfConstraints(19);
  

  // set the initial data of the QP solver
  if(!solver.data()->setHessianMatrix(H_s)) return false;
  if(!solver.data()->setGradient(gradient)) return false;
  if(!solver.data()->setLinearConstraintsMatrix(A_s))return false;

  if(!solver.data()->setLowerBound(lowerBound)) return false;
  if(!solver.data()->setUpperBound(upperBound)) return false;
  if(!solver.initSolver()) return false;

  return true;
}

inline bool updateQP(OsqpEigen::Solver& solver, const vector<cbf> &cbfs, double gamma, double k, int taskset){

  
  Eigen::SparseMatrix<c_float> A_s(19,12);
  // Eigen::MatrixXd A;
  // A.resize(19,12);

  // A <<  MatrixXd::Identity(7,7), MatrixXd::Zero(7, 5),
  //       -cbfs[0].dhdq(0),               0,               0,                0,                 0,                0,                0,  0,  0,  0, 0, 0,
  //                     0, -cbfs[1].dhdq(1),               0,                0,                 0,                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0, -cbfs[2].dhdq(2),               0,                 0,                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0, -cbfs[3].dhdq(3),                0,                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0,                0, -cbfs[4].dhdq(4),                0,                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0,                0,                0, -cbfs[5].dhdq(5),                0,  0,  0,  0, 0, 0,
  //                     0,                0,                0,                0,                0,                0, -cbfs[6].dhdq(6),  0,  0,  0, 0, 0,
  //      -cbfs[7].dhdq(0), -cbfs[7].dhdq(1), -cbfs[7].dhdq(2), -cbfs[7].dhdq(3), -cbfs[7].dhdq(4), -cbfs[7].dhdq(5), -cbfs[7].dhdq(6), -1,  0,  0, 0, 0,
  //      -cbfs[8].dhdq(0), -cbfs[8].dhdq(1), -cbfs[8].dhdq(2), -cbfs[8].dhdq(3), -cbfs[8].dhdq(4), -cbfs[8].dhdq(5), -cbfs[8].dhdq(6),  0, -1,  0, 0, 0,
  //      -cbfs[9].dhdq(0), -cbfs[9].dhdq(1), -cbfs[9].dhdq(2), -cbfs[9].dhdq(3), -cbfs[9].dhdq(4), -cbfs[9].dhdq(5), -cbfs[9].dhdq(6),  0,  0, -1, 0, 0,
  //                     0,                0,                0,                0,                0,                0,                0,  k, -1,  0, 1, 0,
  //                     0,                0,                0,                0,                0,                0,                0,  0,  k, -1, 0, k;
  
  // A_s = dense2sparse(A);	


  // Input bounds
  
  for(int i=0; i<7; i++){
    A_s.insert(i,i) = 1;
  }

  // Update joint limits constraints

  for(int i=0; i<7; i++){
    A_s.insert(7+i,i) = -cbfs[i].dhdq(i);
  }

  // Update task constraints

  for(int i=0; i < 3; i++){
    for(int j=0; j<7; j++){
      A_s.insert(i+14,j) =  -cbfs[i+7].dhdq(j);
    }
    A_s.insert(i+14,i+7) = -1;
  }

  // Update priority matrix

  if(taskset == 1){
    A_s.insert(17,7) = k;
    A_s.insert(17,8) = -1;
    A_s.insert(17,10) = 0*1;

    A_s.insert(18,8) = 0;
    A_s.insert(18,9) = -1;
    A_s.insert(18,11) = 0*k; //0*k;

    // A_s.insert(17,7) = 1;
    // A_s.insert(17,8) = -1/k;
    // A_s.insert(17,10) = 1/k;

    // A_s.insert(18,8) = 0;
    // A_s.insert(18,9) = -0;
    // A_s.insert(18,11) = 1;

  }else{
    A_s.insert(17,7) = k;
    A_s.insert(17,8) = -1;
    A_s.insert(17,10) = 0*1;

    A_s.insert(18,8) = k;
    A_s.insert(18,9) = -1;
    A_s.insert(18,11) = 0*k; //0*k;

    // A_s.insert(17,7) = 1;
    // A_s.insert(17,8) = -1/k;
    // A_s.insert(17,10) = 1/k;

    // A_s.insert(18,8) = 1;
    // A_s.insert(18,9) = -1/k;
    // A_s.insert(18,11) = 1;
  }

  Eigen::Matrix<c_float, 7, 1> u_lowerBound; 
  Eigen::Matrix<c_float, 7, 1> u_upperBound;
  u_lowerBound <<-1.2, -1.2, -1.5, -1.2, -1.5, -1.5, -1.5;
  u_upperBound << 1.2,  1.2,  1.5,  1.2,  1.5,  1.5,  1.5;


  Eigen::Matrix<c_float, 19, 1> lowerBound;
  lowerBound <<  u_lowerBound,-OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY;

                
  Eigen::Matrix<c_float, 19, 1> upperBound;
  upperBound <<  u_upperBound, gamma*cbfs[0].h, gamma*cbfs[1].h,gamma*cbfs[2].h,gamma*cbfs[3].h,gamma*cbfs[4].h,gamma*cbfs[5].h,gamma*cbfs[6].h,gamma*cbfs[7].h,gamma*cbfs[8].h,gamma*cbfs[9].h, 0,0;
  
  // solver.data()->clearLinearConstraintsMatrix();
  // if(!solver.data()->setLinearConstraintsMatrix(A_s))return false;
  if(!solver.updateLinearConstraintsMatrix(A_s)) return false;
  if(!solver.updateBounds(lowerBound, upperBound)) return false;
  return true;
}
