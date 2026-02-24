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


inline Eigen::MatrixXd pinv(const Eigen::MatrixXd &vector, double tolerance =0  ){
    
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

inline bool initQP(OsqpEigen::Solver& solver, double lu, double ld, double gamma, const vector<cbf> &cbfs){
  
  constexpr double tolerance = 1e-6;

  // Hessian matrix cost function
  Eigen::SparseMatrix<c_float> H_s(10,10);
  for(int i=0;i<7;i++){
    H_s.insert(i,i) = lu;
  }
  H_s.insert(7,7) = ld;
  H_s.insert(8,8) = ld;
  H_s.insert(9,9) = ld;


  // Constraint matrix A

  Eigen::SparseMatrix<c_float> A_s(17,10);
  // Eigen::MatrixXd A;
  // A.resize(17,12);

  // A <<  MatrixXd::Identity(7,7), MatrixXd::Zero(7, 5),
  //       -cbfs[0].dhdq(0),               0,               0,                0,                 0,                0,                0,  0,  0,  0,
  //                     0, -cbfs[1].dhdq(1),               0,                0,                 0,                0,                0,  0,  0,  0,
  //                     0,                0, -cbfs[2].dhdq(2),               0,                 0,                0,                0,  0,  0,  0,
  //                     0,                0,                0, -cbfs[3].dhdq(3),                0,                0,                0,  0,  0,  0,
  //                     0,                0,                0,                0, -cbfs[4].dhdq(4),                0,                0,  0,  0,  0,
  //                     0,                0,                0,                0,                0, -cbfs[5].dhdq(5),                0,  0,  0,  0,
  //                     0,                0,                0,                0,                0,                0, -cbfs[6].dhdq(6),  0,  0,  0,
  //      -cbfs[7].dhdq(0), -cbfs[7].dhdq(1), -cbfs[7].dhdq(2), -cbfs[7].dhdq(3), -cbfs[7].dhdq(4), -cbfs[7].dhdq(5), -cbfs[7].dhdq(6), -1,  0,  0,
  //      -cbfs[8].dhdq(0), -cbfs[8].dhdq(1), -cbfs[8].dhdq(2), -cbfs[8].dhdq(3), -cbfs[8].dhdq(4), -cbfs[8].dhdq(5), -cbfs[8].dhdq(6),  0, -1,  0,
  //      -cbfs[9].dhdq(0), -cbfs[9].dhdq(1), -cbfs[9].dhdq(2), -cbfs[9].dhdq(3), -cbfs[9].dhdq(4), -cbfs[9].dhdq(5), -cbfs[9].dhdq(6),  0,  0, -1,

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



  // Gradient cost function

  Eigen::Matrix<c_float, 10, 1> gradient;
  gradient.setZero();

  // Constraints rhs terms / Lower and Upper Bounds;

  Eigen::Matrix<c_float, 7, 1> u_lowerBound; 
  Eigen::Matrix<c_float, 7, 1> u_upperBound;
  u_lowerBound <<  -1.2, -1.2, -1.5, -1.2, -1.5, -1.5, -1.5;
  u_upperBound <<   1.2,  1.2,  1.5,  1.2,  1.5,  1.5,  1.5;


  Eigen::Matrix<c_float, 17, 1> lowerBound;
  lowerBound <<  u_lowerBound,
                -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, 
                -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY;

                
  Eigen::Matrix<c_float, 17, 1> upperBound;
  upperBound <<  u_upperBound, gamma*cbfs[0].h, gamma*cbfs[1].h,gamma*cbfs[2].h,gamma*cbfs[3].h,gamma*cbfs[4].h,gamma*cbfs[5].h,gamma*cbfs[6].h, 
                                                                                            gamma*cbfs[7].h,gamma*cbfs[8].h,gamma*cbfs[9].h;

  // Solver settings

  solver.settings()->setVerbosity(false);
  solver.settings()->setAlpha(1.0);
  solver.settings()->setWarmStart(true);
  
  solver.data()->setNumberOfVariables(10);
  solver.data()->setNumberOfConstraints(17);
  

  // set the initial data of the QP solver
  if(!solver.data()->setHessianMatrix(H_s)) return false;
  if(!solver.data()->setGradient(gradient)) return false;
  if(!solver.data()->setLinearConstraintsMatrix(A_s))return false;

  if(!solver.data()->setLowerBound(lowerBound)) return false;
  if(!solver.data()->setUpperBound(upperBound)) return false;
  if(!solver.initSolver()) return false;

  return true;
}

inline bool updateQP(OsqpEigen::Solver& solver, const vector<cbf> &cbfs, double gamma, double lu, double ld, int taskset){

  Eigen::SparseMatrix<c_float> H_s(10,10);
  for(int i=0;i<7;i++){
    H_s.insert(i,i) = lu;
  }

  Eigen::SparseMatrix<c_float> A_s(17,10);
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

  Eigen::MatrixXd N(7,7);
  N.setIdentity();

  // Update task constraints

  for(int i=0; i < 3; i++){
    Eigen::VectorXd task_gradient = cbfs[i+7].dhdq;
    Eigen::VectorXd task_gradient_proj = task_gradient.transpose()*N;

    H_s.insert(i+7,i+7) = 1.0/(1.0/ld+task_gradient_proj.squaredNorm()-task_gradient_proj.transpose()*N*task_gradient_proj);
   
    if(task_gradient.norm()> 0){
      N = N - N*task_gradient*(task_gradient.transpose()*N*task_gradient).inverse()*task_gradient.transpose()*N;
   
      // std::cout <<"N: \n"<< N <<std::endl;
    }
   // H_s.insert(i+7,i+7) = ld;
    
    for(int j=0; j<7; j++){
      A_s.insert(i+14,j) =  -task_gradient_proj(j);
    }
    A_s.insert(i+14,i+7) = -1;
  }



  Eigen::Matrix<c_float, 7, 1> u_lowerBound; 
  Eigen::Matrix<c_float, 7, 1> u_upperBound;
  u_lowerBound <<-1.2, -1.2, -1.5, -1.2, -1.5, -1.5, -1.5;
  u_upperBound << 1.2,  1.2,  1.5,  1.2,  1.5,  1.5,  1.5;


  Eigen::Matrix<c_float, 17, 1> lowerBound;
  lowerBound <<  u_lowerBound,-OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, 
                              -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY;

                
  Eigen::Matrix<c_float, 17, 1> upperBound;
  upperBound <<  u_upperBound, gamma*cbfs[0].h, gamma*cbfs[1].h,gamma*cbfs[2].h,gamma*cbfs[3].h,gamma*cbfs[4].h,gamma*cbfs[5].h,gamma*cbfs[6].h,gamma*cbfs[7].h,gamma*cbfs[8].h,gamma*cbfs[9].h;
  

  // solver.data()->clearLinearConstraintsMatrix();
  // if(!solver.data()->setLinearConstraintsMatrix(A_s))return false;
  if(!solver.updateHessianMatrix(H_s)) return false;
  if(!solver.updateLinearConstraintsMatrix(A_s)) return false;
  if(!solver.updateBounds(lowerBound, upperBound)) return false;
  return true;
}
