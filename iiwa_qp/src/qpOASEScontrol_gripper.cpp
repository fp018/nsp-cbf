#include "ros/ros.h"
#include <ros/package.h>
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/Pose.h"
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Int8.h>
#include <gazebo_msgs/ModelStates.h>

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>
#include <kdl/chaindynparam.hpp>

#include <eigen3/Eigen/Dense>
#include <tf/tf.h>
#include <tf_conversions/tf_eigen.h>
#include <eigen_conversions/eigen_kdl.h>
#include <joint_limits_interface/joint_limits.h>

#include <iiwa_qp/lbr_iiwa_7_robot.h>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <fstream> 
#include <iomanip>
#include "qpOASES_utils_gripper.hpp"
#include <boeing_gazebo_model_attachment_plugin/Attach.h>

#define STACK 1
#define DQPINTRP 1
//#define NOINTERP 1
#define SOL_TIME_EVAL 0

#define dt 0.001


using namespace std;
class KUKA_CONTROL {

	public:

		KUKA_CONTROL();
		~KUKA_CONTROL();
		void run();
		bool init_model();
		bool set_qp();
		void joint_states_cb(const sensor_msgs::JointState & );
		void artag_cb(const geometry_msgs::Pose &);
		void taskset_cb(const std_msgs::Int8 &);
		void modelStatesCallback(const gazebo_msgs::ModelStates::ConstPtr& msg);

		void update_dirkin(Eigen::VectorXd &q_in); 
		bool computeJacobian(KDL::ChainJntToJacSolver& solver,const KDL::JntArray& q, KDL::Jacobian jac, Eigen::MatrixXd& out, const std::string& name);
		
		void ctrl_loop();
		void task_set_input();
		void grasp_sm(int &taskSet);
		bool request_attach_plate();

		cbf cbf_goto_point(const Eigen::Vector3d &dpe);
		cbf cbf_goto_z_3link(const float &d3z);
		cbf cbf_vision();
		cbf cbf_keep_ee_z_horizontal();
		cbf cbf_keep_ee_x_horizontal();
		cbf cbf_approach_grasp_plate(Affine3d plateTf, float r);
		cbf cbf_avoid_plate(const Affine3d linkTf, const MatrixXd J, const Affine3d plateTf, const Vector3d lengths);
		void cbf_joint_limits(vector<cbf> &cbfjl);
		void taskstack(int & taskSet, vector<cbf> & cbfs, std_msgs::Float64MultiArray & h);
		void solveQP(qpOASES::SQProblem* solver, const vector<cbf> &cbfs, int taskset, Vector7d& u);

	private:

		ros::NodeHandle _nh;
		ros::Subscriber _js_sub, _tag_sub, _task_sub, _plate_sub;
		ros::Publisher _js_pub;
		ros::Publisher _u_pub, _h_pub, _jv_pub;
		ros::Publisher _gripper_lf_pub, _gripper_rf_pub;
		vector<ros::Publisher> _j_pub;

		bool _first_js;
		bool _first_fk;
		bool _sync;
		bool _got_tag;

		geometry_msgs::Pose _ptag;

		KDL::JntArray *_q_in;
		KDL::JntArray *_q_3in;
		KDL::JntArray *_q_gripper;
		KDL::JntArray *_q_gripper_lf;
		KDL::JntArray *_q_gripper_rf;


		KDL::Tree iiwa_tree;
		KDL::ChainFkSolverPos_recursive *_fksolver; //Forward position solver
		KDL::ChainFkSolverPos_recursive *_fk3solver; //Forward 3link position solver
		KDL::ChainFkSolverPos_recursive *_fk_gripper_solver;
		KDL::ChainFkSolverPos_recursive *_fk_gripper_lf_solver;
		KDL::ChainFkSolverPos_recursive *_fk_gripper_rf_solver;

		KDL::ChainJntToJacSolver *_J_solver;
		KDL::ChainJntToJacSolver *_J3_solver;
		KDL::ChainJntToJacSolver *_J_gripper_solver;
		KDL::ChainJntToJacSolver *_J_gripper_lf_solver;
		KDL::ChainJntToJacSolver *_J_gripper_rf_solver;

		KDL::ChainDynParam *_dyn_param;
		KDL::Chain _k_chain;
		KDL::Chain _k3_chain;
		KDL::Chain _k_gripper_chain;
		KDL::Chain _k_gripper_lf_chain;
		KDL::Chain _k_gripper_rf_chain;

		KDL::Frame _Te;
		KDL::Frame _T3;
		KDL::Frame _T_gripper;
		KDL::Frame _T_gripper_lf;
		KDL::Frame _T_gripper_rf;

		KDL::FrameVel _dirkin_out;
		KDL::FrameVel _dirkin_out3;

		Eigen::VectorXd _vel; 
		Eigen::MatrixXd _J;
		Eigen::MatrixXd _J3;
		Eigen::MatrixXd _J_gripper;
		Eigen::MatrixXd _J_gripper_lf;
		Eigen::MatrixXd _J_gripper_rf;

		qpOASES::SQProblem* _solver; 
		qpOASES::SQProblem* _solverprv;     

		int _taskSet;
		std::vector<float> _gamma;
		float _lu;
		float _ldelta;
	

		int _Nj; // number of joints
		vector<float> _q_min;
		vector<float> _q_max;
		Eigen::VectorXd _q;
		Eigen::VectorXd _qfb;
		Eigen::Vector3d _xe;
		
		ros::Time _latest_tag;

		double _kp;
		Eigen::Vector3d _p1;
		Eigen::Vector3d _p2;
		Eigen::Vector3d _plate_position;
		Eigen::Vector3d _plate_rpy;
		float _plate_radius;

		int _state;
		bool _start;

		std::vector<double> _tsol;
		std::ofstream _outFile;

};


KUKA_CONTROL::KUKA_CONTROL() {
	
	if (!init_model()) exit(1);

	ROS_INFO("Robot tree correctly loaded from parameter server!");

	if(set_qp()) ROS_WARN("Using default QP values");

	_j_pub.resize(7);
	_js_sub = _nh.subscribe("/kuka_iiwa/joint_states", 1, &KUKA_CONTROL::joint_states_cb, this);
	_jv_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/jointsVelCommand", 1);
	_h_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/h", 1);
	_tag_sub = _nh.subscribe("/kuka/camera/tag_pose", 0, &KUKA_CONTROL::artag_cb, this);  
	_plate_sub = _nh.subscribe("/gazebo/model_states", 10, &KUKA_CONTROL::modelStatesCallback, this);

	for(int k=0; k < _Nj; k++){
		//_j_pub[k] = _nh.advertise<std_msgs::Float64>("/kuka_iiwa/joint" + to_string(k+1) + "_velocity_controller/command", 0); //rostopic pub /kuka/task_set std_msgs/Int8 '2'
		_j_pub[k] = _nh.advertise<std_msgs::Float64>("/kuka_iiwa/joint" + to_string(k+1) + "_position_controller/command", 1);
	}

	_gripper_lf_pub = _nh.advertise<std_msgs::Float64>("/kuka_iiwa/wsg_50_gl/command", 10);
	_gripper_rf_pub = _nh.advertise<std_msgs::Float64>("/kuka_iiwa/wsg_50_gr/command", 10);

	_q_in = new KDL::JntArray(_Nj);
  	_q_3in = new KDL::JntArray(3);
	_q_gripper = new KDL::JntArray(_k_gripper_chain.getNrOfJoints());
	_q_gripper_lf = new KDL::JntArray(_k_gripper_lf_chain.getNrOfJoints());
	_q_gripper_rf = new KDL::JntArray(_k_gripper_rf_chain.getNrOfJoints());


	_fksolver = new KDL::ChainFkSolverPos_recursive( _k_chain);
	_J_solver = new KDL::ChainJntToJacSolver( _k_chain);

	_fk3solver = new KDL::ChainFkSolverPos_recursive( _k3_chain);
	_J3_solver = new KDL::ChainJntToJacSolver(_k3_chain);

	_fk_gripper_solver = new KDL::ChainFkSolverPos_recursive( _k_gripper_chain);
	_J_gripper_solver = new KDL::ChainJntToJacSolver(_k_gripper_chain);


	_fk_gripper_lf_solver = new KDL::ChainFkSolverPos_recursive( _k_gripper_lf_chain);
	_J_gripper_lf_solver = new KDL::ChainJntToJacSolver(_k_gripper_lf_chain);

	_fk_gripper_rf_solver = new KDL::ChainFkSolverPos_recursive( _k_gripper_rf_chain);
	_J_gripper_rf_solver = new KDL::ChainJntToJacSolver(_k_gripper_rf_chain);

	_dyn_param = new KDL::ChainDynParam(_k_chain,KDL::Vector(0,0,-9.81));

	_first_js = false;
	_sync = false;
	_got_tag = false;
	_taskSet = 0;
	_start = false;

	_outFile.open("qp_time.csv");
}

void  KUKA_CONTROL::modelStatesCallback(const gazebo_msgs::ModelStates::ConstPtr& msg)
{
	const std::string target_model = "grasp_plate";

	for (std::size_t i = 0; i < msg->name.size(); ++i)
	{
		if (msg->name[i] == target_model)
		{
			const geometry_msgs::Pose& pose = msg->pose[i];
			return;
		}
	}
}

bool KUKA_CONTROL::init_model(){

	std::string robot_desc_string;
	_nh.param("robot_description", robot_desc_string, std::string());
	if (!kdl_parser::treeFromString(robot_desc_string, iiwa_tree)){
		ROS_ERROR("Failed to construct kdl tree");
		return false;
	}

	std::string base_link = "world";
	std::string tip_link = "iiwa_link_7";
	std::string link3 = "iiwa_link_3";
	std::string gripper_link = "iiwa_gripper_midpoint";
	std::string gripper_rf_link = "wsg_50_tip_finger_right";
	std::string gripper_lf_link = "wsg_50_tip_finger_left";

	if ( !iiwa_tree.getChain(base_link, tip_link, _k_chain) ) return false;
	if ( !iiwa_tree.getChain(base_link, link3, _k3_chain) ) return false; 
	if ( !iiwa_tree.getChain(base_link, gripper_link, _k_gripper_chain) ) return false; 
	if ( !iiwa_tree.getChain(base_link, gripper_lf_link, _k_gripper_lf_chain)) return false; 
	if ( !iiwa_tree.getChain(base_link, gripper_rf_link, _k_gripper_rf_chain)) return false; 


	_Nj = _k_chain.getNrOfJoints();

	_J.resize(6,_k_chain.getNrOfJoints());
	_J3.resize(6,_k3_chain.getNrOfJoints());
	_J_gripper.resize(6,_k_gripper_chain.getNrOfJoints());
	_J_gripper_lf.resize(6,_k_gripper_lf_chain.getNrOfJoints());
	_J_gripper_rf.resize(6,_k_gripper_rf_chain.getNrOfJoints());

	_J = Eigen::MatrixXd::Zero(6,_k_chain.getNrOfJoints());
	_J3 = Eigen::MatrixXd::Zero(6,_k3_chain.getNrOfJoints());
	_J_gripper = Eigen::MatrixXd::Zero(6,_k_gripper_chain.getNrOfJoints());
	_J_gripper_lf = Eigen::MatrixXd::Zero(6,_k_gripper_lf_chain.getNrOfJoints());
	_J_gripper_rf = Eigen::MatrixXd::Zero(6,_k_gripper_rf_chain.getNrOfJoints());

	_vel.resize(6);

	_q.resize(_Nj);
	_qfb.resize(_Nj);

	_nh.getParam("q_max", _q_max);
	_nh.getParam("q_min", _q_min);
	_nh.getParam("kp", _kp);

	vector<float> p1, p2;
	_nh.getParam("p1", p1);
	_nh.getParam("p2", p2);
	vector<float> plate_position;
	vector<float> plate_rpy;
	_nh.getParam("plate_position", plate_position);
	_nh.getParam("plate_rpy", plate_rpy);
	_nh.getParam("plate_radius", _plate_radius);
	std::cout << "Plate position: " << plate_position[0] << ", " << plate_position[1] << ", " << plate_position[2] << std::endl;
	std::cout << "Plate rpy: " << plate_rpy[0] << ", " << plate_rpy[1] << ", " << plate_rpy[2] << std::endl;
	std::cout << "Plate radius: " << _plate_radius << std::endl;

	for(int i=0; i<3; i++){
		_p1(i) = p1[i];
		_p2(i) = p2[i];
		_plate_position(i) = plate_position[i];
		_plate_rpy(i) = plate_rpy[i];
	}
	
	_state = 0;

	return true;

}

bool KUKA_CONTROL::set_qp(){

	bool dflt = false;
	float gamma_hard;
	float gamma;

	if(!_nh.getParam("gamma", gamma) || !_nh.getParam("gamma_hard", gamma_hard) ||  !_nh.getParam("lu", _lu) ||  !_nh.getParam("ldelta", _ldelta)){
		ROS_WARN("Using default QP values");	

		_lu = 1;
		_gamma.push_back(1.0);
		_gamma.push_back(1.0);
		_ldelta = 1000;

		dflt = true;
	}else{
		_gamma.push_back(gamma);
		_gamma.push_back(gamma_hard);

		ROS_INFO("Parameters loaded");
		std::cout << "gamma: " << gamma << std::endl;
		std::cout << "gamma_hard: " << gamma_hard << std::endl;
		std::cout << "lu: " << _lu << std::endl;
		std::cout << "ldelta: " << _ldelta << std::endl;
	}
	return dflt;
}

void KUKA_CONTROL::joint_states_cb(const sensor_msgs::JointState& js){

	for(int i=0; i<7; i++ ) {
		_qfb(i) = js.position[i];
		_q_gripper->data[i] = js.position[i];
	}
	
	// Gripper joints
	_q_gripper_lf->data[7] = js.position[7];
	_q_gripper_rf->data[7] = js.position[8];

	_first_js = true;
	_sync = true;
}

void KUKA_CONTROL::artag_cb(const geometry_msgs::Pose& tag_pose){

	_ptag = tag_pose;
	_latest_tag = ros::Time::now();
	_got_tag = true;

}

void KUKA_CONTROL::taskset_cb(const std_msgs::Int8 & i){
	
	_taskSet = i.data;
	
}

void KUKA_CONTROL::update_dirkin(Eigen::VectorXd & q_in){

	_q_in->data = q_in;
	_q_3in->data = q_in.head(3);
	_q_gripper->data.head(7) = q_in.head(7);

	_q_gripper_lf->data.head(7) =  q_in.head(7);
	_q_gripper_rf->data.head(7) =  q_in.head(7);

	// std::cout << "q_in: " << q_in.transpose() << std::endl;

	_fksolver->JntToCart(*_q_in, _Te);
	
	_xe(0) = _Te.p.x();
	_xe(1) = _Te.p.y();
	_xe(2) = _Te.p.z();

	_fk3solver->JntToCart(*_q_3in, _T3);
	_fk_gripper_solver->JntToCart(*_q_in, _T_gripper);
	_fk_gripper_lf_solver->JntToCart(*_q_gripper_lf, _T_gripper_lf);
	_fk_gripper_rf_solver->JntToCart(*_q_gripper_rf, _T_gripper_rf);

	KDL::Jacobian Jac(_k_chain.getNrOfJoints());
	KDL::Jacobian Jac3(_k3_chain.getNrOfJoints());
	KDL::Jacobian Jac_gripper(_k_gripper_chain.getNrOfJoints());
	KDL::Jacobian Jac_gripper_lf(_k_gripper_lf_chain.getNrOfJoints());
	KDL::Jacobian Jac_gripper_rf(_k_gripper_rf_chain.getNrOfJoints());
	
	computeJacobian(*_J_solver, *_q_in, Jac, _J, "Jacobian");
	computeJacobian(*_J3_solver, *_q_3in, Jac3, _J3, "Jacobian3");
	computeJacobian(*_J_gripper_solver, *_q_in, Jac_gripper, _J_gripper, "gripper Jacobian");
	computeJacobian(*_J_gripper_lf_solver, *_q_gripper_lf, Jac_gripper_lf, _J_gripper_lf, "left finger Jacobian");
	computeJacobian(*_J_gripper_rf_solver, *_q_gripper_rf, Jac_gripper_rf, _J_gripper_rf, "right finger Jacobian");

}

bool KUKA_CONTROL::computeJacobian(KDL::ChainJntToJacSolver& solver,const KDL::JntArray& q, KDL::Jacobian jac, Eigen::MatrixXd& out, const std::string& name)
{	
    if (solver.JntToJac(q, jac) != KDL::ChainJntToJacSolver::E_NOERROR) {
        std::cout << "failing in " << name << " computation!" << std::endl;
		return false;
    }
    out = jac.data;
    return true;
};

cbf KUKA_CONTROL::cbf_goto_point(const Eigen::Vector3d &dpe){
	
 	Eigen::Vector3d perr(_xe(0)- dpe(0), _xe(1)-dpe(1), _xe(2)-dpe(2));
	cbf cbfpos;
	cbfpos.dhdq.resize(_Nj);
	cbfpos.h = -perr.squaredNorm();
	cbfpos.dhdq = -0.5*perr.transpose()*_J.block(0,0,3,_Nj);
	return cbfpos;
}

cbf KUKA_CONTROL::cbf_goto_z_3link(const float &d3z){

	cbf cbf3z;
	cbf3z.dhdq.resize(7);
	cbf3z.h = -50*(_T3.p.z()-d3z)*(_T3.p.z()-d3z);  

	cbf3z.dhdq << -50*2*(_T3.p.z()-d3z)*Eigen::Vector3d(_J3.row(2)),0,0,0,0;

	return cbf3z;
}

cbf KUKA_CONTROL::cbf_keep_ee_z_horizontal(){

	cbf cbfEEhor;
	cbfEEhor.dhdq.resize(_Nj);
	Eigen::Vector3d z(_Te.M.UnitZ()(0), _Te.M.UnitZ()(1), _Te.M.UnitZ()(2));

	cbfEEhor.h = -(z.squaredNorm() + 1 -2* _Te.M.UnitZ()(2));   // -(z - [0,0,1]')'(z - [0,0,1]')
	
	Eigen::MatrixXd omegaJ = _J.block(3,0,3,_Nj);


	Eigen::Vector3d ez;
	ez = z - Eigen::Vector3d(0,0,1);

	cbfEEhor.dhdq =  -2* Eigen::Vector3d(-ez(1)*z(2) + ez(2)*z(1), 
								          ez(0)*z(2) - ez(2)*z(0),
										 -ez(0)*z(1) + ez(1)*z(0)).transpose()*omegaJ;

	// ---------------------------------------------------------------------------------

	return cbfEEhor;
}


cbf KUKA_CONTROL::cbf_keep_ee_x_horizontal(){

	cbf cbfEEhor;
	cbfEEhor.dhdq.resize(_Nj);
	Eigen::Vector3d x(_Te.M.UnitX()(0), _Te.M.UnitX()(1), _Te.M.UnitX()(2));

	cbfEEhor.h = -(x.squaredNorm() + 1 -2* _Te.M.UnitX()(2));   // -(z - [0,0,1]')'(z - [0,0,1]')
	
	Eigen::MatrixXd omegaJ = _J.block(3,0,3,_Nj);


	Eigen::Vector3d ex;
	ex = x - Eigen::Vector3d(0,0,1);

	cbfEEhor.dhdq =  -2* Eigen::Vector3d(-ex(1)*x(2) + ex(2)*x(1), 
								          ex(0)*x(2) - ex(2)*x(0),
										 -ex(0)*x(1) + ex(1)*x(0)).transpose()*omegaJ;

	// ---------------------------------------------------------------------------------

	return cbfEEhor;
}

cbf KUKA_CONTROL::cbf_approach_grasp_plate(Affine3d plateTf, float r){

	cbf cbfEEplate;
	cbfEEplate.dhdq.resize(_Nj);

	Eigen::Vector3d y(_T_gripper.M.UnitY()(0), _T_gripper.M.UnitY()(1), _T_gripper.M.UnitY()(2)); // Gripper sliding direction (y-axis of the end-effector frame)
	Eigen::Vector3d z(_T_gripper.M.UnitZ()(0), _T_gripper.M.UnitZ()(1), _T_gripper.M.UnitZ()(2));
	Eigen::Vector3d p(_T_gripper.p.x(), _T_gripper.p.y(), _T_gripper.p.z()); 
    
	Eigen::Vector3d p_plate = plateTf.inverse() * p;
	Eigen::Vector3d plate_center = plateTf.translation();

	Eigen::Vector3d plate_n = plateTf.rotation().col(2); // plate normal vector (z-axis of the plate frame)

	const Eigen::Matrix3d W = (Eigen::Vector3d(1.0, 1.0, 1.0)).asDiagonal(); // wx, wy, wz

	
	Eigen::MatrixXd omegaJ = _J_gripper.block(3,0,3,_Nj);
	Eigen::MatrixXd Jp = plateTf.rotation()* _J_gripper.block(0,0,3,_Nj);
	Eigen::Vector2d p_plate_xy(p_plate(0), p_plate(1));

	Eigen::Vector3d o_ey, o_ez;
	o_ey = y - plate_n;


	Eigen::Matrix<double,3,3> P = Eigen::Matrix<double,3,3>::Identity(3,3) - plate_n*plate_n.transpose()/ plate_n.squaredNorm();

	Eigen::Vector3d z_d, delta_p;
	delta_p = plateTf.translation() - p;

	z_d = P*delta_p/delta_p.norm();

	o_ez = z - z_d;

	// ------------------------------------------------------------------------------
	cbfEEplate.h = -(p_plate_xy.norm() - r)*(p_plate_xy.norm() - r)
				   -(p_plate(2)*p_plate(2))
				   -o_ey.squaredNorm()
				   -o_ez.squaredNorm(); //(z-z_d).dot((z-z_d));

	cbfEEplate.dhdq = -(p_plate_xy.norm() - r)*(1/p_plate_xy.norm())*(p - plateTf.translation()).head<2>().transpose()*_J_gripper.block(0,0,2,_Nj) 
					  -2*(p(2) -  plate_center(2))*_J_gripper.block(2,0,1,_Nj) 
					  -2*Eigen::Vector3d(-o_ey(1)*y(2) + o_ey(2)*y(1), // -y z + z y
				    				      o_ey(0)*y(2) - o_ey(2)*y(0), // x z - z x
										 -o_ey(0)*y(1) + o_ey(1)*y(0)).transpose()*omegaJ // -x y + y x
					  -2* Eigen::Vector3d(-o_ez(1)*z(2) + o_ez(2)*z(1),                    // gradient of the term (z-z_d).dot((z-z_d)) 
				     				       o_ez(0)*z(2) - o_ez(2)*z(0),
										  -o_ez(0)*z(1) + o_ez(1)*z(0)).transpose()*omegaJ
				      -2*(z - z_d).transpose()*P*(1/delta_p.norm()) * _J_gripper.block(0,0,3,_Nj)
					  +2*(z - z_d).transpose()*P*(delta_p*delta_p.transpose()/pow(delta_p.norm(), 3))* _J_gripper.block(0,0,3,_Nj); 

	// ------------------------------------------------------------------------------
	// cbfEEplate.h = -(p_plate_xy.squaredNorm() - r*r)*(p_plate_xy.squaredNorm() - r*r)
	// 			   -(p_plate(2)*p_plate(2))
	// 			   -o_ey.squaredNorm()
	// 			   -o_ez.squaredNorm(); //(z-z_d).dot((z-z_d));

	// cbfEEplate.dhdq = -4*(p_plate_xy.squaredNorm() - r*r)*(p - plateTf.translation()).head<2>().transpose()*_J_gripper.block(0,0,2,_Nj) 
	// 				  -2*(p(2) -  plate_center(2))*_J_gripper.block(2,0,1,_Nj) 
	// 				  -2*Eigen::Vector3d(-o_ey(1)*y(2) + o_ey(2)*y(1), 
	// 			    				      o_ey(0)*y(2) - o_ey(2)*y(0),
	// 									 -o_ey(0)*y(1) + o_ey(1)*y(0)).transpose()*omegaJ
	// 				  -2* Eigen::Vector3d(-o_ez(1)*z(2) + o_ez(2)*z(1), 
	// 			     				       o_ez(0)*z(2) - o_ez(2)*z(0),
	// 									  -o_ez(0)*z(1) + o_ez(1)*z(0)).transpose()*omegaJ
	// 			      -2*(z -z_d).transpose()*P*(1/(plateTf.translation() - p).norm())*_J_gripper.block(0,0,3,_Nj); // gradient of the term (z-z_d).dot((z-z_d)) w.r.t. q								   

	return cbfEEplate;
}

cbf KUKA_CONTROL::cbf_avoid_plate(const Affine3d linkTf, const MatrixXd J, const Affine3d plateTf, const Vector3d lengths)
{
    cbf cbf_avoid;
    cbf_avoid.dhdq.resize(_Nj);

    Eigen::Vector3d p = linkTf.translation();
    Eigen::Vector3d c = plateTf.translation();
    Eigen::Matrix3d R = plateTf.rotation();

    // point in plate frame, centered at plate origin
    Eigen::Vector3d plate_p = R.transpose() * (p - c);

    Eigen::Vector3d inv_lengths = lengths.cwiseInverse();
    Eigen::Matrix3d L2 = inv_lengths.cwiseProduct(inv_lengths).asDiagonal();

    // h >= 0 means outside ellipsoid
    cbf_avoid.h = plate_p.dot(L2 * plate_p) - 1.0;

    Eigen::MatrixXd Jp = J.block(0,0,3,_Nj);
	Eigen::MatrixXd Jo = J.block(3,0,3,_Nj);


	cbf_avoid.dhdq = 2.0 * plate_p.transpose() * L2 *(R.transpose() * Jp );
	
	// cbf_avoid.h = (p-c).transpose()*(p-c) - (0.15*0.15); // alternative: sphere of radius 15cm around plate center
	// cbf_avoid.dhdq = 2.0 * (p-c).transpose() * Jp;

    return cbf_avoid;
}


cbf KUKA_CONTROL::cbf_vision(){

	// camera
	KDL::Jacobian J_ee;
	KDL::Frame T_ee;
	// camera in end-effector transform
	KDL::Frame ee_T_cam;

	ee_T_cam.M = KDL::Rotation::Identity();
	ee_T_cam.M = KDL::Rotation::RotZ(-1.57);
	
	ee_T_cam.p = KDL::Vector(0,0,0);
	
	// object
	KDL::Frame cam_T_object, T_object;

	// current forward kinematics
	T_ee = _Te;

	KDL::Frame T_cam = T_ee*ee_T_cam;
	T_cam.p = T_ee.p;
	// current jacobians
	J_ee.data = _J;
	KDL::Jacobian J_cam_sp, J_cam;
	KDL::Vector p_ee_cam = T_ee.M*ee_T_cam.p;  
	J_cam.resize(7);
	J_cam_sp.resize(7);

	//changeRefPoint(J_ee, p_ee_cam, J_cam);            // camera jacobian
	J_cam = J_ee;                                       // as matlab
	changeBase(J_cam, T_cam.M.Inverse(), J_cam_sp);     // end-effector spatial jacobian // T_M := rotation part 
																											//Equivalent to J_cam_sp = blkdiag(T_cam(1:3,1:3)', T_cam(1:3,1:3)')*J_cam;
	
	if(ros::Time::now().toSec() - _latest_tag.toSec() > 10 ) _got_tag = false;

	if(_got_tag){
		cam_T_object.p = KDL::Vector(_ptag.position.x,_ptag.position.y,_ptag.position.z);
		Eigen::Quaterniond cam_quat_object(_ptag.orientation.w,_ptag.orientation.x,_ptag.orientation.y,_ptag.orientation.z);
		//Eigen::Quaterniond cam_quat_object(0,1,0,0);
		tf::quaternionEigenToKDL(cam_quat_object.normalized(), cam_T_object.M);

	}else{
		cam_T_object.p = T_cam.M.Inverse()*(KDL::Vector(2, 0.5, 0.0) - T_cam.p);   //cam_p_object = T_cam.M.transpose()*(p_object - p_cam);
	}
	

	Eigen::Matrix<double,3,1> s,s_d;
	tf::vectorKDLToEigen(cam_T_object.p/(cam_T_object.p.Norm()),s);

	Eigen::Matrix<double,3,3> skew_s = utilities::skew(s);
	Eigen::Matrix<double,3,3> P = Eigen::Matrix<double,3,3>::Identity(3,3) - s*s.transpose();
	Eigen::Matrix<double,3,6> L = Eigen::Matrix<double,3,6>::Zero(3,6);
	L.block(0,0,3,3) = -(1/(cam_T_object.p.Norm()))*P;
	L.block(0,3,3,3) = skew_s;
	s_d << 0,0,1;
	
	cbf visioncbf;
	visioncbf.h = -4*0.5*(s-s_d).dot((s-s_d));

	visioncbf.dhdq = -4*(s-s_d).transpose()*L*J_cam_sp.data;

	return visioncbf;

}

void KUKA_CONTROL::cbf_joint_limits(vector<cbf> &cbfjl){

	cbfjl.resize(_Nj);

	int gamma = 4;

	for(int i=0; i<_Nj; i++){
		cbfjl[i].dhdq.resize(_Nj);
		cbfjl[i].dhdq.setZero();
		cbfjl[i].h = gamma*(_q_max[i]*M_PI/180 - _q(i))*(_q(i)-_q_min[i]*M_PI/180)/((_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180 )*(_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180 ));
   		cbfjl[i].dhdq(i) = gamma*(_q_max[i]*M_PI/180  + _q_min[i]*M_PI/180  - 2*_q(i))/((_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180)*(_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180 ));
	}

}

void KUKA_CONTROL::task_set_input(){
	int i=0;
	ros::Rate r(1/dt);
	
	while(ros::ok()){
		//grasp_sm(_taskSet);
		if(scanf("%d",&i)) _start = true; // _taskSet = i;
	}
}

void KUKA_CONTROL::solveQP(qpOASES::SQProblem* solver, const vector<cbf> &cbfs, int taskset, Vector7d& u){
	
	double sec;
	double sec2;
	sec = ros::WallTime::now().toNSec();

	if(!updateQP(solver, cbfs, _gamma, _lu, _ldelta, taskset)){
		ROS_ERROR("Problem solution ERROR");
		
	}else{
		Eigen::Matrix<double, -1, 1> usol;
		
		usol.resize(10);
        qpOASES::real_t xOpt[10];
		
        solver->getPrimalSolution(xOpt);
        for (int i = 0; i < 10; ++i) {
            usol(i) = xOpt[i];
        }
		u << usol(0), usol(1), usol(2), usol(3), usol(4), usol(5), usol(6);
	
	}
	sec2 = ros::WallTime::now().toNSec();
	#if SOL_TIME_EVAL
	sec = sec/(1000*1000); // ms
	sec2 = sec2/(1000*1000); // ms
	_tsol.push_back(sec2-sec);
	//cout << "Solution time [ms]: " << std::fixed << std::setprecision(8) << sec2-sec << std::endl;
	_outFile << std::fixed << std::setprecision(8) << sec2-sec << ",";
	#endif
}

bool KUKA_CONTROL::request_attach_plate(){

  	ros::ServiceClient attach_client = _nh.serviceClient<boeing_gazebo_model_attachment_plugin::Attach>("/gazebo/attach");

	attach_client.waitForExistence();

	boeing_gazebo_model_attachment_plugin::Attach srv;
	srv.request.model_name_1 = "iiwa";
	srv.request.link_name_1  = "wsg_50_gripper_right";
	srv.request.model_name_2 = "plate";
	srv.request.link_name_2  = "plate_grasp_link";

	if (attach_client.call(srv))
		ROS_INFO("Attach call sent successfully.");
	else
		ROS_ERROR("Failed to call /gazebo/attach");

	return srv.response.success;
}

void KUKA_CONTROL::grasp_sm(int &taskSet){

	Eigen::Affine3d plateTf = Eigen::Affine3d::Identity();
	plateTf.translate(_plate_position);  
	plateTf.rotate(Eigen::AngleAxisd(_plate_rpy(0), Eigen::Vector3d::UnitX()));
	plateTf.rotate(Eigen::AngleAxisd(_plate_rpy(1), Eigen::Vector3d::UnitY()));
	plateTf.rotate(Eigen::AngleAxisd(_plate_rpy(2), Eigen::Vector3d::UnitZ()));

	cbf cbf_approach = cbf_approach_grasp_plate(plateTf, _plate_radius + 0.1);
	cbf cbf_grasp = cbf_approach_grasp_plate(plateTf, _plate_radius - 0.07);
	cbf cbf_go2p1 = cbf_goto_point(_p1);
	cbf cbf_go2p2 = cbf_goto_point(_p2);

	if(!_start) return;
	
	switch(_state){
		case 0: {
			taskSet = 6;
			ROS_INFO("Approaching plate...");
			std_msgs::Float64 msg;
			msg.data = -0.5; 
			_gripper_lf_pub.publish(msg);
			msg.data = 0.5;
			_gripper_rf_pub.publish(msg);

			if(cbf_approach.h > -0.025){
				_state = 1;
			}
		}break;
		case 1: {
			ROS_INFO("Grasping plate...");
			taskSet = 7;
			if(cbf_grasp.h > -0.003){
				_state = 2;
			}
		}break;
		case 2: {
			taskSet = 7;
			std_msgs::Float64 msg;
			msg.data = 0.0; 
			_gripper_lf_pub.publish(msg);
			_gripper_rf_pub.publish(msg);
			if(request_attach_plate()){
				ROS_INFO("Attach!");
				_state = 3;
			}
		}break;
		case 3: {
			ROS_INFO( "Moving to P1!");
			taskSet = 3; 
			if(cbf_go2p1.h > -0.015){
				ROS_INFO( "Target point P1 reached, moving to P2!");
				_state = 4;
			}
		}break;
		case 4: {
			ROS_INFO( "Moving to P2!");
			taskSet = 5; 
			if(cbf_go2p2.h > -0.015){
				ROS_INFO( "Target point P2 reached!");
				_state = 5;
			}
		}break;
		case 5: {
			ROS_INFO( "Task completed!");
			taskSet = 5;
		}break;
	}

}

void KUKA_CONTROL::taskstack(int & taskSet, vector<cbf> & cbfs, std_msgs::Float64MultiArray & h){
	
	cbf c01, c02;
	cbf c1;
	cbf c2;
	cbf c3;

	cbfs = {};
	cbf_joint_limits(cbfs); 

	c01.h = 0;
	c01.dhdq.resize(_Nj);
	c01.dhdq.setZero();
	c02.h = 0;
	c02.dhdq.resize(_Nj);
	c02.dhdq.setZero();

	Eigen::Affine3d plateTf = Eigen::Affine3d::Identity();
	plateTf.translate(_plate_position);  
	plateTf.rotate(Eigen::AngleAxisd(_plate_rpy(0), Eigen::Vector3d::UnitX()));
	plateTf.rotate(Eigen::AngleAxisd(_plate_rpy(1), Eigen::Vector3d::UnitY()));
	plateTf.rotate(Eigen::AngleAxisd(_plate_rpy(2), Eigen::Vector3d::UnitZ()));

	Affine3d _Tf_gripper_lf = kdlFrameToEigenAffine(_T_gripper_lf);
	Affine3d _Tf_gripper_rf = kdlFrameToEigenAffine(_T_gripper_rf);

	switch(taskSet){
		case 1:    //  eePoint < eeHoriz  //  eePoint < vision   
		{	
			cbfs.push_back(c01);
			cbfs.push_back(c02);

			c1 = cbf_goto_point(_p1);
			c2 = cbf_keep_ee_x_horizontal();//cbf_vision();
			c3.h = 0;
			c3.dhdq.resize(_Nj);
			c3.dhdq.setZero();
			cbfs.push_back(c1);
			cbfs.push_back(c2);    
			cbfs.push_back(c3);  
			
			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		}break;
		case 2:   //  eePoint < 3linkz < vision
		{
			
			cbfs.push_back(c01);
			cbfs.push_back(c02);

			c1 = cbf_goto_point(_p1);
			c2 = cbf_goto_z_3link(0.45);
			c3 = cbf_vision();
			cbfs.push_back(c1);
			cbfs.push_back(c2);
			cbfs.push_back(c3);

			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		}break;
		case 3:   //  eePoint1 < eePoint2 < vision
		{
			cbfs.push_back(c01);
			cbfs.push_back(c02);

			c1 = cbf_goto_point(_p1);
			c2 = cbf_keep_ee_x_horizontal();
			c3 = cbf_goto_point(_p2);
			
			cbfs.push_back(c1);
			cbfs.push_back(c2);
			cbfs.push_back(c3);

			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		}break;
		case 4:   // eePoint2 < eePoint1 < vision
		{
			cbfs.push_back(c01);
			cbfs.push_back(c02);
			
			c1 = cbf_goto_point(_p1);
			c2 = cbf_goto_point(_p2);
			c3 = cbf_keep_ee_x_horizontal();

			cbfs.push_back(c1);
			cbfs.push_back(c2);
			cbfs.push_back(c3);					

			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		}break;
		case 5:   // eeHoriz < eePoint2
		{	
			cbfs.push_back(c01);
			cbfs.push_back(c02);

			c1 = cbf_keep_ee_x_horizontal();
			c2 = cbf_goto_point(_p2);
			c3 = cbf_goto_point(_p1);
			// c3.dhdq.resize(_Nj);
			// c3.dhdq.setZero();

			cbfs.push_back(c1);
			cbfs.push_back(c2);    
			cbfs.push_back(c3);  
			
			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		}break;
		case 6:
		{	
			
			c01 = cbf_avoid_plate(_Tf_gripper_lf, _J_gripper_lf, plateTf, Eigen::Vector3d(_plate_radius, _plate_radius, 0.01));
			c02 = cbf_avoid_plate(_Tf_gripper_rf, _J_gripper_rf, plateTf, Eigen::Vector3d(_plate_radius, _plate_radius, 0.01));
			cbfs.push_back(c01);
			cbfs.push_back(c02);

			c1 = cbf_approach_grasp_plate(plateTf, _plate_radius + 0.1);

			c2.h = 0;
			c2.dhdq.resize(_Nj);
			c2.dhdq.setZero();

			c3.h = 0;
			c3.dhdq.resize(_Nj);
			c3.dhdq.setZero();

			cbfs.push_back(c1);
			cbfs.push_back(c2);    
			cbfs.push_back(c3);  
			
			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		}break;
		case 7:
		{
			c01 = cbf_avoid_plate(_Tf_gripper_lf, _J_gripper_lf, plateTf, Eigen::Vector3d(_plate_radius, _plate_radius,  0.01));
			c02 = cbf_avoid_plate(_Tf_gripper_rf, _J_gripper_rf, plateTf, Eigen::Vector3d(_plate_radius, _plate_radius,  0.01));
			cbfs.push_back(c01);
			cbfs.push_back(c02);

			c1 = cbf_approach_grasp_plate(plateTf, _plate_radius - 0.07);
			c2 = cbf_goto_point(_p1);

			c3.h = 0;
			c3.dhdq.resize(_Nj);
			c3.dhdq.setZero();

			cbfs.push_back(c1);
			cbfs.push_back(c2);    
			cbfs.push_back(c3);  
			
			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		}break;
		default:
		
			cbfs.push_back(c01);
			cbfs.push_back(c02);
			c1.h = 0;
			c1.dhdq.resize(_Nj);
			c1.dhdq.setZero();
			cbfs.push_back(c1);
			cbfs.push_back(c1);
			cbfs.push_back(c1);	
			
			
			h.data[0] = c01.h;
			h.data[1] = c02.h;
			h.data[2] = c1.h;
			h.data[3] = c2.h;
			h.data[4] = c3.h;
		break;
	}
}

void KUKA_CONTROL::ctrl_loop(){
	ros::Rate r(1/dt);

	vector<std_msgs::Float64> jcmd;
	vector<std_msgs::Float64> vcmd;
	std_msgs::Float64MultiArray jv_command;
	std_msgs::Float64MultiArray jp_command;

	jcmd.resize(7);
	vcmd.resize(7);
  	jv_command.data.resize(14);
	jp_command.data.resize(14);
	
	while(!_first_js){
		ros::spinOnce();
		r.sleep();
	}

	//Init q
	_q = _qfb;
	

	double mean_st = 0;
	int iter = 0;
	int taskSet = 0;
	int prvtaskSet = 0;
	float sigma = 0;
	bool stackChanged = false;
	bool vision_on = false;

	Vector7d u  = Vector7d::Zero(7,1);
	Vector7d u1 = Vector7d::Zero(7,1);
	Vector7d u2 = Vector7d::Zero(7,1);
	Vector7d qd = Vector7d::Zero(7,1);

	
	vector<cbf> cbfs_old = {};
  	vector<cbf> cbfs = {};

	// CBFs msg for logging
	std_msgs::Float64MultiArray h;
	h.data.resize(8);
	
	taskstack(taskSet,cbfs,h);

	// Initialize solvers
	if(!initQP(_solver,_lu, _ldelta, _gamma, cbfs)){
		ROS_ERROR("Problem init failed");
		exit(1);
	}
	if(!initQP(_solverprv,_lu, _ldelta, _gamma, cbfs)){
		ROS_ERROR("Problem init failed");
		exit(1);
	}
	qd = _qfb;

	while( ros::ok() ) {

		//_q = _qfb;                                    
		_q = qd;
		update_dirkin(_q);

		// Desired tasks set
		grasp_sm(_taskSet);


		if(taskSet != _taskSet){
			stackChanged = true;
			prvtaskSet = taskSet;
			taskSet = _taskSet;
			sigma = 0;
		}
		
		
		taskstack(taskSet,cbfs,h);
 
//************************************************************************
		if(stackChanged){
			taskstack(prvtaskSet,cbfs_old,h);
			solveQP(_solverprv, cbfs_old, prvtaskSet, u1);
			solveQP(_solver, cbfs, taskSet, u2);

			// Interpolation
			
			u = u1*(1-sigma) + sigma*u2; 
			if(sigma < 1) sigma = sigma + 0.002f; //0.002f
			if(sigma > 1) sigma = 1;

		}else{
			solveQP(_solver, cbfs, taskSet, u);
		}
		if(sigma == 1) stackChanged = false;
//************************************************************************

		// // Velocity commands
		for(int k=0; k <_Nj; k++){

			// // Velocity commands
			//vcmd[k].data = u(k);

			// if(taskSet == 10){
			// 	ROS_INFO("Home position");
			// 	vcmd[k].data = -_kp*_qfb[k];
			// }
			
			//_j_pub[k].publish(vcmd[k]);
			jv_command.data[k] = qd[k];
			jv_command.data[7+k] = u(k);
		}

		Eigen::Vector2d bounds = computeBounds(cbfs, _gamma, _ldelta);
		std::cout<<"Bounds: "<<bounds.transpose()<<std::endl;
		//Update qd

		qd = qd + u*dt;  

		// Position commands
		for(int k=0; k <_Nj; k++){
	
			jcmd[k].data = qd(k);
			if(taskSet == 10){
				ROS_INFO("Home position");
				jcmd[k].data = 0;
			}
			_j_pub[k].publish(jcmd[k]);
		}

		h.data[5] = taskSet;
		h.data[6] = bounds(0);
		h.data[7] = bounds(1);

		_h_pub.publish(h);
		_jv_pub.publish(jv_command);

		#if SOL_TIME_EVAL
		auto count = static_cast<float>(_tsol.size());
		double mean = std::accumulate(_tsol.begin(), _tsol.end(),0.0) / count;

		double variance = 0.0;
		for (double value : _tsol) {
			variance += (value - mean) * (value - mean);
		}
		variance /= _tsol.size(); 
		double stddev = std::sqrt(variance);

		// std::cout << "Mean solution time [ms]: " << mean << std::endl;
		// std::cout << "Variance: " << stddev << std::endl;

		#endif
		iter++;

		r.sleep();
		ros::spinOnce();
	}

}

KUKA_CONTROL::~KUKA_CONTROL(){
	delete []	_q_in;
	delete []	_q_3in;

	delete []	_fksolver;
	delete []	_fk3solver;
	delete []	_fk_gripper_solver;

	delete []	_J_solver;
	delete []	_J3_solver;
	delete []	_J_gripper_solver;
	delete []	_J_gripper_lf_solver;
	delete []	_J_gripper_rf_solver;
	delete []	_q_gripper;
	delete []	_q_gripper_lf;
	delete []	_q_gripper_rf;
	delete []	_dyn_param;

	_outFile.close();
}

void KUKA_CONTROL::run() {
	cout << "Insert task set (1-5) or 0 for no task, 10 to reset to home position:\n"
	     << "1: eePoint < eeHoriz  < vision \n"
		 << "2: eePoint < 3linkz < vision \n"
		 << "3: eePoint1 < eePoint2 < vision \n"
		 << "4: eePoint2 < eePoint1 < vision \n"
		 << "5: eeHoriz < eePoint2 \n";
	//boost::thread ctrl_loop_t ( &KUKA_CONTROL::ctrl_loop, this);
	boost::thread task_set_t( &KUKA_CONTROL::task_set_input, this);

	//ros::spin();
}


int main(int argc, char** argv) {
	
	ros::init(argc, argv, "iiwa_control");

	KUKA_CONTROL kc;
	kc.run();
	kc.ctrl_loop();

	return 0;
}
