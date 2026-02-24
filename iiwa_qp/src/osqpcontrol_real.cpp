#include "ros/ros.h"
#include <ros/package.h>
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/Pose.h"
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Int8.h>

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

//#include "gurobi_c++.h"
#include <iiwa_qp/lbr_iiwa_7_robot.h>
#include <vector>
#include <cstdlib>

#include "qp_utils.hpp"


#define STACK 1
#define DQPINTRP 1
//#define NOINTERP 1

#define dt 0.003

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

		void update_dirkin(Eigen::VectorXd &q_in); //, Eigen::Matrix<double,7,1> &u);
		KDL::Jacobian Jacobian(KDL::JntArray* q_in);
		KDL::Jacobian Jacobian3(KDL::JntArray* q_in);
		
		void ctrl_loop();
		void task_set_input();

		cbf  cbf_goto_point(const Eigen::Vector3d &dpe);
		cbf  cbf_goto_z_3link(const float &d3z);
		cbf  cbf_vision();
		void cbf_joint_limits(vector<cbf> &cbfjl);
		void taskstack(int & taskSet, vector<cbf> & cbfs, std_msgs::Float64MultiArray & h);
		void taskstack(int & taskSet, vector<cbf> & cbfs);
		void solveQP(OsqpEigen::Solver& solver, const vector<cbf> &cbfs, int taskset, Vector7d& u);

		void log_data();

	private:

		ros::NodeHandle _nh;
		ros::Subscriber _js_sub, _tag_sub, _task_sub;
		ros::Publisher _js_pub;
		ros::Publisher _u_pub, _h_pub, _jv_pub;
		vector<ros::Publisher> _j_pub;

		bool _first_js;
		bool _first_fk;
		bool _sync;
		bool _got_tag;

		geometry_msgs::Pose _ptag;

		KDL::JntArray *_q_in;
		KDL::JntArray *_q_3in;
		KDL::JntArray *_dq_in;
		KDL::JntArray *_dq_3in;
		KDL::JntArray *_q_kdl;
		KDL::JntArray *_dq_kdl;

		KDL::Tree iiwa_tree;
		KDL::ChainFkSolverPos_recursive *_fksolver; //Forward position solver
		KDL::ChainFkSolverPos_recursive *_fk3solver; //Forward 3link position solver

		KDL::ChainJntToJacSolver *_J_solver;
		KDL::ChainJntToJacSolver *_J3_solver;
		KDL::ChainDynParam *_dyn_param;
		KDL::Chain _k_chain;
		KDL::Chain _k3_chain;

		KDL::Frame _Te;
		KDL::Frame _T3;
		KDL::FrameVel _dirkin_out;
		KDL::FrameVel _dirkin_out3;


		Eigen::VectorXd _vel; 
		Eigen::MatrixXd _J;
		Eigen::MatrixXd _J3;

		OsqpEigen::Solver _solver;
		OsqpEigen::Solver _solverprv;

		int _taskSet;
		float _gamma;
		float _lu;
		float _ldelta;
		float _lk;
		float _k;

		int _Nj; // number of joints
		vector<float> _q_min;
		vector<float> _q_max;
		Eigen::VectorXd _q;
		Eigen::VectorXd _qfb;
		Eigen::Vector3d _xe;
		
		ros::Time _latest_tag;

		double _kp;

		vector<Vector7d> q_log;
		vector<Vector7d> qd_log;
		vector<Vector7d> u_log;
		vector<Vector7d> hj_log;
		vector<Vector3d> h_log;
		vector<double> t_log;
		vector<Vector3d> delta_log;
		vector<Vector2d> k_log;
		vector<int> stack_log;
};


KUKA_CONTROL::KUKA_CONTROL() {
	
	if (!init_model()) exit(1);

	ROS_INFO("Robot tree correctly loaded from parameter server!");

	if(set_qp()) ROS_WARN("Using default QP values");

	_js_sub = _nh.subscribe("/iiwa/joint_states", 1, &KUKA_CONTROL::joint_states_cb, this); // REAL robot
	_j_pub.resize(7);
	//_js_sub = _nh.subscribe("/kuka_iiwa/joint_states", 1, &KUKA_CONTROL::joint_states_cb, this); //GAZEBO
	//_js_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/jointsCommand", 0);	
	//_u_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/quad_u", 1);
	_jv_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/jointsVelCommand", 1);
	//_h_pub = _nh.advertise<std_msgs::Float64MultiArray>("/iiwa/h", 1);
	_tag_sub = _nh.subscribe("/kuka/camera/tag_pose", 0, &KUKA_CONTROL::artag_cb, this);  
	//_task_sub = _nh.subscribe("/kuka/task_set", 0, &KUKA_CONTROL::taskset_cb, this);  

	// for(int k=0; k < _Nj; k++){
	// 	_j_pub[k] = _nh.advertise<std_msgs::Float64>("/kuka_iiwa/joint" + to_string(k+1) + "_velocity_controller/command", 0); //rostopic pub /kuka/task_set std_msgs/Int8 '2'
	// }
	
	_q_in = new KDL::JntArray(_Nj);
	_dq_in = new KDL::JntArray(_Nj);
    _q_3in = new KDL::JntArray(3);
	_dq_3in = new KDL::JntArray(3);

	_fksolver = new KDL::ChainFkSolverPos_recursive( _k_chain );
	_J_solver = new KDL::ChainJntToJacSolver( _k_chain );

	_fk3solver = new KDL::ChainFkSolverPos_recursive( _k3_chain );
	_J3_solver = new KDL::ChainJntToJacSolver(_k3_chain);

	_dyn_param = new KDL::ChainDynParam(_k_chain,KDL::Vector(0,0,-9.81));

	_first_js = false;
	_sync = false;
	_got_tag = false;
	_taskSet = 0;
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

	if ( !iiwa_tree.getChain(base_link, tip_link, _k_chain) ) return false;
	if ( !iiwa_tree.getChain(base_link, link3, _k3_chain) ) return false; 

	_Nj = _k_chain.getNrOfJoints();
	cout<<"Number of joints: "<<_Nj<<endl;
	_J.resize(6,_k_chain.getNrOfJoints());
	_J3.resize(6,_k3_chain.getNrOfJoints());
	_J = Eigen::MatrixXd::Zero(6,_k_chain.getNrOfJoints());
	_J3 = Eigen::MatrixXd::Zero(6,_k3_chain.getNrOfJoints());
	_vel.resize(6);

	_q.resize(_Nj);
	_qfb.resize(_Nj);

	_nh.getParam("q_max", _q_max);
	_nh.getParam("q_min", _q_min);
	_nh.getParam("kp", _kp);

	return true;

}

bool KUKA_CONTROL::set_qp(){

	bool dflt = false;

	if(!_nh.getParam("gamma", _gamma) || !_nh.getParam("lu", _lu) ||  !_nh.getParam("ldelta", _ldelta) || !_nh.getParam("lkappa", _lk) || !_nh.getParam("k", _k)){
		ROS_WARN("Using default QP values");	
		// _lu = 1;
		// _gamma = 0.1;
		// _ldelta = 10;
		// _lk = 1;
		// _k = 1;
		_lu = 1;
		_gamma = 0.8;
		_ldelta = 1;
		_lk = 50;
		_k = 100;
		dflt = true;
	}else{
		ROS_INFO("Parameters loaded");
	}
	return dflt;
}

void KUKA_CONTROL::joint_states_cb(const sensor_msgs::JointState& js){

	for(int i=0; i<7; i++ ) {
		_qfb(i) = js.position[i];
		_dq_in->data[i] = js.velocity[i];
	}

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

	_fksolver->JntToCart(*_q_in, _Te);
	
	_xe(0) = _Te.p.x();
	_xe(1) = _Te.p.y();
	_xe(2) = _Te.p.z();



	_fk3solver->JntToCart(*_q_3in, _T3);

	// _J = Jacobian(_q_in).data;
	// _J3 = Jacobian3(_q_3in).data;
	KDL::Jacobian Jac(_k_chain.getNrOfJoints());
	KDL::Jacobian Jac3(_k3_chain.getNrOfJoints());

	if( _J_solver->JntToJac(*_q_in, Jac) != KDL::ChainJntToJacSolver::E_NOERROR )
		cout << "failing in Jacobian computation!" << endl;
	if( _J3_solver->JntToJac(*_q_3in, Jac3) != KDL::ChainJntToJacSolver::E_NOERROR )
		cout << "failing in Jacobian3 computation!" << endl;

	_J = Jac.data;
	_J3 = Jac3.data;

}

KDL::Jacobian KUKA_CONTROL::Jacobian(KDL::JntArray *q_in){

	KDL::Jacobian Jac(_k_chain.getNrOfJoints());

	if( _J_solver->JntToJac(*q_in, Jac) != KDL::ChainJntToJacSolver::E_NOERROR )
		cout << "failing in Jacobian computation!" << endl;

	return Jac;

}

KDL::Jacobian KUKA_CONTROL::Jacobian3(KDL::JntArray *q_in){

	KDL::Jacobian Jac(_k3_chain.getNrOfJoints());

	if( _J3_solver->JntToJac(*q_in, Jac) != KDL::ChainJntToJacSolver::E_NOERROR )
		cout << "failing in Jacobian 3 computation!" << endl;

	return Jac;

}

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

	cbf3z.dhdq << -50*0.5*(_T3.p.z()-d3z)*Eigen::Vector3d(_J3.row(2)),0,0,0,0;

	return cbf3z;
}

cbf KUKA_CONTROL::cbf_vision(){

	// camera
	KDL::Jacobian J_ee;
	KDL::Frame T_ee;
	// camera in end-effector transform
	KDL::Frame ee_T_cam, Ttag;

	ee_T_cam.M = KDL::Rotation::Identity();
	ee_T_cam.M = KDL::Rotation::RotZ(-1.57);

	ee_T_cam.p = KDL::Vector(0.0,0.025,0.04);
	//ee_T_cam.p = KDL::Vector(0.0,0,0.0);	

	// object
	KDL::Frame cam_T_object, T_object;

	// current forward kinematics
	T_ee = _Te;

	KDL::Frame T_cam = T_ee*ee_T_cam;
	//T_cam.p = T_ee.p;
	
	// current jacobians
	J_ee.data = _J;
	KDL::Jacobian J_cam_sp, J_cam;
	KDL::Vector p_ee_cam = T_ee.M*ee_T_cam.p;  
	J_cam.resize(7);
	J_cam_sp.resize(7);

	changeRefPoint(J_ee, p_ee_cam, J_cam);              // camera jacobian
	J_cam = J_ee;                                       // as matlab
	changeBase(J_cam, T_cam.M.Inverse(), J_cam_sp);     // end-effector spatial jacobian // T_M := rotation part 
														// Equivalent to J_cam_sp = blkdiag(T_cam(1:3,1:3)', T_cam(1:3,1:3)')*J_cam;
	
	if(ros::Time::now().toSec() - _latest_tag.toSec() > 10 ) _got_tag = false;

	if(_got_tag){
		cam_T_object.p = (KDL::Vector(_ptag.position.x,_ptag.position.y,_ptag.position.z));
		Eigen::Quaterniond cam_quat_object(_ptag.orientation.w,_ptag.orientation.x,_ptag.orientation.y,_ptag.orientation.z);
		//Eigen::Quaterniond cam_quat_object(0,1,0,0);
		tf::quaternionEigenToKDL(cam_quat_object.normalized(), cam_T_object.M);

	}else{
		cam_T_object.p = T_cam.M.Inverse()*(KDL::Vector(1,0.15,0.5) - T_cam.p);   //cam_p_object = T_cam.M.transpose()*(p_object - p_cam);
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
	visioncbf.h = -20*0.5*(s-s_d).dot((s-s_d));

	visioncbf.dhdq = -20*(s-s_d).transpose()*L*J_cam_sp.data;

	return visioncbf;

}

void KUKA_CONTROL::cbf_joint_limits(vector<cbf> &cbfjl){

	cbfjl.resize(_Nj);

	int gamma = 10;
	
	for(int i=0; i<_Nj; i++){
		cbfjl[i].dhdq.resize(_Nj);
		cbfjl[i].dhdq.setZero();
		cbfjl[i].h = gamma*(_q_max[i]*M_PI/180 - _q(i))*(_q(i)-_q_min[i]*M_PI/180)/((_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180 )*(_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180 ));
   		cbfjl[i].dhdq(i) = gamma*(_q_max[i]*M_PI/180  + _q_min[i]*M_PI/180  - 2*_q(i))/((_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180)*(_q_max[i]*M_PI/180  - _q_min[i]*M_PI/180 ));
	}

}

void KUKA_CONTROL::task_set_input(){
	int i=0;
	while(ros::ok()){
		if( scanf("%d",&i) ) _taskSet = i;
	}
}

void KUKA_CONTROL::solveQP(OsqpEigen::Solver& solver, const vector<cbf> &cbfs, int taskset, Vector7d& u){
	
	if(!updateQP(solver, cbfs, _gamma, _k, taskset)){
		ROS_ERROR("Update failed");
		
	}else{
		Eigen::Matrix<c_float, -1, 1> usol;

		if(solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError){
			ROS_ERROR("Problem solution ERROR");
			
		}else{
			usol = solver.getSolution();
			u << usol(0), usol(1), usol(2), usol(3), usol(4), usol(5), usol(6);
		}
	}
}

void KUKA_CONTROL::taskstack(int & taskSet, vector<cbf> & cbfs){
	
	cbf c1;
	cbf c2;
	cbf c3;

	cbfs = {};
	cbf_joint_limits(cbfs); 

	switch(taskSet){
				case 1:     //  eePoint < vision
				{	
					Eigen::Vector3d dpe(0.3, 0.2, 0.8);
					c1 = cbf_goto_point(dpe);
					c2 = cbf_vision();
					c3.h = 0;
					c3.dhdq.resize(_Nj);
					c3.dhdq.setZero();
					cbfs.push_back(c1);
					cbfs.push_back(c2);    
					cbfs.push_back(c3);  
					
				}break;
				case 2:   //  eePoint < 3linkz < vision
				{
					Eigen::Vector3d dpe(0.3, 0.2, 0.8);
					c1 = cbf_goto_point(dpe);
					c2 = cbf_goto_z_3link(0.45);
					c3 = cbf_vision();
					cbfs.push_back(c1);
					cbfs.push_back(c2);
					cbfs.push_back(c3);

				}break;
				case 3:   //  eePoint1 < eePoint2 < vision
				{
					Eigen::Vector3d dpe(0.3, 0.2, 0.8);
					Eigen::Vector3d dpe2(0.3, -0.2, 0.7);
					c1 = cbf_goto_point(dpe);
					c2 = cbf_goto_point(dpe2);
					c3 = cbf_vision();
					cbfs.push_back(c1);
					cbfs.push_back(c2);
					cbfs.push_back(c3);

				}break;
				case 4:   // eePoint2 < eePoint1 < vision
				{
					Eigen::Vector3d dpe2(0.3, 0.2, 0.8);
					Eigen::Vector3d dpe(0.3, -0.2, 0.7); // 0.2
					c1 = cbf_goto_point(dpe);
					c2 = cbf_goto_point(dpe2);
					c3 = cbf_vision();
					cbfs.push_back(c1);
					cbfs.push_back(c2);
					cbfs.push_back(c3);					

				}break;
				case 5:  //only vision
				{
					c1 = cbf_vision();
					c3.h = 0;
					c3.dhdq.resize(_Nj);
					c3.dhdq.setZero();
					cbfs.push_back(c1);
					cbfs.push_back(c3);    
					cbfs.push_back(c3); 
				}break;
				default:
					c1.h = 0;
					c1.dhdq.resize(_Nj);
					c1.dhdq.setZero();
					cbfs.push_back(c1);
					cbfs.push_back(c1);
					cbfs.push_back(c1);	
				break;
			}
}
void KUKA_CONTROL::taskstack(int & taskSet, vector<cbf> & cbfs, std_msgs::Float64MultiArray & h){
	
	cbf c1;
	cbf c2;
	cbf c3;

	cbfs = {};
	cbf_joint_limits(cbfs); 

	switch(taskSet){
				case 1:     //  eePoint < vision
				{	
					Eigen::Vector3d dpe(0.3, 0.2, 0.8);
					c1 = cbf_goto_point(dpe);
					c2 = cbf_vision();
					c3.h = 0;
					c3.dhdq.resize(_Nj);
					c3.dhdq.setZero();
					cbfs.push_back(c1);
					cbfs.push_back(c2);    
					cbfs.push_back(c3);  
					
					h.data[0] = c1.h;
					h.data[1] = c2.h;
					h.data[2] = 0;
				}break;
				case 2:   //  eePoint < 3linkz < vision
				{
					Eigen::Vector3d dpe(0.3, 0.2, 0.8);
					c1 = cbf_goto_point(dpe);
					c2 = cbf_goto_z_3link(0.45);
					c3 = cbf_vision();
					cbfs.push_back(c1);
					cbfs.push_back(c2);
					cbfs.push_back(c3);
	
					h.data[0] = c1.h;
					h.data[1] = c2.h;
					h.data[2] = c3.h;
				}break;
				case 3:   //  eePoint1 < eePoint2 < vision
				{
					Eigen::Vector3d dpe(0.3, 0.2, 0.8);
					Eigen::Vector3d dpe2(0.3, -0.2, 0.7);
					c1 = cbf_goto_point(dpe);
					c2 = cbf_goto_point(dpe2);
					c3 = cbf_vision();
					cbfs.push_back(c1);
					cbfs.push_back(c2);
					cbfs.push_back(c3);

					h.data[0] = c1.h;
					h.data[1] = c2.h;
					h.data[2] = c3.h;
				}break;
				case 4:   // eePoint2 < eePoint1 < vision
				{
					Eigen::Vector3d dpe2(0.3, 0.2, 0.8);
					Eigen::Vector3d dpe(0.3, -0.2, 0.7); // 0.2
					c1 = cbf_goto_point(dpe);
					c2 = cbf_goto_point(dpe2);
					c3 = cbf_vision();
					cbfs.push_back(c1);
					cbfs.push_back(c2);
					cbfs.push_back(c3);					

					h.data[0] = c1.h;
					h.data[1] = c2.h;
					h.data[2] = c3.h;
				}break;
				default:
					c1.h = 0;
					c1.dhdq.resize(_Nj);
					c1.dhdq.setZero();
					cbfs.push_back(c1);
					cbfs.push_back(c1);
					cbfs.push_back(c1);	
					
					h.data.resize(3);
					h.data[0] = 0;
					h.data[1] = 0;
					h.data[2] = 0;
				break;
			}
}

void KUKA_CONTROL::ctrl_loop(){
	ros::Rate r(1/dt);

	vector<std_msgs::Float64> jcmd;
	vector<std_msgs::Float64> vcmd;
	std_msgs::Float64MultiArray jv_command;
	// std_msgs::Float64MultiArray u_in;
	// u_in.data.resize(7);
	jcmd.resize(7);
	vcmd.resize(7);
  jv_command.data.resize(14);
	
	while(!_first_js){
		ros::spinOnce();
		r.sleep();
	}
	// for(int time = 0; time<5/dt; time++){ 
		
	// 	for(int i=0; i<_Nj; i++){
		
	// 		jcmd[i].data =  -0.8*(_qfb(i));
	// 		_j_pub[i].publish(jcmd[i]);

	// 	}
	// 	r.sleep();
	// 	ros::spinOnce();
	// }
	//Init q
	_q = _qfb;

	double mean_st = 0;
	int iter = 0;
	int taskSet = 0;
	int prvtaskSet = 0;
	float sigma = 0.0f;
	bool stackChanged = false;
	bool vision_on = false;

	Vector7d u  = Vector7d::Zero(7,1);
	Vector7d u1 = Vector7d::Zero(7,1);
	Vector7d u2 = Vector7d::Zero(7,1);
	Vector7d qd = Vector7d::Zero(7,1);

	
	vector<cbf> cbfs_old = {};
  	vector<cbf> cbfs = {};

	std_msgs::Float64MultiArray h;
	h.data.resize(3);
	
	taskstack(taskSet,cbfs);

	double secs;
	// Initialize solvers

	if(!initQP(_solver,_lu, _ldelta, _lk, _gamma, _k, cbfs)){
		ROS_ERROR("Problem init failed");
		exit(1);
	}
	if(!initQP(_solverprv,_lu, _ldelta, _lk, _gamma, _k, cbfs)){
		ROS_ERROR("Problem init failed");
		exit(1);
	}

	qd = _qfb;

	ros::Time start = ros::Time::now();
	int count = 0;

	while( ros::ok()) {

		_q = _qfb;                                    

		update_dirkin(_q);


		if(_taskSet != 0){
			prvtaskSet = taskSet;
			taskSet = _taskSet;
			_taskSet = 0;
			sigma = 0;
			stackChanged = true;		
		}

		// if((count == 1500 || count == 11500|| count == 21500 || count == 28500) && !stackChanged){
		// 	prvtaskSet = taskSet;
		// 	taskSet++;
		// 	sigma = 0;
		// 	stackChanged = true;
		// }

		
		taskstack(prvtaskSet,cbfs_old);
		taskstack(taskSet,cbfs);
	
//************************************************************************
		
		if(stackChanged){

			solveQP(_solverprv, cbfs_old, prvtaskSet, u1);
			solveQP(_solver, cbfs, taskSet, u2);

			// Interpolation
			
			u = u1*(1-sigma) + sigma*u2; 
			if(sigma < 1) sigma = sigma + 0.002f;
			if(sigma > 1) sigma = 1;

		}else{
			solveQP(_solver, cbfs, taskSet, u);
		}
		if(sigma == 1){stackChanged = false;}

//************************************************************************
		for(int k=0; k <_Nj; k++){
			//u_in.data[k] = u(k);
			
			// vcmd[k].data = u(k) + _kp*(qd[k] - _qfb[k]);
			// _j_pub[k].publish(vcmd[k]);
			jv_command.data[k] = qd[k];
			jv_command.data[7+k] = u(k);
		}
		// Save log Data

		secs = (ros::Time::now()-start).toSec();
		t_log.push_back(secs);
		q_log.push_back(_qfb);
		qd_log.push_back(qd);
		u_log.push_back(u);
		Vector3d hcbf = {cbfs[7].h, cbfs[8].h,cbfs[9].h};
		Vector7d hjcbf;
		hjcbf << cbfs[0].h, cbfs[1].h, cbfs[2].h, cbfs[3].h, cbfs[4].h, cbfs[5].h, cbfs[6].h ;
		hj_log.push_back(hjcbf);
		h_log.push_back(hcbf);
		Eigen::Matrix<c_float, -1, 1> usol;
		usol = _solver.getSolution();
		delta_log.push_back({usol(7),usol(8),usol(9)});
		k_log.push_back({usol(10),usol(11)});
		stack_log.push_back(taskSet);

		//Update qd

		qd = qd + u*dt;  

		count++;
		//_h_pub.publish(h);
		_jv_pub.publish(jv_command);
	
		r.sleep();
		ros::spinOnce();
	}

}
void KUKA_CONTROL::log_data(){
	
	std::ofstream file_q, file_qd, file_u, file_h, file_hj, file_d, file_k;

	std::string log_loc;
	log_loc = "/catkin_ws/src/iiwa_qp/log/";

	file_q.open(log_loc + "log_q.txt", std::ios::out);
	file_qd.open(log_loc + "log_qd.txt", std::ios::out);
	file_u.open(log_loc + "log_u.txt", std::ios::out);
	file_h.open(log_loc + "log_h.txt", std::ios::out);
	file_hj.open(log_loc + "log_hj.txt", std::ios::out);
	file_d.open(log_loc + "log_delta.txt", std::ios::out);
	file_k.open(log_loc + "log_kappa.txt", std::ios::out);
  
	for (int i=0; i<t_log.size(); ++i){
    file_q  <<t_log[i]<<','<<q_log[i](0)<<','<<q_log[i](1)<<','<<q_log[i](2)<<','<<q_log[i](3)<<','<<q_log[i](4)<<','<<q_log[i](5)<<','<<q_log[i](6)<<'\n';
		file_qd <<t_log[i]<<','<<qd_log[i](0)<<','<<qd_log[i](1)<<','<<qd_log[i](2)<<','<<qd_log[i](3)<<','<<qd_log[i](4)<<','<<qd_log[i](5)<<','<<qd_log[i](6)<<'\n';
		file_u  <<t_log[i]<<','<<u_log[i](0)<<','<<u_log[i](1)<<','<<u_log[i](2)<<','<<u_log[i](3)<<','<<u_log[i](4)<<','<<u_log[i](5)<<','<<u_log[i](6)<<'\n';
		file_h  <<t_log[i]<<','<<stack_log[i] <<','<< h_log[i](0)<<','<<h_log[i](1)<<','<<h_log[i](2)<<'\n';
		file_hj <<t_log[i]<<','<<hj_log[i](0)<<','<<hj_log[i](1)<<','<<hj_log[i](2)<<','<<hj_log[i](3)<<','<<hj_log[i](4)<<','<<hj_log[i](5)<<','<<hj_log[i](6)<<'\n';
		file_d  <<t_log[i]<<','<<delta_log[i](0)<<','<<delta_log[i](1)<<','<<delta_log[i](2)<<'\n';
		file_k  <<t_log[i]<<','<<k_log[i](0)<<','<<k_log[i](1)<<'\n';
	}
	
	file_q.close();
	file_qd.close();
	file_u.close();
	file_h.close();
	file_hj.close();
	file_d.close();
	file_k.close();
}

KUKA_CONTROL::~KUKA_CONTROL(){
	delete []	_q_in;
	delete []	_dq_in;

	delete []	_q_3in;
	delete []	_dq_3in;

	delete []	_fksolver;
	delete []	_fk3solver;

	delete []	_J_solver;
	delete []	_J3_solver;

}

void KUKA_CONTROL::run() {
	
	//boost::thread ctrl_loop_t ( &KUKA_CONTROL::ctrl_loop, this);
	boost::thread task_set_t( &KUKA_CONTROL::task_set_input, this);

	//ros::spin();
}


int main(int argc, char** argv) {
	
	ros::init(argc, argv, "iiwa_control");

	KUKA_CONTROL kc;
	kc.run();
	kc.ctrl_loop();
	kc.log_data();

	return 0;
}
