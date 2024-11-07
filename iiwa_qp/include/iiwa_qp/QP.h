#include "gurobi_c++.h"
#include "ros/ros.h"
#include <ros/package.h>
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/Pose.h"
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Int8.h>

#include <iiwa_qp/cbf.h>

using namespace std;

class QP{

  public:
    QP();
    ~QP();
    GRBModel init_quadprog(vector<cbf> &c);
		bool solve_quadprog(GRBModel & model, vector<cbf> &, Eigen::Matrix<double,7,1> &u, int taskSet, Eigen::Vector3d &d, Eigen::Vector3d &vk);
		bool update_task_constraints(GRBModel& model,GRBVar* v, GRBConstr* cnstr ,const vector<cbf> &c, const int taskset);
    bool update_priority_constraint(GRBModel& model, GRBVar* v, vector<cbf> &c, int taskset);
    void update_jl_constraints(GRBModel& model, GRBVar* v,GRBConstr* cnstr, const vector<cbf> &c);
		void empty_taskStack(GRBModel& model);
    void emptyStack(GRBModel& model,GRBVar* v);
    void set_opt(float gamma, float lu, float ld, float lk, float k, float e, float k1, float k2);
    
    GRBModel stack2(GRBModel & model, const vector<cbf> &c);
    GRBModel stack3(GRBModel & model,const vector<cbf> &c);

    bool updatestack2(GRBModel & model, const vector<cbf> &c,GRBVar* v);
    bool updatestack3(GRBModel & model, const vector<cbf> &c,GRBVar* v);

  private:
 
    GRBEnv* _env;
    GRBEnv* _env2;
    float _gamma;
    float _lu;
    float _ldelta;
    float _lk;
    float _k;
    float _k1;
    float _k2;
    double _epsi;
    int _taskSet;
    int _NumCon;
    int _NumVar;

};